//! CDR3 input parsing.
//!
//! Accepts GLIPH2-style tab-separated files. If the first line is a header,
//! columns are located by name (`CDR3b`, `TRBV`, `TRBJ`, `subject:condition`,
//! `count`); otherwise the first column is taken as CDR3β and the rest follow
//! the GLIPH2 column order.

use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader};
use std::path::Path;

/// Index into the deduplicated sequence table.
pub type SequenceId = u32;

/// One input row.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TcrRecord {
    pub cdr3: Vec<u8>,
    pub v_gene: Option<String>,
    pub j_gene: Option<String>,
    pub subject: Option<String>,
    pub count: u32,
}

/// Unique CDR3 sequences, in first-seen order, plus the rows that map to each.
#[derive(Debug, Clone, Default)]
pub struct Repertoire {
    pub sequences: Vec<Vec<u8>>,
    /// `records[i]` are the input rows whose CDR3 is `sequences[i]`.
    pub records: Vec<Vec<TcrRecord>>,
}

impl Repertoire {
    pub fn from_records(records: impl IntoIterator<Item = TcrRecord>) -> Self {
        let mut index: HashMap<Vec<u8>, SequenceId> = HashMap::new();
        let mut rep = Repertoire::default();
        for rec in records {
            let id = *index.entry(rec.cdr3.clone()).or_insert_with(|| {
                rep.sequences.push(rec.cdr3.clone());
                rep.records.push(Vec::new());
                (rep.sequences.len() - 1) as SequenceId
            });
            rep.records[id as usize].push(rec);
        }
        rep
    }

    pub fn from_sequences<S: AsRef<[u8]>>(seqs: impl IntoIterator<Item = S>) -> Self {
        Self::from_records(seqs.into_iter().filter_map(|s| {
            normalize_cdr3(s.as_ref()).map(|cdr3| TcrRecord {
                cdr3,
                v_gene: None,
                j_gene: None,
                subject: None,
                count: 1,
            })
        }))
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }
}

/// Uppercases and validates a CDR3. Returns `None` for empty sequences or
/// anything containing characters outside the 20 standard amino acids.
pub fn normalize_cdr3(raw: &[u8]) -> Option<Vec<u8>> {
    let s: Vec<u8> = raw.trim_ascii().to_ascii_uppercase();
    if s.is_empty() || !s.iter().all(|&c| is_amino_acid(c)) {
        return None;
    }
    Some(s)
}

#[inline]
pub fn is_amino_acid(c: u8) -> bool {
    matches!(
        c,
        b'A' | b'C' | b'D' | b'E' | b'F' | b'G' | b'H' | b'I' | b'K' | b'L'
            | b'M' | b'N' | b'P' | b'Q' | b'R' | b'S' | b'T' | b'V' | b'W' | b'Y'
    )
}

struct Columns {
    cdr3: usize,
    v: Option<usize>,
    j: Option<usize>,
    subject: Option<usize>,
    count: Option<usize>,
}

impl Columns {
    const POSITIONAL: Columns = Columns {
        cdr3: 0,
        v: Some(1),
        j: Some(2),
        subject: Some(4),
        count: Some(5),
    };

    fn from_header(fields: &[&str]) -> Option<Columns> {
        let find = |names: &[&str]| {
            fields
                .iter()
                .position(|f| names.iter().any(|n| f.trim().eq_ignore_ascii_case(n)))
        };
        let cdr3 = find(&["CDR3b", "cdr3", "cdr3_b_aa", "junction_aa", "CDR3"])?;
        Some(Columns {
            cdr3,
            v: find(&["TRBV", "v_gene", "v_call", "v_b_gene"]),
            j: find(&["TRBJ", "j_gene", "j_call", "j_b_gene"]),
            subject: find(&["subject:condition", "subject", "patient", "donor"]),
            count: find(&["count", "Freq", "frequency", "duplicate_count"]),
        })
    }
}

/// Reads a TCR table. Rows with invalid CDR3s are skipped; the number skipped
/// is returned alongside the records.
pub fn read_tcr_table(path: impl AsRef<Path>) -> io::Result<(Vec<TcrRecord>, usize)> {
    parse_tcr_table(BufReader::new(File::open(path)?))
}

pub fn parse_tcr_table(reader: impl BufRead) -> io::Result<(Vec<TcrRecord>, usize)> {
    let mut records = Vec::new();
    let mut skipped = 0;
    let mut cols: Option<Columns> = None;

    for (lineno, line) in reader.lines().enumerate() {
        let line = line?;
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if cols.is_none() {
            if let Some(c) = Columns::from_header(&fields) {
                cols = Some(c);
                continue;
            }
            if lineno == 0 && normalize_cdr3(fields[0].as_bytes()).is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "header has no recognisable CDR3 column",
                ));
            }
            cols = Some(Columns::POSITIONAL);
        }
        let c = cols.as_ref().unwrap();
        let get = |i: Option<usize>| {
            i.and_then(|i| fields.get(i))
                .map(|s| s.trim())
                .filter(|s| !s.is_empty() && *s != "NA")
                .map(str::to_owned)
        };
        let Some(cdr3) = fields.get(c.cdr3).and_then(|s| normalize_cdr3(s.as_bytes())) else {
            skipped += 1;
            continue;
        };
        records.push(TcrRecord {
            cdr3,
            v_gene: get(c.v),
            j_gene: get(c.j),
            subject: get(c.subject),
            count: get(c.count).and_then(|s| s.parse().ok()).unwrap_or(1),
        });
    }
    Ok((records, skipped))
}

/// Reads a reference repertoire: any table `read_tcr_table` accepts, or a
/// plain list with one CDR3 per line.
pub fn read_reference(path: impl AsRef<Path>) -> io::Result<Repertoire> {
    let (records, _) = read_tcr_table(path)?;
    Ok(Repertoire::from_records(records))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_header_and_dedups() {
        let input = "CDR3b\tTRBV\tTRBJ\tCDR3a\tsubject:condition\tcount\n\
                     CASSLGQETQYF\tTRBV7-2\tTRBJ2-5\tNA\tP1:A\t3\n\
                     cassLGQETQYF\tTRBV7-2\tTRBJ2-5\tNA\tP2:A\t1\n\
                     CASS*X\tTRBV1\tTRBJ1\tNA\tP1:A\t1\n";
        let (recs, skipped) = parse_tcr_table(input.as_bytes()).unwrap();
        assert_eq!(skipped, 1);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].count, 3);
        assert_eq!(recs[0].subject.as_deref(), Some("P1:A"));
        let rep = Repertoire::from_records(recs);
        assert_eq!(rep.len(), 1);
        assert_eq!(rep.records[0].len(), 2);
    }

    #[test]
    fn parses_headerless_list() {
        let (recs, _) = parse_tcr_table("CASSLGQETQYF\nCASSPGQETQYF\n".as_bytes()).unwrap();
        assert_eq!(recs.len(), 2);
    }
}
