//! CDR3 input parsing.
//!
//! Accepts GLIPH2-style tab-separated files. If the first line is a header,
//! columns are located by name (`CDR3b`/`cdr3aa`, `TRBV`/`v`, `TRBJ`/`j`,
//! `subject:condition`, `count`); otherwise the first column is taken as
//! CDR3β and the rest follow the GLIPH2 column order.
//!
//! Files are parsed in parallel and stored column-wise (interned gene and
//! subject names, one compact [`Row`] per input line) so pooled cohorts with
//! 10⁸ rows fit comfortably in memory.

use std::borrow::Cow;
use std::collections::HashMap;
use std::io;
use std::path::Path;

use rayon::prelude::*;

/// Index into the deduplicated sequence table.
pub type SequenceId = u32;

/// Marks a missing interned field (V, J, subject).
pub const NA: u32 = u32::MAX;

/// Interned strings, ids assigned in first-seen order.
#[derive(Debug, Clone, Default)]
pub struct StringPool {
    pub names: Vec<String>,
    index: HashMap<String, u32>,
}

impl StringPool {
    pub fn intern(&mut self, s: &str) -> u32 {
        if let Some(&id) = self.index.get(s) {
            return id;
        }
        let id = self.names.len() as u32;
        self.names.push(s.to_owned());
        self.index.insert(s.to_owned(), id);
        id
    }

    pub fn get(&self, id: u32) -> Option<&str> {
        self.names.get(id as usize).map(String::as_str)
    }
}

/// One input line, after deduplication of its CDR3.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Row {
    pub seq: SequenceId,
    pub v_gene: u32,
    pub j_gene: u32,
    pub subject: u32,
    pub count: u32,
}

/// Convenience input record for building a repertoire in code.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct TcrRecord {
    pub cdr3: Vec<u8>,
    pub v_gene: Option<String>,
    pub j_gene: Option<String>,
    pub subject: Option<String>,
    pub count: u32,
}

/// Unique CDR3 sequences (in first-seen order) plus every input row.
#[derive(Debug, Clone, Default)]
pub struct Repertoire {
    pub sequences: Vec<Vec<u8>>,
    /// Input rows grouped by sequence: rows of `sequences[i]` are
    /// `rows[row_offsets[i]..row_offsets[i + 1]]`, in input order.
    pub rows: Vec<Row>,
    pub row_offsets: Vec<usize>,
    pub v_genes: StringPool,
    pub j_genes: StringPool,
    pub subjects: StringPool,
    /// Input lines dropped because their CDR3 was empty or not amino acids.
    pub skipped: usize,
}

impl Repertoire {
    pub fn from_records(records: impl IntoIterator<Item = TcrRecord>) -> Self {
        let mut rep = Repertoire::default();
        let mut parsed = Vec::new();
        for r in records {
            let intern = |pool: &mut StringPool, s: &Option<String>| s.as_deref().map_or(NA, |s| pool.intern(s));
            let fields = RowFields {
                v: intern(&mut rep.v_genes, &r.v_gene),
                j: intern(&mut rep.j_genes, &r.j_gene),
                subject: intern(&mut rep.subjects, &r.subject),
                count: r.count,
            };
            parsed.push((Cow::Owned(r.cdr3), fields));
        }
        rep.assign_sequences(parsed);
        rep
    }

    pub fn from_sequences<S: AsRef<[u8]>>(seqs: impl IntoIterator<Item = S>) -> Self {
        let mut skipped = 0;
        let mut rep = Self::from_records(seqs.into_iter().filter_map(|s| {
            let cdr3 = normalize_cdr3(s.as_ref());
            if cdr3.is_none() {
                skipped += 1;
            }
            Some(TcrRecord { cdr3: cdr3?.into_owned(), v_gene: None, j_gene: None, subject: None, count: 1 })
        }));
        rep.skipped = skipped;
        rep
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    pub fn rows_of(&self, id: SequenceId) -> &[Row] {
        &self.rows[self.row_offsets[id as usize]..self.row_offsets[id as usize + 1]]
    }

    /// Deduplicates CDR3s, numbering them in first-seen order, and groups rows
    /// by sequence. Sort-based so it parallelises and is deterministic.
    fn assign_sequences(&mut self, parsed: Vec<(Cow<'_, [u8]>, RowFields)>) {
        let mut order: Vec<u32> = (0..parsed.len() as u32).collect();
        order.par_sort_unstable_by(|&a, &b| parsed[a as usize].0.cmp(&parsed[b as usize].0).then(a.cmp(&b)));

        // Each run of equal CDR3s starts with its first occurrence.
        let mut seq_of_row = vec![0u32; parsed.len()];
        let mut firsts: Vec<u32> = Vec::new();
        let mut i = 0;
        while i < order.len() {
            let head = order[i];
            let mut j = i;
            while j < order.len() && parsed[order[j] as usize].0 == parsed[head as usize].0 {
                seq_of_row[order[j] as usize] = firsts.len() as u32;
                j += 1;
            }
            firsts.push(head);
            i = j;
        }
        // Renumber groups by first occurrence.
        let mut by_first: Vec<u32> = (0..firsts.len() as u32).collect();
        by_first.par_sort_unstable_by_key(|&g| firsts[g as usize]);
        let mut rank = vec![0u32; firsts.len()];
        for (new, &g) in by_first.iter().enumerate() {
            rank[g as usize] = new as u32;
        }
        self.sequences = by_first.par_iter().map(|&g| parsed[firsts[g as usize] as usize].0.to_vec()).collect();

        let mut rows: Vec<Row> = parsed
            .par_iter()
            .enumerate()
            .map(|(r, (_, f))| Row {
                seq: rank[seq_of_row[r] as usize],
                v_gene: f.v,
                j_gene: f.j,
                subject: f.subject,
                count: f.count,
            })
            .collect();
        rows.par_sort_by_key(|r| r.seq); // stable: keeps input order within a sequence
        let mut offsets = vec![0usize; self.sequences.len() + 1];
        for r in &rows {
            offsets[r.seq as usize + 1] += 1;
        }
        for k in 1..offsets.len() {
            offsets[k] += offsets[k - 1];
        }
        self.rows = rows;
        self.row_offsets = offsets;
    }
}

#[derive(Debug, Clone, Copy)]
struct RowFields {
    v: u32,
    j: u32,
    subject: u32,
    count: u32,
}

/// Uppercases and validates a CDR3, borrowing when no change is needed.
/// Returns `None` for empty sequences or anything containing characters
/// outside the 20 standard amino acids.
pub fn normalize_cdr3(raw: &[u8]) -> Option<Cow<'_, [u8]>> {
    let s = raw.trim_ascii();
    if s.is_empty() {
        return None;
    }
    if s.iter().all(|&c| is_amino_acid(c)) {
        return Some(Cow::Borrowed(s));
    }
    let up = s.to_ascii_uppercase();
    up.iter().all(|&c| is_amino_acid(c)).then_some(Cow::Owned(up))
}

#[inline]
pub fn is_amino_acid(c: u8) -> bool {
    matches!(
        c,
        b'A' | b'C' | b'D' | b'E' | b'F' | b'G' | b'H' | b'I' | b'K' | b'L'
            | b'M' | b'N' | b'P' | b'Q' | b'R' | b'S' | b'T' | b'V' | b'W' | b'Y'
    )
}

#[derive(Debug, Clone, Copy)]
struct Columns {
    cdr3: usize,
    v: Option<usize>,
    j: Option<usize>,
    subject: Option<usize>,
    count: Option<usize>,
}

impl Columns {
    const POSITIONAL: Columns = Columns { cdr3: 0, v: Some(1), j: Some(2), subject: Some(4), count: Some(5) };

    fn from_header(fields: &[&[u8]]) -> Option<Columns> {
        let find = |names: &[&str]| {
            fields
                .iter()
                .position(|f| names.iter().any(|n| f.trim_ascii().eq_ignore_ascii_case(n.as_bytes())))
        };
        let cdr3 = find(&["CDR3b", "cdr3aa", "cdr3", "cdr3_b_aa", "junction_aa", "aminoAcid", "amino_acid"])?;
        Some(Columns {
            cdr3,
            v: find(&["TRBV", "v", "v_gene", "v_call", "v_b_gene", "vGeneName"]),
            j: find(&["TRBJ", "j", "j_gene", "j_call", "j_b_gene", "jGeneName"]),
            subject: find(&["subject:condition", "subject", "patient", "donor", "sample_id"]),
            count: find(&["count", "Freq", "frequency", "duplicate_count", "templates"]),
        })
    }
}

/// Reads a TCR table into a deduplicated [`Repertoire`].
pub fn read_repertoire(path: impl AsRef<Path>) -> io::Result<Repertoire> {
    parse_repertoire(&std::fs::read(path)?)
}

/// Reads only the unique CDR3s of a table (e.g. a reference repertoire).
pub fn read_sequences(path: impl AsRef<Path>) -> io::Result<Vec<Vec<u8>>> {
    Ok(read_repertoire(path)?.sequences)
}

/// A parsed chunk; gene/subject ids index the chunk-local pools.
struct ChunkOut<'a> {
    rows: Vec<(Cow<'a, [u8]>, RowFields)>,
    pools: [StringPool; 3],
    skipped: usize,
}

pub fn parse_repertoire(buf: &[u8]) -> io::Result<Repertoire> {
    let mut body = buf;
    let mut cols = Columns::POSITIONAL;

    // Header detection on the first non-empty, non-comment line.
    loop {
        let (line, rest) = split_line(body);
        let trimmed = line.trim_ascii();
        if trimmed.is_empty() || trimmed.starts_with(b"#") {
            if rest.is_empty() {
                return Ok(Repertoire::default());
            }
            body = rest;
            continue;
        }
        let fields: Vec<&[u8]> = line.split(|&c| c == b'\t').collect();
        if let Some(c) = Columns::from_header(&fields) {
            cols = c;
            body = rest;
        } else if normalize_cdr3(fields[0]).is_none() {
            return Err(io::Error::new(io::ErrorKind::InvalidData, "header has no recognisable CDR3 column"));
        }
        break;
    }

    // Split into ~4 chunks per thread at line boundaries; results are
    // concatenated in chunk order, so the outcome is independent of chunking.
    let n_chunks = (rayon::current_num_threads() * 4).max(1);
    let target = body.len() / n_chunks + 1;
    let mut chunks: Vec<&[u8]> = Vec::with_capacity(n_chunks);
    let mut rest = body;
    while !rest.is_empty() {
        let mut end = target.min(rest.len());
        while end < rest.len() && rest[end - 1] != b'\n' {
            end += 1;
        }
        chunks.push(&rest[..end]);
        rest = &rest[end..];
    }

    let outs: Vec<ChunkOut> = chunks.par_iter().map(|chunk| parse_chunk(chunk, cols)).collect::<io::Result<_>>()?;

    // Merge chunk-local pools in chunk order; since each local pool is in
    // first-seen order, global ids end up in global first-seen order.
    let mut rep = Repertoire::default();
    let mut parsed = Vec::with_capacity(outs.iter().map(|o| o.rows.len()).sum());
    for mut out in outs {
        rep.skipped += out.skipped;
        let globals = [&mut rep.v_genes, &mut rep.j_genes, &mut rep.subjects];
        let maps: Vec<Vec<u32>> = out
            .pools
            .iter()
            .zip(globals)
            .map(|(local, global)| local.names.iter().map(|n| global.intern(n)).collect())
            .collect();
        let remap = |map: &[u32], id: u32| if id == NA { NA } else { map[id as usize] };
        out.rows.par_iter_mut().for_each(|(_, f)| {
            f.v = remap(&maps[0], f.v);
            f.j = remap(&maps[1], f.j);
            f.subject = remap(&maps[2], f.subject);
        });
        parsed.append(&mut out.rows);
    }
    rep.assign_sequences(parsed);
    Ok(rep)
}

fn parse_chunk(chunk: &[u8], c: Columns) -> io::Result<ChunkOut<'_>> {
    let mut out = ChunkOut { rows: Vec::new(), pools: Default::default(), skipped: 0 };
    let mut last = [(&b""[..], NA); 3]; // consecutive rows often repeat a name
    for line in chunk.split(|&b| b == b'\n') {
        let line = line.strip_suffix(b"\r").unwrap_or(line);
        if line.trim_ascii().is_empty() || line.starts_with(b"#") {
            continue;
        }
        let mut fields: [&[u8]; 16] = [&[]; 16];
        for (i, f) in line.split(|&b| b == b'\t').take(16).enumerate() {
            fields[i] = f;
        }
        let Some(cdr3) = normalize_cdr3(fields[c.cdr3.min(15)]) else {
            out.skipped += 1;
            continue;
        };
        let field = |i: Option<usize>| i.map_or(&b""[..], |i| fields[i.min(15)].trim_ascii());
        let mut ids = [NA; 3];
        for (k, col) in [c.v, c.j, c.subject].into_iter().enumerate() {
            let raw = field(col);
            if raw.is_empty() || raw == b"NA" {
                continue;
            }
            if last[k].0 != raw {
                let s = std::str::from_utf8(raw).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
                last[k] = (raw, out.pools[k].intern(s));
            }
            ids[k] = last[k].1;
        }
        let count = std::str::from_utf8(field(c.count))
            .ok()
            .and_then(|s| s.parse::<f64>().ok())
            .filter(|v| *v >= 1.0)
            .map_or(1, |v| v as u32);
        out.rows.push((cdr3, RowFields { v: ids[0], j: ids[1], subject: ids[2], count }));
    }
    Ok(out)
}

fn split_line(buf: &[u8]) -> (&[u8], &[u8]) {
    match buf.iter().position(|&b| b == b'\n') {
        Some(i) => (&buf[..i], &buf[i + 1..]),
        None => (buf, &[]),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_header_and_dedups() {
        let input = "CDR3b\tTRBV\tTRBJ\tCDR3a\tsubject:condition\tcount\n\
                     CASSLGQETQYF\tTRBV7-2\tTRBJ2-5\tNA\tP1:A\t3\n\
                     CASSPGQETQYF\tTRBV7-2\tTRBJ2-5\tNA\tP1:A\t1\n\
                     cassLGQETQYF\tTRBV7-2\tTRBJ2-5\tNA\tP2:A\t1\n\
                     CASS*X\tTRBV1\tTRBJ1\tNA\tP1:A\t1\n";
        let rep = parse_repertoire(input.as_bytes()).unwrap();
        assert_eq!(rep.skipped, 1);
        assert_eq!(rep.sequences, vec![b"CASSLGQETQYF".to_vec(), b"CASSPGQETQYF".to_vec()]);
        let rows = rep.rows_of(0);
        assert_eq!(rows.len(), 2);
        assert_eq!(rows[0].count, 3);
        assert_eq!(rep.subjects.get(rows[0].subject), Some("P1:A"));
        assert_eq!(rep.subjects.get(rows[1].subject), Some("P2:A"));
        assert_eq!(rep.rows_of(1).len(), 1);
    }

    #[test]
    fn parses_vdjtools_columns() {
        let input = "count\tfreq\tcdr3nt\tcdr3aa\tv\td\tj\n\
                     5\t0.1\tTGT\tCASSLGQETQYF\tTRBV7-2\t.\tTRBJ2-5\n";
        let rep = parse_repertoire(input.as_bytes()).unwrap();
        assert_eq!(rep.len(), 1);
        let r = rep.rows_of(0)[0];
        assert_eq!((r.count, rep.v_genes.get(r.v_gene)), (5, Some("TRBV7-2")));
        assert_eq!(r.subject, NA);
    }

    #[test]
    fn parses_headerless_list() {
        let rep = parse_repertoire(b"CASSLGQETQYF\nCASSPGQETQYF\r\n").unwrap();
        assert_eq!(rep.len(), 2);
    }

    #[test]
    fn chunking_does_not_change_result() {
        let mut text = String::from("CDR3b\tTRBV\n");
        for i in 0..5000 {
            text.push_str(&format!("CASS{}F\tTRBV{}\n", ["A", "G", "L", "Q", "W"][i % 5].repeat(1 + i % 7), i % 3));
        }
        let fp = |threads| {
            rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap().install(|| {
                let r = parse_repertoire(text.as_bytes()).unwrap();
                (r.sequences, r.rows, r.v_genes.names)
            })
        };
        assert_eq!(fp(1), fp(7));
    }
}
