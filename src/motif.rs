//! k-mer (motif) extraction and reference frequency tables.
//!
//! Motifs are taken from the CDR3 interior: GLIPH2 ignores the first and last
//! three residues, which are germline-encoded and shared by almost every
//! sequence. Counts are *sequence* counts — a motif occurring twice in one
//! CDR3 counts once.

use std::collections::HashMap;
use std::fmt;

use rayon::prelude::*;

use crate::io::SequenceId;

/// A packed amino-acid k-mer (k ≤ 6): 5 bits per residue plus the length in
/// the top 3 bits. Packing keeps tables small and hashing cheap.
#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct Kmer(u32);

const AA: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

#[inline]
fn aa_code(c: u8) -> u32 {
    match c {
        b'A' => 0, b'C' => 1, b'D' => 2, b'E' => 3, b'F' => 4, b'G' => 5, b'H' => 6,
        b'I' => 7, b'K' => 8, b'L' => 9, b'M' => 10, b'N' => 11, b'P' => 12, b'Q' => 13,
        b'R' => 14, b'S' => 15, b'T' => 16, b'V' => 17, b'W' => 18, b'Y' => 19,
        _ => 31,
    }
}

impl Kmer {
    pub const MAX_K: usize = 6;

    pub fn new(s: &[u8]) -> Kmer {
        assert!(!s.is_empty() && s.len() <= Self::MAX_K, "k must be 1..=6");
        let mut v = 0u32;
        for &c in s {
            v = (v << 5) | aa_code(c);
        }
        Kmer(v | ((s.len() as u32) << 29))
    }

    pub fn len(self) -> usize {
        (self.0 >> 29) as usize
    }

    pub fn is_empty(self) -> bool {
        self.len() == 0
    }

    pub fn to_bytes(self) -> Vec<u8> {
        let k = self.len();
        (0..k)
            .map(|i| {
                let code = (self.0 >> (5 * (k - 1 - i))) & 31;
                AA.get(code as usize).copied().unwrap_or(b'X')
            })
            .collect()
    }
}

impl fmt::Display for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(std::str::from_utf8(&self.to_bytes()).unwrap())
    }
}

impl fmt::Debug for Kmer {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Kmer({self})")
    }
}

#[derive(Debug, Clone)]
pub struct MotifParams {
    pub k_values: Vec<usize>,
    /// Residues ignored at the N-terminal end.
    pub trim_start: usize,
    /// Residues ignored at the C-terminal end.
    pub trim_end: usize,
}

impl Default for MotifParams {
    fn default() -> Self {
        MotifParams { k_values: vec![2, 3, 4], trim_start: 3, trim_end: 3 }
    }
}

/// Distinct motifs in one sequence, sorted.
pub fn extract_motifs(seq: &[u8], params: &MotifParams) -> Vec<Kmer> {
    let mut out = Vec::new();
    if seq.len() <= params.trim_start + params.trim_end {
        return out;
    }
    let core = &seq[params.trim_start..seq.len() - params.trim_end];
    for &k in &params.k_values {
        if k == 0 || k > Kmer::MAX_K || k > core.len() {
            continue;
        }
        out.extend(core.windows(k).map(Kmer::new));
    }
    out.sort_unstable();
    out.dedup();
    out
}

/// Motif → number of sequences containing it, for a whole repertoire.
#[derive(Debug, Clone, Default)]
pub struct MotifTable {
    pub counts: HashMap<Kmer, u32>,
    /// Number of sequences the table was built from.
    pub n_sequences: u32,
}

impl MotifTable {
    pub fn build(sequences: &[Vec<u8>], params: &MotifParams) -> MotifTable {
        let counts = sequences
            .par_iter()
            .fold(HashMap::new, |mut m: HashMap<Kmer, u32>, s| {
                for k in extract_motifs(s, params) {
                    *m.entry(k).or_insert(0) += 1;
                }
                m
            })
            .reduce(HashMap::new, |a, b| {
                let (mut big, small) = if a.len() >= b.len() { (a, b) } else { (b, a) };
                for (k, v) in small {
                    *big.entry(k).or_insert(0) += v;
                }
                big
            });
        MotifTable { counts, n_sequences: sequences.len() as u32 }
    }

    pub fn get(&self, k: Kmer) -> u32 {
        self.counts.get(&k).copied().unwrap_or(0)
    }
}

/// Motif → sorted ids of the input sequences that contain it.
pub fn motif_members(sequences: &[Vec<u8>], params: &MotifParams) -> HashMap<Kmer, Vec<SequenceId>> {
    let mut pairs: Vec<(Kmer, SequenceId)> = sequences
        .par_iter()
        .enumerate()
        .flat_map_iter(|(i, s)| extract_motifs(s, params).into_iter().map(move |k| (k, i as SequenceId)))
        .collect();
    pairs.par_sort_unstable();
    let mut out: HashMap<Kmer, Vec<SequenceId>> = HashMap::new();
    for (k, id) in pairs {
        out.entry(k).or_default().push(id);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn kmer_roundtrip() {
        for s in ["A", "SL", "GQE", "WYVT", "CASSLG"] {
            let k = Kmer::new(s.as_bytes());
            assert_eq!(k.to_string(), s);
            assert_eq!(k.len(), s.len());
        }
        assert_ne!(Kmer::new(b"A"), Kmer::new(b"AA"));
    }

    #[test]
    fn extracts_interior_only() {
        let p = MotifParams { k_values: vec![2], trim_start: 3, trim_end: 3 };
        let m: Vec<String> = extract_motifs(b"CASSLGQYF", &p).iter().map(|k| k.to_string()).collect();
        // core = "SLG"
        assert_eq!(m, vec!["LG", "SL"]);
    }

    #[test]
    fn counts_once_per_sequence() {
        let p = MotifParams { k_values: vec![2], trim_start: 0, trim_end: 0 };
        let t = MotifTable::build(&[b"SLSL".to_vec(), b"SLAA".to_vec()], &p);
        assert_eq!(t.get(Kmer::new(b"SL")), 2);
        assert_eq!(t.get(Kmer::new(b"LS")), 1);
        assert_eq!(t.n_sequences, 2);
    }
}
