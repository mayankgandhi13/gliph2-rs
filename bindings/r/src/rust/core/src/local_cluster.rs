//! Length-bucketed Hamming-distance-1 clustering.
//!
//! Instead of comparing every pair in a length bucket (O(n²)), each position
//! `p` is masked in turn and the bucket is sorted by the masked sequence. Two
//! distinct same-length sequences differ at exactly one position iff they
//! compare equal under exactly one mask, so every edge is found once, with no
//! false positives, in O(n·L·log n) per bucket.

use std::cmp::Ordering;
use std::collections::BTreeMap;

use rayon::prelude::*;

use crate::io::SequenceId;

#[derive(Debug, Clone)]
pub struct LocalParams {
    /// Positions this many residues from either end may not be the mismatch
    /// site (the conserved CDR3 stems). 0 allows a mismatch anywhere.
    pub mismatch_trim: usize,
    /// Sequences shorter than this are ignored.
    pub min_length: usize,
}

impl Default for LocalParams {
    fn default() -> Self {
        LocalParams { mismatch_trim: 0, min_length: 1 }
    }
}

/// An edge between two sequences at Hamming distance 1. Always `a < b`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct LocalEdge {
    pub a: SequenceId,
    pub b: SequenceId,
    /// 0-based position of the substitution.
    pub position: u16,
}

/// Groups sequence ids by length. A `BTreeMap` keeps iteration order fixed.
pub fn bucket_by_length(sequences: &[Vec<u8>], min_length: usize) -> BTreeMap<usize, Vec<SequenceId>> {
    let mut buckets: BTreeMap<usize, Vec<SequenceId>> = BTreeMap::new();
    for (i, s) in sequences.iter().enumerate() {
        if s.len() >= min_length {
            buckets.entry(s.len()).or_default().push(i as SequenceId);
        }
    }
    buckets
}

#[inline]
fn cmp_masked(x: &[u8], y: &[u8], p: usize) -> Ordering {
    x[..p].cmp(&y[..p]).then_with(|| x[p + 1..].cmp(&y[p + 1..]))
}

/// Returns all pairs of (unique) sequences at Hamming distance exactly 1,
/// sorted by `(a, b)`. The input should already be deduplicated.
pub fn local_cluster(sequences: &[Vec<u8>], params: &LocalParams) -> Vec<LocalEdge> {
    let buckets = bucket_by_length(sequences, params.min_length.max(1));

    // One task per (length, mask position) so long buckets spread across threads.
    let tasks: Vec<(&[SequenceId], usize)> = buckets
        .iter()
        .filter(|(_, ids)| ids.len() > 1)
        .flat_map(|(&len, ids)| {
            let lo = params.mismatch_trim.min(len);
            let hi = len.saturating_sub(params.mismatch_trim).max(lo);
            (lo..hi).map(move |p| (ids.as_slice(), p))
        })
        .collect();

    let mut edges: Vec<LocalEdge> = tasks
        .par_iter()
        .flat_map_iter(|&(ids, p)| {
            let mut order = ids.to_vec();
            order.sort_unstable_by(|&i, &j| {
                cmp_masked(&sequences[i as usize], &sequences[j as usize], p).then(i.cmp(&j))
            });
            let mut out = Vec::new();
            let mut start = 0;
            while start < order.len() {
                let head = &sequences[order[start] as usize];
                let mut end = start + 1;
                while end < order.len()
                    && cmp_masked(head, &sequences[order[end] as usize], p) == Ordering::Equal
                {
                    end += 1;
                }
                let group = &order[start..end];
                for (k, &a) in group.iter().enumerate() {
                    for &b in &group[k + 1..] {
                        let (a, b) = if a < b { (a, b) } else { (b, a) };
                        out.push(LocalEdge { a, b, position: p as u16 });
                    }
                }
                start = end;
            }
            out
        })
        .collect();

    edges.par_sort_unstable();
    edges
}

#[inline]
pub fn hamming(a: &[u8], b: &[u8]) -> Option<usize> {
    (a.len() == b.len()).then(|| a.iter().zip(b).filter(|(x, y)| x != y).count())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn seqs(v: &[&str]) -> Vec<Vec<u8>> {
        v.iter().map(|s| s.as_bytes().to_vec()).collect()
    }

    /// Brute-force reference for cross-checking.
    fn naive(s: &[Vec<u8>], trim: usize) -> Vec<(u32, u32)> {
        let mut out = vec![];
        for i in 0..s.len() {
            for j in i + 1..s.len() {
                if hamming(&s[i], &s[j]) == Some(1) {
                    let p = s[i].iter().zip(&s[j]).position(|(x, y)| x != y).unwrap();
                    if p >= trim && p < s[i].len().saturating_sub(trim) {
                        out.push((i as u32, j as u32));
                    }
                }
            }
        }
        out
    }

    #[test]
    fn finds_single_substitutions() {
        let s = seqs(&["CASSLGQETQYF", "CASSPGQETQYF", "CASSLGQETQYW", "CASSLGQETQY", "CATTLGQETQYF"]);
        let e = local_cluster(&s, &LocalParams::default());
        let pairs: Vec<_> = e.iter().map(|e| (e.a, e.b, e.position)).collect();
        assert_eq!(pairs, vec![(0, 1, 4), (0, 2, 11)]);
    }

    #[test]
    fn matches_brute_force() {
        // Deterministic pseudo-random sequences over a tiny alphabet so many
        // collide at distance 1.
        let mut x: u64 = 0x9E3779B97F4A7C15;
        let mut s = vec![];
        let mut seen = std::collections::HashSet::new();
        while s.len() < 400 {
            let len = 6 + (x % 3) as usize;
            let mut v = Vec::with_capacity(len);
            for _ in 0..len {
                x ^= x << 13;
                x ^= x >> 7;
                x ^= x << 17;
                v.push(b"ACDG"[(x % 4) as usize]);
            }
            if seen.insert(v.clone()) {
                s.push(v);
            }
        }
        for trim in [0, 2] {
            let got: Vec<_> = local_cluster(&s, &LocalParams { mismatch_trim: trim, min_length: 1 })
                .iter()
                .map(|e| (e.a, e.b))
                .collect();
            assert_eq!(got, naive(&s, trim), "trim={trim}");
        }
    }
}
