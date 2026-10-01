//! Deterministic motif enrichment: one-sided Fisher's exact test of input vs.
//! reference motif frequency.
//!
//! GLIPH2 estimates enrichment by repeatedly subsampling the reference, which
//! is a source of run-to-run variation. Here the test is computed exactly
//! against the full reference table, summing hypergeometric terms in a fixed
//! order, so results are bit-identical across runs and thread counts.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::io::SequenceId;
use crate::motif::{motif_members, Kmer, MotifParams, MotifTable};

/// `ln(n!)` for `n ≤ max`, built once by cumulative summation.
#[derive(Debug, Clone)]
pub struct LogFactorial(Vec<f64>);

impl LogFactorial {
    pub fn new(max: usize) -> Self {
        let mut v = Vec::with_capacity(max + 1);
        let mut acc = 0.0f64;
        v.push(0.0);
        for i in 1..=max {
            acc += (i as f64).ln();
            v.push(acc);
        }
        LogFactorial(v)
    }

    #[inline]
    pub fn get(&self, n: u64) -> f64 {
        self.0[n as usize]
    }

    #[inline]
    fn ln_choose(&self, n: u64, k: u64) -> f64 {
        self.get(n) - self.get(k) - self.get(n - k)
    }
}

/// One-sided (greater) Fisher's exact test for the 2×2 table
/// `[[a, b], [c, d]]`: P(X ≥ a) where X is hypergeometric with the table's
/// margins.
pub fn fisher_greater(a: u64, b: u64, c: u64, d: u64, lf: &LogFactorial) -> f64 {
    if a == 0 {
        return 1.0;
    }
    let n_draw = a + b; // input size
    let k_succ = a + c; // sequences carrying the motif
    let total = a + b + c + d;
    let x_max = n_draw.min(k_succ);
    let ln_denom = lf.ln_choose(total, n_draw);
    let ln_pmf = |x: u64| lf.ln_choose(k_succ, x) + lf.ln_choose(total - k_succ, n_draw - x) - ln_denom;

    // Upper tail from the mode outward is monotone decreasing, so we can stop
    // once terms are negligible. If `a` is below the mode we sum the (shorter,
    // also monotone) lower tail instead and take the complement.
    let mode = ((n_draw + 1) as f64 * (k_succ + 1) as f64 / (total + 2) as f64).floor() as u64;
    if a >= mode {
        let mut sum = 0.0;
        for x in a..=x_max {
            let t = ln_pmf(x).exp();
            sum += t;
            if t < sum * 1e-16 {
                break;
            }
        }
        sum.min(1.0)
    } else {
        let x_min = n_draw.saturating_sub(total - k_succ);
        let mut lower = 0.0;
        let mut x = a - 1;
        loop {
            let t = ln_pmf(x).exp();
            lower += t;
            if x == x_min || t < lower * 1e-16 {
                break;
            }
            x -= 1;
        }
        (1.0 - lower).clamp(0.0, 1.0)
    }
}

#[derive(Debug, Clone)]
pub struct EnrichmentParams {
    /// Minimum number of distinct input sequences carrying the motif.
    pub min_count: u32,
    /// Maximum Fisher p-value.
    pub max_p: f64,
    /// Minimum fold enrichment, per motif length (`(k, fold)`); lengths not
    /// listed use `default_min_fold`. GLIPH2 defaults: 1000/100/10 for k=2/3/4.
    pub min_fold_by_k: Vec<(usize, f64)>,
    pub default_min_fold: f64,
}

impl Default for EnrichmentParams {
    fn default() -> Self {
        EnrichmentParams {
            min_count: 3,
            max_p: 0.001,
            min_fold_by_k: vec![(2, 1000.0), (3, 100.0), (4, 10.0)],
            default_min_fold: 10.0,
        }
    }
}

impl EnrichmentParams {
    fn min_fold(&self, k: usize) -> f64 {
        self.min_fold_by_k
            .iter()
            .find(|(kk, _)| *kk == k)
            .map(|(_, f)| *f)
            .unwrap_or(self.default_min_fold)
    }
}

#[derive(Debug, Clone)]
pub struct MotifResult {
    pub motif: Kmer,
    pub input_count: u32,
    pub reference_count: u32,
    /// `(input_count / n_input) / (reference_count / n_reference)`; infinite
    /// when the motif is absent from the reference.
    pub fold: f64,
    pub p_value: f64,
    /// Sorted ids of input sequences carrying the motif.
    pub members: Vec<SequenceId>,
}

/// Tests every input motif with at least `min_count` carriers and returns the
/// ones passing all thresholds, sorted by motif.
pub fn motif_enrich(
    input: &[Vec<u8>],
    reference: &MotifTable,
    motif_params: &MotifParams,
    params: &EnrichmentParams,
) -> Vec<MotifResult> {
    let members = motif_members(input, motif_params);
    enrich_from_members(members, input.len() as u32, reference, params)
}

pub fn enrich_from_members(
    members: HashMap<Kmer, Vec<SequenceId>>,
    n_input: u32,
    reference: &MotifTable,
    params: &EnrichmentParams,
) -> Vec<MotifResult> {
    let n_ref = reference.n_sequences;
    let lf = LogFactorial::new((n_input as usize) + (n_ref as usize));

    let mut candidates: Vec<(Kmer, Vec<SequenceId>)> = members
        .into_iter()
        .filter(|(_, m)| m.len() as u32 >= params.min_count)
        .collect();
    candidates.sort_unstable_by_key(|(k, _)| *k);

    candidates
        .into_par_iter()
        .filter_map(|(motif, members)| {
            let a = members.len() as u32;
            // A reference count can't exceed the reference size; clamp so a
            // mismatched table never produces an invalid contingency table.
            let c = reference.get(motif).min(n_ref);
            let fold = if c == 0 {
                f64::INFINITY
            } else {
                (a as f64 / n_input as f64) / (c as f64 / n_ref as f64)
            };
            if fold < params.min_fold(motif.len()) {
                return None;
            }
            let p = fisher_greater(a as u64, (n_input - a) as u64, c as u64, (n_ref - c) as u64, &lf);
            (p <= params.max_p).then_some(MotifResult {
                motif,
                input_count: a,
                reference_count: c,
                fold,
                p_value: p,
                members,
            })
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64) -> bool {
        (a - b).abs() <= 1e-10 * b.abs().max(1e-300)
    }

    #[test]
    fn fisher_matches_r() {
        let lf = LogFactorial::new(10_000);
        // fisher.test(matrix(c(a,c,b,d),2), alternative="greater")$p.value
        assert!(close(fisher_greater(3, 1, 1, 3, &lf), 0.24285714285714288));
        assert!(close(fisher_greater(10, 2, 3, 15, &lf), 0.00046518094336290487));
        assert!(close(fisher_greater(1, 9, 11, 3, &lf), 0.99996634809530227));
        assert!(close(fisher_greater(0, 5, 5, 0, &lf), 1.0));
    }

    #[test]
    fn fisher_is_valid_probability_and_monotone() {
        let lf = LogFactorial::new(2_000);
        let mut prev = 1.0;
        for a in 0..=50 {
            let p = fisher_greater(a, 50 - a, 20, 980, &lf);
            assert!((0.0..=1.0).contains(&p));
            assert!(p <= prev + 1e-12);
            prev = p;
        }
    }

    #[test]
    fn enriched_motif_is_found() {
        let mp = MotifParams { k_values: vec![3], trim_start: 0, trim_end: 0 };
        let input: Vec<Vec<u8>> = ["WWWAAA", "WWWCCC", "WWWDDD", "EEEFFF"].iter().map(|s| s.as_bytes().to_vec()).collect();
        let reference: Vec<Vec<u8>> = (0..1000).map(|i| if i % 500 == 0 { b"WWWGGG".to_vec() } else { b"GGGHHH".to_vec() }).collect();
        let table = MotifTable::build(&reference, &mp);
        let ep = EnrichmentParams { min_count: 3, max_p: 0.01, min_fold_by_k: vec![], default_min_fold: 10.0 };
        let res = motif_enrich(&input, &table, &mp, &ep);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].motif.to_string(), "WWW");
        assert_eq!(res[0].members, vec![0, 1, 2]);
        assert_eq!(res[0].reference_count, 2);
    }
}
