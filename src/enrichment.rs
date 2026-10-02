//! Deterministic motif enrichment: one-sided Fisher's exact test of input vs.
//! reference motif frequency.
//!
//! GLIPH2 estimates enrichment by repeatedly subsampling the reference, which
//! is a source of run-to-run variation. Here the test is computed exactly
//! against the full reference table, summing hypergeometric terms in a fixed
//! order, so results are bit-identical across runs and thread counts.

use std::collections::HashSet;

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
    pub fn ln_choose(&self, n: u64, k: u64) -> f64 {
        self.get(n) - self.get(k) - self.get(n - k)
    }
}

/// One-sided (greater) Fisher's exact test for the 2×2 table
/// `[[a, b], [c, d]]`: P(X ≥ a) where X is hypergeometric with the table's
/// margins.
pub fn fisher_greater(a: u64, b: u64, c: u64, d: u64, lf: &LogFactorial) -> f64 {
    hypergeom_upper(a, a + c, b + d, a + b, lf)
}

/// P(X ≥ x) for X ~ Hypergeometric(`m` successes, `n` failures, `k` draws),
/// i.e. R's `phyper(x - 1, m, n, k, lower.tail = FALSE)`. `lf` must cover
/// `m + n`.
pub fn hypergeom_upper(x: u64, m: u64, n: u64, k: u64, lf: &LogFactorial) -> f64 {
    let total = m + n;
    if x == 0 {
        return 1.0;
    }
    let x_max = k.min(m);
    if k > total || x > x_max {
        return 0.0;
    }
    let x_min = k.saturating_sub(n);
    if x <= x_min {
        return 1.0;
    }
    let ln_denom = lf.ln_choose(total, k);
    let ln_pmf = |i: u64| lf.ln_choose(m, i) + lf.ln_choose(n, k - i) - ln_denom;

    // Upper tail from the mode outward is monotone decreasing, so we can stop
    // once terms are negligible. If `x` is below the mode we sum the (also
    // monotone) lower tail instead and take the complement.
    let mode = ((k + 1) as f64 * (m + 1) as f64 / (total + 2) as f64).floor() as u64;
    if x >= mode {
        let mut sum = 0.0;
        for i in x..=x_max {
            let t = ln_pmf(i).exp();
            sum += t;
            if t < sum * 1e-16 {
                break;
            }
        }
        sum.min(1.0)
    } else {
        let mut lower = 0.0;
        let mut i = x - 1;
        loop {
            let t = ln_pmf(i).exp();
            lower += t;
            if i == x_min || t < lower * 1e-16 {
                break;
            }
            i -= 1;
        }
        (1.0 - lower).clamp(0.0, 1.0)
    }
}

#[derive(Debug, Clone)]
pub struct EnrichmentParams {
    /// Minimum number of distinct input sequences carrying the motif.
    pub min_count: u32,
    /// Significance threshold: applied to the raw p-value, or to the
    /// Benjamini–Hochberg q-value when `fdr` is set.
    pub max_p: f64,
    /// Filter on BH-adjusted q-values instead of raw p-values. GLIPH2 uses
    /// raw p-values, which admits many chance motifs on large inputs.
    pub fdr: bool,
    /// Minimum fold enrichment, per number of motif residues (`(k, fold)`);
    /// counts not listed use `default_min_fold`. GLIPH2 defaults:
    /// 1000/100/10 for k=2/3/4.
    pub min_fold_by_k: Vec<(usize, f64)>,
    pub default_min_fold: f64,
}

impl Default for EnrichmentParams {
    fn default() -> Self {
        EnrichmentParams {
            min_count: 3,
            max_p: 0.001,
            fdr: false,
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
    /// Benjamini–Hochberg q-value over all tested motifs (those with at
    /// least `min_count` carriers).
    pub q_value: f64,
    /// Sorted ids of input sequences carrying the motif.
    pub members: Vec<SequenceId>,
}

/// Benjamini–Hochberg adjustment. Ties are broken by input order, so the
/// result is deterministic.
pub fn benjamini_hochberg(p: &[f64]) -> Vec<f64> {
    let m = p.len();
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| p[a].total_cmp(&p[b]).then(a.cmp(&b)));
    let mut q = vec![0.0; m];
    let mut running = 1.0f64;
    for (rank, &i) in order.iter().enumerate().rev() {
        running = running.min(p[i] * m as f64 / (rank + 1) as f64);
        q[i] = running;
    }
    q
}

/// Tests every input motif with at least `min_count` carriers and returns the
/// ones passing all thresholds, sorted by motif.
///
/// Motifs are counted first; member lists are collected in a second pass for
/// the significant motifs only, so memory does not scale with the total
/// number of k-mer occurrences.
pub fn motif_enrich(
    input: &[Vec<u8>],
    reference: &MotifTable,
    motif_params: &MotifParams,
    params: &EnrichmentParams,
) -> Vec<MotifResult> {
    let counts = MotifTable::build(input, motif_params);
    let mut results = test_motifs(&counts, reference, params);
    let wanted: HashSet<Kmer> = results.iter().map(|r| r.motif).collect();
    let mut members = motif_members(input, motif_params, Some(&wanted));
    for r in &mut results {
        r.members = members.remove(&r.motif).unwrap_or_default();
    }
    results
}

/// Enrichment statistics for every motif in `input` (no member lists),
/// filtered by the thresholds in `params` and sorted by motif.
pub fn test_motifs(input: &MotifTable, reference: &MotifTable, params: &EnrichmentParams) -> Vec<MotifResult> {
    let n_input = input.n_sequences;
    let n_ref = reference.n_sequences;
    let lf = LogFactorial::new((n_input as usize) + (n_ref as usize));

    let mut candidates: Vec<(Kmer, u32)> = input
        .counts
        .iter()
        .filter(|(_, &a)| a >= params.min_count)
        .map(|(&k, &a)| (k, a))
        .collect();
    candidates.sort_unstable_by_key(|(k, _)| *k);

    let mut tested: Vec<MotifResult> = candidates
        .into_par_iter()
        .map(|(motif, a)| {
            // A reference count can't exceed the reference size; clamp so a
            // mismatched table never produces an invalid contingency table.
            let c = reference.get(motif).min(n_ref);
            let fold = if c == 0 {
                f64::INFINITY
            } else {
                (a as f64 / n_input as f64) / (c as f64 / n_ref as f64)
            };
            let p = fisher_greater(a as u64, (n_input - a) as u64, c as u64, (n_ref - c) as u64, &lf);
            MotifResult { motif, input_count: a, reference_count: c, fold, p_value: p, q_value: f64::NAN, members: Vec::new() }
        })
        .collect();

    let p: Vec<f64> = tested.iter().map(|r| r.p_value).collect();
    for (r, q) in tested.iter_mut().zip(benjamini_hochberg(&p)) {
        r.q_value = q;
    }
    tested.retain(|r| {
        let sig = if params.fdr { r.q_value } else { r.p_value };
        r.fold >= params.min_fold(r.motif.residues()) && sig <= params.max_p
    });
    tested
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
        assert!(close(fisher_greater(1, 9, 11, 3, &lf), 0.999_966_348_095_302_3));
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
        let mp = MotifParams { k_values: vec![3], trim_start: 0, trim_end: 0, gapped: false };
        let input: Vec<Vec<u8>> = ["WWWAAA", "WWWCCC", "WWWDDD", "EEEFFF"].iter().map(|s| s.as_bytes().to_vec()).collect();
        let reference: Vec<Vec<u8>> = (0..1000).map(|i| if i % 500 == 0 { b"WWWGGG".to_vec() } else { b"GGGHHH".to_vec() }).collect();
        let table = MotifTable::build(&reference, &mp);
        let ep = EnrichmentParams { min_count: 3, max_p: 0.01, fdr: false, min_fold_by_k: vec![], default_min_fold: 10.0 };
        let res = motif_enrich(&input, &table, &mp, &ep);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].motif.to_string(), "WWW");
        assert_eq!(res[0].members, vec![0, 1, 2]);
        assert_eq!(res[0].reference_count, 2);
    }

    #[test]
    fn bh_matches_r() {
        // p.adjust(c(0.01, 0.04, 0.03, 0.005, 0.5), "BH")
        let q = benjamini_hochberg(&[0.01, 0.04, 0.03, 0.005, 0.5]);
        let want = [0.025, 0.05, 0.05, 0.025, 0.5];
        for (a, b) in q.iter().zip(want) {
            assert!(close(*a, b), "{q:?}");
        }
    }
}
