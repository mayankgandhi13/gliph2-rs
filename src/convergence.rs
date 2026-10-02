//! GLIPH2 convergence groups, matching turboGliph's `gliph2()`.
//!
//! GLIPH2 does not merge local and global similarities into connected
//! components. It reports *convergence groups*, which may overlap:
//!
//! * **local** — one group per (enriched motif, position range), tagged
//!   `MOTIF_start_stop`. Motifs come from the CDR3 interior (the first and
//!   last `boundary_size` residues are excluded), are counted per occurrence
//!   in sample and reference, and are selected by a hypergeometric test with
//!   p-values rounded to two significant figures and fold change (OvE) to one
//!   decimal, exactly as turboGliph does. Occurrence start positions are
//!   chained while consecutive distinct starts differ by at most
//!   `motif_distance_cutoff - 1`; each chain is one group.
//! * **global** — sequences whose interiors are identical except at one
//!   position, tagged by the masked interior plus the residues seen there
//!   (e.g. `SLG%ET_QR`). With `all_aa_interchangeable = false` only
//!   BLOSUM62-compatible substitutions connect sequences.
//!
//! Every step is deterministic and parallel. Groups are returned local first
//! (by motif, then position range), then global (by CDR3 length, masked
//! position and masked interior), independent of thread count. Global groups
//! are found one (length, position) slice at a time, with the reference
//! sorted into the same slice, so memory stays proportional to the output
//! rather than to every candidate tag.

use std::collections::{BTreeMap, HashMap};

use rayon::prelude::*;

use crate::enrichment::{hypergeom_upper, LogFactorial};
use crate::io::{Repertoire, SequenceId, NA};
use crate::motif::Kmer;

/// Parameters of turboGliph's `gliph2()`; defaults are its code defaults.
#[derive(Debug, Clone)]
pub struct Gliph2Params {
    /// `lcminp`: maximum (rounded) motif p-value.
    pub lcminp: f64,
    /// `lcminove`: minimum fold change, one value or one per `motif_length`.
    pub lcminove: Vec<f64>,
    /// `kmer_mindepth`: minimum motif occurrences in the sample.
    pub kmer_mindepth: u64,
    pub motif_distance_cutoff: usize,
    /// Residues excluded at each end (`structboundaries = TRUE`).
    pub boundary_size: usize,
    pub motif_length: Vec<usize>,
    /// `discontinuous_motifs`: also test motifs with one interior wildcard.
    pub discontinuous: bool,
    /// Raised to `2 * boundary_size + 1` if smaller, as turboGliph does.
    pub min_seq_length: usize,
    /// `accept_sequences_with_C_F_start_end`: keep only `C…F` CDR3s.
    pub require_cf: bool,
    pub all_aa_interchangeable: bool,
    /// `global_vgene`: global edges also require the same TRBV.
    pub global_vgene: bool,
    pub cluster_min_size: usize,
    pub local: bool,
    pub global: bool,
}

impl Default for Gliph2Params {
    fn default() -> Self {
        Gliph2Params {
            lcminp: 0.01,
            lcminove: vec![1000.0, 100.0, 10.0],
            kmer_mindepth: 3,
            motif_distance_cutoff: 3,
            boundary_size: 3,
            motif_length: vec![2, 3, 4],
            discontinuous: false,
            min_seq_length: 0,
            require_cf: true,
            all_aa_interchangeable: false,
            global_vgene: false,
            cluster_min_size: 2,
            local: true,
            global: true,
        }
    }
}

impl Gliph2Params {
    /// The parameter set from the GLIPH2 paper's distributed parameter files
    /// (`local_min_pvalue=0.001`, `local_min_OVE=10`, `kmer_min_depth=3`,
    /// `cdr3_length_cutoff=8`, `all_aa_interchangeable=1`).
    pub fn gliph2_paper() -> Self {
        Gliph2Params {
            lcminp: 0.001,
            lcminove: vec![10.0],
            min_seq_length: 8,
            all_aa_interchangeable: true,
            ..Default::default()
        }
    }

    fn effective_min_length(&self) -> usize {
        self.min_seq_length.max(2 * self.boundary_size + 1)
    }

    fn min_ove(&self, residues: usize) -> f64 {
        if self.lcminove.len() == 1 {
            return self.lcminove[0];
        }
        self.motif_length
            .iter()
            .position(|&k| k == residues)
            .and_then(|i| self.lcminove.get(i).copied())
            .unwrap_or(0.0)
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum GroupType {
    Local,
    Global,
}

impl GroupType {
    pub fn as_str(self) -> &'static str {
        match self {
            GroupType::Local => "local",
            GroupType::Global => "global",
        }
    }
}

/// A group's tag in compact form; [`ConvergenceGroup::tag`] renders it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GroupTag {
    /// `MOTIF_start_stop` (1-based CDR3 positions).
    Local { motif: Kmer, start: u32, stop: u32 },
    /// Masked interior of any member, plus residues at the masked position
    /// (bitmask over `ACDEFGHIKLMNPQRSTVWY`) and, with `global_vgene`, the V gene.
    Global { boundary: u8, pos: u16, residues: u32, v_gene: Option<u32> },
}

#[derive(Debug, Clone)]
pub struct ConvergenceGroup {
    pub kind: GroupType,
    pub tag: GroupTag,
    /// Nodes in the group (unique CDR3s, or CDR3+TRBV with `global_vgene`).
    pub cluster_size: u32,
    pub unique_cdr3_sample: u32,
    pub unique_cdr3_ref: u64,
    /// Whole-motif fold change for local groups, 0 for global (as turboGliph).
    pub ove: f64,
    pub fisher_score: f64,
    /// Member sequence ids, sorted by CDR3 string.
    pub members: Vec<SequenceId>,
}

const AA_ORDER: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";

fn aa_bit(a: u8) -> u32 {
    AA_ORDER.iter().position(|&c| c == a).map_or(0, |i| 1 << i)
}

impl ConvergenceGroup {
    /// The turboGliph tag, e.g. `SLG_4_17` or `SLG%ET_QR`.
    pub fn tag(&self, rep: &Repertoire) -> String {
        match self.tag {
            GroupTag::Local { motif, start, stop } => format!("{}_{start}_{stop}", motif_tag(motif)),
            GroupTag::Global { boundary, pos, residues, v_gene } => {
                let mut t = interior(&rep.sequences[self.members[0] as usize], boundary as usize).to_vec();
                t[pos as usize] = b'%';
                t.push(b'_');
                if let Some(v) = v_gene {
                    t.extend_from_slice(rep.v_genes.get(v).unwrap_or("NA").as_bytes());
                    t.push(b'_');
                }
                t.extend(AA_ORDER.iter().enumerate().filter(|(i, _)| residues & (1 << i) != 0).map(|(_, &c)| c));
                String::from_utf8(t).unwrap()
            }
        }
    }
}

/// One enriched motif before position splitting, for reporting.
#[derive(Debug, Clone)]
pub struct SelectedMotif {
    pub motif: Kmer,
    pub num_in_sample: u64,
    pub num_in_ref: u64,
    pub fisher_score: f64,
    pub num_fold: f64,
}

#[derive(Debug, Clone, Default)]
pub struct Gliph2Result {
    pub groups: Vec<ConvergenceGroup>,
    pub selected_motifs: Vec<SelectedMotif>,
    /// Size of the filtered sample (unique C…F CDR3s of sufficient length).
    pub n_sample: usize,
    pub n_reference: usize,
}

/// `as.numeric(formatC(p, digits = 1, format = "e"))`: two significant figures.
pub fn round_2sig(p: f64) -> f64 {
    if !p.is_finite() || p == 0.0 {
        return p;
    }
    format!("{p:.1e}").parse().unwrap_or(p)
}

/// `round(x, 1)`.
fn round_1dp(x: f64) -> f64 {
    (x * 10.0).round() / 10.0
}

/// BLOSUM62 pairs with score ≥ 0 (turboGliph's `BlosumVec`), as a 20×20
/// table indexed by residue code.
fn blosum_compatible(a: u8, b: u8) -> bool {
    const PAIRS: &str = "AA CA GA SA TA VA RR NR QR ER HR KR RN NN DN QN EN GN HN KN SN TN ND DD QD ED SD AC CC RQ NQ DQ \
        QQ EQ HQ KQ MQ SQ RE NE DE QE EE HE KE SE AG NG GG SG RH NH QH EH HH YH II LI MI FI VI IL LL ML FL VL RK NK QK EK \
        KK SK QM IM LM MM FM VM IF LF MF FF WF YF PP AS NS DS QS ES GS KS SS TS AT NT ST TT VT FW WW YW HY FY WY YY AV IV \
        LV MV TV VV";
    PAIRS.split_ascii_whitespace().any(|p| p.as_bytes() == [a, b])
}

fn is_working(seq: &[u8], p: &Gliph2Params) -> bool {
    seq.len() >= p.effective_min_length() && (!p.require_cf || (seq.first() == Some(&b'C') && seq.last() == Some(&b'F')))
}

fn interior(seq: &[u8], b: usize) -> &[u8] {
    &seq[b..seq.len() - b]
}

/// Motif counts per occurrence (overlapping), like `stringdist::qgrams`.
fn count_motifs(regions: &[&[u8]], p: &Gliph2Params) -> HashMap<Kmer, u64> {
    regions
        .par_iter()
        .fold(HashMap::new, |mut m: HashMap<Kmer, u64>, r| {
            for_each_motif(r, p, |_, k| *m.entry(k).or_insert(0) += 1);
            m
        })
        .reduce(HashMap::new, |a, b| {
            let (mut big, small) = if a.len() >= b.len() { (a, b) } else { (b, a) };
            for (k, v) in small {
                *big.entry(k).or_insert(0) += v;
            }
            big
        })
}

/// Calls `f(start, motif)` for every motif occurrence in `region`, in
/// increasing start order per motif.
fn for_each_motif(region: &[u8], p: &Gliph2Params, mut f: impl FnMut(usize, Kmer)) {
    let mut buf = [0u8; Kmer::MAX_K];
    for &k in &p.motif_length {
        if k == 0 || k > Kmer::MAX_K || k > region.len() {
            continue;
        }
        for (i, w) in region.windows(k).enumerate() {
            f(i, Kmer::new(w));
        }
        if p.discontinuous && (2..Kmer::MAX_K).contains(&k) && k < region.len() {
            for (i, w) in region.windows(k + 1).enumerate() {
                for gap in 1..k {
                    buf[..=k].copy_from_slice(w);
                    buf[gap] = b'%';
                    f(i, Kmer::new(&buf[..=k]));
                }
            }
        }
    }
}

fn motif_tag(m: Kmer) -> String {
    m.to_string().replace('%', ".")
}

/// Runs GLIPH2 convergence-group detection on `rep` against `reference`
/// (unique CDR3s).
pub fn gliph2(rep: &Repertoire, reference: &[Vec<u8>], p: &Gliph2Params) -> Gliph2Result {
    // turboGliph's `seqs`: unique, C…F, long enough. Repertoire sequences
    // are already unique and amino-acid only.
    let working: Vec<SequenceId> = (0..rep.len() as SequenceId)
        .filter(|&i| is_working(&rep.sequences[i as usize], p))
        .collect();
    let mut ref_unique: Vec<&[u8]> = reference.iter().map(Vec::as_slice).filter(|s| is_working(s, p)).collect();
    ref_unique.par_sort_unstable();
    ref_unique.dedup();

    let n_samp = working.len() as u64;
    let n_ref = ref_unique.len() as u64;
    let lf = LogFactorial::new((n_samp + n_ref) as usize + 1);
    let mut result = Gliph2Result { n_sample: working.len(), n_reference: ref_unique.len(), ..Default::default() };

    if p.local && n_samp > 0 {
        let (groups, motifs) = local_groups(rep, &working, &ref_unique, p, &lf);
        result.groups.extend(groups);
        result.selected_motifs = motifs;
    }
    if p.global && n_samp > 0 {
        result.groups.extend(global_groups(rep, &working, &ref_unique, p, &lf));
    }
    result.groups.retain(|g| g.cluster_size as usize >= p.cluster_min_size);
    result
}

fn local_groups(
    rep: &Repertoire,
    working: &[SequenceId],
    ref_unique: &[&[u8]],
    p: &Gliph2Params,
    lf: &LogFactorial,
) -> (Vec<ConvergenceGroup>, Vec<SelectedMotif>) {
    let b = p.boundary_size;
    let regions: Vec<&[u8]> = working.iter().map(|&i| interior(&rep.sequences[i as usize], b)).collect();
    let ref_regions: Vec<&[u8]> = ref_unique.iter().map(|s| interior(s, b)).collect();
    let n_samp = regions.len() as u64;
    let n_ref = ref_regions.len() as u64;

    let sample_counts = count_motifs(&regions, p);
    let ref_counts = count_motifs(&ref_regions, p);
    let fisher = |n_s: u64, n_r: u64| {
        // phyper(n_s - 1, m = n_r + n_s, n = N_ref + N_samp - n_r - n_s, k = N_samp, lower = FALSE)
        match (n_ref + n_samp).checked_sub(n_r + n_s) {
            Some(n) => hypergeom_upper(n_s, n_r + n_s, n, n_samp, lf),
            None => f64::NAN,
        }
    };

    let mut selected: Vec<SelectedMotif> = sample_counts
        .par_iter()
        .filter_map(|(&motif, &n_s)| {
            let n_r = ref_counts.get(&motif).copied().unwrap_or(0);
            let fisher_score = round_2sig(fisher(n_s, n_r));
            let num_fold = round_1dp(n_s as f64 / (n_r as f64 + 0.01) / n_samp as f64 * n_ref as f64);
            (num_fold >= p.min_ove(motif.residues()) && fisher_score <= p.lcminp && n_s >= p.kmer_mindepth)
                .then_some(SelectedMotif { motif, num_in_sample: n_s, num_in_ref: n_r, fisher_score, num_fold })
        })
        .collect();
    selected.sort_unstable_by_key(|m| m.motif);
    if selected.is_empty() {
        return (Vec::new(), selected);
    }

    // Non-overlapping, leftmost occurrence starts (str_locate_all semantics)
    // of each selected motif in each working sequence, in one pass.
    let index: HashMap<Kmer, u32> = selected.iter().enumerate().map(|(i, m)| (m.motif, i as u32)).collect();
    let mut occ: Vec<(u32, u32, u32)> = regions
        .par_iter()
        .enumerate()
        .fold(Vec::new, |mut out: Vec<(u32, u32, u32)>, (w, r)| {
            let mut hits: Vec<(u32, usize)> = Vec::new();
            for_each_motif(r, p, |start, k| {
                if let Some(&mi) = index.get(&k) {
                    hits.push((mi, start));
                }
            });
            hits.sort_unstable();
            let (mut cur, mut free) = (u32::MAX, 0usize);
            for (mi, start) in hits {
                if mi != cur {
                    (cur, free) = (mi, 0);
                }
                if start >= free {
                    free = start + selected[mi as usize].motif.len();
                    out.push((mi, start as u32 + 1, w as u32)); // 1-based like R
                }
            }
            out
        })
        .reduce(Vec::new, |mut a, mut b| {
            a.append(&mut b);
            a
        });
    occ.par_sort_unstable();

    let max_region = regions.iter().map(|r| r.len()).max().unwrap_or(0) as u32;
    let motif_diffs = p.motif_distance_cutoff.saturating_sub(1) as u32;
    let by_motif: Vec<&[(u32, u32, u32)]> = occ.chunk_by(|a, b| a.0 == b.0).collect();

    let groups: Vec<ConvergenceGroup> = by_motif
        .par_iter()
        .flat_map_iter(|hits| {
            let m = &selected[hits[0].0 as usize];
            let mut starts: Vec<u32> = hits.iter().map(|h| h.1).collect();
            starts.dedup(); // hits are sorted by start
            let breaks: Vec<usize> =
                (0..starts.len().saturating_sub(1)).filter(|&j| starts[j + 1] - starts[j] > motif_diffs).collect();

            let mut ranges: Vec<(u32, u32, u32, u32)> = Vec::new(); // (lo, hi, tag_start, tag_stop)
            if p.motif_distance_cutoff < 1 || breaks.is_empty() {
                ranges.push((0, u32::MAX, 1, max_region));
            } else {
                let mut lo = starts[0];
                for j in breaks.into_iter().chain(std::iter::once(starts.len() - 1)) {
                    ranges.push((lo, starts[j], lo, starts[j]));
                    lo = starts.get(j + 1).copied().unwrap_or(starts[j]);
                }
            }
            ranges.into_iter().map(move |(lo, hi, t0, t1)| {
                let mut members: Vec<SequenceId> =
                    hits.iter().filter(|h| h.1 >= lo && h.1 <= hi).map(|h| working[h.2 as usize]).collect();
                members.sort_unstable_by(|&x, &y| rep.sequences[x as usize].cmp(&rep.sequences[y as usize]));
                members.dedup();
                let n = members.len() as u64;
                ConvergenceGroup {
                    kind: GroupType::Local,
                    tag: GroupTag::Local { motif: m.motif, start: t0 + b as u32, stop: t1 + b as u32 },
                    cluster_size: members.len() as u32,
                    unique_cdr3_sample: members.len() as u32,
                    unique_cdr3_ref: m.num_in_ref,
                    ove: m.num_fold,
                    fisher_score: round_2sig(fisher(n, m.num_in_ref)),
                    members,
                }
            })
        })
        .collect();
    (groups, selected)
}

fn global_groups(
    rep: &Repertoire,
    working: &[SequenceId],
    ref_unique: &[&[u8]],
    p: &Gliph2Params,
    lf: &LogFactorial,
) -> Vec<ConvergenceGroup> {
    let b = p.boundary_size;
    let seq = |i: SequenceId| rep.sequences[i as usize].as_slice();
    // Sample and reference interiors bucketed by length. Within a slice,
    // entries are sample ids, or `REF | index` for reference sequences.
    const REF: u32 = 1 << 31;
    type Bucket = (Vec<u32>, Vec<u32>); // (sample ids, REF-flagged reference indices)
    let mut buckets: BTreeMap<usize, Bucket> = BTreeMap::new();
    for &i in working {
        buckets.entry(seq(i).len() - 2 * b).or_default().0.push(i);
    }
    for (j, r) in ref_unique.iter().enumerate() {
        if let Some(bucket) = buckets.get_mut(&(r.len() - 2 * b)) {
            bucket.1.push(j as u32 | REF);
        }
    }
    let tasks: Vec<(&Bucket, usize)> = buckets
        .iter()
        .filter(|(_, (s, _))| s.len() > 1)
        .flat_map(|(&len, bucket)| (0..len).map(move |pos| (bucket, pos)))
        .collect();

    let n_samp = working.len() as u64;
    let n_ref = ref_unique.len() as u64;
    let body = |e: u32| -> &[u8] {
        if e & REF != 0 {
            interior(ref_unique[(e & !REF) as usize], b)
        } else {
            interior(seq(e), b)
        }
    };

    tasks
        .par_iter()
        .flat_map_iter(|&((sample, refs), pos)| {
            let masked_cmp = |x: u32, y: u32| {
                let (x, y) = (body(x), body(y));
                x[..pos].cmp(&y[..pos]).then_with(|| x[pos + 1..].cmp(&y[pos + 1..]))
            };
            let mut order: Vec<u32> = sample.iter().chain(refs.iter()).copied().collect();
            order.sort_unstable_by(|&x, &y| masked_cmp(x, y).then(x.cmp(&y)));
            let mut out = Vec::new();
            for run in order.chunk_by(|&x, &y| masked_cmp(x, y).is_eq()) {
                // Sample ids sort before REF-flagged entries within a run.
                let n_s = run.partition_point(|&e| e & REF == 0);
                if n_s >= 2 {
                    let num_in_ref = (run.len() - n_s) as u64;
                    global_components(rep, &run[..n_s], pos, num_in_ref, n_samp, n_ref, p, lf, &mut out);
                }
            }
            out
        })
        .collect()
}

/// Splits one tag's sample members into turboGliph's connected components
/// (same V gene if required; residues joined by BLOSUM62 ≥ 0 unless all are
/// interchangeable) and appends a group per component of ≥ 2 nodes.
#[allow(clippy::too_many_arguments)]
fn global_components(
    rep: &Repertoire,
    members: &[SequenceId],
    pos: usize,
    num_in_ref: u64,
    n_samp: u64,
    n_ref: u64,
    p: &Gliph2Params,
    lf: &LogFactorial,
    out: &mut Vec<ConvergenceGroup>,
) {
    let b = p.boundary_size;
    let residue = |m: SequenceId| interior(&rep.sequences[m as usize], b)[pos];
    // Nodes: CDR3s, or (CDR3, TRBV) pairs when V genes must match.
    let mut nodes: Vec<(SequenceId, u32)> = Vec::with_capacity(members.len());
    for &m in members {
        if p.global_vgene {
            let mut vs: Vec<u32> = rep.rows_of(m).iter().map(|r| r.v_gene).collect();
            vs.sort_unstable();
            vs.dedup();
            nodes.extend(vs.into_iter().map(|v| (m, v)));
        } else {
            nodes.push((m, NA));
        }
    }
    // Component id per residue: one component unless BLOSUM restricts edges,
    // then a tiny union-find over the ≤ 20 residues present (identity pairs
    // are BLOSUM-compatible, so equal residues always connect).
    let mut comp_of = [0u8; 256];
    if !p.all_aa_interchangeable {
        let mut present: Vec<u8> = members.iter().map(|&m| residue(m)).collect();
        present.sort_unstable();
        present.dedup();
        let mut parent: Vec<usize> = (0..present.len()).collect();
        fn root(parent: &mut [usize], mut x: usize) -> usize {
            while parent[x] != x {
                x = parent[x];
            }
            x
        }
        for i in 0..present.len() {
            for j in i + 1..present.len() {
                if blosum_compatible(present[i], present[j]) {
                    let (ri, rj) = (root(&mut parent, i), root(&mut parent, j));
                    parent[ri.max(rj)] = ri.min(rj);
                }
            }
        }
        for (i, &a) in present.iter().enumerate() {
            comp_of[a as usize] = root(&mut parent, i) as u8;
        }
    }
    nodes.sort_unstable_by_key(|&(m, v)| (v, comp_of[residue(m) as usize], m));
    for c in nodes.chunk_by(|x, y| x.1 == y.1 && comp_of[residue(x.0) as usize] == comp_of[residue(y.0) as usize]) {
        if c.len() < 2 {
            continue;
        }
        let residues = c.iter().fold(0u32, |acc, &(m, _)| acc | aa_bit(residue(m)));
        let mut ms: Vec<SequenceId> = c.iter().map(|n| n.0).collect();
        ms.sort_unstable_by(|&x, &y| rep.sequences[x as usize].cmp(&rep.sequences[y as usize]));
        ms.dedup();
        let u = ms.len() as u64;
        let fisher_score = match (n_ref + n_samp).checked_sub(num_in_ref + u) {
            Some(n) => hypergeom_upper(u, num_in_ref + u, n, n_samp, lf),
            None => f64::NAN,
        };
        out.push(ConvergenceGroup {
            kind: GroupType::Global,
            tag: GroupTag::Global {
                boundary: b as u8,
                pos: pos as u16,
                residues,
                v_gene: p.global_vgene.then_some(c[0].1),
            },
            cluster_size: c.len() as u32,
            unique_cdr3_sample: ms.len() as u32,
            unique_cdr3_ref: num_in_ref,
            ove: 0.0,
            fisher_score,
            members: ms,
        });
    }
}

/// Writes groups in turboGliph's `cluster_properties` layout (scores
/// omitted); `members` is space-separated, sorted CDR3s.
pub fn write_groups(w: &mut impl std::io::Write, rep: &Repertoire, groups: &[ConvergenceGroup]) -> std::io::Result<()> {
    writeln!(w, "type\ttag\tcluster_size\tunique_cdr3_sample\tunique_cdr3_ref\tOvE\tfisher.score\tmembers")?;
    for g in groups {
        write!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:e}\t",
            g.kind.as_str(),
            g.tag(rep),
            g.cluster_size,
            g.unique_cdr3_sample,
            g.unique_cdr3_ref,
            g.ove,
            g.fisher_score
        )?;
        for (i, &m) in g.members.iter().enumerate() {
            if i > 0 {
                w.write_all(b" ")?;
            }
            w.write_all(&rep.sequences[m as usize])?;
        }
        w.write_all(b"\n")?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rep(seqs: &[&str]) -> Repertoire {
        Repertoire::from_sequences(seqs.iter().map(|s| s.as_bytes()))
    }

    fn members(r: &Repertoire, g: &ConvergenceGroup) -> Vec<String> {
        g.members.iter().map(|&m| String::from_utf8(r.sequences[m as usize].clone()).unwrap()).collect()
    }

    #[test]
    fn rounding_matches_formatc() {
        assert_eq!(round_2sig(0.012345), 0.012);
        assert_eq!(round_2sig(0.000987), 0.00099);
        assert_eq!(round_2sig(1.0), 1.0);
        assert_eq!(round_1dp(12.345), 12.3);
    }

    #[test]
    fn blosum_table_has_112_pairs() {
        let aa = b"ACDEFGHIKLMNPQRSTVWY";
        let n = aa.iter().flat_map(|&a| aa.iter().map(move |&b| (a, b))).filter(|&(a, b)| blosum_compatible(a, b)).count();
        assert_eq!(n, 112);
        assert!(blosum_compatible(b'I', b'V') && blosum_compatible(b'V', b'I') && !blosum_compatible(b'W', b'A'));
    }

    #[test]
    fn global_groups_use_interior_only() {
        // Interiors (boundary 3): GQETQ / GPETQ / GQETQ, and a different length.
        let r = rep(&["CASGQETQYQF", "CASGPETQYQF", "CSRGQETQYQF", "CASGQETQYQAF", "CASGWWWWYQF"]);
        let p = Gliph2Params { local: false, all_aa_interchangeable: true, ..Default::default() };
        let res = gliph2(&r, &[], &p);
        let tags: Vec<String> = res.groups.iter().map(|g| g.tag(&r)).collect();
        // Identical interiors match at every position; the substitution site
        // yields one group with residues P and Q.
        assert!(tags.iter().any(|t| t == "G%ETQ_PQ"));
        assert!(tags.iter().any(|t| t == "%QETQ_G"));
        let g = res.groups.iter().find(|g| g.tag(&r) == "G%ETQ_PQ").unwrap();
        assert_eq!(members(&r, g), vec!["CASGPETQYQF", "CASGQETQYQF", "CSRGQETQYQF"]);
    }

    #[test]
    fn blosum_restricts_global_edges() {
        let r = rep(&["CASGIETQYQF", "CASGVETQYQF", "CASGWETQYQF"]);
        let strict = Gliph2Params { local: false, ..Default::default() };
        let res = gliph2(&r, &[], &strict);
        let g: Vec<String> = res.groups.iter().map(|g| g.tag(&r)).filter(|t| t.starts_with("G%ETQ")).collect();
        assert_eq!(g, vec!["G%ETQ_IV"]); // W is isolated: BLOSUM(W,I), BLOSUM(W,V) < 0
    }

    #[test]
    fn local_groups_split_by_position() {
        // Motif WWW at interior starts 1, 2 (chained) and 7 (break > 2).
        let mut seqs: Vec<String> = Vec::new();
        for (i, pre) in ["", "A"].iter().enumerate() {
            for j in 0..3 {
                seqs.push(format!("CAS{pre}WWW{}GQYF", "DEK".repeat(1 + j + i)));
            }
        }
        for j in 0..3 {
            seqs.push(format!("CASDEKDEWWW{}YF", "GHT".repeat(1 + j)));
        }
        let r = rep(&seqs.iter().map(String::as_str).collect::<Vec<_>>());
        let reference: Vec<Vec<u8>> = (0..200).map(|i| format!("CASSL{}GYTF", "AQ".repeat(1 + i % 5)).into_bytes()).collect();
        let p = Gliph2Params { global: false, lcminove: vec![1.0], ..Default::default() };
        let res = gliph2(&r, &reference, &p);
        let www: Vec<_> = res.groups.iter().filter(|g| g.tag(&r).starts_with("WWW_")).collect();
        let tags: Vec<String> = www.iter().map(|g| g.tag(&r)).collect();
        assert_eq!(tags, vec!["WWW_4_5", "WWW_9_9"]);
        assert_eq!(www[0].members.len(), 6);
        assert_eq!(www[1].members.len(), 3);
    }

    #[test]
    fn deterministic_across_threads() {
        let input = crate::synthetic::repertoire(4000, 3);
        let r = Repertoire::from_sequences(&input);
        let reference = crate::synthetic::reference(20000, 5);
        let p = Gliph2Params::default();
        let run = |t| {
            rayon::ThreadPoolBuilder::new().num_threads(t).build().unwrap().install(|| {
                gliph2(&r, &reference, &p).groups.into_iter().map(|g| (g.tag(&r), g.members)).collect::<Vec<_>>()
            })
        };
        let a = run(1);
        assert!(!a.is_empty());
        assert_eq!(a, run(6));
    }
}
