//! gliph2-rs: a fast, deterministic reimplementation of the GLIPH2 core.
//!
//! Pipeline: [`local_cluster`] (Hamming ≤ 1 within length buckets) →
//! [`motif_enrich`] (k-mer enrichment vs. a reference, Fisher's exact test) →
//! [`assemble_clusters`] (union-find over both edge sets).
//!
//! Naming note: the GLIPH2 paper calls Hamming-distance grouping "global" and
//! motif grouping "local". This crate follows the README's naming instead:
//! *local* = Hamming edges, *motif* = enriched k-mer edges.

pub mod cluster_assembly;
pub mod enrichment;
pub mod io;
pub mod local_cluster;
pub mod motif;
pub mod synthetic;

use std::io::Write;
use std::time::{Duration, Instant};

pub use cluster_assembly::{assemble_clusters, AssemblyParams, Cluster, UnionFind};
pub use enrichment::{fisher_greater, motif_enrich, EnrichmentParams, MotifResult};
pub use io::{read_reference, read_tcr_table, Repertoire, SequenceId, TcrRecord};
pub use local_cluster::{local_cluster, LocalEdge, LocalParams};
pub use motif::{Kmer, MotifParams, MotifTable};

#[derive(Debug, Clone, Default)]
pub struct Params {
    pub local: LocalParams,
    pub motif: MotifParams,
    pub enrichment: EnrichmentParams,
    pub assembly: AssemblyParams,
    /// Skip motif enrichment entirely (local edges only).
    pub skip_motifs: bool,
}

#[derive(Debug, Clone, Default)]
pub struct Timings {
    pub local: Duration,
    pub reference_table: Duration,
    pub enrichment: Duration,
    pub assembly: Duration,
}

#[derive(Debug, Clone)]
pub struct GliphResult {
    pub local_edges: Vec<LocalEdge>,
    pub motifs: Vec<MotifResult>,
    pub clusters: Vec<Cluster>,
    pub timings: Timings,
}

/// Runs the full pipeline on deduplicated `sequences`.
pub fn run(sequences: &[Vec<u8>], reference: Option<&[Vec<u8>]>, params: &Params) -> GliphResult {
    let mut t = Timings::default();

    let now = Instant::now();
    let local_edges = local_cluster(sequences, &params.local);
    t.local = now.elapsed();

    let motifs = match reference {
        Some(r) if !params.skip_motifs => {
            let now = Instant::now();
            let table = MotifTable::build(r, &params.motif);
            t.reference_table = now.elapsed();
            let now = Instant::now();
            let m = motif_enrich(sequences, &table, &params.motif, &params.enrichment);
            t.enrichment = now.elapsed();
            m
        }
        _ => Vec::new(),
    };

    let now = Instant::now();
    let clusters = assemble_clusters(sequences.len(), &local_edges, &motifs, &params.assembly);
    t.assembly = now.elapsed();

    GliphResult { local_edges, motifs, clusters, timings: t }
}

/// Writes one row per (cluster, member sequence).
pub fn write_clusters(w: &mut impl Write, rep: &Repertoire, clusters: &[Cluster]) -> std::io::Result<()> {
    writeln!(w, "cluster_id\tcluster_size\tn_local_edges\tmotifs\tCDR3b\tTRBV\tTRBJ\tsubject\tcount")?;
    for c in clusters {
        let motifs = if c.motifs.is_empty() {
            "-".to_string()
        } else {
            c.motifs.iter().map(|m| m.to_string()).collect::<Vec<_>>().join(",")
        };
        for &id in &c.members {
            let seq = std::str::from_utf8(&rep.sequences[id as usize]).unwrap();
            let recs = &rep.records[id as usize];
            if recs.is_empty() {
                writeln!(w, "{}\t{}\t{}\t{}\t{}\tNA\tNA\tNA\t1", c.id, c.members.len(), c.n_local_edges, motifs, seq)?;
            }
            for r in recs {
                writeln!(
                    w,
                    "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    c.id,
                    c.members.len(),
                    c.n_local_edges,
                    motifs,
                    seq,
                    r.v_gene.as_deref().unwrap_or("NA"),
                    r.j_gene.as_deref().unwrap_or("NA"),
                    r.subject.as_deref().unwrap_or("NA"),
                    r.count
                )?;
            }
        }
    }
    Ok(())
}

pub fn write_motifs(w: &mut impl Write, motifs: &[MotifResult]) -> std::io::Result<()> {
    writeln!(w, "motif\tk\tinput_count\treference_count\tfold\tp_value")?;
    for m in motifs {
        writeln!(
            w,
            "{}\t{}\t{}\t{}\t{}\t{:e}",
            m.motif,
            m.motif.len(),
            m.input_count,
            m.reference_count,
            m.fold,
            m.p_value
        )?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pipeline_is_reproducible_across_thread_counts() {
        let input = synthetic::repertoire(5_000, 7);
        let rep = Repertoire::from_sequences(&input);
        let reference = synthetic::reference(20_000, 11);
        let params = Params::default();
        let fingerprint = |threads: usize| {
            let pool = rayon::ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
            pool.install(|| {
                let r = run(&rep.sequences, Some(&reference), &params);
                let clusters: Vec<_> = r.clusters.iter().map(|c| c.members.clone()).collect();
                let motifs: Vec<_> = r.motifs.iter().map(|m| (m.motif, m.p_value.to_bits())).collect();
                (r.local_edges, motifs, clusters)
            })
        };
        let a = fingerprint(1);
        assert!(!a.0.is_empty() && !a.1.is_empty() && !a.2.is_empty());
        assert_eq!(a, fingerprint(4));
        assert_eq!(a, fingerprint(4));
    }
}
