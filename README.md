# gliph2-rs

> 🚧 **Status: In progress.** This project is in the planning and early development stage. No benchmarks yet; this README describes the design and roadmap.

A fast, deterministic Rust reimplementation of the GLIPH2 core engine for clustering T-cell receptors (TCRs) by predicted shared antigen specificity, with bindings for R and Python.

## Why

GLIPH2 groups TCR CDR3 sequences that likely recognize the same antigen. The original is a single-threaded binary with two practical problems:

- **Speed:** pairwise comparisons scale as O(n²), so large repertoires (100K+ sequences) can take days.
- **Reproducibility:** motif enrichment results have been reported to vary between runs.

Faster alternatives such as clusTCR and GIANA exist in Python, but there is no Rust implementation yet.

**Goal:** reproducible, GLIPH2-equivalent clusters with a large speedup at scale.

## How GLIPH2 works

1. **Input:** CDR3β amino-acid sequences, with V/J genes and donor/HLA metadata
2. **Local similarity:** group same-length CDR3s that differ by at most one amino acid (Hamming distance ≤ 1)
3. **Global similarity:** extract 2–4 aa motifs and group sequences sharing motifs that are enriched vs. a reference repertoire
4. **Motif enrichment:** Fisher's exact test, input vs. reference
5. **Cluster assembly:** merge sequences connected by local or global edges
6. **HLA scoring:** test clusters for shared HLA alleles among donors

## What moves to Rust

| Step | Rust approach |
|---|---|
| Local clustering | Bucket sequences by length, then compare pairs within each bucket in parallel (`rayon`) |
| Motif extraction | Parallel k-mer extraction per sequence |
| Enrichment testing | Reference k-mer frequency table built **once**, then deterministic Fisher's exact tests per motif |
| Cluster assembly | Union-find with a fixed edge order, giving **reproducible** clusters |

**Stays in R/Python:** input parsing and QC, V/J annotation, HLA enrichment statistics, and visualization.

## Design

- **CDR3 sequences:** `Vec<u8>` amino-acid bytes, grouped in a `HashMap<usize, Vec<SequenceId>>` keyed by length. Length bucketing alone removes most of the pairwise work.
- **Reference k-mer table:** `HashMap<Kmer, u32>`, built once and reused
- **Cluster assembly:** union-find (parent array with path compression), single-threaded, which is fast and deterministic

### Planned crate structure

```
gliph2-rs/
├── Cargo.toml
├── src/
│   ├── lib.rs               # public API
│   ├── local_cluster.rs     # length-bucketed Hamming clustering
│   ├── motif.rs             # k-mer extraction + reference table
│   ├── enrichment.rs        # deterministic Fisher's exact test
│   ├── cluster_assembly.rs  # union-find merge
│   └── io.rs                # CDR3 input parsing
├── bindings/
│   ├── r/                   # extendr wrapper
│   └── python/              # PyO3 wrapper (optional)
└── benches/
    └── vs_gliph2.rs
```

### Planned API

The R binding (via `extendr`) is meant to slot into existing GLIPH2 / turboGliph workflows:

```r
local    <- local_cluster(cdr3)
motifs   <- motif_enrich(cdr3, reference)
clusters <- assemble_clusters(local, motifs)
```

A PyO3 binding is planned for comparison against Python tools like clusTCR and GIANA.

## Validation plan

1. Run original GLIPH2 and gliph2-rs on identical TCR repertoire data with the same parameters
2. Compare:
   - **Cluster membership overlap** with the original
   - **Reproducibility** across repeated runs
   - **Wall-clock time** at the current dataset size
   - **Scaling** at 10× and 100× size, using synthetically expanded repertoires

## Roadmap

- [ ] Input parsing + length bucketing
- [ ] Local similarity clustering
- [ ] Union-find cluster assembly
- [ ] Motif extraction + reference frequency table
- [ ] Deterministic enrichment testing
- [ ] R binding (extendr)
- [ ] Python binding (PyO3)
- [ ] Benchmark suite vs. GLIPH2
- [ ] Write-up of results

## Acknowledgements

Based on the GLIPH2 method:

> Huang H, Wang C, Rubelt F, Scriba TJ, Davis MM (2020). *Analyzing the Mycobacterium tuberculosis immune response by T-cell receptor clustering with GLIPH2 and genome-wide antigen screening.* Nature Biotechnology.

This is an independent reimplementation and is not affiliated with the original authors.

## Author

**Mayank Gandhi**, MS Bioinformatics, Northeastern University
[LinkedIn](https://www.linkedin.com/in/mayankgandhi0713) · [GitHub](https://github.com/mayankgandhi13)
