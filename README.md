# gliph2-rs

> 🚧 **Status: In progress.** The core engine works and matches turboGliph's GLIPH2 convergence groups exactly on real data (see [Validation](#validation)). Bindings and scoring are not done yet.

A fast, deterministic Rust reimplementation of the GLIPH2 core engine for clustering T-cell receptors (TCRs) by predicted shared antigen specificity, with bindings for R and Python.

## Why

GLIPH2 groups TCR CDR3 sequences that likely recognize the same antigen. The original is a single-threaded binary with two practical problems:

- **Speed:** pairwise comparisons scale as O(n²), so large repertoires (100K+ sequences) can take days.
- **Reproducibility:** motif enrichment results have been reported to vary between runs.

Faster alternatives such as clusTCR and GIANA exist in Python, but there is no Rust implementation yet.

**Goal:** reproducible, GLIPH2-equivalent clusters with a large speedup at scale.

## How GLIPH2 works

1. **Input:** CDR3β amino-acid sequences, with V/J genes and donor/HLA metadata
2. **Local similarity:** extract 2–4 aa motifs from the CDR3 interior and group sequences sharing motifs that are enriched vs. a reference repertoire
3. **Global similarity:** group same-length CDR3s whose interiors differ by at most one amino acid
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
│   ├── convergence.rs       # GLIPH2 convergence groups (turboGliph-compatible)
│   ├── io.rs                # parallel CDR3 input parsing
│   ├── synthetic.rs         # deterministic synthetic repertoires
│   └── main.rs              # `gliph2-rs` CLI
├── hpc/                     # Slurm jobs, data setup, validation vs. turboGliph
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

## Usage

```bash
cargo build --release

# GLIPH2 convergence groups, matching turboGliph::gliph2()
./target/release/gliph2-rs gliph2 --input tcrs.tsv --reference ref_CD48_v2.0.tsv --paper-params --out results
# -> results_groups.tsv (type, tag, sizes, fisher.score, members) and a JSON timing line

# Merged clusters (connected components of Hamming-1 and motif edges)
./target/release/gliph2-rs cluster --input tcrs.tsv --reference ref_CD48_v2.0.tsv --fdr --out results
# -> results_clusters.tsv, results_motifs.tsv
```

Input is a GLIPH2-style TSV (`CDR3b  TRBV  TRBJ  CDR3a  subject:condition  count`), a table with VDJtools or AIRR column names, or one CDR3 per line. `--paper-params` uses the values from GLIPH2's distributed parameter files (motif p ≤ 0.001, fold ≥ 10, CDR3 length ≥ 8, all substitutions allowed); without it, turboGliph's code defaults apply. `--gapped` / `--discontinuous` add motifs with one wildcard position, and `--fdr` thresholds Benjamini–Hochberg q-values instead of raw p-values.

GLIPH2 reports every global group without a significance filter, which on a pooled cohort means tens of millions of groups. The `gliph2` command adds a q-value (`fdr.q`, Benjamini–Hochberg within each group type) and a donor count (`n_subjects`) to every group, and three opt-in filters: `--max-global-p P`, `--max-global-q Q` and `--min-subjects N` (groups spanning at least N donors). All are off by default, so the output still matches turboGliph.

The global-group p-value compares the sample against the reference, so it only has power when the sample is smaller than the reference (a single donor, say). On the pooled Emerson cohort, 73.6M CDR3s against a 1.19M reference, no global group can reach significance and `--max-global-q 0.05` keeps none; `--min-subjects 3` keeps 43.6M of 72.0M. Cohort-scale filtering needs a test against donor labels instead (see Roadmap).

**Local clustering without O(n²):** within each length bucket, each position is masked in turn and the bucket is sorted by the masked sequence. Two distinct sequences differ at exactly one position iff they collide under exactly one mask, so every Hamming-1 pair is found exactly once in O(n·L·log n).

> Naming: the `gliph2` command follows the GLIPH2 paper (*local* = motif groups, *global* = Hamming groups). The older `cluster` command calls Hamming edges *local*.

### Running on an HPC cluster (Slurm)

```bash
sbatch hpc/bench_scaling.sbatch    # synthetic size + thread scaling, reproducibility hashes
sbatch hpc/setup_data.sbatch       # download Emerson 2017 cohort + GLIPH2 references, install turboGliph
hpc/prep_hip.sh DATA_DIR           # pool the 786 repertoires into one GLIPH2-format table
sbatch hpc/validate.sbatch         # gliph2-rs vs. turboGliph on 2K–100K subsets
sbatch hpc/run_full_cohort.sbatch  # gliph2-rs on the full 151M-row cohort
```

## Validation

The original GLIPH2 binary was unavailable (its download server is offline), so the reference implementation is [turboGliph](https://github.com/HetzDra/turboGliph)'s `gliph2()`. Both tools ran on identical subsets of the Emerson et al. 2017 cohort (786 TCRβ repertoires), against the GLIPH2 v2.0 CD4+CD8 reference (1.19M CDR3s), on 32 cores of Northeastern's Explorer cluster. Every convergence group was matched by tag and compared member by member.

| Input CDR3s | Groups identical (paper settings) | Groups identical (turboGliph defaults) | gliph2-rs | turboGliph |
|---|---|---|---|---|
| 2,000 | 17 / 17 | 13 / 13 | 0.99 s | 18.3 s |
| 10,000 | 210 / 210 | 109 / 109 | 0.78 s | 23.8 s |
| 50,000 | 4,226 / 4,226 | 1,720 / 1,720 | 0.99 s | 46.3 s |
| 100,000 | 14,869 / 14,869 | 6,164 / 6,164 | 1.03 s | 89.7 s |

Times are for the paper settings; gliph2-rs times include reading the input and reference, turboGliph times cover only its `gliph2()` call. Raw results: [`hpc/results/`](hpc/results/).

**Full cohort:** all 151,020,646 rows (74.2M unique CDR3s) in 4.4 minutes on 64 cores, with about 25 GB peak memory per `/usr/bin/time` (71 GB per Slurm accounting): 190 enriched motifs, 194 local groups and 72.0M global groups.

**Reproducibility:** repeated runs produce byte-identical output at any thread count.

### CMV association on the full cohort

`gliph2 --case CMV+ --control CMV-` tests every convergence group for enrichment in CMV-seropositive donors (one-sided Fisher at the donor level, Benjamini–Hochberg over groups seen in at least 10 labelled donors). On the pooled Emerson cohort (340 CMV+, 421 CMV− donors; [`hpc/cmv_association.sbatch`](hpc/cmv_association.sbatch)):

- 12.1M groups tested, **57 CMV-associated at q ≤ 0.05**, forming 18 independent CDR3 families. With CMV labels shuffled among donors, **0** groups pass.
- The strongest families are carried almost only by CMV+ donors, e.g. `SIGPLEHNE` in 42 CMV+ and 0 CMV− donors (p = 4.5×10⁻¹⁶).
- **7 of 18 families are HLA-restricted** (q ≤ 0.05, 626 HLA-typed donors). The one family VDJdb links to a CMV epitope, `SLIGVSSYNE` (pp65 TPRVTGGGAM, presented by HLA-B*07:02), is carried by HLA-B*07 donors 83% of the time against 21% overall (q = 9×10⁻¹³). Others are restricted by HLA-A*01, A*24 and B*08, alleles that present well-known CMV epitopes.
- VDJdb alone cannot confirm specificity: associated groups are equally enriched for CMV and other-virus entries (22.5× vs 22.6×, size-matched), and VDJdb covers few CMV epitopes.

Results: [`hpc/results/cmv-10881270_*`](hpc/results/).

Not yet compared: turboGliph's simulation-based cluster scores (network size, CDR3 length, V gene, clonal expansion, HLA).

## Roadmap

- [x] Input parsing + length bucketing
- [x] Local similarity clustering
- [x] Union-find cluster assembly
- [x] Motif extraction + reference frequency table
- [x] Deterministic enrichment testing
- [x] Synthetic scaling benchmark (Slurm, `hpc/bench_scaling.sbatch`)
- [x] Parallel, column-wise input parsing for 10⁸-row cohorts
- [x] Gapped motifs and optional FDR control
- [x] GLIPH2 convergence groups, validated against turboGliph on real data
- [x] Full-cohort run (Emerson 2017, 151M rows)
- [x] Per-group q-values, donor counts and opt-in global-group filters
- [x] Donor-label association (CMV+ vs CMV−), with permutation, VDJdb and HLA checks
- [ ] Cluster scoring (network size, CDR3 length, V gene, clonal expansion, HLA)
- [ ] R binding (extendr)
- [ ] Python binding (PyO3)
- [ ] Write-up of results

## Acknowledgements

Based on the GLIPH2 method:

> Huang H, Wang C, Rubelt F, Scriba TJ, Davis MM (2020). *Analyzing the Mycobacterium tuberculosis immune response by T-cell receptor clustering with GLIPH2 and genome-wide antigen screening.* Nature Biotechnology.

This is an independent reimplementation and is not affiliated with the original authors.

## Author

**Mayank Gandhi**, MS Bioinformatics, Northeastern University
[LinkedIn](https://www.linkedin.com/in/mayankgandhi0713) · [GitHub](https://github.com/mayankgandhi13)
