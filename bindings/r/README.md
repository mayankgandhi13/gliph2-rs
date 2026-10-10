# gliph2rs

R bindings for [gliph2-rs](../../), the Rust reimplementation of the GLIPH2 core
engine, built with [extendr](https://extendr.github.io/).

Group membership is identical to `turboGliph::gliph2()` — verified on subsets of
the Emerson 2017 cohort from 2,000 to 100,000 CDR3s, under two parameter sets —
while the computation runs in parallel Rust.

## Install

The package compiles the Rust core, so [Rust](https://rustup.rs) must be
installed (`cargo` on your `PATH`).

```r
# from a checkout of the repository
install.packages("bindings/r", repos = NULL, type = "source")

# or straight from GitHub
remotes::install_github("mayankgandhi13/gliph2-rs", subdir = "bindings/r")
```

The core crate is vendored into `src/rust/core` at build time, so the package is
self-contained once built.

## Use

```r
library(gliph2rs)

res <- gliph2(cdr3_sequences = my_tcrs, refdb_beta = ref_CD48)
head(res$cluster_properties)
res$cluster_list[["SLIGVSSYNE_4_17"]]
```

`cdr3_sequences` is a character vector, or a data frame with a `CDR3b` column
and optionally `TRBV`, `patient` and `counts`. The result is a list:

| Element | Contents |
|---|---|
| `cluster_properties` | one row per convergence group: `type`, `tag`, sizes, `OvE`, `fisher.score`, `fdr.q`, `n_subjects`, `members` |
| `selected_motifs` | the enriched motifs behind the local groups |
| `cluster_list` | member CDR3s per group, named by tag |
| `association` | donor counts and test totals, when `case`/`control` are given |

### GLIPH2's published parameters

`gliph2()` defaults follow turboGliph's code. The settings distributed with
GLIPH2 itself are stricter:

```r
res <- do.call(gliph2, c(list(cdr3_sequences = my_tcrs, refdb_beta = ref_CD48),
                         gliph2_paper_params()))
```

### Testing groups against a donor condition

With donors labelled `subject:condition` (as GLIPH2 expects, e.g.
`HIP00110:CMV-`), every group can be tested for enrichment in one condition:

```r
res <- gliph2(cohort, ref_CD48, case = "CMV+", control = "CMV-",
              assoc_min_donors = 10, assoc_max_q = 0.05)
```

This adds `case_donors`, `control_donors`, `assoc.p` and `assoc.q` columns.
`assoc_permute_seed` shuffles the labels among donors as a calibration control.

## Differences from turboGliph

- Group membership, tags and Fisher scores match exactly.
- The simulation-based cluster scores (network size, CDR3 length, V gene, clonal
  expansion, HLA) are **not** computed, so `total.score` and friends are absent.
- Results are deterministic: repeated calls, and any `n_cores`, give identical
  groups. turboGliph's scores vary between runs because they are simulated.
