# gliph2rs (Python)

Python bindings for [gliph2-rs](../../), the Rust reimplementation of the GLIPH2
core engine, built with [PyO3](https://pyo3.rs).

Group membership is identical to `turboGliph::gliph2()`, verified on subsets of
the Emerson 2017 cohort, while the computation runs in parallel Rust and scales
to pooled cohorts of 100 million sequences.

## Install

The package compiles the Rust core, so [Rust](https://rustup.rs) must be
installed (`cargo` on your `PATH`).

```bash
# from a checkout of the repository
pip install ./bindings/python
```

## Use

```python
import gliph2rs
import pandas as pd

res = gliph2rs.gliph2(cdr3, reference)
groups = pd.DataFrame(res["groups"])
groups.sort_values("fisher_score").head()
```

`groups` has one row per convergence group: `type`, `tag`, `cluster_size`,
`unique_cdr3_sample`, `unique_cdr3_ref`, `OvE`, `fisher_score`, `fdr_q`,
`n_subjects` and `members` (a list of CDR3s). `res["motifs"]` holds the enriched
motifs behind the local groups.

### GLIPH2's published parameters

The defaults follow turboGliph's code. The settings distributed with GLIPH2
itself are stricter:

```python
res = gliph2rs.gliph2(cdr3, reference, **gliph2rs.PAPER_PARAMS)
```

### Testing groups against a donor condition

With donors labelled `subject:condition` (as GLIPH2 expects, e.g.
`HIP00110:CMV-`), every group can be tested for enrichment in one condition:

```python
res = gliph2rs.gliph2(cdr3, reference, subject=donors,
                      case="CMV+", control="CMV-",
                      assoc_min_donors=10, assoc_max_q=0.05)
```

This adds `case_donors`, `control_donors`, `assoc_p` and `assoc_q` columns.
`assoc_permute_seed=` shuffles the labels among donors as a calibration control.

### HLA restriction

```python
res = gliph2rs.gliph2(cdr3, reference, subject=donors,
                      hla={"HIP00110:CMV-": ["HLA-A*03", "HLA-B*07"], ...})
```

Each allele among a group's typed carriers is tested against the typed cohort
(one-sided hypergeometric); the best allele is reported with a q-value
(Bonferroni within the group, then Benjamini-Hochberg across groups).

## Differences from turboGliph

- Group membership, tags and Fisher scores match exactly.
- HLA restriction is scored deterministically. The four simulation-based
  cluster scores (network size, CDR3 length, V gene, clonal expansion) are
  **not** computed: turboGliph never seeds that RNG, so its values differ
  between runs and cannot be matched.
- Results are deterministic: repeated calls, and any `threads`, give identical
  groups.
