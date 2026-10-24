"""Fast, deterministic GLIPH2 convergence groups for TCR repertoires.

Groups CDR3 beta sequences likely to share antigen specificity, following the
GLIPH2 method of Huang et al. (2020). Group membership matches
``turboGliph::gliph2()`` exactly, while the computation runs in parallel Rust.

    >>> import gliph2rs
    >>> res = gliph2rs.gliph2(cdr3, reference)
    >>> import pandas as pd
    >>> pd.DataFrame(res["groups"])
"""

from __future__ import annotations

import re
from typing import Any, Iterable, Mapping, Sequence

from ._gliph2rs import run as _run

__all__ = ["gliph2", "PAPER_PARAMS"]
__version__ = "0.1.0"

#: GLIPH2's published parameter-file settings. The defaults of :func:`gliph2`
#: follow turboGliph's code instead, which is more permissive.
PAPER_PARAMS: dict[str, Any] = {
    "lcminp": 0.001,
    "lcminove": [10.0],
    "min_seq_length": 8,
    "all_aa_interchangeable": True,
}

_AMINO_ACID = re.compile(r"^[ACDEFGHIKLMNPQRSTVWY]+$")


def _as_str_list(values: Iterable[Any] | None, n: int, name: str) -> list[str]:
    """Normalises an optional per-sequence column to a list of ``n`` strings."""
    if values is None:
        return [""] * n
    out = ["" if v is None else str(v) for v in values]
    if len(out) != n:
        raise ValueError(f"{name} has {len(out)} entries but cdr3 has {n}")
    return ["" if v.lower() in ("nan", "na", "none") else v for v in out]


def _hla_pairs(hla: Mapping[str, Iterable[str]] | Iterable[tuple[str, str]] | None):
    """Flattens donor typing into parallel donor/allele lists."""
    if hla is None:
        return [], []
    donors: list[str] = []
    alleles: list[str] = []
    items = hla.items() if isinstance(hla, Mapping) else hla
    for donor, value in items:
        found = [value] if isinstance(value, str) else list(value)
        for allele in found:
            if donor and allele and str(allele).upper() not in ("NA", "NAN", "NONE"):
                donors.append(str(donor))
                alleles.append(str(allele))
    return donors, alleles


def gliph2(
    cdr3: Sequence[str],
    reference: Sequence[str],
    *,
    v_gene: Sequence[str] | None = None,
    subject: Sequence[str] | None = None,
    count: Sequence[int] | None = None,
    lcminp: float = 0.01,
    lcminove: Sequence[float] = (1000.0, 100.0, 10.0),
    kmer_mindepth: int = 3,
    motif_length: Sequence[int] = (2, 3, 4),
    discontinuous_motifs: bool = False,
    motif_distance_cutoff: int = 3,
    boundary_size: int = 3,
    min_seq_length: int = 0,
    accept_sequences_with_c_f_start_end: bool = True,
    all_aa_interchangeable: bool = False,
    global_vgene: bool = False,
    cluster_min_size: int = 2,
    local_similarities: bool = True,
    global_similarities: bool = True,
    min_subjects: int = 0,
    global_max_p: float | None = None,
    global_max_q: float | None = None,
    case: str | None = None,
    control: str | None = None,
    assoc_min_donors: int = 10,
    assoc_max_p: float | None = None,
    assoc_max_q: float | None = None,
    assoc_permute_seed: int | None = None,
    hla: Mapping[str, Iterable[str]] | None = None,
    threads: int = 0,
) -> dict[str, Any]:
    """Finds GLIPH2 convergence groups.

    Local groups come from motifs enriched against ``reference``; global groups
    from CDR3s whose interiors differ at one position. They are reported
    separately, not merged, as GLIPH2 does.

    Results are deterministic: repeated calls, and any ``threads``, give
    identical groups.

    Args:
        cdr3: CDR3 beta amino-acid sequences. Anything outside the 20 standard
            amino acids is dropped.
        reference: reference repertoire CDR3s, for motif enrichment.
        v_gene, subject, count: optional per-sequence columns, each the same
            length as ``cdr3``. ``subject`` holds GLIPH2's
            ``subject:condition`` strings, e.g. ``"HIP00110:CMV-"``.
        case, control: donor conditions to compare, read from the text after
            the last ``:`` of ``subject``. Given both, every group gets a
            donor-level Fisher p-value and q-value.
        assoc_permute_seed: shuffle the case/control labels among donors, as a
            calibration control.
        hla: donor typing as ``{donor: [alleles]}`` or ``(donor, allele)``
            pairs. Each group is then scored for HLA restriction.
        threads: worker threads; 0 lets Rust choose.

    Returns:
        A dict with ``groups`` and ``motifs`` (each a dict of equal-length
        columns, ready for ``pandas.DataFrame``), ``association`` and ``hla``
        summaries (or ``None``), and the input sizes. ``groups["members"]``
        holds one list of CDR3s per group.

    Raises:
        ValueError: if no valid CDR3s remain, if a reference is needed but
            empty, or if only one of ``case``/``control`` is given.
    """
    cdr3 = [str(c) for c in cdr3]
    n = len(cdr3)
    v_gene_l = _as_str_list(v_gene, n, "v_gene")
    subject_l = _as_str_list(subject, n, "subject")
    if count is None:
        count_l = [1] * n
    else:
        count_l = [int(c) if c and int(c) >= 1 else 1 for c in count]
        if len(count_l) != n:
            raise ValueError(f"count has {len(count_l)} entries but cdr3 has {n}")

    keep = [i for i, c in enumerate(cdr3) if _AMINO_ACID.match(c.upper())]
    if not keep:
        raise ValueError("no valid CDR3 amino-acid sequences in cdr3")
    reference = [str(r).upper() for r in reference]
    reference = [r for r in reference if _AMINO_ACID.match(r)]
    if not reference and local_similarities:
        raise ValueError("reference is empty; motif enrichment needs a reference repertoire")
    if (case is None) != (control is None):
        raise ValueError("case and control must be given together")

    donors, alleles = _hla_pairs(hla)
    params = {
        "lcminp": float(lcminp),
        "lcminove": [float(x) for x in lcminove],
        "kmer_mindepth": int(kmer_mindepth),
        "motif_length": [int(x) for x in motif_length],
        "discontinuous_motifs": bool(discontinuous_motifs),
        "motif_distance_cutoff": int(motif_distance_cutoff),
        "boundary_size": int(boundary_size),
        "min_seq_length": int(min_seq_length),
        "accept_sequences_with_c_f_start_end": bool(accept_sequences_with_c_f_start_end),
        "all_aa_interchangeable": bool(all_aa_interchangeable),
        "global_vgene": bool(global_vgene),
        "cluster_min_size": int(cluster_min_size),
        "local_similarities": bool(local_similarities),
        "global_similarities": bool(global_similarities),
        "min_subjects": int(min_subjects),
        "global_max_p": global_max_p,
        "global_max_q": global_max_q,
        "case": case,
        "control": control,
        "assoc_min_donors": int(assoc_min_donors),
        "assoc_max_p": assoc_max_p,
        "assoc_max_q": assoc_max_q,
        "assoc_permute_seed": assoc_permute_seed,
        "hla_donor": donors,
        "hla_allele": alleles,
        "threads": int(threads),
    }
    return _run(
        [cdr3[i].upper() for i in keep],
        [v_gene_l[i] for i in keep],
        [subject_l[i] for i in keep],
        [count_l[i] for i in keep],
        reference,
        params,
    )
