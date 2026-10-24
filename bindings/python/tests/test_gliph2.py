import math

import pytest

import gliph2rs

AA = "ACDEFGHIKLMNPQRSTVWY"

# 400 distinct reference CDR3s; the first 300 carry the G%ETQ interior.
REFERENCE = [
    f"C{AA[i % 20]}{AA[(i // 20) % 20]}G{AA[i % 20]}ETQYQF"
    if i < 300
    else f"CASSL{AA[i % 20]}{AA[(i // 20) % 20]}GYEQYF"
    for i in range(400)
]


def tags(res):
    return list(res["groups"]["tag"])


def row(res, tag):
    i = tags(res).index(tag)
    return {k: v[i] for k, v in res["groups"].items()}


def test_global_groups_use_interior_and_blosum():
    seqs = ["CASGIETQYQF", "CASGVETQYQF", "CASGWETQYQF"]
    strict = gliph2rs.gliph2(seqs, REFERENCE, local_similarities=False)
    # I and V are BLOSUM62-compatible; W is not, so it stays out.
    assert "G%ETQ_IV" in tags(strict)
    assert not any("W" in t for t in tags(strict))

    loose = gliph2rs.gliph2(seqs, REFERENCE, local_similarities=False, all_aa_interchangeable=True)
    assert "G%ETQ_IVW" in tags(loose)


def test_members_match_cluster_size():
    res = gliph2rs.gliph2(["CASGIETQYQF", "CASGVETQYQF"], REFERENCE, local_similarities=False)
    g = row(res, "G%ETQ_IV")
    assert sorted(g["members"]) == ["CASGIETQYQF", "CASGVETQYQF"]
    assert g["cluster_size"] == len(g["members"])


def test_enriched_motifs_give_local_groups():
    # WWW sits in the interior of four sample CDR3s and never in the reference.
    seqs = [
        "CASWWWDEKGQYF",
        "CASWWWDEKDEKGQYF",
        "CASWWWDEKDEKDEKGQYF",
        "CASWWWQQQGQYF",
        "CASSLAAGYEQYF",
    ]
    res = gliph2rs.gliph2(seqs, REFERENCE, lcminove=[1.0], global_similarities=False)
    assert "WWW" in res["motifs"]["motif"]
    assert all(t == "local" for t in res["groups"]["type"])
    assert any(t.startswith("WWW_") for t in tags(res))


def test_paper_params_are_accepted():
    seqs = ["CASGIETQYQF", "CASGVETQYQF", "CASSLNTEAFF"]
    res = gliph2rs.gliph2(seqs, REFERENCE, **gliph2rs.PAPER_PARAMS)
    # all_aa_interchangeable is on, so I/V still group.
    assert "G%ETQ_IV" in tags(res)


def test_donor_association():
    donors = [f"D{i}:CMV+" for i in range(1, 7)] + [f"D{i}:CMV-" for i in range(7, 13)]
    cdr3 = ["CASGQETQYQF"] * 3 + ["CASGPETQYQF"] * 3 + ["CASWAWWWYQF"] * 3 + ["CASWCWWWYQF"] * 3 + ["CASSLKKKYQF"] * 3
    subject = donors[:6] + donors[0:3] + donors[6:9] + donors[9:12]

    res = gliph2rs.gliph2(
        cdr3, REFERENCE, subject=subject, local_similarities=False,
        all_aa_interchangeable=True, case="CMV+", control="CMV-", assoc_min_donors=2,
    )
    assert res["association"]["n_case"] == 6
    assert res["association"]["n_control"] == 6

    case_only = row(res, "G%ETQ_PQ")
    assert (case_only["case_donors"], case_only["control_donors"]) == (6, 0)
    # The only table at least this extreme: 1 / choose(12, 6).
    assert case_only["assoc_p"] == pytest.approx(1 / math.comb(12, 6), rel=1e-12)

    shared = row(res, "W%WWW_AC")
    assert (shared["case_donors"], shared["control_donors"]) == (3, 3)
    assert shared["assoc_p"] > 0.5

    filtered = gliph2rs.gliph2(
        cdr3, REFERENCE, subject=subject, local_similarities=False,
        all_aa_interchangeable=True, case="CMV+", control="CMV-",
        assoc_min_donors=2, assoc_max_p=0.01,
    )
    assert tags(filtered) == ["G%ETQ_PQ"]


def test_hla_scoring():
    # Carriers D1-D5 all carry A*01; 5 of the other 15 donors do too.
    cdr3 = ["CASGIETQYQF", "CASGVETQYQF"] * 5 + [f"CASSL{a}KKYEQYF" for a in "ACDEFGHIKLMNPQR"]
    subject = [f"D{i}" for i in range(1, 6) for _ in range(2)] + [f"D{i}" for i in range(6, 21)]
    hla = {
        f"D{i}": (["HLA-B*07"] if i % 2 == 0 else ["HLA-B*08"]) + (["HLA-A*01"] if i <= 10 else [])
        for i in range(1, 21)
    }

    res = gliph2rs.gliph2(cdr3, REFERENCE, subject=subject, local_similarities=False, hla=hla)
    assert res["hla"]["typed_donors"] == 20

    g = row(res, "G%ETQ_IV")
    assert g["hla_allele"] == "HLA-A*01"
    assert (g["hla_donors"], g["hla_with_allele"]) == (5, 5)
    # phyper(4, 10, 10, 5, lower.tail = FALSE) = choose(10, 5) / choose(20, 5)
    assert g["hla_p"] == pytest.approx(math.comb(10, 5) / math.comb(20, 5), rel=1e-12)
    assert g["hla_q"] >= g["hla_p"]


def test_results_do_not_depend_on_thread_count():
    seqs = list(dict.fromkeys(
        [f"CASSL{AA[i % 20]}GQ{AA[i % 7]}ETQYF" for i in range(200)] + ["CASGIETQYQF", "CASGVETQYQF"]
    ))
    one = gliph2rs.gliph2(seqs, REFERENCE, threads=1)["groups"]
    many = gliph2rs.gliph2(seqs, REFERENCE, threads=4)["groups"]
    assert one == many


def test_bad_input_is_rejected():
    with pytest.raises(ValueError, match="no valid CDR3"):
        gliph2rs.gliph2([], REFERENCE)
    with pytest.raises(ValueError, match="no valid CDR3"):
        gliph2rs.gliph2(["XXXX", "1234"], REFERENCE)
    with pytest.raises(ValueError, match="reference is empty"):
        gliph2rs.gliph2(["CASSLGQETQYF"], [])
    with pytest.raises(ValueError, match="together"):
        gliph2rs.gliph2(["CASSLGQETQYF"], REFERENCE, case="CMV+")
    with pytest.raises(ValueError, match="entries but cdr3 has"):
        gliph2rs.gliph2(["CASSLGQETQYF"], REFERENCE, subject=["a", "b"])
