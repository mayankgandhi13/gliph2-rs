#!/usr/bin/env python3
"""Checks CMV-associated convergence groups against known antigen-specific TCRs.

Usage: assoc_vdjdb.py GROUPS.tsv VDJDB_SLIM.txt TOP_OUT.tsv [Q_CUTOFF]

GROUPS.tsv is `gliph2-rs gliph2 --case ... --control ...` output. Among tested
groups (those with an association p-value), it compares groups significant at
q <= Q_CUTOFF (default 0.05) with the rest: how often each contains a CDR3
that VDJdb lists as CMV-specific, and, as a specificity control, how often it
contains one specific to another virus. Prints one JSON summary line and
writes the top associated groups to TOP_OUT.tsv.

Larger groups are more likely to contain any VDJdb CDR3, and associated
groups tend to be large, so the summary also gives size-matched expectations:
for each significant group, the hit rate of non-significant tested groups
with the same member count. Overlapping significant groups (sharing a CDR3,
e.g. one family seen at several masked positions) are merged into families.
"""
import csv
import json
import math
import sys
from collections import defaultdict

csv.field_size_limit(sys.maxsize)


def log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_greater(a, b, c, d):
    """One-sided Fisher p-value P(X >= a) for the 2x2 table [[a, b], [c, d]]."""
    n1, k, n = a + b, a + c, a + b + c + d
    denom = log_choose(n, n1)
    return min(1.0, sum(math.exp(log_choose(k, x) + log_choose(n - k, n1 - x) - denom) for x in range(a, min(k, n1) + 1)))


def poisson_upper(k, lam):
    """P(X >= k) for X ~ Poisson(lam)."""
    if k <= 0:
        return 1.0
    term = math.exp(-lam)
    cdf = term
    for i in range(1, k):
        term *= lam / i
        cdf += term
    return max(0.0, 1.0 - cdf)


def families(sig):
    """Union-find over significant groups that share a member CDR3."""
    parent = list(range(len(sig)))

    def root(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    first = {}
    for i, g in enumerate(sig):
        for m in g["members"]:
            if m in first:
                a, b = root(i), root(first[m])
                parent[max(a, b)] = min(a, b)
            else:
                first[m] = i
    fams = defaultdict(list)
    for i in range(len(sig)):
        fams[root(i)].append(sig[i])
    return list(fams.values())


def table(sig_hit, sig_miss, rest_hit, rest_miss):
    odds = (sig_hit + 0.5) * (rest_miss + 0.5) / ((sig_miss + 0.5) * (rest_hit + 0.5))
    return {
        "significant": {"with_hit": sig_hit, "total": sig_hit + sig_miss,
                        "fraction": round(sig_hit / max(1, sig_hit + sig_miss), 6)},
        "other_tested": {"with_hit": rest_hit, "total": rest_hit + rest_miss,
                         "fraction": round(rest_hit / max(1, rest_hit + rest_miss), 6)},
        "odds_ratio": round(odds, 3),
        "fisher_p": fisher_greater(sig_hit, sig_miss, rest_hit, rest_miss),
    }


def main():
    groups_path, vdjdb_path, top_path = sys.argv[1:4]
    q_cut = float(sys.argv[4]) if len(sys.argv) > 4 else 0.05

    # Human TRB CDR3 -> antigen species and epitopes.
    species = defaultdict(set)
    epitopes = defaultdict(set)
    with open(vdjdb_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            if row["species"] == "HomoSapiens" and row["gene"] == "TRB":
                species[row["cdr3"]].add(row["antigen.species"])
                epitopes[row["cdr3"]].add(f'{row["antigen.species"]}:{row["antigen.epitope"]}')
    viral_other = {"EBV", "InfluenzaA", "SARS-CoV-2", "HIV-1", "YFV", "HCV", "HSV-2", "DENV1", "DENV2", "DENV3/4", "HPV"}
    n_cmv = sum("CMV" in s for s in species.values())
    n_other = sum(bool(s & viral_other) for s in species.values())

    counts = {k: [0, 0, 0, 0] for k in ("cmv", "other_virus")}  # sig_hit, sig_miss, rest_hit, rest_miss
    tested = significant = 0
    top = []
    sig_groups = []
    cap = 40  # member-count bins: 1..39, 40+
    bin_total = [0] * (cap + 1)
    bin_hits = {k: [0] * (cap + 1) for k in ("cmv", "other_virus")}
    with open(groups_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            p = float(row["assoc.p"])
            if math.isnan(p):
                continue
            tested += 1
            q = float(row["assoc.q"])
            sig = q <= q_cut
            significant += sig
            members = row["members"].split()
            hit_cmv = any("CMV" in species.get(m, ()) for m in members)
            hit_other = any(species.get(m, set()) & viral_other for m in members)
            for key, hit in (("cmv", hit_cmv), ("other_virus", hit_other)):
                counts[key][(0 if sig else 2) + (0 if hit else 1)] += 1
            size = min(len(members), cap)
            if not sig:
                bin_total[size] += 1
                bin_hits["cmv"][size] += hit_cmv
                bin_hits["other_virus"][size] += hit_other
            else:
                sig_groups.append({"tag": row["tag"], "p": p, "size": size, "members": set(members),
                                   "cmv": hit_cmv, "other_virus": hit_other,
                                   "case": int(row["case_donors"]), "control": int(row["control_donors"])})
                hits = sorted({e for m in members for e in epitopes.get(m, ()) if e.startswith("CMV:")})
                top.append((p, row["type"], row["tag"], int(row["case_donors"]), int(row["control_donors"]),
                            q, len(members), ";".join(hits)))

    top.sort()
    with open(top_path, "w") as out:
        out.write("type\ttag\tcase_donors\tcontrol_donors\tassoc.p\tassoc.q\tmembers\tvdjdb_cmv_epitopes\n")
        for p, kind, tag, case, ctrl, q, n, hits in top:
            out.write(f"{kind}\t{tag}\t{case}\t{ctrl}\t{p:.3e}\t{q:.3e}\t{n}\t{hits or '-'}\n")

    matched = {}
    for key in ("cmv", "other_virus"):
        expected = sum(bin_hits[key][g["size"]] / max(1, bin_total[g["size"]]) for g in sig_groups)
        observed = sum(g[key] for g in sig_groups)
        matched[key] = {"observed": observed, "expected_size_matched": round(expected, 3),
                        "obs_over_exp": round(observed / expected, 3) if expected else None,
                        "poisson_p": poisson_upper(observed, expected)}

    fams = families(sig_groups)
    fam_rows = []
    for fam in fams:
        best = min(fam, key=lambda g: g["p"])
        cdr3s = set().union(*(g["members"] for g in fam))
        hits = sorted({e for m in cdr3s for e in epitopes.get(m, ())})
        fam_rows.append((best["p"], best["tag"], len(fam), len(cdr3s), best["case"], best["control"],
                         any(g["cmv"] for g in fam), any(g["other_virus"] for g in fam), ";".join(hits) or "-",
                         " ".join(sorted(cdr3s))))
    fam_rows.sort()
    with open(top_path.replace(".tsv", "_families.tsv"), "w") as out:
        out.write("best_p\tbest_tag\tgroups\tcdr3s\tcase_donors\tcontrol_donors\tvdjdb_cmv\tvdjdb_other_virus\tvdjdb_epitopes\tmembers\n")
        for r in fam_rows:
            out.write("\t".join(f"{x:.3e}" if isinstance(x, float) else str(x) for x in r) + "\n")

    print(json.dumps({
        "q_cutoff": q_cut,
        "tested_groups": tested,
        "significant_groups": significant,
        "vdjdb_human_trb_cdr3s": len(species),
        "vdjdb_cmv_cdr3s": n_cmv,
        "vdjdb_other_virus_cdr3s": n_other,
        "cmv": table(*counts["cmv"]),
        "other_virus": table(*counts["other_virus"]),
        "size_matched": matched,
        "families": len(fams),
        "families_with_cmv_hit": sum(r[6] for r in fam_rows),
        "families_with_other_virus_hit": sum(r[7] for r in fam_rows),
    }))


if __name__ == "__main__":
    main()
