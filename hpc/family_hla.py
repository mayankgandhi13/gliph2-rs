#!/usr/bin/env python3
"""HLA restriction of CMV-associated CDR3 families.

Usage: family_hla.py FAMILIES.tsv POOLED.tsv HLA.tsv OUT.tsv

FAMILIES.tsv is the *_families.tsv written by assoc_vdjdb.py (one family per
row, member CDR3s space-separated). POOLED.tsv is the GLIPH2-format cohort
table, HLA.tsv lists `subject:condition` then the subject's HLA alleles.

For each family, the carriers are the HLA-typed donors with any member CDR3.
Every allele seen among carriers is tested (one-sided Fisher) for being more
common in carriers than in non-carriers; Benjamini-Hochberg q-values are taken
over all family x allele tests. Antigen-specific public TCRs are expected to
be restricted by one HLA allele. Writes the best allele per family to OUT.tsv
and prints a JSON summary.
"""
import csv
import json
import math
import sys
from collections import defaultdict


def log_choose(n, k):
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_greater(a, b, c, d):
    n1, k, n = a + b, a + c, a + b + c + d
    denom = log_choose(n, n1)
    return min(1.0, sum(math.exp(log_choose(k, x) + log_choose(n - k, n1 - x) - denom) for x in range(a, min(k, n1) + 1)))


def bh(p):
    order = sorted(range(len(p)), key=lambda i: (p[i], i))
    q = [0.0] * len(p)
    running = 1.0
    for rank in range(len(order) - 1, -1, -1):
        i = order[rank]
        running = min(running, p[i] * len(p) / (rank + 1))
        q[i] = running
    return q


def main():
    fam_path, pooled_path, hla_path, out_path = sys.argv[1:5]

    families = []
    with open(fam_path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            families.append({"tag": row["best_tag"], "cdr3s": set(row["members"].split()),
                             "case": int(row["case_donors"]), "control": int(row["control_donors"])})
    cdr3_fams = defaultdict(list)
    for i, fam in enumerate(families):
        for c in fam["cdr3s"]:
            cdr3_fams[c].append(i)

    hla = {}
    with open(hla_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            hla[parts[0]] = {a for a in parts[1:] if a}

    carriers = [set() for _ in families]
    with open(pooled_path) as f:
        next(f)
        for line in f:
            cdr3, rest = line.split("\t", 1)
            for i in cdr3_fams.get(cdr3, ()):
                carriers[i].add(rest.split("\t")[3])

    typed = set(hla)
    tests = []  # (family index, allele, a, b, c, d, p)
    for i, fam in enumerate(families):
        car = carriers[i] & typed
        non = typed - car
        alleles = set().union(*(hla[s] for s in car)) if car else set()
        for allele in sorted(alleles):
            a = sum(allele in hla[s] for s in car)
            c = sum(allele in hla[s] for s in non)
            tests.append((i, allele, a, len(car) - a, c, len(non) - c, fisher_greater(a, len(car) - a, c, len(non) - c)))
    qs = bh([t[6] for t in tests])

    best = {}
    for t, q in zip(tests, qs):
        i = t[0]
        if i not in best or t[6] < best[i][0][6]:
            best[i] = (t, q)

    with open(out_path, "w") as out:
        out.write("family\tcase_donors\tcontrol_donors\ttyped_carriers\tbest_allele\tcarriers_with_allele\t"
                  "carrier_frac\tbackground_frac\tp\tq\n")
        for i, fam in enumerate(families):
            if i not in best:
                continue
            (_, allele, a, b, c, d, p), q = best[i]
            out.write(f"{fam['tag']}\t{fam['case']}\t{fam['control']}\t{a + b}\t{allele}\t{a}\t"
                      f"{a / max(1, a + b):.3f}\t{c / max(1, c + d):.3f}\t{p:.3e}\t{q:.3e}\n")

    print(json.dumps({
        "families": len(families),
        "hla_typed_donors": len(typed),
        "tests": len(tests),
        "families_with_hla_q05": sum(1 for i in best if best[i][1] <= 0.05),
    }))


if __name__ == "__main__":
    main()
