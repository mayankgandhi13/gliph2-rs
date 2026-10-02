#!/usr/bin/env python3
"""Compares two convergence-group tables (type, tag, ..., members) by tag.

Usage: compare_groups.py REFERENCE.tsv CANDIDATE.tsv [LABEL]
Prints one JSON line: tag overlap, exact member-set agreement, mean Jaccard
over shared tags, and per-type breakdown. Mismatch examples go to stderr.
"""
import csv
import json
import sys

csv.field_size_limit(sys.maxsize)


def load(path):
    groups = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            groups[(row["type"], row["tag"])] = frozenset(row["members"].split())
    return groups


def summarize(ref, cand, kinds):
    ref = {k: v for k, v in ref.items() if k[0] in kinds}
    cand = {k: v for k, v in cand.items() if k[0] in kinds}
    shared = ref.keys() & cand.keys()
    exact = sum(ref[k] == cand[k] for k in shared)
    jac = [len(ref[k] & cand[k]) / len(ref[k] | cand[k]) for k in shared]
    # Sequence-pair view: do the same pairs of CDR3s end up grouped together?
    return {
        "ref_groups": len(ref),
        "cand_groups": len(cand),
        "shared_tags": len(shared),
        "only_ref": len(ref.keys() - cand.keys()),
        "only_cand": len(cand.keys() - ref.keys()),
        "identical_members": exact,
        "mean_jaccard_shared": round(sum(jac) / len(jac), 6) if jac else None,
        "exact_match_rate": round(exact / len(ref), 6) if ref else None,
    }


def main():
    ref, cand = load(sys.argv[1]), load(sys.argv[2])
    out = {"label": sys.argv[3] if len(sys.argv) > 3 else ""}
    out["all"] = summarize(ref, cand, {"local", "global"})
    out["local"] = summarize(ref, cand, {"local"})
    out["global"] = summarize(ref, cand, {"global"})
    print(json.dumps(out))

    for name, keys in (("only in reference", ref.keys() - cand.keys()), ("only in candidate", cand.keys() - ref.keys())):
        for k in sorted(keys)[:5]:
            src = ref if k in ref else cand
            print(f"{name}: {k[0]} {k[1]} n={len(src[k])}", file=sys.stderr)
    diff = [k for k in sorted(ref.keys() & cand.keys()) if ref[k] != cand[k]]
    for k in diff[:5]:
        print(f"members differ: {k[1]}: -{sorted(ref[k] - cand[k])[:3]} +{sorted(cand[k] - ref[k])[:3]}", file=sys.stderr)


if __name__ == "__main__":
    main()
