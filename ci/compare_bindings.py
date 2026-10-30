#!/usr/bin/env python3
"""Checks that the CLI, the R package and the Python package agree.

All three drive the same engine, so on one input they must report exactly the
same convergence groups with the same members. Run from the repository root:

    python ci/compare_bindings.py

Exits non-zero, naming the groups that differ, if any pair disagrees.
"""

from __future__ import annotations

import csv
import json
import pathlib
import subprocess
import sys
import tempfile

ROOT = pathlib.Path(__file__).resolve().parent.parent
CLI = ROOT / "target" / "release" / "gliph2-rs"


def run(cmd, **kw):
    return subprocess.run(cmd, check=True, capture_output=True, text=True, **kw)


def read_cdr3(path: pathlib.Path) -> list[str]:
    with path.open() as f:
        return [row["CDR3b"] for row in csv.DictReader(f, delimiter="\t")]


def groups_from_cli(out_prefix: pathlib.Path) -> set[tuple[str, str, str]]:
    csv.field_size_limit(sys.maxsize)
    with pathlib.Path(f"{out_prefix}_groups.tsv").open() as f:
        return {
            (r["type"], r["tag"], ",".join(sorted(r["members"].split())))
            for r in csv.DictReader(f, delimiter="\t")
        }


def groups_from_python(cdr3, reference) -> set[tuple[str, str, str]]:
    import gliph2rs

    g = gliph2rs.gliph2(cdr3, reference)["groups"]
    return {
        (t, tag, ",".join(sorted(m)))
        for t, tag, m in zip(g["type"], g["tag"], g["members"])
    }


def groups_from_r(input_tsv: pathlib.Path, ref_tsv: pathlib.Path, tmp: pathlib.Path):
    script = tmp / "run.R"
    out = tmp / "r_groups.json"
    script.write_text(f"""
        suppressPackageStartupMessages(library(gliph2rs))
        read_col <- function(p) read.delim(p, colClasses = "character")$CDR3b
        res <- gliph2rs::gliph2(read_col({str(input_tsv)!r}), read_col({str(ref_tsv)!r}))
        cp <- res$cluster_properties
        out <- mapply(function(type, tag, members) {{
            list(type = type, tag = tag, members = paste(sort(members), collapse = ","))
        }}, cp$type, cp$tag, res$cluster_list, SIMPLIFY = FALSE, USE.NAMES = FALSE)
        writeLines(jsonlite::toJSON(out, auto_unbox = TRUE), {str(out)!r})
    """)
    run(["Rscript", str(script)])
    return {(g["type"], g["tag"], g["members"]) for g in json.loads(out.read_text())}


def report(name_a, a, name_b, b) -> bool:
    if a == b:
        print(f"  {name_a} vs {name_b}: {len(a)} groups, identical")
        return True
    only_a, only_b = sorted(a - b), sorted(b - a)
    print(f"  {name_a} vs {name_b}: MISMATCH — {len(only_a)} only in {name_a}, {len(only_b)} only in {name_b}")
    for kind, rows in ((name_a, only_a), (name_b, only_b)):
        for t, tag, members in rows[:5]:
            print(f"     only in {kind}: {t} {tag} ({members.count(',') + 1} members)")
    return False


def main() -> int:
    if not CLI.exists():
        print(f"build the CLI first: cargo build --release ({CLI} missing)", file=sys.stderr)
        return 2

    with tempfile.TemporaryDirectory() as td:
        tmp = pathlib.Path(td)
        input_tsv, ref_tsv = tmp / "in.tsv", tmp / "ref.tsv"
        run([str(CLI), "synth", "--n", "4000", "--seed", "11", "--out", str(input_tsv)])
        run([str(CLI), "synth", "--reference", "--n", "60000", "--seed", "3", "--out", str(ref_tsv)])
        run([str(CLI), "gliph2", "--input", str(input_tsv), "--reference", str(ref_tsv),
             "--out", str(tmp / "cli")])

        cdr3, reference = read_cdr3(input_tsv), read_cdr3(ref_tsv)
        cli = groups_from_cli(tmp / "cli")
        python = groups_from_python(cdr3, reference)
        r = groups_from_r(input_tsv, ref_tsv, tmp)

        print(f"comparing {len(cli)} convergence groups across three interfaces")
        ok = report("CLI", cli, "Python", python)
        ok &= report("CLI", cli, "R", r)
        ok &= report("Python", python, "R", r)

    if not ok:
        print("\nthe interfaces disagree", file=sys.stderr)
        return 1
    print("\nall three interfaces agree")
    return 0


if __name__ == "__main__":
    sys.exit(main())
