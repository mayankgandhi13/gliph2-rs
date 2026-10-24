#!/bin/sh
# Writes inst/AUTHORS from a vendored crate directory, so the bundled
# dependencies' authors and licences are listed as CRAN requires.
#
# Usage: tools/authors.sh VENDOR_DIR   (called by tools/vendor.sh)

set -eu
pkg_root=$(cd "$(dirname "$0")/.." && pwd)
mkdir -p "$pkg_root/inst"
VENDOR_DIR="$1" OUT="$pkg_root/inst/AUTHORS" python3 - <<'PY'
import os, re
root, out = os.environ["VENDOR_DIR"], os.environ["OUT"]
rows = []
for d in sorted(os.listdir(root)):
    toml = os.path.join(root, d, "Cargo.toml")
    if not os.path.exists(toml):
        continue
    text = open(toml, encoding="utf-8", errors="replace").read()
    def field(name):
        m = re.search(rf'^{name}\s*=\s*"([^"]*)"', text, re.M)
        return m.group(1) if m else ""
    m = re.search(r'^authors\s*=\s*\[(.*?)\]', text, re.M | re.S)
    auth = "; ".join(n.split(" <")[0] for n in re.findall(r'"([^"]+)"', m.group(1))) if m else ""
    rows.append((field("name") or d, field("version"), field("license"), auth))

lines = [
    "The R code and the gliph2-rs core crate in this package are",
    "Copyright (c) 2026 Mayank Gandhi, released under the MIT licence",
    "(see the LICENSE file).",
    "",
    "The package bundles its Rust dependencies in src/rust/vendor.tar.xz so",
    "that it builds without network access. Those crates remain the work and",
    "copyright of their own authors, under the licences shown below. Each",
    "crate's full licence text is included in its directory inside the",
    "bundle.",
    "",
]
w = max(len(r[0]) for r in rows) + 1
for name, version, lic, auth in rows:
    lines.append(f"{name:<{w}} {version:<10} {lic}")
    if auth:
        lines.append(f"{'':<{w}} {auth}")
open(out, "w").write("\n".join(lines) + "\n")
print(f"{len(rows)} crates listed in inst/AUTHORS")
PY
