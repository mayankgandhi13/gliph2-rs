#!/bin/sh
# Regenerates src/rust/vendor.tar.xz, the bundled Cargo dependencies.
#
# CRAN requires that a package's Rust dependencies ship with it and that the
# build does not download anything, so the crates listed in Cargo.lock are
# vendored and compressed here. Re-run after changing Cargo.toml or updating
# Cargo.lock, then commit the result.
#
# Usage: tools/vendor.sh   (from the package root)

set -eu

pkg_root=$(cd "$(dirname "$0")/.." && pwd)
rust_dir="$pkg_root/src/rust"
staging=$(mktemp -d)
trap 'rm -rf "$staging"' EXIT

# The core crate is a path dependency and ships as src/rust/core, so make sure
# it is present before resolving the lock file.
"$pkg_root/tools/sync-core.sh"

cd "$rust_dir"
cargo vendor --versioned-dirs --locked "$staging/vendor" > /dev/null

# Ship the bundle exactly as cargo produced it, so its checksums stay valid.
(cd "$staging" && tar cJf "$rust_dir/vendor.tar.xz" vendor)
# Keep the licence listing in step with what is actually bundled.
"$pkg_root/tools/authors.sh" "$staging/vendor"

echo "wrote $rust_dir/vendor.tar.xz ($(du -h "$rust_dir/vendor.tar.xz" | cut -f1), $(ls "$staging/vendor" | wc -l | tr -d ' ') crates)"
