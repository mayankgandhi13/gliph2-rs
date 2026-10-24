#!/bin/sh
# Copies the gliph2-rs core crate into src/rust/core, so the package carries
# the Rust sources it compiles.
#
# In a repository checkout the crate sits three levels above the package and
# is copied in (refreshing any older copy). A source tarball has no repository
# above it and ships the copy instead, so this is then a no-op.
#
# Usage: tools/sync-core.sh   (from anywhere)

set -eu

pkg_root=$(cd "$(dirname "$0")/.." && pwd)
vendor="$pkg_root/src/rust/core"
core="$pkg_root/../.."

if [ -f "$core/Cargo.toml" ]; then
    mkdir -p "$vendor"
    cp "$core/Cargo.toml" "$vendor/Cargo.toml"
    rm -rf "$vendor/src"
    cp -R "$core/src" "$vendor/src"
elif [ ! -f "$vendor/Cargo.toml" ]; then
    echo "ERROR: core crate missing at $core and not vendored in $vendor" >&2
    exit 1
fi
