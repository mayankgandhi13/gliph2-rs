#!/bin/sh
# Builds the PyPI artefacts and validates their metadata.
#
# Usage: tools/release.sh [OUTPUT_DIR]   (from the package root)
#
# Wheels are per-interpreter and per-platform; this builds one for whichever
# Python is active. Use maturin-action in CI for the full matrix.

set -eu
out=${1:-dist}

# maturin picks the interpreter from VIRTUAL_ENV, which a stray CONDA_PREFIX
# silently overrides — build against the wrong Python and the module will not
# import.
if [ -n "${VIRTUAL_ENV:-}" ] && [ -n "${CONDA_PREFIX:-}" ]; then
    echo "Both VIRTUAL_ENV and CONDA_PREFIX are set; unset one so maturin picks the right interpreter." >&2
    exit 1
fi

maturin sdist -o "$out"
maturin build --release -o "$out"
python -m twine check "$out"/*
echo "artefacts in $out:"
ls -lh "$out"
