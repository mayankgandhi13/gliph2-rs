# Sourced by the Slurm scripts: points at the Rust toolchain on /scratch and
# repairs a gutted crate cache.
#
# Explorer purges /scratch, which can leave a directory tree in place with its
# files removed. For cargo that means an extracted crate whose Cargo.toml is
# gone, and every build then fails with "unable to get packages from source".
# Dropping registry/src makes cargo re-extract from registry/cache, or
# re-download if those are gone too.

export RUSTUP_HOME="${RUSTUP_HOME:-/scratch/$USER/rustup}"
export CARGO_HOME="${CARGO_HOME:-/scratch/$USER/cargo}"
export PATH="$CARGO_HOME/bin:$PATH"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-/scratch/$USER/gliph2-rs-target}"

for crate in "$CARGO_HOME"/registry/src/*/*/; do
  [ -d "$crate" ] || continue
  if [ ! -f "$crate/Cargo.toml" ]; then
    echo "crate cache is incomplete ($crate); clearing $CARGO_HOME/registry/src" >&2
    rm -rf "${CARGO_HOME:?}/registry/src"
    break
  fi
done
