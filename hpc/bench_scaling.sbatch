#!/bin/bash
#SBATCH --job-name=gliph2rs-scaling
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=hpc/out/%x-%j.log

# Scaling benchmark on Explorer. Compute nodes have no system cargo and the
# login node kills compiles, so use a rustup toolchain on /scratch (shared):
#   RUSTUP_HOME=/scratch/$USER/rustup CARGO_HOME=/scratch/$USER/cargo \
#     sh <(curl -sSf https://sh.rustup.rs) -y --no-modify-path --profile minimal
# then submit from the repo root:
#   sbatch hpc/bench_scaling.sbatch
# Results: hpc/out/scaling-<jobid>.jsonl (one JSON object per run).

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
mkdir -p hpc/out
export RUSTUP_HOME="${RUSTUP_HOME:-/scratch/$USER/rustup}"
export CARGO_HOME="${CARGO_HOME:-/scratch/$USER/cargo}"
export PATH="$CARGO_HOME/bin:$PATH"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-/scratch/$USER/gliph2-rs-target}"

cargo build --release --locked
BIN="$CARGO_TARGET_DIR/release/gliph2-rs"
DATA="/scratch/$USER/gliph2-rs-data"
WORK="/scratch/$USER/gliph2-rs-work/$SLURM_JOB_ID"
mkdir -p "$DATA" "$WORK"
RESULTS="hpc/out/scaling-${SLURM_JOB_ID}.jsonl"

echo "node=$(hostname) cpus=${SLURM_CPUS_PER_TASK} $(lscpu | grep 'Model name' | sed 's/  */ /g')"

REF="$DATA/ref_1M.tsv"
[ -f "$REF" ] || "$BIN" synth --reference --n 1000000 --seed 42 --out "$REF"

run() { # label threads input
  local out
  out=$(/usr/bin/time -f '%M' -o "$WORK/mem" "$BIN" cluster --input "$3" --reference "$REF" \
        --threads "$2" --out "$WORK/$1_t$2")
  echo "${out%\}},\"label\":\"$1\",\"max_rss_kb\":$(cat "$WORK/mem")}" | tee -a "$RESULTS"
}

# Size scaling at full thread count.
for n in 10000 100000 1000000 10000000; do
  IN="$DATA/rep_${n}.tsv"
  [ -f "$IN" ] || "$BIN" synth --n "$n" --seed 7 --out "$IN"
  run "n$n" "$SLURM_CPUS_PER_TASK" "$IN"
done

# Thread scaling at 1M sequences.
for t in 1 2 4 8 16 32; do
  run "n1000000" "$t" "$DATA/rep_1000000.tsv"
done

# Reproducibility: repeated runs must produce byte-identical output.
for i in 1 2 3; do
  "$BIN" cluster --input "$DATA/rep_1000000.tsv" --reference "$REF" \
    --threads "$SLURM_CPUS_PER_TASK" --out "$WORK/repro$i" > /dev/null
done
sha256sum "$WORK"/repro*_clusters.tsv "$WORK"/repro*_motifs.tsv | tee -a "hpc/out/repro-${SLURM_JOB_ID}.txt"
