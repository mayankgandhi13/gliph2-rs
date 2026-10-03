#!/bin/bash
#SBATCH --job-name=gliph2rs-filters
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=hpc/out/%x-%j.log

# Effect of the global-group filters on the full pooled HIP / Emerson 2017
# cohort (paper settings). Submit from the repo root:
#   sbatch hpc/full_cohort_filters.sbatch
# Results: hpc/out/filters-<jobid>.jsonl

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export RUSTUP_HOME=/scratch/$USER/rustup CARGO_HOME=/scratch/$USER/cargo
export PATH="$CARGO_HOME/bin:$PATH"
export CARGO_TARGET_DIR=/scratch/$USER/gliph2-rs-target
cargo build --release --locked
BIN="$CARGO_TARGET_DIR/release/gliph2-rs"

DATA=/scratch/$USER/gliph2-rs-data
REF="$DATA/reference/human_v2.0_CD48.tsv"
OUT=/scratch/$USER/gliph2-rs-work/filters-$SLURM_JOB_ID
RESULTS="hpc/out/filters-${SLURM_JOB_ID}.jsonl"
mkdir -p "$OUT" hpc/out
CPUS=${SLURM_CPUS_PER_TASK:-64}

run() { # label, extra flags...
  local label=$1; shift
  local out
  out=$(/usr/bin/time -f '%M %e' -o "$OUT/time_$label" "$BIN" gliph2 --input "$DATA/hip_pooled.tsv" \
        --reference "$REF" --paper-params --threads "$CPUS" --out "$OUT/$label" "$@")
  read -r rss wall < "$OUT/time_$label"
  local bytes
  bytes=$(stat -c %s "$OUT/${label}_groups.tsv")
  echo "{\"label\":\"$label\",\"flags\":\"$*\",\"max_rss_kb\":$rss,\"wall_s\":$wall,\"output_bytes\":$bytes,\"run\":$out}" | tee -a "$RESULTS"
}

run fdr05 --max-global-q 0.05
run subjects3 --min-subjects 3
run fdr05_subjects3 --max-global-q 0.05 --min-subjects 3
ls -la "$OUT"
