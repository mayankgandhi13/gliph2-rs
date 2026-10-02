#!/bin/bash
#SBATCH --job-name=gliph2rs-full
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=240G
#SBATCH --time=12:00:00
#SBATCH --output=hpc/out/%x-%j.log

# Runs gliph2-rs on the entire pooled HIP / Emerson 2017 cohort (786 subjects,
# ~151M rows) — the scale at which the original GLIPH2 is impractical.
# Requires setup_data.sbatch and prep_hip.sh. Submit from the repo root:
#   sbatch hpc/run_full_cohort.sbatch
# Results: hpc/out/full-<jobid>.jsonl

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export RUSTUP_HOME=/scratch/$USER/rustup CARGO_HOME=/scratch/$USER/cargo
export PATH="$CARGO_HOME/bin:$PATH"
export CARGO_TARGET_DIR=/scratch/$USER/gliph2-rs-target
cargo build --release --locked
BIN="$CARGO_TARGET_DIR/release/gliph2-rs"

DATA=/scratch/$USER/gliph2-rs-data
REF="$DATA/reference/human_v2.0_CD48.tsv"
OUT=/scratch/$USER/gliph2-rs-work/full-$SLURM_JOB_ID
RESULTS="hpc/out/full-${SLURM_JOB_ID}.jsonl"
mkdir -p "$OUT" hpc/out
CPUS=${SLURM_CPUS_PER_TASK:-64}
echo "node=$(hostname) cpus=$CPUS $(lscpu | grep 'Model name' | sed 's/  */ /g')"

timed() { # label, command...
  local label=$1; shift
  local out
  out=$(/usr/bin/time -f '%M %e' -o "$OUT/time_$label" "$@")
  read -r rss wall < "$OUT/time_$label"
  echo "{\"label\":\"$label\",\"max_rss_kb\":$rss,\"wall_s\":$wall,\"run\":$out}" | tee -a "$RESULTS"
}

# GLIPH2 convergence groups with the GLIPH2 parameter-file settings.
timed gliph2_paper "$BIN" gliph2 --input "$DATA/hip_pooled.tsv" --reference "$REF" \
  --paper-params --threads "$CPUS" --out "$OUT/gliph2_paper"

# Merged-cluster mode with FDR-controlled motifs.
timed cluster_fdr "$BIN" cluster --input "$DATA/hip_pooled.tsv" --reference "$REF" \
  --fdr --mismatch-trim 3 --threads "$CPUS" --out "$OUT/cluster_fdr"

ls -la "$OUT"
