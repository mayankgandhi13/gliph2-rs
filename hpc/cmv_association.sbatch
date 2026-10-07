#!/bin/bash
#SBATCH --job-name=gliph2rs-cmv
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=04:00:00
#SBATCH --output=hpc/out/%x-%j.log

# CMV association on the full pooled HIP / Emerson 2017 cohort: every GLIPH2
# convergence group (paper settings) is tested for enrichment in CMV+ over
# CMV- donors, then checked against VDJdb's CMV-specific TCRs (raw and
# size-matched) and for HLA restriction. A second run with shuffled CMV labels
# is the calibration control.
# Requires setup_data.sbatch and prep_hip.sh. Submit from the repo root:
#   sbatch hpc/cmv_association.sbatch
# Results: hpc/out/cmv-<jobid>.jsonl, top groups in $OUT/top_cmv_groups.tsv

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export RUSTUP_HOME=/scratch/$USER/rustup CARGO_HOME=/scratch/$USER/cargo
export PATH="$CARGO_HOME/bin:$PATH"
export CARGO_TARGET_DIR=/scratch/$USER/gliph2-rs-target
cargo build --release --locked
BIN="$CARGO_TARGET_DIR/release/gliph2-rs"

DATA=/scratch/$USER/gliph2-rs-data
REF="$DATA/reference/human_v2.0_CD48.tsv"
OUT=/scratch/$USER/gliph2-rs-work/cmv-$SLURM_JOB_ID
RESULTS="hpc/out/cmv-${SLURM_JOB_ID}.jsonl"
mkdir -p "$OUT" hpc/out "$DATA/vdjdb"
CPUS=${SLURM_CPUS_PER_TASK:-64}

# VDJdb (Goncharov et al. 2022), release 2026-06-03.
VDJDB="$DATA/vdjdb/vdjdb.slim.txt"
if [ ! -f "$VDJDB" ]; then
  curl -sSfL --retry 5 -o "$DATA/vdjdb/vdjdb.zip" \
    https://github.com/antigenomics/vdjdb-db/releases/download/2026-06-03-ZENODO/vdjdb-2026-06-03.zip
  unzip -o -j "$DATA/vdjdb/vdjdb.zip" '*vdjdb.slim.txt' -d "$DATA/vdjdb"
fi

run() { # label, extra flags...
  local label=$1; shift
  local out
  out=$(/usr/bin/time -f '%M %e' -o "$OUT/time_$label" "$BIN" gliph2 --input "$DATA/hip_pooled.tsv" \
        --reference "$REF" --paper-params --threads "$CPUS" --min-subjects 10 \
        --case CMV+ --control CMV- --assoc-min-donors 10 --out "$OUT/$label" "$@")
  read -r rss wall < "$OUT/time_$label"
  echo "{\"label\":\"$label\",\"flags\":\"$*\",\"max_rss_kb\":$rss,\"wall_s\":$wall,\"run\":$out}" | tee -a "$RESULTS"
}

run cmv
run cmv_permuted --assoc-permute-seed 20261006

for label in cmv cmv_permuted; do
  summary=$(python3 hpc/assoc_vdjdb.py "$OUT/${label}_groups.tsv" "$VDJDB" "$OUT/top_${label}_groups.tsv" 0.05)
  echo "{\"label\":\"vdjdb_$label\",\"summary\":$summary}" | tee -a "$RESULTS"
done

# HLA restriction of the CMV-associated families (626 HLA-typed donors).
hla=$(python3 hpc/family_hla.py "$OUT/top_cmv_groups_families.tsv" "$DATA/hip_pooled.tsv" "$DATA/hip_hla.tsv" "$OUT/family_hla.tsv")
echo "{\"label\":\"family_hla\",\"summary\":$hla}" | tee -a "$RESULTS"
column -t -s $'\t' "$OUT/family_hla.tsv"
