#!/bin/bash
#SBATCH --job-name=gliph2rs-validate
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --output=hpc/out/%x-%j.log

# Validates `gliph2-rs gliph2` against turboGliph::gliph2 on subsets of the
# HIP / Emerson 2017 cohort, with the GLIPH2 v2.0 CD4+CD8 reference.
# Requires setup_data.sbatch and prep_hip.sh to have run. Submit from repo root:
#   sbatch hpc/validate.sbatch            # sizes default to "10000 50000 100000"
#   SIZES="10000" sbatch hpc/validate.sbatch
# Results: hpc/out/validate-<jobid>.jsonl

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export RUSTUP_HOME=/scratch/$USER/rustup CARGO_HOME=/scratch/$USER/cargo
export PATH="$CARGO_HOME/bin:/scratch/$USER/envs/gliph-r/bin:$PATH"
export CARGO_TARGET_DIR=/scratch/$USER/gliph2-rs-target
cargo build --release --locked
BIN="$CARGO_TARGET_DIR/release/gliph2-rs"

DATA=/scratch/$USER/gliph2-rs-data
REF="$DATA/reference/human_v2.0_CD48.tsv"
WORK=/scratch/$USER/gliph2-rs-work/$SLURM_JOB_ID
RESULTS="hpc/out/validate-${SLURM_JOB_ID}.jsonl"
mkdir -p "$WORK" hpc/out
CPUS=${SLURM_CPUS_PER_TASK:-32}

# Subsets: unique CDR3s (first TRBV seen) from the first 40 subjects,
# shuffled with a fixed seed, so every tool sees identical input.
POOL="$WORK/pool_unique.tsv"
python3 - "$DATA/hip_pooled.tsv" "$POOL" <<'EOF'
import random, sys
seen, rows, subjects = set(), [], []
with open(sys.argv[1]) as f:
    next(f)
    for line in f:
        cdr3, v, _, _, subj, _ = line.rstrip("\n").split("\t")
        if not subjects or subjects[-1] != subj:
            subjects.append(subj)
            if len(subjects) > 40:
                break
        if cdr3 not in seen:
            seen.add(cdr3)
            rows.append((cdr3, v))
random.Random(20260930).shuffle(rows)
with open(sys.argv[2], "w") as out:
    out.write("CDR3b\tTRBV\n")
    out.writelines(f"{c}\t{v}\n" for c, v in rows)
print(f"pool: {len(rows)} unique CDR3s from {min(len(subjects), 40)} subjects", file=sys.stderr)
EOF

for n in ${SIZES:-10000 50000 100000}; do
  IN="$WORK/in_$n.tsv"
  head -n $((n + 1)) "$POOL" > "$IN"
  for params in paper default; do
    flag=""; [ "$params" = paper ] && flag="--paper-params"
    rs=$("$BIN" gliph2 --input "$IN" --reference "$REF" --out "$WORK/rs_${n}_$params" --threads "$CPUS" $flag)
    echo "{\"tool\":\"gliph2-rs\",\"params\":\"$params\",\"n\":$n,\"run\":$rs}" | tee -a "$RESULTS"
    # turboGliph's progress goes to the job log (timestamped); its last line is the JSON summary.
    tg=$(timeout "${TG_TIMEOUT:-20h}" Rscript hpc/run_turbogliph.R "$IN" "$REF" "$WORK/tg_${n}_$params.tsv" "$CPUS" "$params" \
         | while IFS= read -r l; do echo "[$(date +%T)] turboGliph: $l" >&2; echo "$l"; done | tail -1) || tg='{"error":"turboGliph failed or timed out"}'
    echo "{\"tool\":\"turboGliph\",\"params\":\"$params\",\"n\":$n,\"run\":$tg}" | tee -a "$RESULTS"
    if [ -f "$WORK/tg_${n}_$params.tsv" ]; then
      python3 hpc/compare_groups.py "$WORK/tg_${n}_$params.tsv" "$WORK/rs_${n}_${params}_groups.tsv" "n=$n params=$params" | tee -a "$RESULTS"
    fi
  done
done
