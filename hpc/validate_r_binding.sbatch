#!/bin/bash
#SBATCH --job-name=gliph2rs-rbinding
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=hpc/out/%x-%j.log

# Installs the gliph2rs R package on Linux and compares it with
# turboGliph::gliph2() inside one R session, on the same Emerson 2017 subsets.
# Requires setup_data.sbatch and prep_hip.sh. Submit from the repo root:
#   sbatch hpc/validate_r_binding.sbatch
# Results: hpc/out/rbinding-<jobid>.jsonl

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
export RUSTUP_HOME=/scratch/$USER/rustup CARGO_HOME=/scratch/$USER/cargo
module load miniconda3/25.9.1
ENV=/scratch/$USER/envs/gliph-r
export PATH="$CARGO_HOME/bin:$ENV/bin:$PATH"
export R_LIBS_USER="$ENV/lib/R/library"
export NOT_CRAN=true

DATA=/scratch/$USER/gliph2-rs-data
REF="$DATA/reference/human_v2.0_CD48.tsv"
WORK=/scratch/$USER/gliph2-rs-work/rbinding-$SLURM_JOB_ID
RESULTS="hpc/out/rbinding-${SLURM_JOB_ID}.jsonl"
mkdir -p "$WORK" hpc/out
CPUS=${SLURM_CPUS_PER_TASK:-32}

R CMD INSTALL --no-multiarch --preclean bindings/r
Rscript -e 'cat("gliph2rs", as.character(packageVersion("gliph2rs")), "\n")'

# Same shuffled pool of unique CDR3s the CLI validation uses.
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
EOF

for n in ${SIZES:-2000 10000 50000}; do
  head -n $((n + 1)) "$POOL" > "$WORK/in_$n.tsv"
  Rscript hpc/compare_r_binding.R "$WORK/in_$n.tsv" "$REF" "$CPUS" | tee -a "$RESULTS"
done
