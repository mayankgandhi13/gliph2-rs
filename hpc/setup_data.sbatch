#!/bin/bash
#SBATCH --job-name=gliph2rs-setup
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=hpc/out/%x-%j.log

# One-time setup on Explorer: downloads public data and the turboGliph
# comparator. Submit from the repo root:  sbatch hpc/setup_data.sbatch
#
# Data:
#   HIP / Emerson et al. 2017 (Nat Genet, doi:10.1038/ng.3822): 786 TCRβ
#     repertoires with CMV and HLA labels, via the huggingface.co/datasets/isalgo/airr_hip
#     mirror of immuneACCESS. Observe the immuneACCESS data-use terms.
#   GLIPH2 reference repertoires (Huang et al. 2020): zenodo.org/records/18925758

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-.}"
DATA="/scratch/$USER/gliph2-rs-data"
mkdir -p hpc/out "$DATA/hip/corr" "$DATA/reference"

HF=https://huggingface.co/datasets/isalgo/airr_hip/resolve/main
for f in metadata.txt metadata.cmv.txt metadata_A02.txt pool.aa.summary.txt; do
  curl -sSfL --retry 5 -o "$DATA/hip/$f" "$HF/$f"
done
curl -sSfL --retry 5 "https://huggingface.co/api/datasets/isalgo/airr_hip/tree/main/corr" \
  | python3 -c 'import json,sys; [print(x["path"]) for x in json.load(sys.stdin)]' > "$DATA/hip/files.txt"
echo "downloading $(wc -l < "$DATA/hip/files.txt") repertoires"
# -C - resumes partial files, so the job is safe to rerun.
xargs -P 8 -I{} curl -sSfL --retry 5 -C - -o "$DATA/hip/{}" "$HF/{}" < "$DATA/hip/files.txt"
echo "hip: $(ls "$DATA/hip/corr" | wc -l) files, $(du -sh "$DATA/hip/corr" | cut -f1)"

[ -s "$DATA/reference/reference_list.RData" ] || curl -sSfL --retry 5 -o "$DATA/reference/reference_list.RData" \
  "https://zenodo.org/records/18925758/files/reference_list.RData?download=1"

# VDJdb (Goncharov et al. 2022), used by the CMV association checks.
mkdir -p "$DATA/vdjdb"
if [ ! -s "$DATA/vdjdb/vdjdb.slim.txt" ]; then
  curl -sSfL --retry 5 -o "$DATA/vdjdb/vdjdb.zip" \
    https://github.com/antigenomics/vdjdb-db/releases/download/2026-06-03-ZENODO/vdjdb-2026-06-03.zip
  unzip -o -j "$DATA/vdjdb/vdjdb.zip" '*vdjdb.slim.txt' -d "$DATA/vdjdb"
fi

# turboGliph (R implementation of GLIPH2) as the comparator, in a conda-forge
# R so binary deps (igraph etc.) don't have to compile against the module R.
# immGLIPH, its successor, needs R >= 4.5 and Bioconductor devel.
ENV="/scratch/$USER/envs/gliph-r"
module load miniconda3/25.9.1
# /scratch is purged, which can leave the env present but gutted (empty
# lib/R/bin/exec), so test that R actually runs rather than that it exists.
if ! "$ENV/bin/Rscript" -e 'invisible(1)' >/dev/null 2>&1; then
  echo "creating conda R environment at $ENV"
  rm -rf "$ENV"
  conda create -y -p "$ENV" -c conda-forge --override-channels \
    r-base=4.4 r-igraph r-stringdist r-dplyr r-visnetwork r-viridis r-stringr \
    r-doparallel r-foreach r-remotes r-data.table
fi
export PATH="$ENV/bin:$PATH"
export R_LIBS_USER="$ENV/lib/R/library"
Rscript -e '
  options(repos = c(CRAN = "https://cloud.r-project.org"), Ncpus = 8)
  # grr was archived on CRAN; install its last release from the archive.
  if (!requireNamespace("grr", quietly = TRUE))
    install.packages("https://cran.r-project.org/src/contrib/Archive/grr/grr_0.9.5.tar.gz", repos = NULL, type = "source")
  if (!requireNamespace("plotfunctions", quietly = TRUE)) install.packages("plotfunctions")
  if (!requireNamespace("turboGliph", quietly = TRUE)) remotes::install_github("HetzDra/turboGliph", upgrade = "never")
  library(turboGliph); cat("turboGliph", as.character(packageVersion("turboGliph")), "\n")'

Rscript hpc/export_reference.R "$DATA/reference/reference_list.RData" "$DATA/reference"
ls -la "$DATA/reference"
