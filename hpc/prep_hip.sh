#!/bin/bash
# Converts the HIP / Emerson 2017 per-subject VDJtools tables into one pooled
# GLIPH2-format TSV (CDR3b TRBV TRBJ CDR3a subject:condition count), plus a
# GLIPH2-style HLA table. Out-of-frame CDR3s (containing _ or *) are dropped.
#
# Usage: hpc/prep_hip.sh DATA_DIR      (DATA_DIR contains hip/ from setup_data.sbatch)
# Output: DATA_DIR/hip_pooled.tsv, DATA_DIR/hip_hla.tsv

set -euo pipefail
DATA="$1"
HIP="$DATA/hip"
POOLED="$DATA/hip_pooled.tsv"
PARTS="$DATA/hip_parts"
mkdir -p "$PARTS"

# sample_id -> CMV label (+, -, NA), keyed by file name.
tail -n +2 "$HIP/metadata.txt" | awk -F'\t' '{print $1 "\t" $2 ":CMV" $6}' > "$PARTS/labels.tsv"

convert() { # file_name subject:condition
  zcat "$HIP/$1" | awk -F'\t' -v subj="$2" 'NR > 1 && $4 !~ /[_*]/ && $4 != "" {
    split($5, v, ","); split($7, j, ",");
    print $4 "\t" v[1] "\t" j[1] "\tNA\t" subj "\t" $1 }' > "$PARTS/$(basename "$1" .txt.gz).tsv"
}
export -f convert
export HIP PARTS
xargs -P "${SLURM_CPUS_PER_TASK:-8}" -L 1 bash -c 'convert "$0" "$1"' < "$PARTS/labels.tsv"

# Concatenate in metadata order so the pooled file is reproducible.
{
  printf 'CDR3b\tTRBV\tTRBJ\tCDR3a\tsubject:condition\tcount\n'
  cut -f1 "$PARTS/labels.tsv" | while read -r f; do cat "$PARTS/$(basename "$f" .txt.gz).tsv"; done
} > "$POOLED"

# HLA table: subject:condition followed by alleles, one per column.
tail -n +2 "$HIP/metadata.txt" | awk -F'\t' '$7 != "NA" && $7 != "" {
  gsub(",", "\t", $7); print $2 ":CMV" $6 "\t" $7 }' > "$DATA/hip_hla.tsv"

rm -rf "$PARTS"
echo "pooled: $(($(wc -l < "$POOLED") - 1)) rows -> $POOLED"
echo "hla: $(wc -l < "$DATA/hip_hla.tsv") subjects -> $DATA/hip_hla.tsv"
