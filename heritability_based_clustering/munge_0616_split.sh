#!/bin/bash -l
#SBATCH --account=stats
#SBATCH --job-name=LDSC_0616
#SBATCH --cpus-per-task=12
#SBATCH --time=1-00:00
#SBATCH --mem-per-cpu=10gb
#SBATCH --array=1-22%12
set -euo pipefail


# PYTHON="${PYTHON:-$HOME/.conda/envs/ldsc/bin/python}"
# activate existing conda env
  source /burg-archive/opt/anaconda3-2023.09/etc/profile.d/conda.sh
  conda activate ldsc


PHENO="0616"
N="31968"

BASE="$HOME/heritability"
SUMDIR="$BASE/UKB0616"
LDSC_DIR="$BASE/ldsc"
REFDIR="$BASE/eur_w_ld_chr"
MERGE="$REFDIR/w_hm3.snplist"


cd "$SUMDIR"

CHR=$(printf "%02d" "${SLURM_ARRAY_TASK_ID}")
FILE="${PHENO}_${CHR}.txt"
PREFIX="${PHENO}_${CHR}"

echo "Using python: $(which python)"
python --version
echo "Host: $(hostname)"
echo "Start time: $(date)"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}  CHR=${CHR}  FILE=${FILE}"

if [ ! -f "$FILE" ]; then
  echo "Missing file: $FILE" >&2
  exit 1
fi

~/.conda/envs/ldsc/bin/python "$LDSC_DIR/munge_sumstats.py" \
  --sumstats "$FILE" \
  --N "$N" \
  --out "${PREFIX}_ldsc" \
  --merge-alleles "$MERGE"

~/.conda/envs/ldsc/bin/python "$LDSC_DIR/ldsc.py" \
  --h2 "${PREFIX}_ldsc.sumstats.gz" \
  --ref-ld-chr "$REFDIR/" \
  --w-ld-chr  "$REFDIR/" \
  --out "${PREFIX}_h2"

echo "End time: $(date)"
