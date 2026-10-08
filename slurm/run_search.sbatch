#!/bin/bash
#SBATCH --job-name=search
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=logs/search_%j.out
#SBATCH --error=logs/search_%j.err

# Template. Set BASE_DIR to the project folder holding the cohort manifest and the data.
set -euo pipefail
BASE_DIR="${BASE_DIR:-$PWD}"
ENV_PREFIX="${ENV_PREFIX:-$BASE_DIR/envs/search}"
MANIFEST="${MANIFEST:-$BASE_DIR/cohorts.yaml}"
PRESET="${PRESET:-buffa2010}"
OUT="${OUT:-$BASE_DIR/runs}"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
fi
conda activate "$ENV_PREFIX"
cd "$BASE_DIR"
mkdir -p logs

# one BLAS thread per worker process; Search parallelises over replicates
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

search run --manifest "$MANIFEST" --preset "$PRESET" --out "$OUT" \
  --n-jobs "${SLURM_CPUS_PER_TASK:-1}" --resume "$@"
