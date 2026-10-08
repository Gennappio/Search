#!/bin/bash
#SBATCH --job-name=gbm_buffa_test
#SBATCH --partition=defq
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/logs/gbm_buffa_test_%j.out
#SBATCH --error=/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/logs/gbm_buffa_test_%j.err

set -euo pipefail

BASE_DIR="/home/MoscatoG/GBM_HYPOXIA_SIGNATURE"
ENV_PREFIX="$BASE_DIR/envs/gbm_buffa_r"
mkdir -p "$BASE_DIR/logs"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

conda activate "$ENV_PREFIX"
cd "$BASE_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1

Rscript "$BASE_DIR/CODE/run_buffa2010_hpc.R" \
  --base_dir "$BASE_DIR" \
  --b_boot 20 \
  --n_mc 50 \
  --cores "${SLURM_CPUS_PER_TASK:-8}" \
  --boot_chunk_size 10 \
  --mc_chunk_size 10 \
  --resume FALSE \
  --test_mode TRUE
