#!/bin/bash
#SBATCH --job-name=gbm_buffa_consensus
#SBATCH --partition=compute
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --output=/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/logs/gbm_buffa_consensus_%j.out
#SBATCH --error=/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/logs/gbm_buffa_consensus_%j.err

set -euo pipefail

BASE_DIR="/home/MoscatoG/GBM_HYPOXIA_SIGNATURE"
ENV_PREFIX="$BASE_DIR/envs/gbm_buffa_r"

ANNOTATION_FILE="$BASE_DIR/DATA/master_annotation_ITFR_Ensembl104_TCGA_GDC_consensus.csv"

mkdir -p "$BASE_DIR/logs" "$BASE_DIR/RESULTS" "$BASE_DIR/CODE" "$BASE_DIR/SLURM"

if [ -f "$HOME/miniconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/miniconda3/etc/profile.d/conda.sh"
elif [ -f "$HOME/anaconda3/etc/profile.d/conda.sh" ]; then
  source "$HOME/anaconda3/etc/profile.d/conda.sh"
fi

if command -v conda >/dev/null 2>&1; then
  conda activate "$ENV_PREFIX"
else
  source "$ENV_PREFIX/bin/activate"
fi

cd "$BASE_DIR"

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

CORES=${SLURM_CPUS_PER_TASK:-24}

B_BOOT=${B_BOOT:-1000}
N_MC=${N_MC:-5000}
BOOT_CHUNK_SIZE=${BOOT_CHUNK_SIZE:-100}
MC_CHUNK_SIZE=${MC_CHUNK_SIZE:-100}

if [ ! -f "$ANNOTATION_FILE" ]; then
  echo "ERROR: missing annotation file: $ANNOTATION_FILE" >&2
  exit 1
fi

if [ ! -f "$BASE_DIR/CODE/run_buffa2010_hpc.R" ]; then
  echo "ERROR: missing R script: $BASE_DIR/CODE/run_buffa2010_hpc.R" >&2
  exit 1
fi

Rscript "$BASE_DIR/CODE/run_buffa2010_hpc.R" \
  --base_dir "$BASE_DIR" \
  --annotation_file "$ANNOTATION_FILE" \
  --b_boot "$B_BOOT" \
  --n_mc "$N_MC" \
  --cores "$CORES" \
  --alpha 0.05 \
  --mc_p_cutoff 0.01 \
  --product_cutoff 0.5 \
  --boot_chunk_size "$BOOT_CHUNK_SIZE" \
  --mc_chunk_size "$MC_CHUNK_SIZE" \
  --resume FALSE \
  --test_mode FALSE