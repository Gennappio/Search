#!/usr/bin/env bash
# Parity harness, R side: the lab's GBM script (reference/r_sources/gbm), run unmodified.
#
# The script hard-codes the layout of its base directory, so a staging folder mirrors it:
#   stage/DATA/NORMALIZED                              -> the real matrices (symbolic link)
#   stage/DATA/gencode.v36.annotation.gtf.gene.probemap   stand-in built from the lab annotation
#                                                         (the script uses it for symbols only)
# The script writes its results under stage/RESULTS/BUFFA2010_EXACT_HPC/run_<timestamp>; a second
# R script then exports them, with the random draws, as text tables for the Python comparison.
#
# Usage: SEARCH_GBM_DATA=<GBM project folder> run_gbm_script.sh <output dir> [b_boot] [n_mc] [cores]
set -euo pipefail

repo=$(cd "$(dirname "$0")/../.." && pwd)
out=${1:?usage: run_gbm_script.sh <output dir> [b_boot] [n_mc] [cores]}
b_boot=${2:-20}
n_mc=${3:-50}
cores=${4:-8}
: "${SEARCH_GBM_DATA:?set SEARCH_GBM_DATA to the GBM project folder}"

mkdir -p "$out"
out=$(cd "$out" && pwd)
stage="$out/stage"
rm -rf "$stage/RESULTS"
mkdir -p "$stage/DATA"
ln -sfn "$SEARCH_GBM_DATA/DATA/NORMALIZED" "$stage/DATA/NORMALIZED"

Rscript -e 'suppressPackageStartupMessages(library(data.table)); a <- commandArgs(TRUE); ann <- fread(a[1]); fwrite(ann[, .(id = ensg, gene = symbol_consensus)], a[2], sep = "\t")' \
  "$SEARCH_GBM_DATA/DATA/master_annotation_ITFR_Ensembl104_TCGA_GDC_consensus.csv" \
  "$stage/DATA/gencode.v36.annotation.gtf.gene.probemap"

start=$(date +%s)
Rscript "$repo/reference/r_sources/gbm/run_buffa2010_hpc_wo_consensus.R" \
  --base_dir "$stage" --b_boot "$b_boot" --n_mc "$n_mc" \
  --boot_chunk_size 10 --mc_chunk_size 10 --cores "$cores" --resume FALSE \
  > "$out/r_script_stdout.log" 2>&1
echo "R script: $(( $(date +%s) - start )) s for $b_boot bootstrap and $n_mc Monte Carlo replicates per cohort" | tee "$out/r_script_seconds.txt"

Rscript "$repo/validation/r_parity/r/gbm_script_export.R" "$repo" "$stage" "$out/r_export"
