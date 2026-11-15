#!/usr/bin/env bash
# Runs the whole R parity harness and leaves the three reports under results/.
#
#   validation/r_parity/run_all.sh [b_boot] [n_mc] [cores]
#
# Needs R (data.table, optparse, matrixStats, jsonlite; RcppArmadillo for the nodes.fastLm
# cases) and the project environment (.venv, with the glm extra). The comparisons on the GBM
# cohorts run only when SEARCH_GBM_DATA points to the GBM project folder; the others use the
# public breast matrices under tests/data.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
out="$here/out"
python="$repo/.venv/bin/python"
b_boot=${1:-20}
n_mc=${2:-50}
cores=${3:-8}
mkdir -p "$here/results"

echo "== SEARCH 0.0.2 package against the search2014 preset"
Rscript "$here/r/search2014_run.R" "$repo" "$out/search2014"
"$python" "$here/py/compare_search2014.py" "$repo" "$out/search2014" > "$out/search2014/comparison.log"
cp "$out/search2014/report_search2014.csv" "$here/results/"

echo "== functions of the GBM script against the buffa2010 preset, public data"
Rscript "$here/r/gbm_functions_public.R" "$repo" "$out/buffa2010_public"
"$python" "$here/py/compare_buffa2010_public.py" "$repo" "$out/buffa2010_public" > "$out/buffa2010_public/comparison.log"
cp "$out/buffa2010_public/report_buffa2010_public.csv" "$here/results/"

if [ -n "${SEARCH_GBM_DATA:-}" ]; then
  echo "== GBM script, unmodified, against the buffa2010 preset: $b_boot bootstrap, $n_mc Monte Carlo"
  "$here/run_gbm_script.sh" "$out/buffa2010" "$b_boot" "$n_mc" "$cores"
  "$python" "$here/py/compare_buffa2010.py" "$repo" "$out/buffa2010" "$cores" > "$out/buffa2010/comparison.log"
  cp "$out/buffa2010/report_buffa2010.csv" "$here/results/report_buffa2010_${b_boot}_${n_mc}.csv"
else
  echo "== SEARCH_GBM_DATA not set: the comparisons on the GBM cohorts are skipped"
fi

"$python" "$here/py/summarise_reports.py" "$here/results"
