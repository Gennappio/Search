#!/usr/bin/env bash
# Verbatim port (search.r2014) against the SEARCH 0.0.2 package, unmodified, run by R.
#
#   validation/r_parity/run_r2014.sh [random cases of SEARCH()] [parallel R processes]
#
# Writes the cases, runs every case file through r/r2014_oracle.R, runs the port on the same
# cases and leaves results/report_r2014.csv. Needs R (no package; RcppArmadillo for the
# nodes.fastLm cases, which are skipped without it) and the project environment. The GBM cases
# run when SEARCH_GBM_DATA points to the GBM project folder.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
out="$here/out/r2014"
python="$repo/.venv/bin/python"
n_fuzz=${1:-20000}
jobs=${2:-8}
gbm=""
if [ -n "${SEARCH_GBM_DATA:-}" ]; then gbm="--gbm"; fi

"$python" "$here/py/compare_r2014.py" cases "$out" --fuzz "$n_fuzz" --chunks "$jobs" $gbm
rm -rf "$out/oracle"
mkdir -p "$out/oracle" "$here/results"
ls "$out/cases" | grep '\.R$' | sed 's/\.R$//' | xargs -P "$jobs" -I{} sh -c \
  'Rscript "$1/r/r2014_oracle.R" "$2" "$3/cases/$0.R" "$3/oracle/$0" > "$3/oracle/$0.log" 2>&1 || echo "R failed on $0, see $3/oracle/$0.log"' \
  {} "$here" "$repo" "$out"
rm -rf "$out"/oracle/*/work
"$python" "$here/py/compare_r2014.py" check "$out" --report "$here/results/report_r2014.csv"
