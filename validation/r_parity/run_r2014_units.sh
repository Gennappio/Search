#!/usr/bin/env bash
# Unit cases of the verbatim port (search.r2014): every function of the SEARCH 0.0.2 package
# called on its own, in R and in the port.
#
#   validation/r_parity/run_r2014_units.sh [scale]
#
# Writes the cases, has r/r2014_units.R answer them with the unmodified package, runs the port
# on the same cases and leaves results/report_r2014_units.md (one row per function) and
# report_r2014_units.csv (one row per case). With scale 1 (the default) it also rewrites
# tests/data/r2014/units_r44_arm64.json.gz, the answers the test suite uses. A larger scale
# multiplies the random cases, leaves the test data alone and writes
# results/report_r2014_units_x<scale>.md only; scale 20 is 25,000 cases and two hours of R.
# Needs R with jsonlite (and RcppArmadillo for the nodes.fastLm cases) and the project
# environment. To compare the port again with answers R already gave:
#   .venv/bin/python validation/r_parity/py/r2014_units.py check validation/r_parity/out/r2014_units
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
python="$repo/.venv/bin/python"
scale=${1:-1}

mkdir -p "$here/results"
if [ "$scale" = "1" ]; then
  out="$here/out/r2014_units"
  "$python" "$here/py/r2014_units.py" run "$out" \
    --report "$here/results/report_r2014_units.csv" --summary "$here/results/report_r2014_units.md"
  "$python" "$here/py/r2014_units.py" golden "$out" "$repo/tests/data/r2014"
else
  out="$here/out/r2014_units_x$scale"
  "$python" "$here/py/r2014_units.py" run "$out" --scale "$scale" \
    --report "$out/report.csv" --summary "$here/results/report_r2014_units_x$scale.md"
fi
