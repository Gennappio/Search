#!/usr/bin/env bash
# The same cases in the R the package was written for: R 3.1.0 (released 2014-04-10) for
# x86-64, from the image rocker/r-ver:3.1.0, emulated where the host is not x86.
#
#   validation/r_parity/run_r2014_r31.sh [parallel containers]
#
# Run run_r2014.sh first: it writes the cases. Under emulation the regression cases on
# Breastdata1 and the GBM cohort take hours and are left out by default (R2014_SKIP, a regular
# expression on the names of the case files); the nodes.fastLm cases are skipped anyway (the
# image has no RcppArmadillo). Leaves results/report_r2014_r31.csv. What this run is for: R 3.1.0
# draws its samples in another way, rounds in another way and, on x86, adds up in 80-bit
# extended precision; the report shows where that changes a result of the package.
set -euo pipefail

here=$(cd "$(dirname "$0")" && pwd)
repo=$(cd "$here/../.." && pwd)
out="$here/out/r2014"
python="$repo/.venv/bin/python"
jobs=${1:-6}
image=${R2014_IMAGE:-rocker/r-ver:3.1.0}
if [ ! -d "$out/cases" ]; then echo "no cases under $out: run run_r2014.sh first"; exit 1; fi
gbm_mount=""
if [ -n "${SEARCH_GBM_DATA:-}" ]; then gbm_mount="-v $SEARCH_GBM_DATA:/gbm:ro -e SEARCH_GBM_DATA=/gbm"; fi

rm -rf "$out/oracle_r31"
mkdir -p "$out/oracle_r31" "$here/results"
skip=${R2014_SKIP:-'named_(breast_(glm|fastlm|manual_glm|unknown_family)|it_)'}
ls "$out/cases" | grep '\.R$' | sed 's/\.R$//' | grep -Ev "$skip" | xargs -P "$jobs" -I{} sh -c \
  'docker run --rm --platform linux/amd64 -v "$1":/repo:ro -v "$2":/out $4 "$3" Rscript /repo/validation/r_parity/r/r2014_oracle.R /repo /out/cases/$0.R /out/oracle_r31/$0 > "$2/oracle_r31/$0.log" 2>&1 || echo "R failed on $0, see $2/oracle_r31/$0.log"' \
  {} "$repo" "$out" "$image" "$gbm_mount"
rm -rf "$out"/oracle_r31/*/work
"$python" "$here/py/compare_r2014.py" check "$out" --oracle oracle_r31 --report "$here/results/report_r2014_r31.csv"
