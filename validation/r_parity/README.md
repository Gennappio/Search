# R parity harness

Runs the two R originals and `search` on the same inputs and compares every number they
produce. The findings are written up in [docs/porting.md](../../docs/porting.md); this file
says how to run the comparison again.

Nothing here is part of the test suite: the tests need no R. The harness is run by hand when
the engine changes in a way that could move results, and before a release.

## What it runs

| comparison | R side | Python side | data |
|---|---|---|---|
| SEARCH 0.0.2 package against the `search2014` preset | `r/search2014_run.R` | `py/compare_search2014.py` | Breastdata1 (in the repository), GBM IT cohort |
| functions of the GBM script against the `buffa2010` preset | `r/gbm_functions_public.R` | `py/compare_buffa2010_public.py` | Breastdata1 |
| GBM script, unmodified, against the `buffa2010` preset | `run_gbm_script.sh`, `r/gbm_script_export.R` | `py/compare_buffa2010.py` | GBM IT, FR and TCGA cohorts |

The R sources are the ones under `reference/r_sources`, read from there and never copied.
The package is not installed: its files are sourced into an environment, with
`parallel::mclapply` in place of the removed `multicore` and the plotting function stubbed
out. The GBM script is run as a whole on a staging folder that mirrors the layout it expects.

Two ideas make the comparison exact instead of approximate:

- **Patched runs.** Where an original has a defect, it is run twice: as it is, and with a
  one-line correction (`PATCHES` in `r/search2014_run.R`). The port must reproduce the
  patched run exactly; the unpatched run shows what the defect does.
- **Same draws.** R and NumPy have different random number generators. The draws of the R run
  (bootstrap sample indices, random seed sets) are exported and handed to the pipeline of
  `search` in place of its own, so every replicate can be compared number by number.

## Requirements

- R 4.x with `data.table`, `optparse`, `matrixStats`, `jsonlite`, and `RcppArmadillo` for the
  `nodes.fastLm` cases (skipped without it). To keep them out of the main library:
  `Rscript -e 'install.packages(c("data.table", "optparse", "RcppArmadillo"), lib = "<dir>")'`
  and `export R_LIBS=<dir>`.
- The project environment with the `glm` extra (`uv sync --all-extras`).
- For the GBM comparisons, `SEARCH_GBM_DATA` pointing to the GBM project folder (the one
  holding `DATA/NORMALIZED`). Without it only the public-data comparisons run.

## Running

```bash
export SEARCH_GBM_DATA=/path/to/GBM_HYPOXIA_SIGNATURE     # optional
validation/r_parity/run_all.sh            # 20 bootstrap and 50 Monte Carlo replicates per cohort
validation/r_parity/run_all.sh 100 500 8  # replicates and cores for the GBM script
```

About ten minutes with the default sizes: five for the growth loop on the IT cohort in each
language, three for the GBM script. The GBM script takes about one second per replicate and
cohort on eight cores, so the sizes of the preset (1000 and 5000) take about four hours.

Single steps:

```bash
Rscript validation/r_parity/r/search2014_run.R . validation/r_parity/out/search2014 [case ...]
.venv/bin/python validation/r_parity/py/compare_search2014.py . validation/r_parity/out/search2014 [case ...]
```

## Output

- `out/` (not in git): everything the two sides write. For the GBM cohorts it holds tables
  derived from patient data and must stay out of the repository.
- `results/report_*.csv` (in git): one row per compared quantity, with the number of values
  compared, how many differ, the largest absolute difference and a status.
- `results/summary.md`: the reports condensed to one row per case.

Status values: `equal` (every value within 1e-12, or within the stated tolerance for
quantities R rounds or fits iteratively), `differs`, `r_error` (the R original stopped with an
error on that case), `info` (a description, not a check).

## Adding a case

Cases of the package comparison are the `CASES` list in `r/search2014_run.R`; each is a call
of `SEARCH()` or `SEARCH.randomise()` with its arguments. The Python side reads the case from
the `case.json` the R side writes and translates the arguments into overrides of the
`search2014` preset (`config_for` in `py/compare_search2014.py`), so a new case needs no
Python change unless it uses an argument not translated yet.
