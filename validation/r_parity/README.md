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

A second set of scripts compares the verbatim port, `search.r2014`, with the same package
(findings in [docs/replication.md](../../docs/replication.md)):

| comparison | R side | Python side | data |
|---|---|---|---|
| SEARCH 0.0.2 package, unmodified, against `search.r2014`: everything returned and written, bit for bit | `r/r2014_oracle.R` | `py/compare_r2014.py`, `py/r2014_lib.py`, cases from `py/r2014_cases.py` | 25,000 random matrices, Breastdata1, GBM IT cohort |
| the functions of R the port reimplements (`sample`, reading and writing numbers, `round`, `sd`, `cor`, `p.adjust`, `rowMeans`) | `r/r2014_primitives.R` | `py/check_r2014_primitives.py` | generated |
| two installations of R with each other | the two above, run twice | `py/compare_r_versions.py`, `py/compare_platforms.py` | the same cases |

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
  `nodes.fastLm` cases (skipped without it). The scripts of the verbatim port (`r2014_*`) need
  no R package at all, so that they also run in R 3.1.0. To keep them out of the main library:
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

The verbatim port:

```bash
export R_LIBS=<dir>                         # optional, for the nodes.fastLm cases
validation/r_parity/run_r2014.sh            # 20000 random SEARCH() cases and the rest, 15 minutes
validation/r_parity/run_r2014.sh 2000 4     # fewer cases, four R processes
validation/r_parity/run_r2014_r31.sh 2      # the same cases in R 3.1.0 (2014), in containers
Rscript validation/r_parity/r/r2014_primitives.R <dir>
.venv/bin/python validation/r_parity/py/check_r2014_primitives.py <dir>/primitives.json
.venv/bin/python validation/r_parity/py/make_r2014_golden.py validation/r_parity/out/r2014 tests/data/r2014
```

`run_r2014_r31.sh` needs Docker and the image `rocker/r-ver:3.1.0`, which exists for x86-64
only: elsewhere it is emulated, about ten times slower, and each container holds a whole R.
Give it a number of parallel containers the machine can spare; with seven of them next to
other containers, on a virtual machine of 4 GB, the container engine restarted. Do not edit
a script of this folder while it runs: bash and Rscript read their files as they go.
`make_r2014_golden.py` refreshes the answers the test suite compares the port with.

## Output

- `out/` (not in git): everything the two sides write. For the GBM cohorts it holds tables
  derived from patient data and must stay out of the repository.
- `results/report_*.csv` (in git): one row per compared quantity, with the number of values
  compared, how many differ, the largest absolute difference and a status.
- `results/summary.md`: the reports condensed to one row per case.
- `results/report_r2014.md` and `report_r2014.csv` (R of the machine), `report_r2014_r31.*`
  (R 3.1.0): the verbatim port, one row per family of cases; `report_r2014_cases.csv` lists
  the named cases and every random case that is not simply equal; `report_r2014_primitives*.csv`
  the building blocks.

Status values: `equal` (every value within 1e-12, or within the stated tolerance for
quantities R rounds or fits iteratively), `differs`, `r_error` (the R original stopped with an
error on that case), `info` (a description, not a check).

## Adding a case

Cases of the package comparison are the `CASES` list in `r/search2014_run.R`; each is a call
of `SEARCH()` or `SEARCH.randomise()` with its arguments. The Python side reads the case from
the `case.json` the R side writes and translates the arguments into overrides of the
`search2014` preset (`config_for` in `py/compare_search2014.py`), so a new case needs no
Python change unless it uses an argument not translated yet.
