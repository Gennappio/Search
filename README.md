# Search

Python engine for seed-based gene co-expression clustering: **S**eed **A**gglomerative and
**R**ecursive **C**lustering with **H**ypothesis-oriented initialisation (SEARCH).

Starting from a set of *seed* genes known to drive a biological process, Search finds the genes
consistently co-expressed with the seeds across patients, scores every gene by its connectivity to
the seed network, estimates the stability of that score by bootstrap, its significance by Monte
Carlo with random seed sets, combines several cohorts, and derives a compact gene signature plus a
per-sample score.

One engine, two named presets:

- `buffa2010` reproduces the method of Buffa et al., *Br J Cancer* 2010 (single pass, Bonferroni
  edges, cluster-size membership, bootstrap ranks, pooled Monte Carlo null, cross-cohort rank
  product, hypoxia score).
- `search2014` reproduces the behaviour of the original SEARCH R package (hard correlation
  threshold, binary membership, recursive network growth with D and O pruning), with the null
  model of `buffa2010`.

Free configurations are allowed; only presets are validated, and the run manifest records which
one was used.

## Status

Version 0.1.1: the single-pass engine, the recursive growth loop, correlation and linear-model
association with covariates, bootstrap, Monte Carlo and external null models, cohort combination,
signature score, checkpoints with resume, run manifest, variant
comparison, plots, command line and documentation site (`uv run mkdocs serve`). See `docs/` for the method, the presets, the decision log and the validation.

The port was compared with both R originals number by number; the result, and every difference
with its cause, is in `docs/porting.md`. The comparison can be repeated with
`validation/r_parity/run_all.sh` (it needs R).

Next to the engine, `search.r2014` is a second, verbatim port of the SEARCH R package: every
function transliterated with its defects, to obtain again, without R, what the 2014 package
returns. It reproduces the unmodified package bit for bit on its correlation path, in R 4.4.3
and in the R 3.1.0 of 2014; `docs/replication.md` has the study and
`validation/r_parity/run_r2014.sh` repeats it.

## Install

```bash
uv sync --all-extras          # development environment (.venv in the project, or UV_PROJECT_ENVIRONMENT)
uv run search --help
uv run pytest                 # tests; SEARCH_GBM_DATA=<folder> enables the real-data tests
```

On a SLURM cluster without `uv`, `environment.yml` creates a conda environment and installs the
package with pip; `slurm/run_search.sbatch` is a job template.

## Inputs

Already-normalised expression matrices (genes x samples, CSV/TSV/parquet, optional h5ad), one
per cohort, listed in a YAML cohort manifest together with a seed list and an optional sample
table per cohort. Normalisation (DESeq2 VST, ComBat) is out of scope and must happen upstream.

## Origin and licence

Search is a Python port and generalisation of the SEARCH R package (Syed Haider and
Francesca M. Buffa, University of Oxford, 2014, GPL-2) and of a 2026 R script implementing
Buffa et al. 2010 on RNA-seq cohorts. The original R sources are kept under `reference/` for
attribution, for consultation and for the comparison harness under `validation/r_parity`.
Search is released under the GPL-2.0-or-later licence.
