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
  threshold, binary membership, recursive network growth with D and O pruning).

Free configurations are allowed; only presets are validated, and the run manifest records which
one was used.

## Status

Version 0.1: the single-pass engine, the recursive growth loop, correlation and linear-model
association with covariates, bootstrap, Monte Carlo and external null models, cohort combination,
signature score (internal, or through pysigscore), checkpoints with resume, run manifest, variant
comparison, plots, command line and documentation site (`uv run mkdocs serve`). See `docs/` for the method, the presets, the decision log and the validation.

## Install

```bash
uv sync --all-extras          # development environment (.venv in the project, or UV_PROJECT_ENVIRONMENT)
uv run search --help
uv run pytest                 # tests; SEARCH_GBM_DATA=<folder> enables the real-data tests
```

pysigscore, the optional scorer backend, is not on PyPI: install it by hand into the same
environment (`uv pip install --python .venv/bin/python <path or git url>`); without it the
internal score of Buffa 2010 is used.

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
attribution and consultation only. Search is released under the GPL-2.0-or-later licence.
