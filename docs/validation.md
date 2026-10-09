# Validation

No comparison with the R originals is used as an acceptance criterion. Each layer must hold by
construction or by common sense.

| layer | where | what |
|---|---|---|
| hand-computed toy | `tests/unit/test_core_toy.py` | 3 seeds, 3 genes: D, C, O, S, F, degrees written by hand |
| independent references | `tests/unit/test_association.py`, `test_stats.py` | Spearman and Pearson equal to scipy pair by pair; p-values equal to the manual t formula; p.adjust known values |
| properties | `tests/property/` | gene and sample order do not matter; ranges and symmetries; self edges; disabled growth equals one pass |
| mechanics | `tests/unit/test_bootstrap_null.py`, `test_resampling_checkpoint.py` | replicates reproducible from derived seeds; chunks match direct replicates; resume reuses valid chunks only; pooled and per-gene counts by hand |
| synthetic structure | `tests/synthetic/` | a planted module is recovered from five of its genes (AUC > 0.95); an unrelated module does not overlap |
| calibration | `tests/statistical/` | on pure noise the null p-values show no excess of small values |
| pipeline | `tests/unit/test_pipeline.py` | end to end on three synthetic cohorts, through the API and the CLI, with resume |
| real data | `tests/realdata/` | GBM cohorts through `SEARCH_GBM_DATA`: the published 51-gene metagene is enriched at the top |
| benchmark | `tests/bench/` | one pass on 16548 x 95 under 1 s |

## Real data, phase 2

With the ten seeds of Buffa 2010 and the `buffa2010` preset:

| cohort | metagene genes in the top 200 | hypergeometric p |
|---|---|---|
| IT (95 samples) | 8 of 40 | 3e-9 |
| FR (279 samples) | 10 of 36 | 8e-12 |

A smoke run (20 bootstrap, 50 Monte Carlo replicates per cohort) on IT + FR takes 20 seconds on a
laptop and selects 26 candidate genes, the same 26 genes the R script selected in June 2026 with
1000 and 5000 replicates. This is reported as information, not as a benchmark: the two runs used
slightly different versions of the input matrices.

## Growth loop, phase 3

| test | what it checks |
|---|---|
| `tests/synthetic/test_growth.py` | five seeds of a planted module of 30 recruit at least 25 module genes and at most 3 background genes; a foreign seed is pruned by D at iteration 1, or by Omax at iteration 2 when seeds are not protected; one iteration without changes equals the single pass; the sd prefilter is keyed by gene name; bootstrap over growth and an external null from random-seed runs |
| `tests/unit/test_breast_search2014.py` | the `search2014` preset grows a network on Breastdata1 and records its stop reason |
| `tests/unit/test_pipeline.py` | the pipeline with growth writes the growth tables and `D_final` |

## Full run on the GBM cohorts, phase 2

1000 bootstrap and 5000 Monte Carlo replicates per cohort, IT + FR training, TCGA scoring,
8 worker processes on a laptop:

| stage | IT | FR |
|---|---|---|
| core pass | 0.1 s | 0.2 s |
| bootstrap | 16 s | 39 s |
| null model | 59 s | 163 s |

Whole run 4 min 39 s (the R script took 1 h 50 min on 32 cluster cores), checkpoints 113 MB
(1.5 GB in R). The 26 selected genes are the same 26 of the June 2026 R run; the bootstrap means
of the two runs correlate at 0.999 and the Monte Carlo p-values at 0.995, with the two runs
using different random draws and slightly different input matrices.

## Linear models, phase 4

| test | what it checks |
|---|---|
| `tests/unit/test_linear.py` | the vectorised coefficients and p-values equal least squares pair by pair with numeric and categorical covariates; without covariates they equal `scipy.stats.linregress`; they equal statsmodels OLS and gaussian GLM; a planted batch effect creates spurious edges without the covariate and none with it, while the real module is still recovered (AUC > 0.95); resampling keeps covariates aligned with the columns |

## Usage audit and variant comparison, 2026-10-09

After the audit every command was exercised on the GBM cohorts through the command line:
`qc`, `run`, `resume` after a deleted chunk, `combine` with other cutoffs, `score`, `compare`,
`plot`, `summarise`, a strict run (refused on the sample names of the IT cohort), a missing seed
file (one-line error), the `search2014` preset on Breastdata1 through `examples/breast`, and the
conda environment of `environment.yml`, which reproduces the uv results exactly.

The three open statistical questions were then run with the full replicate numbers of the
`buffa2010` preset (1000 bootstrap, 5000 Monte Carlo per cohort, IT + FR), six pipelines of
about five minutes each on eight processes:

| membership | null pooling | candidates | Spearman vs script | top-50 overlap | signature Jaccard |
|---|---|---|---|---|---|
| cluster_size (script) | pooled | 26 | 1.000 | 1.00 | 1.00 |
| cluster_size | per gene | 28 | 1.000 | 1.00 | 0.54 |
| binary (SEARCH 2014) | pooled | 216 | 0.998 | 0.88 | 0.10 |
| binary | per gene | 216 | 0.998 | 0.88 | 0.12 |
| gene_normalised (paper as printed) | pooled | 0 | 0.999 | 0.60 | 0.00 |
| gene_normalised | per gene | 0 | 0.999 | 0.60 | 0.00 |

Reading the table: the ranking of the genes barely depends on the membership weighting
(Spearman above 0.998 in every variant), the **selection** does. With the script's weighting the
pooled and the per-gene null share 19 genes; the pooled null adds AKAP12, ENO2, GAPDH, KISS1R,
MAFF, MCHR1 and ZNF395, the per-gene null adds APLN, FAM180A, IGFBP5, KCNE4, OXSR1, PCOLCE2, PI3,
SLC6A8 and ZBTB21. The binary weighting makes the product criterion pass for hundreds of genes
because D takes only eleven distinct values. The paper's formula as printed selects nothing. With
a tenth of the replicates the per-gene variant had given 19 candidates instead of 28: that null
needs the full replicate count to settle. The tables of every variant are kept with the project
results (`RESULTS/SEARCH_PY/variants_full/`), outside the repository.
