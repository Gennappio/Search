# Validation

The test suite uses no comparison with the R originals as an acceptance criterion: each layer
must hold by construction or by common sense, and the tests need no R. The comparison with the
originals is a separate harness, `validation/r_parity`, run by hand; what it found is in
[Porting from R](porting.md).

| layer | where | what |
|---|---|---|
| hand-computed toy | `tests/unit/test_core_toy.py` | 3 seeds, 3 genes: D, C, O, S, F, degrees written by hand; connectivities equal in exact arithmetic stay equal, bit for bit, whatever the order of the seeds; the degree of a seed with an association that is not symmetric |
| independent references | `tests/unit/test_association.py`, `test_stats.py` | Spearman and Pearson equal to scipy pair by pair; p-values equal to the manual t formula; p.adjust known values |
| properties | `tests/property/` | gene, sample and seed order do not matter; ranges and symmetries; self edges; disabled growth equals one pass |
| mechanics | `tests/unit/test_bootstrap_null.py`, `test_resampling_checkpoint.py` | replicates reproducible from derived seeds; chunks match direct replicates; resume reuses valid chunks only; pooled and per-gene counts by hand; a gene without an observed value gets no p-value |
| synthetic structure | `tests/synthetic/` | a planted module is recovered from five of its genes (AUC > 0.95); an unrelated module does not overlap |
| calibration | `tests/statistical/` | on pure noise the null p-values show no excess of small values |
| pipeline | `tests/unit/test_pipeline.py` | end to end on three synthetic cohorts, through the API and the CLI, with resume |
| real data | `tests/realdata/` | GBM cohorts through `SEARCH_GBM_DATA`: the published 51-gene metagene is enriched at the top |
| benchmark | `tests/bench/` | one pass on 16548 x 95 under 1 s |
| R originals | `validation/r_parity` (not part of the tests) | every number of the GBM script and of the SEARCH package against the two presets, on the same inputs and the same random draws |
| verbatim port | `tests/unit/test_r2014_units.py`, `test_r2014_golden.py`, `test_r2014_defects.py` | `search.r2014` against answers recorded from the unmodified R package: every function of the package and every rule of R the port wrote down, called on its own (1,636 calls, R 4.4.3); whole runs in R 4.4.3 and in R 3.1.0 (445 cases) and a sample of every building block of R it reimplements; bit for bit; each defect of the package on an example worked out by hand. [Tests of the verbatim port](r2014_tests.md) describes them |

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
| `tests/synthetic/test_growth.py` | five seeds of a planted module of 30 recruit at least 25 module genes and at most 3 background genes; a foreign seed is pruned by D at iteration 1, or by Omax at iteration 2 when seeds are not protected; one iteration without changes equals the single pass; the sd prefilter is keyed by gene name; bootstrap over growth and an external null from random-seed runs; the internal null counts the random seed sets that lose every seed (19 of 20 here) and still separates the module (p < 0.01) from the background (p > 0.5) |
| `tests/unit/test_breast_search2014.py` | the `search2014` preset grows a network on Breastdata1 and records its stop reason |
| `tests/unit/test_pipeline.py` | the pipeline with growth writes the growth tables and `D_final` |

## Null model of `search2014`, 2026-10-12

The preset now carries the null model of `buffa2010`. Run as it is:

| data | genes with D > 0 | p < 0.05 | p < 0.01 | null stage |
|---|---|---|---|---|
| Breastdata1, the six seeds of `examples/breast` | 15 of 1129 | 15 | 15 | 100 s on one process |
| GBM IT cohort (21254 genes x 95 samples), the ten seeds of Buffa 2010 | 2132 of 21253 | 1355 (D >= 0.2) | 921 (D >= 0.222) | 7 min 29 s on eight processes |

With the pooled null the p-value is a function of D alone, so the cutoffs are thresholds on D.
On Breastdata1 only 0.01% of the probeset pairs pass the correlation threshold, random seed sets
connect nothing, and every connected probeset is significant: that dataset shows that the stage
runs, not that it discriminates. With 500 replicates instead of 5000 the breast example selects
the same 15 probesets in 11 seconds.

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
| `tests/unit/test_linear.py` | the vectorised coefficients and p-values equal least squares pair by pair with numeric and categorical covariates; without covariates they equal `scipy.stats.linregress`; they equal statsmodels OLS and gaussian GLM; a planted batch effect creates spurious edges without the covariate and none with it, while the real module is still recovered (AUC > 0.95); resampling keeps covariates aligned with the columns; the GLM path reads covariates whatever their names, turns a model that cannot be built and a run of failed fits into errors, and, like least squares, gives no association to a gene with a missing value or without variance |

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

## Comparison with the R originals, 2026-11-15 and 2026-11-17

Both originals were run next to Search on the same inputs
([Porting from R](porting.md) has the method, the tables and every difference). The first
run found six defects in Search, corrected on 2026-11-17; the numbers are from the corrected
code:

| comparison | compared quantities | equal | different | R stops |
|---|---|---|---|---|
| GBM script, 20 bootstrap and 50 Monte Carlo replicates per cohort, same draws | 85 | 85 | 0 | 0 |
| GBM script, 100 and 500 replicates, same draws | 88 | 76 | 12 | 0 |
| functions of the GBM script on Breastdata1, same draws | 27 | 27 | 0 | 0 |
| SEARCH package as it is, 28 cases | 317 | 268 | 44 | 5 |
| SEARCH package with its defects corrected, 11 cases | 156 | 155 | 1 | 0 |

The 12 differences of the second row come from one bootstrap replicate of 200, in which two
seed clusters have the same size: the script ranks apart genes whose connectivities are equal
in exact arithmetic, and Search gives them the rank of the exact values. The 26 selected genes
are the same, in the same order. With independent random draws at that size every gene falls
on the same side of the p < 0.01 cutoff in both cohorts and the same 26 genes are selected.
The 44 differences from the package as it is are 34 effects of its defects and 10 deliberate
differences; none comes from a defect of Search.

After the corrections the full run on the GBM cohorts (1000 bootstrap, 5000 Monte Carlo) was
repeated: the 26 selected genes are the same, in the same order; connectivity, observed ranks,
Monte Carlo p-values and scores are unchanged; the bootstrap means move by at most 1.5e-6.
Whole run 4 min 54 s on eight processes.

## Verbatim port of the SEARCH package, 2026-11-19

`search.r2014` is checked against R, not against common sense: for it the R package is the
specification. [Replicating the 2014 package](replication.md) has the method and every table,
[Tests of the verbatim port](r2014_tests.md) the tests that keep it so;
in short, against the unmodified package run by R 4.4.3 on the same machine:

| family | cases | equal | differ |
|---|---|---|---|
| correlation path: `SEARCH()`, `SEARCH.randomise()` and `summarise.iterations()` on Breastdata1, on the GBM IT cohort and on random matrices | 22,535 | 22,535 (37.9 million numbers bit for bit, files character for character, 5,907 cases in which R stops with the same error) | 0 |
| regression functions, to nine significant digits | 2,198 | 2,192 | 6 (standard errors of badly conditioned models, equal to six digits) |
| the building blocks of R (random draws, reading and writing numbers, `round`, `sd`, `cor`, `p.adjust`, `rowMeans`) | 389,832 values | all | 0 |

Against R 3.1.0 on Linux x86-64, with the port set to that platform: 26 of 26 cases on
Breastdata1 and all 389,671 values of the building blocks. The random cases were not run in
that R.
