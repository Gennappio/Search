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
