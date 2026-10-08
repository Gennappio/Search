# Presets

A preset is a named `RunConfig`. `search.presets.matching_preset` tells whether a configuration
reproduces one, ignoring runtime fields (number of replicates, chunk size, random seeds,
parallelism, checkpoint precision). The run manifest records the answer as `preset_validated`.

## buffa2010

Buffa et al. 2010 as implemented by the GBM script.

| stage | values |
|---|---|
| association | Spearman, no covariates |
| edges | `p_threshold`, alpha 0.05 / number of genes, self edge |
| membership | `cluster_size` |
| connectivity | self excluded, fraction, fractional rank |
| growth | off |
| resampling | bootstrap, 1000 replicates, chunks of 100, seeds from 1 000 000, statistic `C_rank` |
| null | internal random seed sets, 5000 replicates, pooled, statistic `C_raw` |
| combine | p < 0.01 in every cohort, rank product > 0.5, positive consistency, seeds excluded, top 51 |
| score | median of |expression| of the signature genes, fractional rank |

## search2014

The SEARCH R package with the defaults of `SEARCH.randomise`.

| stage | values |
|---|---|
| association | Spearman, no covariates |
| edges | `abs_threshold` 0.5, self edge |
| membership | `binary` |
| connectivity | self excluded, fraction, no rank (this is D) |
| growth | on: max_iter 1, D_lower 0.1, D_higher 0.24, O_lower 0.01, seeds kept, sd threshold 0.1 |
| resampling | bootstrap, 10 replicates, statistic `C_raw` |
| null | external runs, per gene |
| combine, score | off |
