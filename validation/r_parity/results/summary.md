# R parity reports, summary

## report_buffa2010_100_500.csv

92 compared quantities: 80 equal, 12 different, 0 where R stopped with an error.

| case | R run | quantities | equal | differ | R error | values compared | largest difference among the equal |
|---|---|---|---|---|---|---|---|
| IT | same draws | 27 | 22 | 5 | 0 | 10602699 | 3.3e-14 |
| FR | same draws | 24 | 24 | 0 | 0 | 10652316 | 5.3e-15 |
| IT x FR | same draws | 29 | 22 | 7 | 0 | 330832 | 6.1e-16 |
| scores | same draws | 8 | 8 | 0 | 0 | 4704 | 1.8e-15 |
| IT | own draws | 2 | 2 | 0 | 0 | 33078 | 2.2e-16 |
| FR | own draws | 2 | 2 | 0 | 0 | 33078 | 2.2e-16 |

Information rows:

- IT, own draws, C_rank_boot_mean: Spearman 0.9650, largest absolute difference 0.1176
- IT, own draws, p_MC_C: Spearman 1.0000, largest absolute difference 0.0034; p < 0.01: R 73, Python 73, same side of the cutoff 16539 of 16539
- FR, own draws, C_rank_boot_mean: Spearman 0.9549, largest absolute difference 0.1197
- FR, own draws, p_MC_C: Spearman 1.0000, largest absolute difference 0.0011; p < 0.01: R 202, Python 202, same side of the cutoff 16539 of 16539
- IT x FR, own draws, selection: genes of the signature: R 26, Python 26, shared 26
- timing, same draws, seconds per stage: {'IT': {'core': 0.1, 'bootstrap': 9.1, 'null': 33.4}, 'FR': {'core': 0.2, 'bootstrap': 23.6, 'null': 94.4}}
- timing, own draws, seconds per stage: {'IT': {'core': 0.1, 'bootstrap': 4.4, 'null': 15.0}, 'FR': {'core': 0.2, 'bootstrap': 7.1, 'null': 30.0}}
- timing, R script, whole run: R script: 1249 s for 100 bootstrap and 500 Monte Carlo replicates per cohort

## report_buffa2010_20_50.csv

89 compared quantities: 89 equal, 0 different, 0 where R stopped with an error.

| case | R run | quantities | equal | differ | R error | values compared | largest difference among the equal |
|---|---|---|---|---|---|---|---|
| IT | same draws | 24 | 24 | 0 | 0 | 1885746 | 3.3e-14 |
| FR | same draws | 24 | 24 | 0 | 0 | 1885746 | 5.3e-15 |
| IT x FR | same draws | 29 | 29 | 0 | 0 | 446605 | 6.1e-16 |
| scores | same draws | 8 | 8 | 0 | 0 | 4704 | 1.8e-15 |
| IT | own draws | 2 | 2 | 0 | 0 | 33078 | 2.2e-16 |
| FR | own draws | 2 | 2 | 0 | 0 | 33078 | 2.2e-16 |

Information rows:

- IT, own draws, C_rank_boot_mean: Spearman 0.8792, largest absolute difference 0.2614
- IT, own draws, p_MC_C: Spearman 1.0000, largest absolute difference 0.0192; p < 0.01: R 72, Python 75, same side of the cutoff 16536 of 16539
- FR, own draws, C_rank_boot_mean: Spearman 0.8600, largest absolute difference 0.2710
- FR, own draws, p_MC_C: Spearman 1.0000, largest absolute difference 0.0168; p < 0.01: R 164, Python 197, same side of the cutoff 16506 of 16539
- IT x FR, own draws, selection: genes of the signature: R 26, Python 26, shared 26
- timing, same draws, seconds per stage: {'IT': {'core': 0.1, 'bootstrap': 2.0, 'null': 3.4}, 'FR': {'core': 0.2, 'bootstrap': 4.9, 'null': 9.7}}
- timing, own draws, seconds per stage: {'IT': {'core': 0.1, 'bootstrap': 1.6, 'null': 1.5}, 'FR': {'core': 0.2, 'bootstrap': 1.5, 'null': 3.0}}
- timing, R script, whole run: R script: 154 s for 20 bootstrap and 50 Monte Carlo replicates per cohort

## report_buffa2010_public.csv

42 compared quantities: 30 equal, 12 different, 0 where R stopped with an error.

| case | R run | quantities | equal | differ | R error | values compared | largest difference among the equal |
|---|---|---|---|---|---|---|---|
| Breast1 | same draws | 27 | 27 | 0 | 0 | 81358 | 3.3e-16 |
| Breast1, one missing value in 30 probesets | by design | 15 | 3 | 12 | 0 | 3693 | 0.0e+00 |

## report_search2014.csv

473 compared quantities: 423 equal, 45 different, 5 where R stopped with an error.

| case | R run | quantities | equal | differ | R error | values compared | largest difference among the equal |
|---|---|---|---|---|---|---|---|
| breast_bootstrap | verbatim | 11 | 10 | 1 | 0 | 17367 | 2.0e-06 |
| breast_bootstrap_6iter | verbatim | 15 | 13 | 2 | 0 | 53189 | 4.3e-06 |
| breast_bootstrap_defaults | verbatim | 7 | 1 | 5 | 1 | 10 | 0.0e+00 |
| breast_bootstrap_defaults | patched | 10 | 9 | 1 | 0 | 15888 | 3.3e-06 |
| breast_doc_example | verbatim | 12 | 12 | 0 | 0 | 2636 | 3.3e-16 |
| breast_doc_example | patched | 12 | 12 | 0 | 0 | 2636 | 3.3e-16 |
| breast_doc_not_immune | verbatim | 1 | 0 | 0 | 1 | 0 | - |
| breast_doc_not_immune | patched | 16 | 16 | 0 | 0 | 3789 | 2.2e-16 |
| breast_fastlm | verbatim | 10 | 10 | 0 | 0 | 1315 | 7.1e-15 |
| breast_fastlm [association.method=glm] | verbatim | 10 | 10 | 0 | 0 | 1315 | 1.4e-14 |
| breast_fastlm_doc | verbatim | 14 | 7 | 7 | 0 | 94 | 5.6e-17 |
| breast_fastlm_doc | patched | 14 | 14 | 0 | 0 | 2849 | 1.2e-14 |
| breast_glm_bh | verbatim | 10 | 10 | 0 | 0 | 1413 | 1.4e-14 |
| breast_glm_bh [association.method=glm] | verbatim | 10 | 10 | 0 | 0 | 1413 | 1.1e-14 |
| breast_glm_covariates | verbatim | 10 | 10 | 0 | 0 | 1597 | 1.3e-14 |
| breast_glm_covariates [association.method=glm] | verbatim | 10 | 10 | 0 | 0 | 1597 | 3.0e-13 |
| breast_glm_covariates_ignored [covariates declared] | verbatim | 10 | 4 | 6 | 0 | 7 | 0.0e+00 |
| breast_glm_covariates_ignored [no covariates] | verbatim | 10 | 10 | 0 | 0 | 1592 | 1.4e-14 |
| breast_glm_covariates_plain | verbatim | 10 | 10 | 0 | 0 | 1597 | 1.3e-14 |
| breast_glm_covariates_plain [association.method=glm] | verbatim | 10 | 10 | 0 | 0 | 1597 | 3.0e-13 |
| breast_glm_doc | verbatim | 14 | 7 | 7 | 0 | 94 | 5.6e-17 |
| breast_glm_doc | patched | 14 | 14 | 0 | 0 | 2849 | 6.9e-15 |
| breast_glm_gamma | verbatim | 10 | 10 | 0 | 0 | 1315 | 8.1e-12 |
| breast_glm_poisson | verbatim | 1 | 0 | 0 | 1 | 0 | - |
| breast_glm_poisson_loose | verbatim | 10 | 10 | 0 | 0 | 3412 | 1.7e-11 |
| breast_growth | verbatim | 12 | 12 | 0 | 0 | 2558 | 3.3e-16 |
| breast_jackknife | verbatim | 10 | 9 | 1 | 0 | 8642 | 3.3e-06 |
| breast_missing_values | verbatim | 1 | 0 | 0 | 1 | 0 | - |
| breast_missing_values | patched | 24 | 24 | 0 | 0 | 6368 | 3.3e-16 |
| breast_omax_pruning | verbatim | 16 | 10 | 6 | 0 | 2499 | 1.1e-16 |
| breast_omax_pruning | patched | 16 | 16 | 0 | 0 | 3775 | 3.3e-16 |
| breast_pearson | verbatim | 8 | 8 | 0 | 0 | 1278 | 6.7e-16 |
| breast_pruned_members | verbatim | 24 | 19 | 5 | 0 | 6196 | 1.1e-16 |
| breast_pruned_members | patched | 24 | 24 | 0 | 0 | 6368 | 3.3e-16 |
| breast_sd_prefilter | verbatim | 8 | 4 | 4 | 0 | 16 | 5.6e-17 |
| breast_sd_prefilter | patched | 8 | 8 | 0 | 0 | 1284 | 3.3e-16 |
| breast_seeds_not_immune | verbatim | 1 | 0 | 0 | 1 | 0 | - |
| breast_seeds_not_immune | patched | 12 | 12 | 0 | 0 | 2558 | 3.3e-16 |
| breast_single_pass | verbatim | 10 | 10 | 0 | 0 | 1316 | 3.3e-16 |
| breast_universe_subset | verbatim | 8 | 8 | 0 | 0 | 652 | 3.3e-16 |
| breast_unscaled | verbatim | 8 | 8 | 0 | 0 | 1288 | 3.3e-16 |
| it_growth | verbatim | 10 | 10 | 0 | 0 | 67368 | 1.1e-16 |
| it_sd_prefilter | verbatim | 6 | 6 | 0 | 0 | 22176 | 1.1e-16 |
| it_sd_prefilter | patched | 6 | 6 | 0 | 0 | 22176 | 1.1e-16 |
