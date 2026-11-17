# Porting from R

Search is a Python port of two R programs. This page records what was ported, how the port
was checked against the originals, where the results are the same, where they are not, and
why. It is written for whoever has to decide whether results obtained with Search can stand
in for results obtained with the R code.

The comparison was first run on 2026-11-15, against commit `5177011` of Search (version
0.1.1). It found defects in the originals and in the port; those of the port were corrected
on 2026-11-17, and every number on this page is from the corrected code. Versions: R 4.4.3,
data.table 1.18.2, matrixStats 1.5.0, RcppArmadillo 15.2.4, Python 3.12.9, NumPy 2.5.3,
SciPy 1.18.1, pandas 3.0.6 and statsmodels 0.15.0, on macOS (arm64).

## Summary

- **The lab's GBM script against the `buffa2010` preset: the same results.** The script was
  run unmodified on the IT, FR and TCGA matrices, and Search was given the same random
  draws. With 20 bootstrap and 50 Monte Carlo replicates per cohort all 85 compared
  quantities are equal: 4.2 million numbers, the largest difference 3.3e-14, and the same 26
  genes selected in the same order. With 100 and 500 replicates 76 of 88 quantities are
  equal; the other 12 all come from one bootstrap replicate out of 200, in which the script
  ranks apart genes whose connectivities are equal in exact arithmetic, while Search gives
  them the rank of the exact values (see [Floating point](#floating-point)). The effect on
  that run is at most 1.2e-5 on a bootstrap mean and none on the selected genes.
- **The SEARCH 0.0.2 package against the `search2014` preset: the same results wherever the
  package works as intended.** In 28 cases that exercise its three association functions,
  the growth thresholds, both resampling schemes and the summary step, 423 of 473 compared
  quantities are equal. The package stops with an error in 5 cases. Of the 45
  quantities that differ, 34 are caused by defects of the package, each of which was
  reproduced and then removed with a one-line correction of the R source, after which R and
  Search agree. The other 11 are deliberate differences.
- **Eight defects of the R package** are involved: three were already in the
  [decision log](decisions.md), one extends an entry of it, four were found by this
  comparison. Run with all its default arguments on its own test data, the package stops
  with an error at the summary step.
- **Apart from that tie effect, random numbers are the only reason a `buffa2010` run
  differs between the two languages.** R and NumPy use different generators, so the same
  seed gives different bootstrap samples. With independent draws and 100 and 500
  replicates, all 16539 genes fall on the same side of the p < 0.01 cutoff in both cohorts
  and the same 26 genes are selected.
- **Six defects were found in the port, and corrected on 2026-11-17**, each with a test.
  Four came out of the comparison: a silent failure of the GLM association with some
  covariate names, a rank that depended on the order of the seeds when two seed clusters
  have the same size (the tie effect seen from the Python side), a transposed node degree
  for non-gaussian GLMs, and a p-value given to genes that have no observed value. Two more
  came out while correcting the first: with the GLM association a constant gene got an edge
  with every seed, and with least squares it got a meaningless statistic. Only the
  correction of the rank touches `buffa2010` results: run again with 1000 and 5000
  replicates, the bootstrap means move by at most 1.5e-6 and the 26 selected genes are the
  same, in the same order. They are listed in
  [Defects found in the port and corrected](#defects-found-in-the-port-and-corrected).

## What was ported

| original | where | size | what it does | preset |
|---|---|---|---|---|
| SEARCH 0.0.2, an R package (Haider and Buffa, 2014) | `reference/r_sources/SEARCH` | 19 files, 1548 lines | grows a network around seed genes: correlation or regression between seeds and genes, connectivity D, overlap O, recursive addition and removal of seeds, bootstrap, comparison with runs on random seeds | `search2014` |
| `run_buffa2010_hpc.R`, the lab's script (2026) | `reference/r_sources/gbm` | 733 lines (794 in the variant with the annotation filter) | one pass of the method of Buffa et al. 2010 on the GBM cohorts: Spearman correlations, Bonferroni edges, membership, connectivity C, overlaps S and F, bootstrap of ranks, Monte Carlo null model, combination of two cohorts, hypoxia score | `buffa2010` |

Search has one engine for both (`search.core.corepass.run_core`); a preset is a named
configuration of it. The two variants of the script differ only in a gene filter on an
annotation table, which removes no gene from the current matrices; the variant without the
filter is the one that produced the lab's results and the one compared here.

Not ported: the plots of the package (igraph and gplots; `search plot` draws its own
figures), `get.program.defaults` (it only finds the test data inside an installed package)
and the hard-coded folder layout of both programs.

## How the comparison was made

The comparison is a set of scripts, `validation/r_parity`, that can be run again at any time
([Repeating the comparison](#repeating-the-comparison)). Four principles:

**The originals run as they are.** The package depends on the removed `multicore` and cannot
be installed, so its source files are read from `reference/r_sources` into an R session,
with `parallel::mclapply` in place of multicore's and the plotting function replaced by an
empty one. The GBM script is executed as a whole, on a folder that mirrors the layout it
expects. No line of either is edited.

**Same inputs.** Both sides read the same files: Breastdata1, the 1231 x 77 microarray
matrix bundled with the package (in the repository under `tests/data/breast`), and the
three GBM matrices (IT 21254 x 95, FR 17062 x 279, TCGA 32987 x 214; 16539 genes are common
to IT and FR).

**Same draws.** A bootstrap replicate draws samples at random, and a null replicate draws
seeds at random. R (Mersenne Twister) and NumPy (PCG64) produce different numbers from the
same seed, so two correct implementations give different replicates. To compare replicate
by replicate, the draws made by R are written to files and Search is handed those instead of
drawing its own. Nothing else in Search is replaced. A second run, with Search's own
generator, shows what the difference in random numbers alone does.

**Corrected runs.** Where the R package has a defect, it is run twice: as it is, and with a
correction that changes one line of its source (the corrections are the `PATCHES` list in
`validation/r_parity/r/search2014_run.R`). If Search differs from the package as it is and
agrees with the corrected package, the difference is explained by that defect and by
nothing else.

Two numbers are called equal when they differ by less than 1e-12. Both languages compute in
double precision (about 16 significant digits) with different libraries, so the last digits
can differ; the largest difference seen among numbers reported as equal is 3.3e-14. Two
exceptions with a wider tolerance are stated where they apply: the R package rounds D to
five decimals in its files, and GLM families other than the gaussian are fitted by an
iterative algorithm that R and statsmodels stop at slightly different points.

## The GBM script and the `buffa2010` preset

### Code

| in the script | in Search |
|---|---|
| `read_expr_csv` | `search.data.read_matrix` |
| restriction to the genes common to IT and FR | `gene_filter.common_genes_across_training` in the manifest |
| `buffa2010_core_one_dataset`: correlations and p-values | `search.core.association.spearman_block`, `search.stats.t_approx_pvalue` (`spearman_p_approx`) |
| edges `delta` (p < alpha / number of genes) | `search.core.edges.compute_edges`, rule `p_threshold`, `bonferroni = n_genes` |
| membership `gamma` (edge divided by the size of the seed's cluster) | `search.core.membership.compute_membership`, mode `cluster_size` |
| connectivity `C_raw`, fractional rank `C_rank` (`rank01`) | `search.core.connectivity.compute_connectivity`, `search.stats.rank01` |
| edge counts, `positive_consistent`, `weighted_seed_rho` | `search.core.connectivity.edge_statistics` |
| overlaps `S` and `F` between seeds | `search.core.overlap.shared_neighbourhood` |
| `bootstrap_checkpointed` | `search.bootstrap.run_bootstrap` |
| `monte_carlo_checkpointed` (random seed sets, pooled p-values) | `search.null.run_internal_null`, `pooled_counts`, `empirical_p` |
| `chunk_indices`, `.rds` checkpoints | `search.checkpoint` (parquet chunks with a metadata file) |
| `parallel_lapply` | `search.parallel.parallel_map` |
| `run_dataset` | `search.run.run_cohort` |
| combination of IT and FR, selection (top-level code) | `search.combine.combine_cohorts` |
| `buffa_hypoxia_score`, `score_all` | `search.signature.median_abs_rank01`, `score_samples` |
| plots, run summary, `sessionInfo` | `search plot`, `manifest.json` |

### Options

| script | Search |
|---|---|
| `--base_dir` and the file names written in the script | the cohort manifest (YAML) |
| the ten seeds written in the script | the seeds file (`examples/gbm/seeds_buffa2010_setA.csv`) |
| `--b_boot` (1000) | `resampling.n_rep` |
| `--n_mc` (5000) | `null.n_rep` |
| `--alpha` (0.05) | `edge.alpha` |
| `--mc_p_cutoff` (0.01) | `combine.mc_p_cutoff` |
| `--product_cutoff` (0.5) | `combine.product_cutoff` |
| the 51 of the top list, the variance floor 1e-8, the minimum of 3 genes for a score | `combine.top_n`, `combine.var_floor`, `score.min_genes` |
| `--boot_chunk_size`, `--mc_chunk_size` (100) | `resampling.chunk_size`, `null.chunk_size` |
| seeds `1000000 + b` and `2000000 + i` of the replicates | `resampling.seed_base`, `null.seed_base` (same bases, different generator) |
| `--cores` | `--n-jobs` |
| `--resume` | `search resume`, `search run --resume` |
| `--test_mode` | `--set resampling.n_rep=20 --set null.n_rep=50` |

### Output files

| script (`csv/`) | Search (`results/`) |
|---|---|
| `Buffa2010_<cohort>_gene_connectivity.csv` | `<cohort>_gene_connectivity.csv`, same columns plus `boot_n_ok`, `MC_count_ge`, `MC_n_null` |
| `Buffa2010_<cohort>_shared_neighbourhood_S.csv`, `..._sign_concordance_F.csv` | `<cohort>_shared_neighbourhood_S.csv`, `<cohort>_sign_concordance_F.csv` |
| `Buffa2010_<cohort>_MC_random_seed_S_summary.csv` | `<cohort>_MC_random_seed_S_summary.csv`, plus the random seeds of every replicate |
| `Buffa2010_exact_ITxFR_full_ranking.csv` | `combined_full_ranking.csv`: neutral column names (`pC_boot`, `C_IVW`, `candidate`, ...) and, for two cohorts, the script's names as well |
| `GBM_Hypoxia_Buffa2010_exact_ITxFR_signature.csv`, `..._top51.csv` | `combined_signature.csv`, `combined_signature_top.csv` |
| `Buffa2010_exact_ITxFR_HS_all_signature_IT_FR_TCGA.csv`, `..._HS_top51_...` | `scores_signature_all.csv`, `scores_signature_top.csv` |
| `Buffa2010_seed_set.csv`, `..._seed_presence_IT_FR.csv`, `run_parameters.csv`, `..._run_summary.csv` | `<cohort>_seeds.csv`, `config.yaml`, `manifest.json` |
| none | `<cohort>_seed_overlap_O.csv`, `<cohort>_seed_clusters.csv`, `gene_symbols.csv` |

### Results with the same draws

Every row is a set of numbers compared one by one between the script's `.rds` files and the
tables of Search. "Values" counts the numbers of the run with 100 and 500 replicates.

| stage | what is compared | values | 20 and 50 replicates | 100 and 500 replicates |
|---|---|---|---|---|
| core | correlation of every gene with every seed, edges, memberships, cluster sizes, both cohorts | 992,360 | equal, largest difference 1.3e-15 | equal |
| core | S and F between seeds | 400 | equal, 3.3e-14 | equal |
| core | C_raw, C_rank, edge counts, sign consistency, weighted correlation of every gene | 264,624 | equal, 6.1e-16 | equal |
| bootstrap | C_rank of every gene in every replicate | 3,307,800 | equal, 2.2e-16 | 199 of 200 replicates equal; in one the script ranks 149 tied genes apart, by up to 1.2e-3, and Search equals the ranks of the exact values |
| bootstrap | mean, variance and standard deviation over the replicates | 99,234 | equal | FR equal; IT differs in the same 149 genes, the mean by up to 1.2e-5 and the standard deviation by up to 4.1e-5 |
| null model | C_raw of every gene for every random seed set | 16,539,000 | equal, 1.0e-16 | equal |
| null model | genes eligible as random seeds, mean and median S of every random set, p-value of every gene | 68,136 | equal | equal |
| combination | per-cohort columns, rank products, inverse-variance score, the three criteria and the candidate flag | 430,014 | equal | equal except the columns derived from the IT bootstrap ranks, up to 7.4e-5 |
| combination | order of the 16539 rows | 16,539 | equal | equal for the first 1080 rows |
| selection | genes of the signature and of the top list, in order | 26 | equal | equal |
| score | median expression, score and gene counts of 588 samples (IT, FR, TCGA), two signatures | 4,704 | equal, 1.8e-15 | equal |

The same comparison on public data, so that it can be repeated without the GBM cohorts: the
script's functions `run_dataset`, `buffa2010_core_one_dataset`, `bootstrap_checkpointed` and
`monte_carlo_checkpointed`, taken from the script file unmodified and run on Breastdata1
with six seeds, 10 bootstrap and 20 Monte Carlo replicates. All 27 quantities are equal
(81,358 numbers, largest difference 3.3e-16). Two of the ten bootstrap replicates have two
seed clusters of the same size; in both, R and Search give the ranks of the exact values.
The combination and the score are top-level code of the script and can only be compared on
the GBM cohorts.

### Results with independent draws

Search with its own random generator against the same R runs:

| | 20 and 50 replicates | 100 and 500 replicates |
|---|---|---|
| core (no random numbers) | equal | equal |
| bootstrap mean of C_rank, Spearman correlation between R and Search, IT / FR | 0.879 / 0.860 | 0.965 / 0.955 |
| Monte Carlo p-value, Spearman correlation, IT / FR | 1.0000 / 1.0000 | 1.0000 / 1.0000 |
| Monte Carlo p-value, largest absolute difference, IT / FR | 0.019 / 0.017 | 0.0034 / 0.0011 |
| genes with p < 0.01, R and Search, IT | 72 and 75 | 73 and 73 |
| genes with p < 0.01, R and Search, FR | 164 and 197 | 202 and 202 |
| genes on the same side of the p < 0.01 cutoff, of 16539, IT / FR | 16536 / 16506 | 16539 / 16539 |
| selected genes, R / Search / shared | 26 / 26 / 26 | 26 / 26 / 26 |

The differences shrink as the replicates grow, as expected of sampling noise. At the size
of the preset (1000 and 5000) the stored run of Search and the lab's R run of June 2026
select the same 26 genes, with correlations of 0.999 between the bootstrap means and 0.995
between the p-values ([Validation](validation.md)); that R run used matrices with 16548
common genes instead of today's 16539, so it is not an exact comparison. The script was not
run again at that size on today's matrices: it needs about four hours on eight cores.

### Time

On the same laptop the script takes 154 s for 20 and 50 replicates per cohort (eight cores)
and 1249 s for 100 and 500 (seven cores), reading and writing included. Search takes 8 s and
57 s for the core pass, the bootstrap and the null model of the two cohorts (eight
processes), and 4 min 54 s for a whole run with 1000 and 5000 replicates. The script spends
its time in `cor(..., use = "pairwise.complete.obs")`, which ranks every pair of rows
separately.

## The SEARCH package and the `search2014` preset

### Code

| in the package | in Search |
|---|---|
| `SEARCH()`, the loop over iterations | `search.growth.run_growth` |
| one iteration: neighbourhoods of the seeds, D | `search.core.corepass.run_core`, `compute_membership` (mode `binary`), `compute_connectivity` |
| `nodes.correlation` | `search.core.association.spearman_block` or `pearson_block`, edge rule `abs_threshold` |
| `nodes.glm` | `search.core.linear.glm_loop` (statsmodels), or `ols_fwl` for the gaussian family; edge rule `p_threshold` |
| `nodes.fastLm` | `search.core.linear.ols_fwl` (least squares for every seed and gene at once) |
| `getO`, the O matrix of `get.matrix.format` | `search.core.overlap.jaccard_binary` |
| `getOmax` | `search.core.overlap.omax` |
| `removeNode`, `removeByName` | list operations inside `run_growth` |
| the correlation matrix of `get.matrix.format` | `CoreResult.neighbourhoods`, table `<cohort>_seed_clusters.csv` |
| `estimate.nodes.degree` | `search.core.connectivity.nodes_degree` |
| `SEARCH.randomise`, `.SEARCH.randomise.inner` | `search.bootstrap.run_bootstrap`, `search.resampling.resample_indices`, `search.run.run_cohort` |
| `write.data.files` | `search.export.tables` |
| `summarise.iterations`: mean over the replicates | `BootstrapSummary` (columns `C_raw_boot_mean`, `_var`, `_sd`) |
| `summarise.iterations`: p-value against runs with random seeds | `search.null.external_null`, command `search summarise` |
| `make.default.plots`, `make.networkplot`, `make.scatter.plot`, `get.adjacency.matrices` | `search plot` (different figures, not compared) |

### Arguments of `SEARCH.randomise`

| package (default) | Search (`search2014` preset) |
|---|---|
| `data` | the matrix of a cohort in the manifest |
| `seed.nodes` | the seeds file; a seed absent from the matrix is an error, the package drops it silently |
| `universe.nodes` | `universe` in the manifest; the package adds the seeds to it, Search asks that they be in it |
| `covariates`, `formula.string` | `association.covariates` (columns of the sample table), `association.formula` for extra terms |
| `scale` (TRUE) | `connectivity.normalise`: `fraction` or `count` |
| `D.lower` (0.1), `D.higher` (0.24), `O.lower` (0.01) | `growth.D_lower`, `growth.D_higher`, `growth.O_lower` |
| `keep.seeds` (TRUE) | `growth.keep_seeds` |
| `max.iter` (1) | `growth.max_iter` |
| `sd.threshold` (0.1) | `growth.sd_threshold` |
| `n.boot` (10), `n.boot.start` | `resampling.n_rep`; replicates are numbered from 1 and stored in chunks |
| `sampling` (`bootstrap`, `jacknife`) | `resampling.scheme` (`bootstrap`, `jackknife`) |
| `sampling.seed` | `resampling.seed_base` |
| `nodes.association.function` | `association.method`: `spearman`, `pearson`, `ols`, `glm` |
| `correlation.method`, `correlation.threshold` (0.5) | `association.method`, `edge.threshold` |
| `correlation.use` | none: a gene with a missing value has no association, as with the default `everything` |
| `significance.threshold` (0.01) | `edge.alpha`, with `edge.rule = p_threshold` and `edge.bonferroni = none` |
| `adjust.p` (any method of `p.adjust`) | `edge.p_adjust`: `none`, `bonferroni`, `bh` |
| `glm.family` (gaussian, Gamma, inverse.gaussian, poisson, quasi, quasipoisson) | `association.glm_family`: gaussian, gamma, inverse_gaussian, poisson |
| `mc.cores` (processes over genes) | `--n-jobs` (processes over replicates) |
| `output.directory` | `--out` |

With these correspondences the preset needs no override to stand for a call of
`SEARCH.randomise` with all its defaults, apart from its null model (see
[Deliberate differences](#deliberate-differences)).

### Output files

| package (`output/`, one file per replicate and iteration) | Search |
|---|---|
| `connectivity__iteration_<replicate>_<iteration>.txt.gz` (D rounded to five decimals) | `<cohort>_growth_connectivity.csv` (D of every iteration), `D_final` in `<cohort>_gene_connectivity.csv`, one column per replicate in `checkpoints/<cohort>/bootstrap_chunk_*.parquet` |
| `network_growth__iteration_*.txt.gz` (seeds added) | `<cohort>_growth_events.csv` (added, removed on D, removed on Omax), `<cohort>_growth_iterations.csv`, `<cohort>_growth_final_seeds.csv` |
| `shared_neighbourhood__iteration_*.txt.gz` (O, upper triangle) | `<cohort>_seed_overlap_O.csv` (symmetric, final seeds) |
| `subnetworks_correlation__iteration_*.txt.gz` | `<cohort>_seed_clusters.csv` |
| `nodes_degree__iteration_*.txt.gz` | `nodes_degree` in `<cohort>_seeds.csv` |
| `other_data__iteration_*.txt.gz` (all GLM coefficients) | `<cohort>_covariate_coefficients.csv` (coefficients of the covariates) |
| the table returned by `summarise.iterations` | `<cohort>_external_null.csv` |

In a `search2014` run the column to set against the package's connectivity file is
`D_final`. The column `C_raw` of the same table is the connectivity to the final seed set,
which includes the genes promoted in the last iteration; the package computes it only if
another iteration follows.

### Results

Each case is one call of `SEARCH()` or of `SEARCH.randomise()` plus
`summarise.iterations()`. For every case the comparison covers what the package returns: D
of every gene at every iteration, the seeds added at every iteration, the seeds removed on
D, the genes that left the universe, the members of every neighbourhood with their
correlation or coefficient, the O matrices and the node degrees. For the bootstrap cases it
also covers every replicate (on the samples drawn by R), the means and the p-values.

"As it is" and "corrected" count the quantities equal to Search's out of those compared.

| case | what it runs | R as it is | R corrected | cause of the differences, remarks |
|---|---|---|---|---|
| `breast_single_pass` | one pass, thresholds of `SEARCH.randomise`, six co-expressed seeds | 10 of 10 | | |
| `breast_growth` | growth to convergence, same thresholds | 12 of 12 | | |
| `breast_doc_example` | the correlation example of the package manual (ten seeds, three iterations) | 12 of 12 | 12 of 12 | |
| `breast_pearson` | Pearson instead of Spearman | 8 of 8 | | |
| `breast_unscaled` | `scale = FALSE` | 8 of 8 | | |
| `breast_universe_subset` | universe of 600 probesets | 8 of 8 | | |
| `it_growth` | GBM IT cohort (21254 genes), ten hypoxia seeds, three iterations | 10 of 10 | | |
| `it_sd_prefilter` | GBM IT cohort, one pass with the sd prefilter | 6 of 6 | 6 of 6 | no IT gene is below the threshold |
| `breast_sd_prefilter` | every default of `SEARCH.randomise`, sd prefilter included | 4 of 8 | 8 of 8 | R1 |
| `breast_pruned_members` | lower `D.higher`: seeds are added, then removed on D | 19 of 24 | 24 of 24 | R4 |
| `breast_omax_pruning` | higher `O.lower`: new seeds are removed on Omax | 10 of 16 | 16 of 16 | R2, R4 |
| `breast_seeds_not_immune` | `keep.seeds = FALSE` | R stops | 12 of 12 | R2, R3 |
| `breast_doc_not_immune` | the manual's example with `keep.seeds = FALSE` | R stops | 16 of 16 | R2, R3 |
| `breast_missing_values` | one missing value in 25 probesets | R stops | 24 of 24 | R7, R4 |
| `breast_glm_doc` | the GLM example of the manual (`nodes.glm`, two iterations) | 7 of 14 | 14 of 14 | R5 |
| `breast_fastlm_doc` | the same with `nodes.fastLm` | 7 of 14 | 14 of 14 | R5 |
| `breast_fastlm` | `nodes.fastLm`, six seeds; Search with `ols` and with `glm` | 10 of 10, both | | |
| `breast_glm_bh` | `nodes.glm` with `adjust.p = "BH"`; Search with `ols` and with `glm` | 10 of 10, both | | |
| `breast_glm_covariates` | `nodes.glm` with covariates `t.rfs` and `e.rfs`; Search with `ols` and with `glm` | 10 of 10, both | | 2 of 10 with `glm` before the correction of P1 |
| `breast_glm_covariates_plain` | the same covariates renamed `t_rfs` and `e_rfs`; Search with `ols` and with `glm` | 10 of 10, both | | |
| `breast_glm_covariates_ignored` | covariates passed without `formula.string`; Search without covariates | 10 of 10 | | |
| | the same; Search with the covariates declared | 4 of 10 | | D3 |
| `breast_glm_gamma` | `nodes.glm`, Gamma family | 10 of 10 | | tolerance 1e-8, largest difference 8e-12 |
| `breast_glm_poisson_loose` | `nodes.glm`, poisson family, threshold 0.9 | 10 of 10 | | tolerance 1e-8, largest difference 2e-11; the node degree differed before the correction of P3 |
| `breast_glm_poisson` | poisson family, threshold 1e-5: no association passes | R stops | | R6 |
| `breast_bootstrap` | `SEARCH.randomise`, five replicates, two iterations, null runs, `summarise.iterations` | 10 of 11 | | D2 |
| `breast_bootstrap_6iter` | eight replicates, up to six iterations | 13 of 15 | | D2, D4 |
| `breast_jackknife` | `sampling = "jacknife"` | 9 of 10 | | D2 |
| `breast_bootstrap_defaults` | `SEARCH.randomise` with every default, ten replicates | 1 of 6, then the summary stops | 9 of 10 | R1, R8; D2 |

On the GBM IT cohort the growth loop gives the same result in both languages to the last
digit: 912, 819 and 297 seeds added in the three iterations, 106 seeds removed on D, 1635
seeds holding a neighbourhood at the end, D equal for 21253, 21197 and 21148 genes. The
package takes 249 s and Search 220 s.

In total: 317 quantities compared with the package as it is (268 equal, 44 different, 5
stops) and 156 with the corrected package (155 equal; the one that differs is D2, a
deliberate difference).

## Differences and their causes

### Defects of the R package

Each was reproduced with the package as it is and, where a correction exists, removed by
changing one line. Search does not have them. R1, R2 and R3 were already in the
[decision log](decisions.md); R4 extends its entry on R3; R5 to R8 were found by this
comparison.

| | where | what happens | measured effect | correction used |
|---|---|---|---|---|
| R1 | `SEARCH.R:26-27` | The standard deviations of the sd prefilter are computed in the row order of the matrix and named in the order of the universe (seeds first), so the filter removes the wrong genes. | Breastdata1: 102 probesets removed on each side, 30 in common; one seed removed by mistake; D differs, or exists on one side only, for 157 probesets. No effect on the IT cohort, where no gene is below 0.1. | names taken from the rows |
| R2 | `getOmax.R:25` | A seed is compared only with the seeds after it in the list. A new seed misses the new seeds before it; with `keep.seeds = FALSE` the last seed always gets Omax 0 and is removed when `O.lower > 0`. | `breast_omax_pruning`: a seed removed on Omax by R is removed on D by Search, and D differs for 7 probesets. | comparison with every other seed |
| R3 | `SEARCH.R:124` | With `keep.seeds = FALSE`, an old seed removed on Omax keeps its neighbourhood list. D then has missing names, missing genes are "added" as seeds, and the next iteration stops with `subscript out of bounds`. | The package cannot complete a run with `keep.seeds = FALSE` once a seed is removed on Omax. | the list of the removed seed is dropped |
| R4 | `SEARCH.R:204`, `SEARCH.R:138` | A gene removed from the universe stays as a member in the neighbourhoods of the other seeds, and keeps counting in every later O. D is not affected, because it is counted over the current universe. | `breast_pruned_members`: O differs for 24 of 36 seed pairs, by up to 0.375 (0.375 in R, 0.75 in Search). | removed genes are dropped from every list |
| R5 | `nodes.glm.R:107-108`, `nodes.fastLm.R:93-94` | A neighbourhood with a single member loses its name (R drops names when one row and one column of a table are selected). A seed connected only to itself is then not counted in its own D, gets D = -1 / (K - 1) and is removed even with `D.lower = 0`. | The GLM example of the manual: 7 of 10 seeds removed on D; at iteration 2 D differs, or exists on one side only, for 25 probesets. No effect with the default `D.lower = 0.1`, which removes such a seed anyway. | the name is restored |
| R6 | `nodes.glm.R:109` | A seed with no significant association at all, itself included, disappears from the result (assigning an empty value to a list element deletes it). If this happens to every seed, `SEARCH` stops with `subscript out of bounds`. | `breast_glm_poisson`: R stops; Search ends the growth with every seed removed on D. | none |
| R7 | `nodes.correlation.R:31` | With a missing value the correlation is missing, and selecting with a condition that is missing adds an entry named `NA` to the neighbourhood. The run stops with `subscript out of bounds` at the second iteration. | `breast_missing_values`: R stops. Corrected, R and Search agree: both give no association to a gene with a missing value. | missing conditions are skipped |
| R8 | `summarise.iterations.R:124` | The table of null values has one row per gene of the real runs. A null run that holds any other gene stops the function with `subscript out of bounds`. | With every default on Breastdata1 this happens, because R1 removes different genes for different seed sets. | none; it does not occur once R1 is corrected |

Three more points of the decision log were confirmed by reading the code and are visible in
the outputs: the O and correlation matrices exist only for iterations followed by another
one, and O is filled above the diagonal only (`get.matrix.format.R:17`); with a single seed D
is divided by zero (`SEARCH.R:185`); the list `removed` returned by `SEARCH` holds only the
seeds removed on D.

### Deliberate differences

| | difference | effect |
|---|---|---|
| D1 | **Null model of `search2014`.** The package compares a run with other runs, made by hand with random seeds, read from folders with fixed names. The preset draws the random seed sets inside the run and pools their values, as `buffa2010` does (decision of 2026-10-12). The package's scheme remains available as `null.scheme = external_runs` with `search summarise`. | Not comparable by construction: the package has no counterpart of the pooled p-value. |
| D2 | **Genes absent from a null replicate.** For the p-value of a gene the package counts the null values at least as large as the observed one and divides by the number of all null replicates plus one, also when the gene is absent from some of them (a seed removed in that run). `search.null.per_gene_counts` divides by the number of replicates in which the gene has a value, plus one. A gene absent from every null replicate gets 1 / (replicates + 1) from the package and no p-value from Search. | Equal for every gene present in all null replicates (1219 of 1231 in `breast_bootstrap`). For the others the package gives the smaller p-value: 0.545 against 1 for a gene present in 5 of 10 null replicates. |
| D3 | **Covariates.** `nodes.glm` and `nodes.fastLm` use the covariates only when `formula.string` names them. Search uses every declared covariate. | `breast_glm_covariates_ignored`: Search with the covariates declared differs from R (D of 37 probesets); Search without covariates is equal to R. |
| D4 | **What is averaged over bootstrap replicates of a growth run.** `summarise.iterations` averages D separately for every iteration number, over the replicates that reached it. Search averages the final D of every replicate, whatever the iteration it stopped at. | Equal when every replicate stops at the same iteration, which is always so with the preset (`max_iter = 1`). In `breast_bootstrap_6iter` two replicates stop at iteration 6 and six at iteration 4: the mean of iteration 6 in R is over two replicates and differs from Search's mean over eight for 560 probesets, by up to 0.31. Iteration by iteration the two agree. |
| D5 | **Rounding in the files.** The package writes D rounded to five decimals and averages what it reads back. Search keeps full precision. | Means differ by at most 4.3e-6. |
| D6 | **Missing values in the `buffa2010` path.** The script correlates over the samples that are present (`pairwise.complete.obs`) and keeps the full number of samples in the p-value. Search gives no association to a gene with a missing value, and the quality check warns about it. | None on the GBM matrices, which have no missing value. On Breastdata1 with one missing value in 30 probesets, five of them inside seed clusters: 8 edges lost, one cluster from 5 to 2 genes, C_raw different for 39 probesets, and the rank of nearly all of them, because ranks are relative. |
| D7 | **Self edge.** A seed always belongs to its own cluster (`edge.self_edge`). In both originals this is a consequence of the seed correlating perfectly with itself, and fails when it does not pass the edge rule (R5, R6). | None beyond R5 and R6. |
| D8 | **Edge decisions on rounded values.** Search compares statistics and p-values rounded to 12 decimals with the threshold. | No edge differs from R in any comparison: the 330,780 gene and seed pairs of the two GBM cohorts, the replicates built on them, and every neighbourhood of the package cases. |
| D9 | **Stricter inputs.** Seeds missing from the matrix, fewer than two seeds and a growth that removes every seed are errors in Search; the package continues with what is left. | Not a numerical difference. |
| D10 | **Not available in Search.** Kendall correlation, the GLM families `quasi` and `quasipoisson`, the `p.adjust` methods other than `bonferroni` and `BH`, and values of `correlation.use` other than the default. | |
| D11 | **Missing values in the GLM path.** R's `glm` leaves out the samples in which the gene has a missing value and fits on the others. Search gives the gene no association, as it does with every other method. | Read in the code, not run: the package was not run with `nodes.glm` on data with missing values. |

### Random numbers

`set.seed(1000000 + b)` in R and `numpy.random.default_rng(1000000 + b)` in Search start
different generators, so replicate *b* resamples different columns in the two languages.
The comparison with the same draws shows that this is the only source of difference in the
bootstrap and in the null model, and the comparison with independent draws shows its size
([Results with independent draws](#results-with-independent-draws)). Search stores the
sample indices of every bootstrap replicate and the random seeds of every null replicate in
its checkpoints, so a run can be re-examined replicate by replicate; the script stores
neither.

### Floating point

Numbers reported as equal differ by at most 3.3e-14 (S between two seeds, a ratio of two
sums over 16539 genes) and typically by 1e-16. Coefficients of the Gamma and poisson
families differ by up to 2e-11, because the iterative fit stops at slightly different
points; no neighbourhood changes because of it.

One effect is larger and was found only with 100 bootstrap replicates. C_raw of a gene is
the sum, over the seeds it is connected to, of one divided by the size of the seed's
cluster. When two seeds have clusters of the same size, genes connected to different
combinations of seeds can have the same connectivity: equal in exact arithmetic, but a
floating-point sum depends on the order of the additions, and the two sums can differ in the
last digit. Genes that should share a rank then get different ranks.

In replicate 22 of the IT cohort two seeds have clusters of 2548 genes, and there are 308
distinct connectivity values in exact arithmetic, computed here with fractions. The script
adds the terms in the order of the seeds and gets 318 distinct numbers: ten ties are split,
149 genes are up to 0.0012 away from the rank of the exact values. Search adds the terms of
each gene from the smallest to the largest, so that equal terms give equal sums: it gets 308
distinct numbers and the ranks of the exact values, whatever the order of the seeds. This is
the only replicate of the 200 in which two clusters have the same size, and the only one in
which the two programs differ. The difference moves the IT bootstrap mean of 149 genes by up
to 1.2e-5 and, through it, the order of the combined ranking from row 1081 on. Equal cluster
sizes do not always split a tie: on Breastdata1 two replicates of ten have them and the
script agrees with the exact ranks.

Before its correction Search added in the order of the seeds too, with NumPy's own grouping
of the additions: 314 distinct numbers in that replicate, 87 genes away from the exact
ranks, and a rank that changed with the order of the seeds (P2 below).

## Defects found in the port and corrected

P1 to P4 were found by the comparison in commit `5177011`; P5 and P6 while correcting P1.
All were corrected on 2026-11-17. Each has a test that failed before the correction, and
the comparison was run again afterwards: nothing that was equal became different.

| | where | defect | correction | checked by |
|---|---|---|---|---|
| P1 | `glm_loop` in `src/search/core/linear.py` | With `association.method = glm`, covariate names were pasted into a statsmodels formula as they are. A name that is not a valid Python name, such as `t.rfs`, made every fit fail, and the failures were caught and ignored: the run ended with an empty network and no error. | Covariate names are quoted in the formula. A model that cannot be built is a configuration error; a fit that fails is counted and logged; a run in which every fit fails stops. The covariates are checked as in the `ols` path. | `test_glm_reads_covariates_whatever_their_names`, `test_glm_failures_are_errors_not_empty_networks`; `breast_glm_covariates` with `glm`, from 2 to 10 of 10 |
| P2 | `compute_connectivity` in `src/search/core/connectivity.py` | C_raw was added in the order of the seeds. When two seed clusters have the same size, connectivities that are equal in exact arithmetic could differ in the last digit ([Floating point](#floating-point)), tied genes got different ranks, and which ones depended on the order of the seeds: in the IT replicate above, listing the seeds in reverse changed C_rank of 173 genes by up to 0.0012. | The terms of each gene are added from the smallest to the largest. | `test_tied_connectivities_stay_tied_whatever_the_seed_order`, `test_seed_order_does_not_matter`; in the IT replicate Search equals the ranks of the exact values, and six other orders of the seeds give the same C_raw bit for bit |
| P3 | `nodes_degree` in `src/search/core/connectivity.py` | It counted in how many seed neighbourhoods a seed sits; `estimate.nodes.degree` counts how many seeds sit in its neighbourhood. The two are the same number when the association is symmetric (correlation, least squares), and not for the other GLM families. | It counts as the package does. | `test_nodes_degree_counts_the_seeds_inside_the_neighbourhood`; `breast_glm_poisson_loose`, the six degrees equal R's |
| P4 | `per_gene_counts` in `src/search/null.py` | With per-gene pooling a gene without an observed value had no null value above it and so the smallest possible p-value, 1 / (replicates + 1). The pooled path already set it to missing; the package has the same defect. It showed through `search summarise`, for a gene whose bootstrap mean is missing. | Such a gene has no p-value. | `test_null_counting_helpers`, `test_a_gene_without_observed_value_gets_no_p_value` |
| P5 | `glm_loop` in `src/search/core/linear.py` | A gene with the same value in every sample is collinear with the intercept; statsmodels still returns a coefficient and a p-value for it, and on expression-like values (mean 8, standard deviation 1) every such gene got an edge with every seed, with p-values down to 8e-42. A gene with a missing value was fitted on the remaining samples. The rule of the engine is that neither has an association. | Both are left without association, as in every other method. | `test_glm_gives_no_association_to_genes_the_engine_cannot_use` |
| P6 | `ols_fwl` in `src/search/core/linear.py` | A constant gene was recognised by residuals summing to exactly zero, which rounding prevents: it got a meaningless statistic, the ratio of two rounding residues. No edge came from it in any check made (0 of 6000 pairs). | A constant gene is recognised by its largest value being its smallest, before the residuals are taken. | the same test |

What the corrections change in results already produced. For `buffa2010` only P2 matters.
The run on the GBM cohorts with 1000 and 5000 replicates was repeated: C_raw, the observed
C_rank, the Monte Carlo p-values and the scores are unchanged; the bootstrap mean of C_rank
moves for 285 genes in IT and 618 in FR, by at most 5.4e-7 and 1.5e-6; the 26 selected genes
are the same, in the same order, and 4 of the 16539 rows of the full ranking change place,
from row 4578 on. For `search2014` runs with correlations nothing changes. Runs with
`association.method = glm` can change through P1 and P5, and their `nodes_degree` through P3.

Not corrected, and already known: the growth loop is slow on a real-size cohort (220 s for
three iterations on IT, against 249 s for the R package) because
`_membership_from_neighbourhoods` (`src/search/growth.py:137`) tests every neighbourhood
against the universe once per seed. A smaller point of the same kind as P6 is left as it is:
with Pearson correlation a constant gene whose value is not exactly representable in binary
gets a correlation of about 1e-17 where a missing value is promised; no edge can come from
it.

## What was not compared

- The figures of both programs.
- The package with one seed, with the `inverse.gaussian` family, with `n.boot.start` other
  than 1, with `nodes.glm` on data with missing values (D11), and with the options that
  Search does not have (D10).
- The GBM script with 1000 and 5000 replicates on today's matrices (about four hours of R).
- Operating systems other than macOS on arm64. The last digits, and with them the tie effect
  described under [Floating point](#floating-point), can differ on another platform or with
  another numerical library.

## Repeating the comparison

```bash
export SEARCH_GBM_DATA=/path/to/GBM_HYPOXIA_SIGNATURE   # optional: without it only the public data
validation/r_parity/run_all.sh            # 20 and 50 replicates for the GBM script, about ten minutes
validation/r_parity/run_all.sh 100 500 8  # about half an hour
```

The reports are in `validation/r_parity/results`: one row per compared quantity, with the
number of values, how many differ, the largest difference and, where something differs, the
first example. `validation/r_parity/README.md` describes the scripts and their requirements.
The reports in the repository were produced on 2026-11-17, after the corrections of the
port, with the versions listed at the top of this page.
