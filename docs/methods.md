# Method

Search computes, from an expression matrix (genes x samples) and a set of seed genes, a
genes x seeds **membership matrix**; everything else derives from it. One pass of the engine
(`search.core.corepass.run_core`) has five steps. Growth, bootstrap and the null model call
the same pass repeatedly with different seeds, columns or random seed sets.

## 1. Association

For every gene *g* and seed *j* an association statistic. Spearman (default) ranks every row
once, centres and scales it to unit norm and takes one matrix product, so the whole genes x seeds
block costs one multiplication. Pearson skips the ranking. Linear models with covariates arrive in
phase 4. Genes with missing values or zero variance get no association (NaN) and never an edge.

The two-sided p-value uses the t approximation with *n - 2* degrees of freedom:
`t = rho sqrt((n - 2) / (1 - rho^2))`, `p = 2 P(T > |t|)`.

## 2. Edges (delta)

`abs_threshold`: an edge where |stat| > threshold (SEARCH 2014, 0.5).
`p_threshold`: an edge where p < alpha, with alpha divided by the number of genes of the universe
when `bonferroni = n_genes` (Buffa 2010, 0.05), or with per-seed adjusted p-values (`p_adjust`).
The seed's own entry is forced to `self_edge`.

## 3. Membership (gamma)

| mode | gamma[g, j] | origin |
|---|---|---|
| `binary` | delta[g, j] | SEARCH 2014 |
| `cluster_size` | delta[g, j] / number of genes in cluster *j* | GBM script |
| `gene_normalised` | delta[g, j] / number of seeds *g* is connected to | Buffa 2010, eq. 2 as printed |
| `sigmoid` | delta[g, j] / (1 + exp(-k (|stat| - x0))) | soft membership |

With step edges `gene_normalised` gives every connected gene the same connectivity (see the
decision log).

## 4. Connectivity

`C_raw[g] = sum_j w[g, j] gamma[g, j] / sum_j w[g, j]`, with `w = 0` on a seed's own column
when `exclude_self`. With binary membership this is the D of SEARCH 2014: seeds divided by
K - 1, other genes by K. `rank = fractional` adds `C_rank`, the fractional rank in [0, 1]
with average ranks for ties.

Per gene the engine also reports the number of edges, of positive and negative edges,
`positive_consistent` (at least one edge and no negative one) and the membership-weighted
association.

## 5. Seed overlaps

`O` is the binary Jaccard index between the clusters of two seeds, all genes counted
(SEARCH 2014). `S` and `F` follow Buffa 2010, equations 3 and 4, with the rows of the two seeds
excluded: `S = sum min(gamma_a, gamma_b) / sum max(gamma_a, gamma_b)`,
`F = sum min(gamma_a, gamma_b) sign(stat_a stat_b) / sum min(gamma_a, gamma_b)`.
With binary membership `S` is the Jaccard index on the non-seed rows. `S` and `F` cost one pass
over the genes per seed pair and are computed automatically up to `overlap.auto_max_seeds` seeds.

## 6. Recursive growth (SEARCH 2014)

`growth.enabled` wraps the pass in a loop over the seed set. Each iteration:

1. the neighbourhoods of the seeds added in the previous iteration are computed with the engine
   (all seeds at iteration 1);
2. from iteration 2, the binary overlap O between every pair of seeds is computed from the stored
   neighbourhoods restricted to the current universe; the seeds under check (the new ones, or all
   of them when `keep_seeds` is false) whose largest overlap with any other seed is below
   `O_lower` are pruned;
3. the connectivity D of every universe gene to the current seeds is computed from the
   neighbourhoods with the configured membership and connectivity rules;
4. seeds with D below `D_lower` are pruned; genes outside the seed set with D above `D_higher`
   become seeds;
5. the loop stops when nothing changed, when fewer than two seeds remain, when the universe is
   empty, or at `max_iter`.

Pruned genes leave the universe and every neighbourhood for good. `sd_threshold` removes flat
genes before the loop, keyed by gene name. The final connectivity, aligned to the initial
universe, is the statistic that bootstrap and the null model resample when growth is enabled.
The tables of a run then also contain `D_final`, the iteration table, the events (added,
removed by D, removed by O) and the connectivity of every iteration.

## 7. Null model and p-values

`internal_random_seeds` repeats the pass (or the growth loop) with random seed sets of the same
size drawn from the genes with standard deviation above `eligible_min_sd`, the real seeds
excluded. The empirical p-value of a gene is `(1 + count_ge) / (1 + n_null)`, where `count_ge`
counts null values greater than or equal to the observed statistic (within 1e-12). With
`pooled`, every gene's null values form one distribution (n_null = genes x replicates); with
`per_gene`, only the gene's own replicates count (n_null = replicates). Both are computed from
the same stored chunks. `external_runs` reads the per-replicate statistics stored by other runs,
typically runs with random seeds, through `search summarise --null-runs`.

## 8. Bootstrap, cohort combination and score

Bootstrap draws the samples with replacement (jackknife: all but one), recomputes the statistic
per replicate from a seed derived as `seed_base + replicate`, stores replicates in chunks and
summarises mean, variance and standard deviation per gene. For several training cohorts the
product of the bootstrap-mean ranks (`pC_boot`), the product of the observed ranks (`pC_obs`)
and the inverse-variance weighted mean (`C_IVW`) are computed; a gene is a candidate when it is
not a seed, is positively consistent in every cohort, has a null p-value below `mc_p_cutoff` in
every cohort and `pC_boot` above `product_cutoff`; candidates are sorted by `pC_boot`, `pC_obs`,
`C_IVW` and capped at `top_n`. The per-sample score is the median of the absolute expression of
the signature genes, ranked between 0 and 1 across the samples of each cohort.
