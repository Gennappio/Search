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
