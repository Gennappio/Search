# Decision log

Deviations from the R originals and open questions, with the default chosen in the meantime.

## Defects of the originals, fixed without a compatibility option

1. **Standard-deviation prefilter named by universe order** (SEARCH.R). The R code computed the
   standard deviations in matrix row order and named them in universe order, so the two were
   misaligned unless the universe followed the row order. Search names them by gene.
2. **Triangular Omax** (getOmax.R). Seed *i* was compared only with seeds after it in the list,
   so later seeds were under-compared and, with seeds not protected, the last one always got 0.
   Search takes the maximum over every other seed.
3. **Matrices stored only in iterations that changed something** (SEARCH.R). Search stores every
   iteration.
4. **Division by zero with one seed**. Search raises an explicit error: connectivity needs at
   least two seeds.
5. **Pruned old seeds left in the neighbourhood lists** (SEARCH.R with keep.Alist = FALSE),
   which produced NA names in D. Search removes them everywhere.
6. **Covariates ignored unless a formula string was given** (nodes.glm.R). Search uses the
   declared covariates always (phase 4).
7. **Correlation with pairwise handling of missing values** (GBM script). The script paid for
   per-pair ranking on data without missing values, 12 times the cost of the rest of the pass.
   Search ranks once per gene; genes with missing values get no edges and the QC warns.

## Open statistical questions, parametrised, default = the GBM script

- **Membership normalisation** (`membership.mode`). `cluster_size` (script), `binary` (SEARCH)
  and `gene_normalised` (paper, equation 2 as printed) give different rankings. With step edges
  `gene_normalised` is degenerate: every connected gene gets 1/K. Default `cluster_size`.
- **Null pooling** (`null.pooling`). `pooled` throws every gene's null values into one
  distribution; `per_gene` compares each gene with its own. Default `pooled`.
- **Selection criteria** (`combine.*`). Rank product > 0.5 and the top-51 cap are not in the
  paper; p < 0.01 and positive consistency are. All free parameters.
- **Self edge** (`edge.self_edge`, `connectivity.exclude_self`). Both presets keep the seed in
  its own cluster and exclude it from its own connectivity, as both originals did.

## Choices fixed by the engine

- p-values of correlations through the t approximation; average ranks for ties; strict
  inequalities at thresholds; the Bonferroni denominator is the size of the universe.
- The fractional rank of tied values is their average rank, so the lowest rank is not 0 when
  many genes share the lowest connectivity.
