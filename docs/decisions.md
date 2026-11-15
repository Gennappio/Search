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
   declared covariates always; they are columns of the sample table, which travels with the
   matrix and follows every resampling.
7. **Correlation with pairwise handling of missing values** (GBM script). The script paid for
   per-pair ranking on data without missing values, 12 times the cost of the rest of the pass.
   Search ranks once per gene; genes with missing values get no edges and the QC warns.

## Open statistical questions, parametrised, default = the GBM script

- **Membership normalisation** (`membership.mode`). `cluster_size` (script), `binary` (SEARCH)
  and `gene_normalised` (paper, equation 2 as printed) give different rankings.
  `gene_normalised` is degenerate by construction: the memberships of a gene sum to 1, so with
  unit weights equation 5 gives 1/K to every non-seed gene with at least one nonzero membership,
  whatever the edge rule (with step edges the connected genes get 1/K and the others 0). Table 2
  of the paper shows graded scores, which the formulas as printed cannot give; the script's
  normalisation within each seed cluster does. Default `cluster_size`.
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
- Edge decisions compare statistics and p-values rounded to 12 decimals with the threshold:
  with few samples a Spearman correlation can sit exactly on 0.5, and without rounding the
  floating-point noise of the matrix product, which depends on the order of samples, decided
  the edge.

## Operational decisions, 2026-10-09

- **One scorer, no external scoring package.** The per-sample score is the built-in
  median-of-absolute-expression rank of Buffa 2010; Search depends on no other scoring project.
- **Strict mode** refuses configurations outside a validated preset and QC warnings about the
  data; a missing sample table is informational, otherwise cohorts without clinical tables could
  never run strictly.
- **A growth loop that prunes every seed** is a failed run with an explicit error, not a result:
  the growth tables are written first so the pruning can be inspected.
- **External null per cohort.** Replicates of different cohorts are never mixed; `summarise`
  matches cohorts by name across runs.
- **Tag v0.1.0** points at the phase 5 commit, two commits before the fixtures were added: a
  fresh clone of that tag fails two tests. It is left in place; v0.1.1 is the first complete tag.

## Method decisions, 2026-10-12

- **`search2014` uses the null model of `buffa2010`** (random seed sets drawn inside the run,
  pooled, 5000 replicates), following the rule that `buffa2010` decides wherever a choice is
  open. The R package compared a run with other runs made by hand with random seeds. Ported as
  it was, every random seed set is a run of its own, and a run whose growth prunes every seed is
  a failed run that leaves no replicates: 22 of 40 random sets of ten seeds on a GBM cohort of
  21254 genes and 95 samples, 40 of 40 random sets of six on Breastdata1. The null was then built
  from the surviving sets only, or could not be built at all. Inside the run a random seed set
  that loses its seeds is a replicate like any other. `null.scheme = external_runs` and
  `search summarise` remain available.

## Found by the comparison with R, 2026-11-15

The comparison described in [Porting from R](porting.md) ran both originals next to Search.
It confirmed entries 1, 2 and 5 above with numbers, and added the following.

Defects of the SEARCH package that were not recorded. Search never had them; they explain
differences from the package, and three of them stop it with an error:

8. **Pruned genes stayed members of the other neighbourhoods** (SEARCH.R, any value of
   keep.Alist). Entry 5 is the special case of an old seed removed on Omax. A gene removed on D
   or on Omax kept counting in every later overlap O; D was not affected.
9. **A neighbourhood of one member lost its name** (nodes.glm.R, nodes.fastLm.R). A seed
   associated only with itself was not counted in its own D, got D = -1 / (K - 1) and was removed
   even with D.lower = 0. Search keeps the self edge whatever the association method.
10. **A seed with no significant association disappeared from the result** (nodes.glm.R), and
    when this happened to every seed the run stopped with an error.
11. **Missing values added an entry named NA to a neighbourhood** (nodes.correlation.R), and the
    run stopped at the second iteration. With that corrected the package and Search agree: a
    gene with a missing value has no association.
12. **summarise.iterations stopped on a gene present only in a null run.** With every default
    argument on Breastdata1 this follows from entry 1.

Choices of the engine that differ from the package and were not recorded:

- **Genes absent from a null replicate** (`search.null.per_gene_counts`). The p-value divides
  by the number of null replicates in which the gene has a value, plus one. The package divided
  by the number of all null replicates plus one, so a gene absent from half of them could get
  0.55 where Search gives 1.
- **Bootstrap of a growth run.** The replicates are summarised on the final D of each. The
  package summarised each iteration number separately, over the replicates that reached it.
  The two coincide when every replicate stops at the same iteration, as with `max_iter = 1`.

Open defects of the port, found by the same comparison and listed with their evidence in
[Porting from R](porting.md#defects-found-in-the-port): covariate names that are not valid
Python names silently empty the network with `association.method = glm`; `C_rank` depends on
the order of the seeds when two seed clusters have the same size, because mathematically tied
connectivities can differ in the last digit; `nodes_degree` is transposed for association
methods that are not symmetric; with per-gene pooling a gene without an observed value gets
the smallest p-value. None is corrected in version 0.1.1.

