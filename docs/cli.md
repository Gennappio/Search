# Command line

```
search run --manifest cohorts.yaml [--preset buffa2010] [--set key.path=value ...]
           [--out runs] [--n-jobs N] [--resume] [--strict] [--name NAME]
search resume [runs/latest] [--n-jobs N]
search qc --manifest cohorts.yaml
search combine [runs/latest] [--mc-p-cutoff P] [--product-cutoff X] [--top-n N]
               [--no-positive-consistent] [--keep-seeds] [--prefix recombined]
search score --manifest cohorts.yaml --signature genes.csv [--out scores.csv] [--method M]
search compare <run_a> <run_b> [--top-n 50] [--out metrics.csv]
search variants --manifest cohorts.yaml --grid membership.mode=cluster_size,binary,gene_normalised
                [--grid null.pooling=pooled,per_gene] [--preset P] [--set ...] [--out runs/variants]
search plot [runs/latest]
search config show [--preset NAME | --file cfg.yaml] [--set key.path=value ...]
search config validate cfg.yaml
```

`variants` runs the pipeline once per combination of the grid, with the first combination as
baseline, and writes `variants_summary.csv`: candidates per variant, Spearman correlation and
top-N overlap of the combined score with the baseline, Jaccard of the signatures, top-10 genes.
This is the table to bring to a domain reviewer for the open statistical questions.
`compare` reports the same agreement measures between any two runs. `plot` draws the seed
networks (S and O), the rank scatter between the first two cohorts, the product-score histogram
and the score distribution per cohort (needs the `plots` extra).

`--set` takes dotted paths into the configuration and YAML values: `--set null.pooling=per_gene`,
`--set resampling.n_rep=200`, `--set association.covariates=[batch]`. Any configuration outside a
validated preset runs with a warning and `preset_validated: false` in the manifest; `--strict`
refuses it. `--n-jobs` defaults to 1, or to `SLURM_CPUS_PER_TASK` when set.

## Cohort manifest

```yaml
name: gbm_hypoxia
seeds: seeds.csv                 # columns: gene[,symbol]
universe: null                   # optional gene list restricting the analysis
gene_filter:
  strip_ensembl_version: true    # ENSG00000000003.15 -> ENSG00000000003
  common_genes_across_training: true
  annotation_table: annotation.csv   # optional; gives symbols and allows a query
  annotation_query: null             # pandas query selecting the universe, e.g. "biotype == 'protein_coding'"
cohorts:
  - name: IT
    role: training               # training cohorts build the network; scoring cohorts only get scores
    matrix: DATA/IT.csv          # genes x samples (csv, tsv, parquet, h5ad), relative to this file
    format: csv                  # optional, inferred from the extension
    orientation: genes_x_samples # or samples_x_genes
    sample_table: DATA/IT_samples.csv   # optional; first column = sample id; batch, patient, tissue, time, event
  - {name: FR, role: training, matrix: DATA/FR.csv}
  - {name: TCGA, role: scoring, matrix: DATA/TCGA.csv}
config:
  preset: buffa2010              # or file: config.yaml
  overrides: {resampling.n_rep: 1000, null.n_rep: 5000}
```

## Run directory

```
runs/<name>_<timestamp>/
  config.yaml                 the full configuration of the run
  manifest.json               versions, input hashes, QC warnings, timings, outputs, status
  log.txt
  checkpoints/<cohort>/       bootstrap_chunk_*.parquet, null_chunk_*.parquet (+ .meta.json)
  results/
    <cohort>_gene_connectivity.csv      C_raw, C_rank, edges, bootstrap summary, p_MC_C
    <cohort>_shared_neighbourhood_S.csv, <cohort>_sign_concordance_F.csv, <cohort>_seed_overlap_O.csv
    <cohort>_seed_clusters.csv          seed, gene, statistic, p, gamma
    <cohort>_MC_random_seed_S_summary.csv
    combined_full_ranking.csv, combined_signature.csv, combined_signature_top.csv
    scores_signature_all.csv, scores_signature_top.csv
runs/latest/run_path.txt      points to the most recent run
```

`search resume` reuses every checkpoint whose metadata (configuration hash, input hash, replicate
range, seed) still matches; anything else is recomputed.
