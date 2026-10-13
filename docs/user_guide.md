# SEARCH user guide

A practical guide to preparing inputs, running an analysis and understanding the results.
Written for the current Python implementation, version 0.1.1, on 13 October 2026.

**For the GBM workflow, start one run with IT and FR listed together. SEARCH analyses each
cohort separately, then combines their gene-level results. TCGA can be supplied in the same
run to receive scores from the resulting signature.**

The `buffa2010` preset follows the selection rules of the GBM R script. The procedure in
`GBM_hypoxia_signature.pdf`, which selects 439 genes and compacts them to 93, has additional
method choices. Section 9 explains the difference.

## Contents

1. [What SEARCH does](#1-what-search-does)
2. [Prepare the input files](#2-prepare-the-input-files)
3. [Describe your cohorts in a manifest](#3-describe-your-cohorts-in-a-manifest)
4. [Run the GBM analysis on this computer](#4-run-the-gbm-analysis-on-this-computer)
5. [Understand the analysis steps](#5-understand-the-analysis-steps)
6. [Find and interpret the results](#6-find-and-interpret-the-results)
7. [Resume an interrupted run](#7-resume-an-interrupted-run)
8. [Reuse results or score another cohort](#8-reuse-results-or-score-another-cohort)
9. [Understand the R script and report differences](#9-understand-the-r-script-and-report-differences)
10. [Troubleshooting](#10-troubleshooting)
11. [Run on another computer or an HPC cluster](#11-run-on-another-computer-or-an-hpc-cluster)
12. [Terms and further reading](#12-terms-and-further-reading)

## 1. What SEARCH does

SEARCH starts with a short list of **seed genes** associated with a biological process.
It measures how other genes follow those seeds across samples, evaluates the stability and
statistical evidence of those relationships, and produces gene-level results.

For the GBM hypoxia workflow, the seeds are ADM, AK4, BNIP3, CA9, ENO1, HK2, LDHA, PGK1,
SLC2A1 and VEGFA. The `buffa2010` preset combines evidence across construction cohorts,
selects a gene signature and calculates a score for each sample.

| Your task | Command | What happens |
|---|---|---|
| Check the inputs | `search qc` | Loads the inputs and reports quality warnings |
| Derive a signature | `search run` | Runs the analysis on all training cohorts, combines results and scores samples |
| Continue an interrupted run | `search resume` | Reuses valid saved checkpoints |
| Change final selection thresholds | `search combine` | Reuses completed gene-level results |
| Apply an existing signature | `search score` | Calculates sample scores without deriving a new signature |
| Generate figures | `search plot` | Writes plots for an existing run |

The cohort roles are important:

| Cohort | Role in the GBM manifest | Contribution |
|---|---|---|
| IT | `training` | Builds gene-level evidence for the signature |
| FR | `training` | Builds independent gene-level evidence for the same signature |
| TCGA | `scoring` | Receives sample scores after the signature is selected |

The samples of IT and FR stay in separate matrices. Listing both in one manifest does not
pool their samples. The software matches the results by gene identifier and combines those
results after analysing each cohort.

A scoring cohort's expression values do not enter the derivation calculations. Changing
its role to `training` would make it part of construction and change the analysis design.

## 2. Prepare the input files

### Expression matrices

Supply **one already-normalised expression matrix per cohort**. For a CSV matrix, the first
column contains gene identifiers, each remaining column is a sample, and the cells contain
numeric expression values. This is called `genes_x_samples` orientation.

Example of the file structure, with illustrative values:

```csv
gene,IT_sample_1,IT_sample_2,IT_sample_3
ENSG00000148926,5.2,6.1,4.8
ENSG00000107159,3.0,4.7,2.9
ENSG00000112715,7.4,8.0,6.5
```

SEARCH accepts CSV, TSV and Parquet; H5AD support requires the optional `h5ad` dependency.
For a table with samples in rows, declare `orientation: samples_x_genes` in the manifest.

Before running:

- Apply the appropriate normalisation and batch correction upstream. SEARCH does not run
  DESeq2, VST or ComBat. The existing GBM files are already transformed: IT and FR use
  DESeq2 VST plus ComBat, and TCGA uses DESeq2 VST.
- Use the same identifier system in every training matrix and in the seed list. An Ensembl
  identifier such as `ENSG00000148926` and a symbol such as `ADM` are different strings;
  supplying an annotation table does not automatically translate matrix identifiers.
- Keep gene identifiers unique. Removing Ensembl version suffixes can expose duplicates;
  resolve those duplicates upstream.
- Check sample identity. Multiple measurements from the same patient are not independent
  patients. SEARCH's default bootstrap resamples sample columns and does not automatically
  group, merge or remove repeated patient samples.

### Seed list

The seed CSV needs a `gene` column. A `symbol` column is optional and helps label results.
Use identifiers that match the expression matrices. At least two usable seeds are required
for the analysis; the GBM workflow uses ten.

Every seed listed in the file must be present in each training cohort after filtering;
SEARCH raises an error if a requested seed is absent. For the intended GBM comparison, keep
all ten reference seeds. Deliberately removing seeds from the file changes the biological
reference used by the analysis.

```csv
gene,symbol
ENSG00000148926,ADM
ENSG00000162433,AK4
ENSG00000176171,BNIP3
ENSG00000107159,CA9
ENSG00000074800,ENO1
ENSG00000159399,HK2
ENSG00000134333,LDHA
ENSG00000102144,PGK1
ENSG00000117394,SLC2A1
ENSG00000112715,VEGFA
```

`strip_ensembl_version: true` removes suffixes from matrix identifiers, for example
`ENSG00000148926.12` becomes `ENSG00000148926`. Use stable, unversioned IDs in the seed CSV.

### Optional sample and annotation tables

A sample table has one row per sample, with the sample identifier in its first column.
Its identifiers must cover the samples in the expression matrix. Columns such as `patient`,
`batch` and `tissue` help the quality checks; covariates can also be used in configured linear
models. Adding these columns alone does not change the `buffa2010` Spearman analysis or
automatically correct batch effects.

An annotation table can provide gene symbols. **A table path alone does not apply a
consensus filter.** An `annotation_query` selects allowed genes; alternatively, `universe`
can point to an explicit gene list. Any such restriction can affect the final signature.

## 3. Describe your cohorts in a manifest

A **manifest** is a YAML file listing the inputs, cohort roles and analysis settings. YAML
uses indentation: use spaces rather than tabs, and preserve the nesting shown below.

For a new project, a minimal GBM-style manifest could be:

```yaml
name: gbm_hypoxia
seeds: seeds.csv
gene_filter:
  strip_ensembl_version: true
  common_genes_across_training: true
cohorts:
  - name: IT
    role: training
    matrix: data/IT.csv
    format: csv
    orientation: genes_x_samples
  - name: FR
    role: training
    matrix: data/FR.csv
    format: csv
  - name: TCGA
    role: scoring
    matrix: data/TCGA.csv
    format: csv
config:
  preset: buffa2010
```

This is a template: replace the filenames with your actual files. Relative input paths are
resolved **from the directory containing the manifest**, regardless of the terminal's current
directory. The default orientation is `genes_x_samples`. Cohort names must be unique, and a
valid manifest must contain at least one `training` cohort.

`common_genes_across_training: true` restricts IT and FR to genes present in both matrices.
It aligns the genes, while leaving each cohort's samples separate. TCGA is not part of this
training-gene intersection.

**For the existing local GBM data, use the ready-to-run manifest:**

```text
/Users/gennaroabbruzzese/Activities/SEARCH/GBM_HYPOXIA_SIGNATURE/RESULTS/SEARCH_PY/cohorts_used.yaml
```

It already points to the three normalised matrices and the seed list. It loads an annotation
table for symbols, but currently has no `annotation_query`. The separate template under
`examples/gbm/` needs its data paths adjusted before use from its current location.

## 4. Run the GBM analysis on this computer

Run the following commands in a terminal. Keep the same terminal open: the variables defined
in step A are reused by the later commands. These examples create new output directories and
leave the archived analyses available for comparison.

### A. Set the paths and check the installation

```bash
SEARCH_PROJECT="/Users/gennaroabbruzzese/Activities/SEARCH"
SEARCH_BIN="$SEARCH_PROJECT/.venv/bin/search"
SEARCH_MANIFEST="$SEARCH_PROJECT/GBM_HYPOXIA_SIGNATURE/RESULTS/SEARCH_PY/cohorts_used.yaml"
SEARCH_SMOKE_OUT="$SEARCH_PROJECT/GBM_HYPOXIA_SIGNATURE/RESULTS/SEARCH_PY/user_smoke"
SEARCH_FULL_OUT="$SEARCH_PROJECT/GBM_HYPOXIA_SIGNATURE/RESULTS/SEARCH_PY/user_full"

"$SEARCH_BIN" --version
"$SEARCH_BIN" run --help
```

`SEARCH_BIN` points directly to the installed command, so environment activation is unnecessary.
If it does not exist, follow the installation instructions in section 11.

### B. Check the inputs

```bash
"$SEARCH_BIN" qc --manifest "$SEARCH_MANIFEST"
```

Check the cohort names, roles, sample counts and seed presence. With the current local
matrices, the check reports IT with 95 samples, FR with 279 and TCGA with 214. IT and FR
share 16,539 loaded genes. This differs from the 16,548-gene universe in the archived R
analysis; record the actual input versions when comparing results.

The current IT file triggers warnings about control-like and replicate-like sample names.
No sample tables are supplied in the ready-to-run manifest. These name-based checks are clues,
not confirmed sample classifications; clarify the sample identities before biological
interpretation. SEARCH reports warnings and continues by default. `run --strict` refuses
non-informational data warnings, so it currently refuses these inputs. A missing sample table
alone is informational.

### C. Make a short trial run

```bash
"$SEARCH_BIN" run \
  --manifest "$SEARCH_MANIFEST" \
  --preset buffa2010 \
  --set resampling.n_rep=20 \
  --set null.n_rep=50 \
  --out "$SEARCH_SMOKE_OUT" \
  --n-jobs 8
```

This checks that loading, analysis and export work together. The reduced repetition counts
make it a **smoke test**, not the final statistical analysis. Similar gene counts in a smoke
test and a full run do not establish that the shorter analysis is sufficient.

### D. Run the full analysis

```bash
"$SEARCH_BIN" run \
  --manifest "$SEARCH_MANIFEST" \
  --preset buffa2010 \
  --out "$SEARCH_FULL_OUT" \
  --n-jobs 8
```

The preset uses **1,000 bootstrap resamples and 5,000 random seed sets per training cohort**.
One command handles both IT and FR, combines their results and scores IT, FR and TCGA.

`--n-jobs 8` means eight worker processes for repetitions. Reduce it for a smaller machine;
it does not mean eight cohorts or eight repetitions. Runtime depends on the data and hardware.
Without `--resume`, running the command again creates a new timestamped run directory.

## 5. Understand the analysis steps

For the `buffa2010` preset, SEARCH performs these steps:

1. **Measure association within IT.** Spearman correlation asks whether a gene and a seed
   tend to have higher expression in the same samples. A significant association forms an
   edge after a Bonferroni correction over the gene universe.
2. **Calculate connectivity.** Connections are weighted by the size of each seed's cluster.
   SEARCH ranks the resulting gene connectivity scores within the cohort.
3. **Estimate stability.** Bootstrap resamples sample columns with replacement, recomputes
   ranks and reports their mean and variance. It checks sensitivity to the sampled columns.
4. **Build a random comparison.** Monte Carlo repeats the analysis with random seed sets of
   the same size. The preset pools null connectivity values across genes and repetitions
   within that cohort. It does not pool IT and FR samples or their null distributions.
5. **Repeat independently in FR.** The same method uses the French samples only.
6. **Select the common signature.** A gene must satisfy every default condition below.
7. **Score samples.** For each sample, take the median of the absolute expression values
   of the available signature genes, then rank the sample scores within its cohort.

Default final selection conditions:

| Condition | Meaning |
|---|---|
| Not a seed | The ten reference genes are excluded from the final candidate list |
| Positive consistency in every training cohort | At least one significant edge, with no significant negative edges, in each cohort |
| Monte Carlo p-value < 0.01 in every training cohort | IT and FR must each pass their own statistical threshold |
| `pC_boot` > 0.5 | The product of the mean bootstrap connectivity ranks must exceed 0.5 |

`pC_boot` is a product of ranked connectivity scores, **not a p-value**. With two cohorts,
mean ranks of 0.90 and 0.80 give a product of 0.72. This condition is separate from the
Monte Carlo p-value test.

Passing genes are ordered by `pC_boot`, then `pC_obs`, then `C_IVW` (a mean weighted by
inverse bootstrap variance). The top signature contains **at most 51** passing genes.

## 6. Find and interpret the results

At completion, the terminal prints the run directory and candidate gene count. Retrieve
that directory as follows:

```bash
cat "$SEARCH_FULL_OUT/latest/run_path.txt"
SEARCH_RUN="$(cat "$SEARCH_FULL_OUT/latest/run_path.txt")"
```

`latest` contains a pointer to the most recently started run. It is not the directory
holding the result tables and does not certify successful completion. Check `manifest.json`
inside the real run directory: a completed run has `"status": "finished"`.

The most useful files are:

| File inside the run directory | What to use it for |
|---|---|
| `results/IT_gene_connectivity.csv` | IT associations, connectivity, bootstrap summaries and Monte Carlo p-values |
| `results/FR_gene_connectivity.csv` | The same evidence calculated independently in FR |
| `results/combined_full_ranking.csv` | All analysed genes, with combined values and selection flags |
| `results/combined_signature.csv` | All genes that pass the final selection conditions |
| `results/combined_signature_top.csv` | The first passing genes, capped at the configured `top_n` |
| `results/scores_signature_all.csv` | Sample scores using all selected genes |
| `results/scores_signature_top.csv` | Sample scores using the capped signature |
| `config.yaml` | The resolved analysis settings |
| `manifest.json` | Completion status, software versions, input hashes, warnings and timings |
| `log.txt` | Progress messages and analysis warnings |
| `checkpoints/IT/` and `checkpoints/FR/` | Saved repetition chunks used for resuming |

Start with `combined_signature.csv` for the gene list. To understand an exclusion, inspect
that gene in `combined_full_ranking.csv`: check `is_seed`, the per-cohort `p_MC_C_*` columns,
`positive_consistent_all`, `product_pass` and `candidate`.

The archived local full analysis selected **26 genes**. This is an observed result for those
inputs, not a fixed target. With 26 passing genes, both signature files contain 26 rows;
the top-51 cap does not create 25 additional genes.

### Reading sample scores

Each score row identifies a cohort and sample. Important columns are:

- `E_median_abs_expression`: the expression summary before ranking.
- `HS_rank01`: the relative position of that summary within the cohort, on a 0-1 scale.
- `n_signature_genes_present`: how many signature genes were available in that cohort.
- `n_signature_genes_total`: the number of distinct genes in the supplied signature.

A score near 1 places a sample high relative to the other samples in its own cohort. It is
not a measured oxygen concentration. Because ranks are calculated separately, an IT score
of 0.9 and a TCGA score of 0.9 express similar relative positions, not equal absolute
expression levels. A new set of samples can change those relative ranks.

By default, at least three signature genes must be present for scoring. Check coverage even
when scoring succeeds: using different subsets of a signature can affect comparability.
The derivation run skips scores when fewer than three genes are selected.

For figures, with the `plots` dependency installed:

```bash
"$SEARCH_BIN" plot "$SEARCH_RUN"
```

The figures are written under the run's `plots/` directory.

## 7. Resume an interrupted run

For a run created using the full-analysis output directory above:

```bash
"$SEARCH_BIN" resume "$SEARCH_FULL_OUT/latest" --n-jobs 8
```

SEARCH reopens the saved configuration and the recorded input manifest. It reuses checkpoint
chunks whose metadata still match the inputs, configuration and repetition settings, and
recomputes missing or incompatible chunks. Keep the input manifest and matrices accessible.

Use a particular timestamped run directory instead of `latest` if a newer run has since
changed the pointer. A copied or moved run can retain old paths in `manifest.json` and
`run_path.txt`; update the paths deliberately before attempting to resume it.

The archived `full_run_1000boot_5000mc` copy records an old temporary manifest path that is
no longer available. For a straightforward new analysis, use `cohorts_used.yaml` and a new
output directory, as in section 4, rather than resuming that archive directly.

## 8. Reuse results or score another cohort

### Change final selection without repeating bootstrap and Monte Carlo

For example, request up to 20 genes under the same selection thresholds:

```bash
"$SEARCH_BIN" combine "$SEARCH_RUN" \
  --top-n 20 \
  --prefix top20
```

This writes `top20_full_ranking.csv`, `top20_signature.csv` and `top20_signature_top.csv`
under `results/`. It reuses the per-cohort gene tables. The full selected list remains
uncapped; only the `_signature_top.csv` file is capped at 20.

Other supported options include `--mc-p-cutoff`, `--product-cutoff` and `--keep-seeds`.
Changing them changes the final selection rule. `combine` does not regenerate sample scores
or change the stored configuration of the original run. Use `score` with the new signature
to calculate corresponding scores, and retain the recombination command with your analysis
record. A repeated prefix replaces the corresponding recombination files.

```bash
"$SEARCH_BIN" score \
  --manifest "$SEARCH_MANIFEST" \
  --signature "$SEARCH_RUN/results/top20_signature_top.csv" \
  --out "$SEARCH_RUN/results/top20_scores.csv"
```

### Apply a fixed signature to another cohort

Make a copy of your cohort manifest and add the new cohort with `role: scoring` and its matrix
path. Keep the existing valid manifest structure and training entries. Then use `search score`
with that copied manifest and the saved signature CSV, following the command above.

`score` computes scores for **all cohorts listed in the manifest**, regardless of role.
It does not rerun correlations, bootstrap, Monte Carlo or signature selection. Adding the
cohort as `scoring` and using `score` therefore preserves the chosen gene list. The input
loading and gene filters in the manifest still apply.

## 9. Understand the R script and report differences

The current GBM preset follows the R script's selection and scoring choices. The report
describes a different combination and compaction procedure:

| Choice | Current `buffa2010` preset / GBM R selection | Procedure described in the report |
|---|---|---|
| Analyse construction cohorts | IT and FR separately | IT and FR separately |
| Final Monte Carlo evidence | p < 0.01 separately in every training cohort | Fisher-combined p < 0.05 |
| Allow seed genes in the selected set | Excluded | Eligible when the other criteria pass |
| Expression summary for sample scores | Median of absolute expression | Mean expression |
| Choose compact size | Fixed cap, default 51 | Smallest size reaching Spearman correlation >= 0.95 with the full score in both construction cohorts |
| Observed gene counts in the compared outputs | 26 selected genes in the archived run | 439 full genes, reduced to 93 |

Fisher combines the two p-values, while leaving the original samples in separate cohorts.
In the report, IT reaches the score-stability threshold at 59 genes and FR at 93, so FR
determines the compact size. TCGA does not determine that size.

**There is currently no Fisher-combination option or automatic score-stability compaction
to 93 genes in SEARCH.** Setting `combine.top_n=93` only changes a cap. Keeping seeds or
changing the per-cohort p-value threshold also does not implement Fisher. Reproducing the
report requires those method changes plus matching input data and gene filtering.

The R files under `reference/` are reference material; the Python command does not execute
them. In the supplied consensus R variant, the expected annotation column
`in_tcga_rowdata` differs from the available `in_tcga_rowranges`, which explains the archived
consensus-run error. Python applies any requested annotation restriction through the YAML
`annotation_query` instead.

### Which preset should you choose?

Use `buffa2010` for the multi-cohort GBM workflow in this guide. `search2014` configures the
original SEARCH growth procedure: it updates the seed network using recruitment and pruning
rules. The current preset uses one growth iteration by default and has cohort combination
and sample scoring disabled. It therefore exports cohort-level network results rather than
the combined signature and scores described above. Both presets currently use an internal
random-seed null model.

Inspect the actual settings before changing a method:

```bash
"$SEARCH_BIN" config show --preset buffa2010
"$SEARCH_BIN" config show --preset search2014
```

Method overrides can produce a configuration that matches neither tested preset; SEARCH
records this as `preset_validated: false`. This field describes a match to the tested
configuration, not biological or clinical validation of the resulting signature.

## 10. Troubleshooting

| Symptom | What to check or do |
|---|---|
| The `search` command is not found | Use the full `SEARCH_BIN` path from section 4, or install the project environment |
| Matrix or seed file not found | Check spelling and manifest-relative paths; the existing example manifest needs its data paths edited |
| A seed is missing or too few are usable | Check ID type, version suffixes, gene filters and whether the seed varies across samples |
| Duplicate gene identifiers | Resolve duplicates upstream, including duplicates created by stripping version suffixes |
| Non-numeric matrix values | Remove annotation text from expression columns; check the delimiter and orientation |
| Sample table lacks samples | Match its first-column IDs to the matrix sample names |
| `annotation_query` fails | Use column names actually present in the annotation table; the local consensus table uses `in_tcga_rowranges` |
| `--strict` refuses the run | Read the recorded warnings and configuration mismatch; the current IT inputs have name-based warnings |
| No genes are selected | Inspect the separate selection flags in the full ranking and check inputs and roles before altering thresholds |
| Fewer than 51 genes in the top file | The file contains at most 51 passing genes, and can contain fewer |
| No sample score files | Check whether at least three genes were selected; also check whether scoring is enabled for the preset |
| `resume` cannot find the manifest | Inspect recorded paths, especially after moving or copying a run |
| `plot` reports a missing library | Install the optional `plots` dependencies, for example with `uv sync --all-extras` |
| Results differ from an earlier run | Compare input hashes, roles, filters, seeds, resolved settings, repetition counts and software versions |

For different runs, keep the timestamped run directories and their `manifest.json` files.
The `compare` command can report ranking agreement and gene-list overlap between two runs;
see the CLI reference for its options. A matching gene count alone does not show that the
gene identities or the methods match.

## 11. Run on another computer or an HPC cluster

### Fresh local installation

The project requires Python 3.12 or newer. From your repository checkout, after installing
`uv`, set `SEARCH_PROJECT` to that checkout's absolute path and run:

```bash
cd "$SEARCH_PROJECT"
uv sync --all-extras
"$SEARCH_PROJECT/.venv/bin/search" --help
```

Update the manifest to point to files on that computer. Absolute paths starting with
`/Users/gennaroabbruzzese/` belong to the current Mac and will not locate data on a cluster.
Relative paths make the data layout easier to move when the manifest and its data folder
travel together.

### SLURM cluster with Conda

The repository provides `environment.yml` and `slurm/run_search.sbatch`. On the cluster,
set `SEARCH_PROJECT` to the repository checkout, then install from that directory:

```bash
cd "$SEARCH_PROJECT"
conda env create --prefix "$SEARCH_PROJECT/envs/search" --file "$SEARCH_PROJECT/environment.yml"
```

Prepare a cluster-specific manifest, for example at `$SEARCH_PROJECT/cohorts.yaml`, with
paths that exist on the cluster. Inspect the SLURM template's partition and resource requests
and adapt them to the cluster. Create the log directory before submission, since SLURM opens
the log files before the job starts:

```bash
mkdir -p "$SEARCH_PROJECT/logs"
sbatch \
  --export="ALL,BASE_DIR=$SEARCH_PROJECT,ENV_PREFIX=$SEARCH_PROJECT/envs/search,MANIFEST=$SEARCH_PROJECT/cohorts.yaml,OUT=$SEARCH_PROJECT/runs" \
  "$SEARCH_PROJECT/slurm/run_search.sbatch"
```

The template uses the allocated CPUs for worker processes and passes `--resume` to reuse
matching checkpoints of the latest run under its output directory. Choose a new output
directory for a deliberately separate analysis.

## 12. Terms and further reading

| Term | Meaning in this workflow |
|---|---|
| Cohort | A group of samples analysed together |
| Expression | The numeric RNA-expression value supplied for a gene in a sample |
| Seed | A reference gene associated with the process being investigated |
| Co-expression | Genes tending to vary together across samples |
| Connectivity | A summary of a gene's connections to the seeds |
| Bootstrap | Repeatedly resampling sample columns to assess stability |
| Null model | A random comparison describing connectivity obtained with random seeds |
| Monte Carlo p-value | How often null values reach or exceed the observed connectivity, with a pseudocount; it is not the probability that a gene truly belongs to the biological process |
| Signature | The selected gene list used to calculate sample scores |
| Manifest | The YAML file listing input files, cohort roles and settings |
| Preset | A named bundle of method settings |
| Checkpoint | Saved repetition results that can be reused when resuming |

For additional commands and mathematical details, read the [CLI reference](cli.md),
[preset definitions](presets.md), [method](methods.md), [decision log](decisions.md) and
[validation notes](validation.md).

The implementation behind this guide is in `src/search/run.py`, `src/search/combine.py`,
`src/search/signature.py` and `src/search/cli.py`. The local GBM manifest is under
`GBM_HYPOXIA_SIGNATURE/RESULTS/SEARCH_PY/cohorts_used.yaml`.
