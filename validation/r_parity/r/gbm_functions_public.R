#!/usr/bin/env Rscript
# Parity harness, R side: the functions of the lab's GBM script on public data.
#
# The script as a whole runs only on the GBM cohorts (file names and seeds are written in it),
# and those are private. Its per-dataset functions are general: run_dataset() and what it calls
# (buffa2010_core_one_dataset, bootstrap_checkpointed, monte_carlo_checkpointed) are taken from
# the script file unmodified and run on Breastdata1, the matrix bundled with the SEARCH package.
# The cross-cohort combination and the score are top-level code of the script and are covered
# only by the GBM comparison.
#
# Two runs: the complete matrix (core, bootstrap and null model), and the matrix with one
# missing value in about thirty probesets (core only), where the two implementations differ
# by design.
#
# Usage: Rscript gbm_functions_public.R <repository root> <output directory> [b_boot] [n_mc]

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(args[1])
out <- args[2]
b_boot <- if (length(args) > 2) as.integer(args[3]) else 10L
n_mc <- if (length(args) > 3) as.integer(args[4]) else 20L
export <- file.path(out, "r_export")
unlink(out, recursive = TRUE)
dir.create(export, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo, "validation", "r_parity", "r", "export_helpers.R"))

fun <- extract_functions(file.path(repo, "reference", "r_sources", "gbm", "run_buffa2010_hpc_wo_consensus.R"))
# the two globals of the script that its functions use: a logger and the csv folder
fun$stamp <- function(...) invisible(NULL)
fun$csv_dir <- file.path(out, "r_run", "csv")
for (d in c("csv", "rds", "checkpoints")) dir.create(file.path(out, "r_run", d), recursive = TRUE, showWarnings = FALSE)

x <- as.matrix(read.table(
  file.path(repo, "tests", "data", "breast", "Breastdata1", "mRNA_abundance.txt.gz"),
  row.names = 1, header = TRUE, sep = "\t"
))
storage.mode(x) <- "double"
seeds <- as.character(read.csv(file.path(repo, "examples", "breast", "seeds_breast_cluster.csv"), stringsAsFactors = FALSE)$gene)
ds <- "Breast1"
alpha <- 0.05

result <- fun$run_dataset(
  expr_mat = x, seeds = seeds, dataset_name = ds, alpha = alpha, B = b_boot, n_mc = n_mc,
  cores = 1L, boot_chunk_size = 5L, mc_chunk_size = 5L,
  ckpt_dir = file.path(out, "r_run", "checkpoints"), rds_dir = file.path(out, "r_run", "rds"),
  csv_dir = fun$csv_dir, resume = FALSE
)
write_run_info(export, list(b_boot = b_boot, n_mc = n_mc, alpha = alpha, dataset = ds))
writeLines(seeds, file.path(export, "seeds.txt"))
export_draws(x, seeds, ds, b_boot, n_mc, fun$safe_row_sds, export)
export_dataset(result, ds, export, file.path(out, "r_run", "checkpoints"))

# one missing value (third sample) in every 50th probeset from the 25th and in the five
# probesets most connected to the seeds, so that some of them sit inside the seed clusters
most_connected <- head(result$table[is_seed == FALSE][order(-C_raw)]$gene, 5)
na_rows <- setdiff(union(rownames(x)[seq(25, nrow(x), by = 50)], most_connected), seeds)
x_na <- x
x_na[na_rows, 3] <- NA
writeLines(na_rows, file.path(export, "na_rows.txt"))
core_na <- fun$buffa2010_core_one_dataset(x_na, seeds = seeds, alpha = alpha)
export_core(core_na, "Breast1NA", export)
write_tsv(core_na$gene_table, file.path(export, "Breast1NA_gene_connectivity.tsv"))
cat("exported to", export, "\n")
