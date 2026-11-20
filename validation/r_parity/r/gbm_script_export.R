#!/usr/bin/env Rscript
# Parity harness, R side: export of a finished run of the lab's GBM script.
#
# Every table of the run is written from its .rds files with 17 significant digits, so the
# comparison is not limited by the 15 digits of the script's csv files, and the random draws of
# the run are written next to them (see export_helpers.R).
#
# Usage: Rscript gbm_script_export.R <repository root> <staging folder> <export folder>

args <- commandArgs(trailingOnly = TRUE)
repo <- normalizePath(args[1])
stage <- normalizePath(args[2])
export <- args[3]
dir.create(export, recursive = TRUE, showWarnings = FALSE)
source(file.path(repo, "validation", "r_parity", "r", "export_helpers.R"))

script <- file.path(repo, "reference", "r_sources", "gbm", "run_buffa2010_hpc_wo_consensus.R")
fun <- extract_functions(script)
run_dir <- readLines(file.path(stage, "RESULTS", "BUFFA2010_EXACT_HPC", "latest", "run_path.txt"))[1]
rds <- function(name) readRDS(file.path(run_dir, "rds", name))

params <- rds("run_parameters.rds")
write_run_info(export, list(
  run_dir = run_dir, b_boot = params$B_BOOT, n_mc = params$N_MC, alpha = params$alpha,
  mc_p_cutoff = params$mc_p_cutoff, product_cutoff = params$product_cutoff, seed = params$seed
))
seeds <- rds("Buffa2010_seed_set.rds")$gene

# the matrices as the script reads and aligns them, needed for the draws
data_dir <- file.path(stage, "DATA", "NORMALIZED")
expr <- list(
  IT = fun$read_expr_csv(file.path(data_dir, "IT_GBM_RNA_Deseq2_vst_Combat.csv")),
  FR = fun$read_expr_csv(file.path(data_dir, "FR_GBM_RNA_Deseq2_vst_Combat.csv"))
)
common <- Reduce(intersect, lapply(expr, rownames))
expr <- lapply(expr, function(x) x[common, , drop = FALSE])

for (ds in names(expr)) {
  export_draws(expr[[ds]], seeds, ds, params$B_BOOT, params$N_MC, fun$safe_row_sds, export)
  export_dataset(rds(paste0("Buffa2010_", ds, "_result.rds")), ds, export, file.path(run_dir, "checkpoints"))
}

write_tsv(rds("Buffa2010_exact_ITxFR_full_ranking.rds"), file.path(export, "combined_full_ranking.tsv"))
write_tsv(rds("GBM_Hypoxia_Buffa2010_exact_ITxFR_signature.rds"), file.path(export, "combined_signature.tsv"))
write_tsv(rds("GBM_Hypoxia_Buffa2010_exact_ITxFR_signature_top51.rds"), file.path(export, "combined_signature_top51.tsv"))
for (suffix in c("all_signature", "top51")) {
  path <- file.path(run_dir, "rds", paste0("Buffa2010_exact_ITxFR_HS_", suffix, "_IT_FR_TCGA.rds"))
  if (file.exists(path)) write_tsv(readRDS(path), file.path(export, paste0("scores_", suffix, ".tsv")))
}
cat("exported to", export, "\n")
