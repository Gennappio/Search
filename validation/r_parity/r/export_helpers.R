# Parity harness, R side: helpers shared by the exports of the GBM script.

suppressPackageStartupMessages({
  library(data.table)
  library(matrixStats)
})

# The script's own functions, taken from its file: only top-level `name <- function(...)`
# expressions are evaluated, so nothing of the script runs.
extract_functions <- function(script) {
  env <- new.env(parent = globalenv())
  for (e in parse(script, keep.source = FALSE)) {
    if (is.call(e) && identical(e[[1]], as.name("<-")) && is.call(e[[3]]) &&
      identical(e[[3]][[1]], as.name("function"))) {
      eval(e, envir = env)
    }
  }
  env
}

# 17 significant digits: a double printed this way reads back as the same double.
num <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))

write_tsv <- function(tab, path) {
  tab <- as.data.frame(tab, stringsAsFactors = FALSE)
  for (col in names(tab)) {
    if (is.double(tab[[col]])) tab[[col]] <- num(tab[[col]])
    if (is.logical(tab[[col]])) tab[[col]] <- ifelse(is.na(tab[[col]]), "NA", ifelse(tab[[col]], "TRUE", "FALSE"))
  }
  write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_gene_matrix <- function(m, path) {
  tab <- data.frame(gene = rownames(m), matrix(num(as.numeric(m)), nrow = nrow(m)), stringsAsFactors = FALSE)
  colnames(tab) <- c("gene", colnames(m))
  write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

# The random draws of a run. The script does not save them: they are reproduced with the calls
# and the seeds it uses, set.seed(1000000 + b) and sample.int for the bootstrap,
# set.seed(2000000 + i) and sample for the random seed sets. The comparison checks the
# reproduction: Python is handed these draws and matches the script replicate by replicate only
# if they are the ones the script made.
export_draws <- function(x, seeds, ds, b_boot, n_mc, safe_row_sds, export) {
  genes <- rownames(x)
  n_samples <- ncol(x)
  k <- length(seeds)
  boot <- t(vapply(seq_len(b_boot), function(b) {
    set.seed(1000000L + b)
    sample.int(n_samples, size = n_samples, replace = TRUE)
  }, integer(n_samples)))
  write.table(boot, file.path(export, paste0(ds, "_bootstrap_indices_1based.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  eligible <- genes[safe_row_sds(x) > 0]
  eligible <- setdiff(eligible, seeds)
  writeLines(eligible, file.path(export, paste0(ds, "_eligible_genes.txt")))
  sets <- t(vapply(seq_len(n_mc), function(i) {
    set.seed(2000000L + i)
    sample(eligible, size = k, replace = FALSE)
  }, character(k)))
  write.table(sets, file.path(export, paste0(ds, "_random_seed_sets.tsv")), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

export_core <- function(core, ds, export) {
  write_gene_matrix(core$rho, file.path(export, paste0(ds, "_core_rho.tsv")))
  write_gene_matrix(core$delta * 1, file.path(export, paste0(ds, "_core_delta.tsv")))
  write_gene_matrix(core$gamma, file.path(export, paste0(ds, "_core_gamma.tsv")))
  write_gene_matrix(core$S, file.path(export, paste0(ds, "_S.tsv")))
  write_gene_matrix(core$F, file.path(export, paste0(ds, "_F.tsv")))
  write_tsv(
    data.frame(seed = names(core$cluster_size), cluster_size = as.integer(core$cluster_size)),
    file.path(export, paste0(ds, "_cluster_sizes.tsv"))
  )
}

# Everything run_dataset() returns or leaves in its checkpoints, for one dataset.
export_dataset <- function(result, ds, export, ckpt_dir) {
  write_tsv(result$table, file.path(export, paste0(ds, "_gene_connectivity.tsv")))
  export_core(result$core, ds, export)
  write_tsv(result$mc$random_seed_summary, file.path(export, paste0(ds, "_mc_random_seed_summary.tsv")))
  ck <- file.path(ckpt_dir, ds)
  boot_files <- sort(list.files(file.path(ck, "bootstrap"), pattern = "^bootstrap_chunk_[0-9]+\\.rds$", full.names = TRUE))
  write_gene_matrix(do.call(cbind, lapply(boot_files, readRDS)), file.path(export, paste0(ds, "_bootstrap_C_rank.tsv")))
  mc_files <- sort(list.files(file.path(ck, "monte_carlo"), pattern = "^mc_C_chunk_[0-9]+\\.rds$", full.names = TRUE))
  write_gene_matrix(do.call(cbind, lapply(mc_files, readRDS)), file.path(export, paste0(ds, "_mc_C_raw.tsv")))
}

write_run_info <- function(export, values) {
  values <- c(values, list(
    r_version = R.version.string,
    data.table = as.character(packageVersion("data.table")),
    matrixStats = as.character(packageVersion("matrixStats"))
  ))
  writeLines(paste0(names(values), "\t", vapply(values, as.character, character(1))), file.path(export, "run_info.tsv"))
}
