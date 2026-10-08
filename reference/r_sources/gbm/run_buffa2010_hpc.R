#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(matrixStats)
  library(ggplot2)
})

option_list <- list(
  make_option("--base_dir", type = "character", default = "/home/MoscatoG/GBM_HYPOXIA_SIGNATURE"),
  make_option("--annotation_file", type = "character", default = NA_character_),
  make_option("--b_boot", type = "integer", default = 1000),
  make_option("--n_mc", type = "integer", default = 5000),
  make_option("--cores", type = "integer", default = NA_integer_),
  make_option("--alpha", type = "double", default = 0.05),
  make_option("--mc_p_cutoff", type = "double", default = 0.01),
  make_option("--product_cutoff", type = "double", default = 0.5),
  make_option("--boot_chunk_size", type = "integer", default = 100),
  make_option("--mc_chunk_size", type = "integer", default = 100),
  make_option("--seed", type = "integer", default = 1234),
  make_option("--resume", type = "character", default = "FALSE"),
  make_option("--test_mode", type = "character", default = "FALSE")
)
opt <- parse_args(OptionParser(option_list = option_list))

str_to_logical <- function(x) {
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
}

resume <- str_to_logical(opt$resume)
test_mode <- str_to_logical(opt$test_mode)

if (test_mode) {
  opt$b_boot <- min(opt$b_boot, 20L)
  opt$n_mc <- min(opt$n_mc, 50L)
  opt$boot_chunk_size <- min(opt$boot_chunk_size, 10L)
  opt$mc_chunk_size <- min(opt$mc_chunk_size, 10L)
}

if (is.na(opt$cores)) {
  slurm_cores <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
  if (!is.na(slurm_cores)) {
    opt$cores <- as.integer(slurm_cores)
  } else {
    opt$cores <- max(1L, parallel::detectCores(logical = FALSE) - 1L)
  }
}
opt$cores <- max(1L, opt$cores)

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

set.seed(opt$seed)

base_dir <- opt$base_dir
data_dir <- file.path(base_dir, "DATA", "NORMALIZED")
ann_file <- if (is.na(opt$annotation_file) || opt$annotation_file == "") {
  file.path(base_dir, "DATA", "master_annotation_ITFR_Ensembl104_TCGA_GDC_consensus.csv")
} else {
  opt$annotation_file
}
code_dir <- file.path(base_dir, "CODE")
run_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
run_dir <- file.path(base_dir, "RESULTS", "BUFFA2010_EXACT_HPC", paste0("run_", run_tag))
latest_dir <- file.path(base_dir, "RESULTS", "BUFFA2010_EXACT_HPC", "latest")

rds_dir <- file.path(run_dir, "rds")
csv_dir <- file.path(run_dir, "csv")
plot_dir <- file.path(run_dir, "plots")
log_dir <- file.path(run_dir, "logs")
ckpt_dir <- file.path(run_dir, "checkpoints")

for (d in c(run_dir, rds_dir, csv_dir, plot_dir, log_dir, ckpt_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# If resuming from latest, use latest target if it exists and has a run_path.txt file.
if (resume && file.exists(file.path(latest_dir, "run_path.txt"))) {
  previous_run <- readLines(file.path(latest_dir, "run_path.txt"), warn = FALSE)[1]
  if (!is.na(previous_run) && dir.exists(previous_run)) {
    message("Resume requested: using previous run directory: ", previous_run)
    run_dir <- previous_run
    rds_dir <- file.path(run_dir, "rds")
    csv_dir <- file.path(run_dir, "csv")
    plot_dir <- file.path(run_dir, "plots")
    log_dir <- file.path(run_dir, "logs")
    ckpt_dir <- file.path(run_dir, "checkpoints")
  }
}

dir.create(latest_dir, recursive = TRUE, showWarnings = FALSE)
writeLines(run_dir, file.path(latest_dir, "run_path.txt"))

log_file <- file.path(log_dir, "run_buffa2010_hpc.log")
log_message_con <- file(log_file, open = "at")
sink(log_file, split = TRUE)
sink(log_message_con, type = "message")

on.exit({
  try(sink(type = "message"), silent = TRUE)
  try(close(log_message_con), silent = TRUE)
  try(sink(), silent = TRUE)
}, add = TRUE)

stamp <- function(...) {
  cat(format(Sys.time(), "[%Y-%m-%d %H:%M:%S]"), ..., "\n")
  flush.console()
}

stamp("Starting Buffa 2010 exact HPC run")
stamp("Base dir:", base_dir)
stamp("Run dir:", run_dir)
stamp("B_BOOT:", opt$b_boot, "N_MC:", opt$n_mc, "CORES:", opt$cores)
stamp("ALPHA:", opt$alpha, "MC cutoff:", opt$mc_p_cutoff, "Product cutoff:", opt$product_cutoff)

params <- list(
  base_dir = base_dir,
  data_dir = data_dir,
  ann_file = ann_file,
  run_dir = run_dir,
  B_BOOT = opt$b_boot,
  N_MC = opt$n_mc,
  cores = opt$cores,
  alpha = opt$alpha,
  mc_p_cutoff = opt$mc_p_cutoff,
  product_cutoff = opt$product_cutoff,
  boot_chunk_size = opt$boot_chunk_size,
  mc_chunk_size = opt$mc_chunk_size,
  seed = opt$seed,
  resume = resume,
  test_mode = test_mode,
  start_time = as.character(Sys.time())
)
saveRDS(params, file.path(rds_dir, "run_parameters.rds"))
fwrite(as.data.table(params), file.path(csv_dir, "run_parameters.csv"))

# -----------------------------
# Read input data
# -----------------------------
read_expr_csv <- function(path) {
  if (!file.exists(path)) stop("Missing expression file: ", path)
  dt <- fread(path)
  gene <- dt[[1]]
  dt[[1]] <- NULL
  mat <- as.matrix(dt)
  storage.mode(mat) <- "double"
  gene <- sub("\\..*$", "", gene)
  if (anyDuplicated(gene)) {
    dup <- unique(gene[duplicated(gene)])
    stop("Duplicated stable Ensembl IDs after stripping version in ", path, ": ", paste(head(dup, 20), collapse = ", "))
  }
  rownames(mat) <- gene
  mat
}

IT_file <- file.path(data_dir, "IT_GBM_RNA_Deseq2_vst_Combat.csv")
FR_file <- file.path(data_dir, "FR_GBM_RNA_Deseq2_vst_Combat.csv")
TCGA_file <- file.path(data_dir, "TCGA_GBM_RNA_Deseq2_vst.csv")

stamp("Reading expression matrices")
IT_GBM_RNA <- read_expr_csv(IT_file)
FR_GBM_RNA <- read_expr_csv(FR_file)
TCGA_GBM_RNA <- read_expr_csv(TCGA_file)

stamp("IT dim:", paste(dim(IT_GBM_RNA), collapse = " x "))
stamp("FR dim:", paste(dim(FR_GBM_RNA), collapse = " x "))
stamp("TCGA dim:", paste(dim(TCGA_GBM_RNA), collapse = " x "))

if (!file.exists(ann_file)) stop("Missing consensus annotation file: ", ann_file)
consensus_ann <- fread(ann_file)
required_ann_cols <- c("ensg", "in_itfr_gtf104", "in_tcga_rowdata", "symbol_consensus", "symbol_itfr", "symbol_tcga", "annotation_status", "gene_biotype_consensus")
missing_ann_cols <- setdiff(required_ann_cols, colnames(consensus_ann))
if (length(missing_ann_cols) > 0) {
  stop("Consensus annotation file is missing required columns: ", paste(missing_ann_cols, collapse = ", "))
}

as_logical_safe <- function(x) {
  if (is.logical(x)) return(x)
  toupper(as.character(x)) %in% c("TRUE", "T", "1", "YES", "Y")
}

consensus_ann[, in_itfr_gtf104_bool := as_logical_safe(in_itfr_gtf104)]
consensus_ann[, in_tcga_rowdata_bool := as_logical_safe(in_tcga_rowdata)]
consensus_primary <- consensus_ann[in_itfr_gtf104_bool == TRUE & in_tcga_rowdata_bool == TRUE]
consensus_primary <- unique(consensus_primary, by = "ensg")
primary_consensus_genes <- consensus_primary$ensg

ann <- unique(consensus_primary[, .(
  gene = ensg,
  symbol = symbol_consensus,
  symbol_itfr = symbol_itfr,
  symbol_tcga = symbol_tcga,
  annotation_status = annotation_status,
  gene_biotype_consensus = gene_biotype_consensus
)], by = "gene")

stamp("Consensus annotation file:", ann_file)
stamp("Primary consensus genes:", length(primary_consensus_genes))

filter_expr_to_consensus <- function(mat, dataset_name, primary_genes) {
  before <- nrow(mat)
  keep <- intersect(rownames(mat), primary_genes)
  removed <- setdiff(rownames(mat), primary_genes)
  stamp(dataset_name, "genes before consensus safety filter:", before)
  stamp(dataset_name, "genes retained by consensus safety filter:", length(keep))
  stamp(dataset_name, "genes removed by consensus safety filter:", length(removed))
  if (length(keep) == 0) stop(dataset_name, ": no genes left after consensus filter.")
  mat[keep, , drop = FALSE]
}

input_consensus_qc <- rbindlist(list(
  data.table(dataset = "IT", n_before = nrow(IT_GBM_RNA), n_in_consensus = sum(rownames(IT_GBM_RNA) %in% primary_consensus_genes), n_removed = sum(!rownames(IT_GBM_RNA) %in% primary_consensus_genes)),
  data.table(dataset = "FR", n_before = nrow(FR_GBM_RNA), n_in_consensus = sum(rownames(FR_GBM_RNA) %in% primary_consensus_genes), n_removed = sum(!rownames(FR_GBM_RNA) %in% primary_consensus_genes)),
  data.table(dataset = "TCGA", n_before = nrow(TCGA_GBM_RNA), n_in_consensus = sum(rownames(TCGA_GBM_RNA) %in% primary_consensus_genes), n_removed = sum(!rownames(TCGA_GBM_RNA) %in% primary_consensus_genes))
))
fwrite(input_consensus_qc, file.path(csv_dir, "Buffa2010_input_consensus_annotation_QC.csv"))

IT_GBM_RNA <- filter_expr_to_consensus(IT_GBM_RNA, "IT", primary_consensus_genes)
FR_GBM_RNA <- filter_expr_to_consensus(FR_GBM_RNA, "FR", primary_consensus_genes)
TCGA_GBM_RNA <- filter_expr_to_consensus(TCGA_GBM_RNA, "TCGA", primary_consensus_genes)

stamp("IT dim after consensus safety filter:", paste(dim(IT_GBM_RNA), collapse = " x "))
stamp("FR dim after consensus safety filter:", paste(dim(FR_GBM_RNA), collapse = " x "))
stamp("TCGA dim after consensus safety filter:", paste(dim(TCGA_GBM_RNA), collapse = " x "))

seed_primary_name <- "Buffa2010_setA_conservative"
seed_primary_symb <- c("ADM", "AK4", "BNIP3", "CA9", "ENO1", "HK2", "LDHA", "PGK1", "SLC2A1", "VEGFA")
seed_primary_ensg <- c(
  "ENSG00000148926",
  "ENSG00000162433",
  "ENSG00000176171",
  "ENSG00000107159",
  "ENSG00000074800",
  "ENSG00000159399",
  "ENSG00000134333",
  "ENSG00000102144",
  "ENSG00000117394",
  "ENSG00000112715"
)
seed_map <- data.table(gene = seed_primary_ensg, symbol = seed_primary_symb)
fwrite(seed_map, file.path(csv_dir, "Buffa2010_seed_set.csv"))
saveRDS(seed_map, file.path(rds_dir, "Buffa2010_seed_set.rds"))

expr_train <- list(IT = IT_GBM_RNA, FR = FR_GBM_RNA)
common_train_genes <- Reduce(intersect, lapply(expr_train, rownames))
stamp("Common genes in IT + FR:", length(common_train_genes))
expr_train <- lapply(expr_train, function(x) x[common_train_genes, , drop = FALSE])

seed_presence <- rbindlist(lapply(names(expr_train), function(ds) {
  data.table(
    dataset = ds,
    gene = seed_map$gene,
    symbol = seed_map$symbol,
    present = seed_map$gene %in% rownames(expr_train[[ds]])
  )
}))
fwrite(seed_presence, file.path(csv_dir, "Buffa2010_seed_presence_IT_FR.csv"))
print(seed_presence)
if (any(!seed_presence$present)) stop("One or more seed genes are missing from IT/FR training matrices.")

# -----------------------------
# Math/helper functions
# -----------------------------
rank01 <- function(x, higher_is_better = TRUE) {
  out <- rep(NA_real_, length(x))
  names(out) <- names(x)
  ok <- is.finite(x) & !is.na(x)
  if (sum(ok) <= 1) return(out)
  xx <- x[ok]
  if (!higher_is_better) xx <- -xx
  rr <- rank(xx, ties.method = "average")
  out[ok] <- (rr - 1) / (length(rr) - 1)
  out
}

spearman_p_approx <- function(rho, n) {
  rho <- pmax(pmin(rho, 0.999999), -0.999999)
  tstat <- rho * sqrt((n - 2) / pmax(1e-12, 1 - rho^2))
  2 * pt(-abs(tstat), df = n - 2)
}

safe_row_sds <- function(x) {
  out <- matrixStats::rowSds(as.matrix(x), na.rm = TRUE)
  out[!is.finite(out)] <- 0
  out
}

buffa2010_core_one_dataset <- function(expr_mat, seeds, alpha = 0.05) {
  expr_mat <- as.matrix(expr_mat)
  genes <- rownames(expr_mat)
  n_samples <- ncol(expr_mat)
  n_genes <- nrow(expr_mat)
  k_seeds <- length(seeds)
  seeds_present <- intersect(seeds, genes)
  if (length(seeds_present) != length(seeds)) {
    missing <- setdiff(seeds, genes)
    stop("Missing seeds: ", paste(missing, collapse = ", "))
  }

  rho <- cor(
    t(expr_mat),
    t(expr_mat[seeds_present, , drop = FALSE]),
    method = "spearman",
    use = "pairwise.complete.obs"
  )
  rho[!is.finite(rho)] <- NA_real_

  pval <- spearman_p_approx(rho, n = n_samples)
  bonf_threshold <- alpha / n_genes
  delta <- pval < bonf_threshold
  delta[is.na(delta)] <- FALSE

  # Membership is normalized within each seed-defined cluster.
  cluster_size <- colSums(delta, na.rm = TRUE)
  gamma <- matrix(0, nrow = nrow(delta), ncol = ncol(delta), dimnames = dimnames(delta))
  for (j in seq_len(ncol(delta))) {
    if (cluster_size[j] > 0) gamma[, j] <- as.numeric(delta[, j]) / cluster_size[j]
  }

  S <- matrix(NA_real_, nrow = k_seeds, ncol = k_seeds, dimnames = list(seeds_present, seeds_present))
  F <- matrix(NA_real_, nrow = k_seeds, ncol = k_seeds, dimnames = list(seeds_present, seeds_present))

  for (a in seq_len(k_seeds)) {
    for (b in seq_len(k_seeds)) {
      if (a == b) {
        S[a, b] <- 1
        F[a, b] <- 1
      } else {
        seed_a <- seeds_present[a]
        seed_b <- seeds_present[b]
        idx <- !(genes %in% c(seed_a, seed_b))
        ga <- gamma[idx, a]
        gb <- gamma[idx, b]
        min_g <- pmin(ga, gb)
        max_g <- pmax(ga, gb)
        denom_S <- sum(max_g, na.rm = TRUE)
        denom_F <- sum(min_g, na.rm = TRUE)
        S[a, b] <- ifelse(denom_S > 0, sum(min_g, na.rm = TRUE) / denom_S, NA_real_)
        sign_prod <- sign(rho[idx, a] * rho[idx, b])
        sign_prod[!is.finite(sign_prod)] <- NA_real_
        F[a, b] <- ifelse(denom_F > 0, sum(min_g * sign_prod, na.rm = TRUE) / denom_F, NA_real_)
      }
    }
  }

  weight_mat <- matrix(1, nrow = n_genes, ncol = k_seeds, dimnames = dimnames(gamma))
  for (j in seq_len(k_seeds)) weight_mat[genes == seeds_present[j], j] <- 0
  denom_w <- rowSums(weight_mat)
  C_raw <- rowSums(weight_mat * gamma, na.rm = TRUE) / denom_w
  C_raw[!is.finite(C_raw)] <- NA_real_
  names(C_raw) <- genes
  C_rank <- rank01(C_raw, higher_is_better = TRUE)

  n_edges <- rowSums(delta)
  n_positive_edges <- rowSums(delta & rho > 0, na.rm = TRUE)
  n_negative_edges <- rowSums(delta & rho < 0, na.rm = TRUE)
  positive_consistent <- n_edges > 0 & n_negative_edges == 0
  weighted_rho <- rowSums(gamma * rho, na.rm = TRUE) / rowSums(gamma)
  weighted_rho[!is.finite(weighted_rho)] <- NA_real_

  gene_table <- data.table(
    gene = genes,
    C_raw = as.numeric(C_raw),
    C_rank = as.numeric(C_rank),
    n_edges = as.integer(n_edges),
    n_positive_edges = as.integer(n_positive_edges),
    n_negative_edges = as.integer(n_negative_edges),
    positive_consistent = positive_consistent,
    weighted_seed_rho = as.numeric(weighted_rho),
    is_seed = genes %in% seeds_present
  )

  list(
    gene_table = gene_table,
    S = S,
    F = F,
    rho = rho,
    delta = delta,
    gamma = gamma,
    cluster_size = cluster_size,
    bonf_threshold = bonf_threshold,
    n_samples = n_samples,
    n_genes = n_genes,
    seeds = seeds_present
  )
}

chunk_indices <- function(n, chunk_size) {
  split(seq_len(n), ceiling(seq_len(n) / chunk_size))
}

parallel_lapply <- function(X, FUN, cores) {
  if (.Platform$OS.type == "windows" || cores <= 1L) {
    lapply(X, FUN)
  } else {
    parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = FALSE)
  }
}

bootstrap_checkpointed <- function(expr_mat, seeds, dataset_name, B, alpha, chunk_size, cores, ckpt_dir, resume = TRUE, seed_base = 1000000L) {
  expr_mat <- as.matrix(expr_mat)
  genes <- rownames(expr_mat)
  n_samples <- ncol(expr_mat)
  chunks <- chunk_indices(B, chunk_size)
  ds_ckpt_dir <- file.path(ckpt_dir, dataset_name, "bootstrap")
  dir.create(ds_ckpt_dir, recursive = TRUE, showWarnings = FALSE)

  for (ci in seq_along(chunks)) {
    b_idx <- chunks[[ci]]
    out_file <- file.path(ds_ckpt_dir, sprintf("bootstrap_chunk_%04d.rds", ci))
    if (resume && file.exists(out_file)) {
      stamp(dataset_name, "bootstrap chunk", ci, "already exists; skipping")
      next
    }
    stamp(dataset_name, "bootstrap chunk", ci, "/", length(chunks), "reps", min(b_idx), "-", max(b_idx))
    boot_list <- parallel_lapply(b_idx, function(b) {
      set.seed(seed_base + b)
      idx <- sample.int(n_samples, size = n_samples, replace = TRUE)
      expr_b <- expr_mat[, idx, drop = FALSE]
      core_b <- buffa2010_core_one_dataset(expr_b, seeds = seeds, alpha = alpha)
      out <- core_b$gene_table$C_rank
      names(out) <- core_b$gene_table$gene
      out[genes]
    }, cores = cores)
    rank_chunk <- do.call(cbind, boot_list)
    rownames(rank_chunk) <- genes
    colnames(rank_chunk) <- paste0("boot_", b_idx)
    saveRDS(rank_chunk, out_file, compress = FALSE)
    rm(boot_list, rank_chunk)
    gc()
  }

  chunk_files <- file.path(ds_ckpt_dir, sprintf("bootstrap_chunk_%04d.rds", seq_along(chunks)))
  if (!all(file.exists(chunk_files))) stop("Missing bootstrap checkpoint chunks for ", dataset_name)
  stamp(dataset_name, "combining bootstrap chunks")
  rank_mat <- do.call(cbind, lapply(chunk_files, readRDS))
  out <- data.table(
    gene = genes,
    C_rank_boot_mean = rowMeans(rank_mat, na.rm = TRUE),
    C_rank_boot_var = matrixStats::rowVars(rank_mat, na.rm = TRUE),
    C_rank_boot_sd = sqrt(matrixStats::rowVars(rank_mat, na.rm = TRUE))
  )
  saveRDS(out, file.path(ds_ckpt_dir, "bootstrap_summary.rds"), compress = TRUE)
  rm(rank_mat)
  gc()
  out
}

monte_carlo_checkpointed <- function(expr_mat, seeds, observed_C_raw, dataset_name, n_mc, alpha, chunk_size, cores, ckpt_dir, resume = TRUE, seed_base = 2000000L) {
  expr_mat <- as.matrix(expr_mat)
  genes <- rownames(expr_mat)
  k <- length(seeds)
  eligible <- genes[safe_row_sds(expr_mat) > 0]
  eligible <- setdiff(eligible, seeds)
  if (length(eligible) < k) stop("Not enough eligible random seed genes for ", dataset_name)

  chunks <- chunk_indices(n_mc, chunk_size)
  ds_ckpt_dir <- file.path(ckpt_dir, dataset_name, "monte_carlo")
  dir.create(ds_ckpt_dir, recursive = TRUE, showWarnings = FALSE)

  for (ci in seq_along(chunks)) {
    i_idx <- chunks[[ci]]
    c_file <- file.path(ds_ckpt_dir, sprintf("mc_C_chunk_%04d.rds", ci))
    s_file <- file.path(ds_ckpt_dir, sprintf("mc_S_summary_chunk_%04d.rds", ci))
    if (resume && file.exists(c_file) && file.exists(s_file)) {
      stamp(dataset_name, "MC chunk", ci, "already exists; skipping")
      next
    }
    stamp(dataset_name, "MC chunk", ci, "/", length(chunks), "reps", min(i_idx), "-", max(i_idx))
    mc_list <- parallel_lapply(i_idx, function(i) {
      set.seed(seed_base + i)
      random_seeds <- sample(eligible, size = k, replace = FALSE)
      core_i <- buffa2010_core_one_dataset(expr_mat, seeds = random_seeds, alpha = alpha)
      C_i <- core_i$gene_table$C_raw
      names(C_i) <- core_i$gene_table$gene
      C_i <- C_i[genes]
      C_i[!is.finite(C_i)] <- NA_real_
      S_i <- core_i$S
      list(
        C = C_i,
        summary = data.table(
          mc = i,
          mean_S = mean(S_i[upper.tri(S_i)], na.rm = TRUE),
          median_S = median(S_i[upper.tri(S_i)], na.rm = TRUE)
        )
      )
    }, cores = cores)
    C_mat <- do.call(cbind, lapply(mc_list, `[[`, "C"))
    rownames(C_mat) <- genes
    colnames(C_mat) <- paste0("mc_", i_idx)
    S_summary <- rbindlist(lapply(mc_list, `[[`, "summary"))
    saveRDS(C_mat, c_file, compress = FALSE)
    saveRDS(S_summary, s_file, compress = TRUE)
    rm(mc_list, C_mat, S_summary)
    gc()
  }

  c_files <- file.path(ds_ckpt_dir, sprintf("mc_C_chunk_%04d.rds", seq_along(chunks)))
  s_files <- file.path(ds_ckpt_dir, sprintf("mc_S_summary_chunk_%04d.rds", seq_along(chunks)))
  if (!all(file.exists(c_files))) stop("Missing MC C checkpoint chunks for ", dataset_name)
  if (!all(file.exists(s_files))) stop("Missing MC S checkpoint chunks for ", dataset_name)

  stamp(dataset_name, "combining MC null values and computing p-values")
  null_C_values <- unlist(lapply(c_files, function(f) as.vector(readRDS(f))), use.names = FALSE)
  null_C_values <- null_C_values[is.finite(null_C_values)]
  sorted_null <- sort(null_C_values)
  n_null <- length(sorted_null)
  obs <- observed_C_raw
  names(obs) <- names(observed_C_raw)
  p_mc <- rep(NA_real_, length(obs))
  ok <- is.finite(obs)
  # Count values >= observed using sorted null distribution and pseudocount.
  count_less <- findInterval(obs[ok] - 1e-12, sorted_null)
  count_ge <- n_null - count_less
  p_mc[ok] <- (1 + count_ge) / (1 + n_null)

  p_table <- data.table(gene = names(obs), p_MC_C = as.numeric(p_mc))
  random_seed_summary <- rbindlist(lapply(s_files, readRDS))
  saveRDS(p_table, file.path(ds_ckpt_dir, "mc_p_table.rds"), compress = TRUE)
  saveRDS(random_seed_summary, file.path(ds_ckpt_dir, "mc_random_seed_summary.rds"), compress = TRUE)
  fwrite(random_seed_summary, file.path(csv_dir, paste0("Buffa2010_", dataset_name, "_MC_random_seed_S_summary.csv")))
  rm(null_C_values, sorted_null)
  gc()
  list(p_table = p_table, random_seed_summary = random_seed_summary)
}

run_dataset <- function(expr_mat, seeds, dataset_name, alpha, B, n_mc, cores, boot_chunk_size, mc_chunk_size, ckpt_dir, rds_dir, csv_dir, resume = TRUE) {
  final_rds <- file.path(rds_dir, paste0("Buffa2010_", dataset_name, "_result.rds"))
  final_csv <- file.path(csv_dir, paste0("Buffa2010_", dataset_name, "_gene_connectivity.csv"))
  if (resume && file.exists(final_rds) && file.exists(final_csv)) {
    stamp(dataset_name, "final result exists; loading")
    return(readRDS(final_rds))
  }

  stamp("Running core for", dataset_name)
  core <- buffa2010_core_one_dataset(expr_mat, seeds = seeds, alpha = alpha)
  saveRDS(core, file.path(rds_dir, paste0("Buffa2010_", dataset_name, "_core.rds")), compress = TRUE)

  observed_C <- core$gene_table$C_raw
  names(observed_C) <- core$gene_table$gene

  stamp("Running bootstrap for", dataset_name)
  boot <- bootstrap_checkpointed(expr_mat, seeds, dataset_name, B, alpha, boot_chunk_size, cores, ckpt_dir, resume = resume)

  stamp("Running Monte Carlo for", dataset_name)
  mc <- monte_carlo_checkpointed(expr_mat, seeds, observed_C, dataset_name, n_mc, alpha, mc_chunk_size, cores, ckpt_dir, resume = resume)

  out_table <- merge(core$gene_table, boot, by = "gene", all.x = TRUE)
  out_table <- merge(out_table, mc$p_table, by = "gene", all.x = TRUE)
  out_table[, dataset := dataset_name]
  setcolorder(out_table, c("dataset", setdiff(names(out_table), "dataset")))

  S_long <- as.data.table(as.table(core$S))
  setnames(S_long, c("seed_x", "seed_y", "S"))
  S_long[, dataset := dataset_name]
  setcolorder(S_long, c("dataset", "seed_x", "seed_y", "S"))

  F_long <- as.data.table(as.table(core$F))
  setnames(F_long, c("seed_x", "seed_y", "F"))
  F_long[, dataset := dataset_name]
  setcolorder(F_long, c("dataset", "seed_x", "seed_y", "F"))

  fwrite(out_table, final_csv)
  fwrite(S_long, file.path(csv_dir, paste0("Buffa2010_", dataset_name, "_shared_neighbourhood_S.csv")))
  fwrite(F_long, file.path(csv_dir, paste0("Buffa2010_", dataset_name, "_sign_concordance_F.csv")))

  result <- list(table = out_table, S = core$S, F = core$F, core = core, boot = boot, mc = mc)
  saveRDS(result, final_rds, compress = TRUE)
  result
}

# -----------------------------
# Run IT and FR
# -----------------------------
res_IT <- run_dataset(
  expr_mat = expr_train$IT,
  seeds = seed_primary_ensg,
  dataset_name = "IT",
  alpha = opt$alpha,
  B = opt$b_boot,
  n_mc = opt$n_mc,
  cores = opt$cores,
  boot_chunk_size = opt$boot_chunk_size,
  mc_chunk_size = opt$mc_chunk_size,
  ckpt_dir = ckpt_dir,
  rds_dir = rds_dir,
  csv_dir = csv_dir,
  resume = resume
)

res_FR <- run_dataset(
  expr_mat = expr_train$FR,
  seeds = seed_primary_ensg,
  dataset_name = "FR",
  alpha = opt$alpha,
  B = opt$b_boot,
  n_mc = opt$n_mc,
  cores = opt$cores,
  boot_chunk_size = opt$boot_chunk_size,
  mc_chunk_size = opt$mc_chunk_size,
  ckpt_dir = ckpt_dir,
  rds_dir = rds_dir,
  csv_dir = csv_dir,
  resume = resume
)

# -----------------------------
# Combine IT and FR
# -----------------------------
stamp("Combining IT and FR")

tab_IT <- copy(res_IT$table)[, .(
  gene,
  C_raw_IT = C_raw,
  C_rank_IT = C_rank,
  C_rank_boot_mean_IT = C_rank_boot_mean,
  C_rank_boot_var_IT = C_rank_boot_var,
  C_rank_boot_sd_IT = C_rank_boot_sd,
  p_MC_C_IT = p_MC_C,
  positive_consistent_IT = positive_consistent,
  weighted_seed_rho_IT = weighted_seed_rho,
  is_seed
)]

tab_FR <- copy(res_FR$table)[, .(
  gene,
  C_raw_FR = C_raw,
  C_rank_FR = C_rank,
  C_rank_boot_mean_FR = C_rank_boot_mean,
  C_rank_boot_var_FR = C_rank_boot_var,
  C_rank_boot_sd_FR = C_rank_boot_sd,
  p_MC_C_FR = p_MC_C,
  positive_consistent_FR = positive_consistent,
  weighted_seed_rho_FR = weighted_seed_rho,
  is_seed_FR = is_seed
)]

buffa_ITFR <- merge(tab_IT, tab_FR, by = "gene", all = TRUE)
buffa_ITFR[, is_seed := fifelse(is.na(is_seed), FALSE, is_seed) | fifelse(is.na(is_seed_FR), FALSE, is_seed_FR)]
buffa_ITFR[, is_seed_FR := NULL]

buffa_ITFR[, var_IT_safe := pmax(C_rank_boot_var_IT, 1e-8, na.rm = TRUE)]
buffa_ITFR[, var_FR_safe := pmax(C_rank_boot_var_FR, 1e-8, na.rm = TRUE)]

buffa_ITFR[, pC_hat_ITxFR_obs := C_rank_IT * C_rank_FR]
buffa_ITFR[, pC_hat_ITxFR_boot := C_rank_boot_mean_IT * C_rank_boot_mean_FR]
buffa_ITFR[, C_hat_IVW_ITFR := (C_rank_boot_mean_IT / var_IT_safe + C_rank_boot_mean_FR / var_FR_safe) / (1 / var_IT_safe + 1 / var_FR_safe)]

buffa_ITFR[, positive_consistent_ITFR := positive_consistent_IT & positive_consistent_FR]
buffa_ITFR[, MC_significant_ITFR := p_MC_C_IT < opt$mc_p_cutoff & p_MC_C_FR < opt$mc_p_cutoff]
buffa_ITFR[, Buffa_product_pass := pC_hat_ITxFR_boot > opt$product_cutoff]
buffa_ITFR[, Buffa_exact_candidate := !is_seed & positive_consistent_ITFR & MC_significant_ITFR & Buffa_product_pass]

buffa_ITFR <- merge(buffa_ITFR, ann, by = "gene", all.x = TRUE)
front_cols <- intersect(
  c("gene", "symbol", "symbol_itfr", "symbol_tcga", "annotation_status", "gene_biotype_consensus"),
  names(buffa_ITFR)
)
setcolorder(buffa_ITFR, c(front_cols, setdiff(names(buffa_ITFR), front_cols)))
setorder(buffa_ITFR, -Buffa_exact_candidate, -pC_hat_ITxFR_boot, -pC_hat_ITxFR_obs, -C_hat_IVW_ITFR)

signature_ITFR <- buffa_ITFR[Buffa_exact_candidate == TRUE]
signature_ITFR_top51 <- head(signature_ITFR, 51)

fwrite(buffa_ITFR, file.path(csv_dir, "Buffa2010_exact_ITxFR_full_ranking.csv"))
fwrite(signature_ITFR, file.path(csv_dir, "GBM_Hypoxia_Buffa2010_exact_ITxFR_signature.csv"))
fwrite(signature_ITFR_top51, file.path(csv_dir, "GBM_Hypoxia_Buffa2010_exact_ITxFR_signature_top51.csv"))
saveRDS(buffa_ITFR, file.path(rds_dir, "Buffa2010_exact_ITxFR_full_ranking.rds"), compress = TRUE)
saveRDS(signature_ITFR, file.path(rds_dir, "GBM_Hypoxia_Buffa2010_exact_ITxFR_signature.rds"), compress = TRUE)
saveRDS(signature_ITFR_top51, file.path(rds_dir, "GBM_Hypoxia_Buffa2010_exact_ITxFR_signature_top51.rds"), compress = TRUE)

stamp("Candidate genes:", nrow(signature_ITFR))
stamp("Top51 genes actually saved:", nrow(signature_ITFR_top51))

# -----------------------------
# Buffa-style hypoxia score
# -----------------------------
buffa_hypoxia_score <- function(expr_mat, signature_genes, dataset_name) {
  expr_mat <- as.matrix(expr_mat)
  genes_present <- intersect(signature_genes, rownames(expr_mat))
  if (length(genes_present) < 3) stop(dataset_name, ": too few signature genes present.")
  E <- matrixStats::colMedians(abs(expr_mat[genes_present, , drop = FALSE]), na.rm = TRUE)
  HS <- rank01(E, higher_is_better = TRUE)
  data.table(
    dataset = dataset_name,
    sample = colnames(expr_mat),
    E_median_abs_expression = as.numeric(E),
    HS_rank01 = as.numeric(HS),
    n_signature_genes_present = length(genes_present),
    n_signature_genes_total = length(signature_genes)
  )
}

score_all <- function(signature_genes, suffix) {
  if (length(signature_genes) < 3) {
    stamp("Skipping score", suffix, "because signature has fewer than 3 genes")
    return(NULL)
  }
  out <- rbindlist(list(
    buffa_hypoxia_score(IT_GBM_RNA, signature_genes, "IT"),
    buffa_hypoxia_score(FR_GBM_RNA, signature_genes, "FR"),
    buffa_hypoxia_score(TCGA_GBM_RNA, signature_genes, "TCGA")
  ))
  fwrite(out, file.path(csv_dir, paste0("Buffa2010_exact_ITxFR_HS_", suffix, "_IT_FR_TCGA.csv")))
  saveRDS(out, file.path(rds_dir, paste0("Buffa2010_exact_ITxFR_HS_", suffix, "_IT_FR_TCGA.rds")), compress = TRUE)
  out
}

HS_all <- score_all(signature_ITFR$gene, "all_signature")
HS_top51 <- score_all(signature_ITFR_top51$gene, "top51")

# -----------------------------
# Basic QC plots
# -----------------------------
try({
  p_rank_product <- ggplot(buffa_ITFR[is_seed == FALSE], aes(x = C_rank_IT, y = C_rank_FR)) +
    geom_point(aes(color = Buffa_exact_candidate), alpha = 0.55, size = 1.0) +
    geom_hline(yintercept = sqrt(opt$product_cutoff), linetype = "dashed") +
    geom_vline(xintercept = sqrt(opt$product_cutoff), linetype = "dashed") +
    theme_classic(base_size = 12) +
    labs(
      title = "Buffa connectivity ranks in IT and FR",
      x = "IT fractional rank of C",
      y = "FR fractional rank of C",
      color = "Selected"
    )
  ggsave(file.path(plot_dir, "Buffa2010_IT_FR_connectivity_rank_scatter.pdf"), p_rank_product, width = 6, height = 5)

  p_product_distribution <- ggplot(buffa_ITFR[is_seed == FALSE], aes(x = pC_hat_ITxFR_boot)) +
    geom_histogram(bins = 60, color = "white") +
    geom_vline(xintercept = opt$product_cutoff, linetype = "dashed") +
    theme_classic(base_size = 12) +
    labs(
      title = "Distribution of IT x FR product score",
      x = "Product score C_hat_IT x C_hat_FR",
      y = "Number of genes"
    )
  ggsave(file.path(plot_dir, "Buffa2010_ITxFR_product_score_distribution.pdf"), p_product_distribution, width = 6, height = 4)

  if (!is.null(HS_top51)) {
    p_HS <- ggplot(HS_top51, aes(x = dataset, y = HS_rank01, fill = dataset)) +
      geom_violin(trim = TRUE, alpha = 0.7) +
      geom_boxplot(width = 0.12, outlier.shape = NA, alpha = 0.8) +
      geom_jitter(width = 0.08, size = 0.7, alpha = 0.45) +
      theme_classic(base_size = 12) +
      labs(
        title = "Buffa-style hypoxia score using ITxFR top signature",
        x = NULL,
        y = "Hypoxia score rank 0-1"
      ) +
      theme(legend.position = "none")
    ggsave(file.path(plot_dir, "Buffa2010_exact_ITxFR_HS_IT_FR_TCGA.pdf"), p_HS, width = 5.5, height = 4.5)
  }
}, silent = FALSE)

# -----------------------------
# Summary and session info
# -----------------------------
summary_dt <- data.table(
  metric = c(
    "common_train_genes",
    "B_BOOT",
    "N_MC",
    "cores",
    "alpha",
    "mc_p_cutoff",
    "product_cutoff",
    "n_candidate_genes",
    "n_top51_saved"
  ),
  value = c(
    length(common_train_genes),
    opt$b_boot,
    opt$n_mc,
    opt$cores,
    opt$alpha,
    opt$mc_p_cutoff,
    opt$product_cutoff,
    nrow(signature_ITFR),
    nrow(signature_ITFR_top51)
  )
)
fwrite(summary_dt, file.path(csv_dir, "Buffa2010_exact_ITxFR_run_summary.csv"))
saveRDS(sessionInfo(), file.path(rds_dir, "sessionInfo.rds"), compress = TRUE)
writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))

stamp("Finished Buffa 2010 exact HPC run")
stamp("Results directory:", run_dir)
