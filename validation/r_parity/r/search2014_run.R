#!/usr/bin/env Rscript
# Parity harness, R side: the SEARCH 0.0.2 package (Haider and Buffa, 2014).
#
# The package is not installed (it depends on the removed `multicore`): its files under
# reference/r_sources/SEARCH/R are sourced as they are into an environment, with
# parallel::mclapply standing in for multicore's and the plotting function stubbed out.
#
# Every case is run with the sources "verbatim" and, when the case lists patches, "patched":
# the same sources with the listed one-line corrections, each of which is a defect recorded in
# docs/decisions.md. The Python comparison must reproduce the patched run exactly, and the
# verbatim run wherever the defect does not come into play.
#
# Usage: Rscript search2014_run.R <repository root> <output directory> [case id ...]
# The GBM cases run only when SEARCH_GBM_DATA points to the GBM project folder.

suppressPackageStartupMessages({
  library(parallel)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("usage: search2014_run.R <repository root> <output directory> [case id ...]")
repo <- normalizePath(args[1])
out_root <- args[2]
only <- args[-(1:2)]
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)
out_root <- normalizePath(out_root)
src_dir <- file.path(repo, "reference", "r_sources", "SEARCH", "R")
cores <- max(1L, min(8L, parallel::detectCores(logical = FALSE)))

have_fastlm <- requireNamespace("RcppArmadillo", quietly = TRUE)
if (have_fastlm) suppressPackageStartupMessages(library(RcppArmadillo))

SRC_FILES <- c(
  "SEARCH.R", "nodes.correlation.R", "nodes.glm.R", "nodes.fastLm.R", "getO.R", "getOmax.R",
  "removeNode.R", "removeByName.R", "get.matrix.format.R", "estimate.nodes.degree.R",
  "write.data.files.R", "SEARCH.randomise.R", "SEARCH-internal.R", "summarise.iterations.R"
)

# One-line corrections of the reference sources. `find` must occur exactly once in `file`.
PATCHES <- list(
  # decisions.md 1: the standard deviations are computed in matrix row order and were named in
  # universe order.
  sd_names = list(
    file = "SEARCH.R",
    find = "names(sds) <- Blist;",
    replace = "names(sds) <- rownames(x)[rownames(x) %in% Blist];"
  ),
  # decisions.md 2: Omax compared a seed only with the seeds after it in the list.
  omax_full = list(
    file = "getOmax.R",
    find = "if(y.index > x.index) {",
    replace = "if(y.index != x.index) {"
  ),
  # decisions.md 5: seeds pruned on Omax stayed in the neighbourhood lists when they were old
  # seeds (keep.Alist = FALSE).
  drop_pruned_lists = list(
    file = "SEARCH.R",
    find = "newCL <- removeByName(to.remove, newCL);",
    replace = "newCL <- removeByName(to.remove, newCL); correlationLists <- removeByName(to.remove, correlationLists);"
  ),
  # decisions.md 5, wider than recorded: a gene pruned from the universe (a seed removed on D or
  # on Omax) stayed as a member in the neighbourhoods of the other seeds, so it kept counting in
  # every later overlap O.
  drop_d_pruned_members = list(
    file = "SEARCH.R",
    find = "correlationLists <- removeByName(to.remove, correlationLists);",
    replace = "correlationLists <- removeByName(to.remove, correlationLists); correlationLists <- lapply(correlationLists, function(cl) cl[!(names(cl) %in% to.remove)]);"
  ),
  drop_o_pruned_members = list(
    file = "SEARCH.R",
    find = "correlationLists <- c(newCL, correlationLists);",
    replace = "correlationLists <- lapply(c(newCL, correlationLists), function(cl) cl[!(names(cl) %in% to.remove)]);"
  ),
  # not in decisions.md, found by this harness: subsetting one row and one column of the result
  # table drops the name, so a seed whose neighbourhood is itself alone gets an unnamed list,
  # is not counted in its own D and ends at D = -1 / (K - 1).
  glm_single_name = list(
    file = "nodes.glm.R",
    find = "glm.res.filtered <- unlist(glm.res.filtered);",
    replace = "glm.res.filtered <- unlist(glm.res.filtered); names(glm.res.filtered) <- x.rownames[which(glm.res.table[, \"P\"] < significance.threshold)];"
  ),
  # not in decisions.md, found by this harness: with missing values the correlation is NA, and
  # subsetting by a logical vector that holds NA adds an entry named NA to the neighbourhood;
  # the entry then counts as a shared member in every overlap O.
  cor_na_names = list(
    file = "nodes.correlation.R",
    find = "function(q) { q[abs(q) > correlation.threshold] }",
    replace = "function(q) { q[which(abs(q) > correlation.threshold)] }"
  ),
  fastlm_single_name = list(
    file = "nodes.fastLm.R",
    find = "glm.res.filtered <- unlist(glm.res.filtered);",
    replace = "glm.res.filtered <- unlist(glm.res.filtered); names(glm.res.filtered) <- x.rownames[which(glm.res.table[, \"P\"] < significance.threshold)];"
  )
)

load_sources <- function(patches = character()) {
  env <- new.env(parent = globalenv())
  for (f in SRC_FILES) {
    original <- readLines(file.path(src_dir, f), warn = FALSE)
    txt <- original
    for (p in patches) {
      spec <- PATCHES[[p]]
      if (is.null(spec)) stop("unknown patch: ", p)
      if (spec$file != f) next
      hit <- grep(spec$find, original, fixed = TRUE) # located on the unpatched text
      if (length(hit) != 1) stop("patch ", p, ": expected one match in ", f, ", found ", length(hit))
      txt[hit] <- sub(spec$find, spec$replace, txt[hit], fixed = TRUE)
    }
    eval(parse(text = txt, keep.source = FALSE), envir = env)
  }
  # the plots need igraph and gplots and are not part of the comparison
  env$make.default.plots <- function(...) invisible(NULL)
  env
}

# ---- data ------------------------------------------------------------------------------------
read_breast <- function(name) {
  # exactly as in the examples of SEARCH.randomise.Rd
  as.matrix(read.table(
    file = file.path(repo, "tests", "data", "breast", name, "mRNA_abundance.txt.gz"),
    row.names = 1, header = TRUE, sep = "\t"
  ))
}

read_breast_annotation <- function(name, samples) {
  ann <- read.table(
    file = file.path(repo, "tests", "data", "breast", name, "patient_annotation.txt"),
    header = TRUE, sep = "\t", stringsAsFactors = FALSE
  )
  rownames(ann) <- make.names(ann$ID)
  ann[samples, , drop = FALSE]
}

gbm_dir <- Sys.getenv("SEARCH_GBM_DATA", unset = "")
read_gbm <- function(cohort) {
  path <- file.path(gbm_dir, "DATA", "NORMALIZED", paste0(cohort, "_GBM_RNA_Deseq2_vst_Combat.csv"))
  tab <- read.csv(path, row.names = 1, check.names = FALSE)
  m <- as.matrix(tab)
  storage.mode(m) <- "double"
  m
}

seed_file <- function(path) {
  as.character(read.csv(file.path(repo, path), stringsAsFactors = FALSE)$gene)
}

DATA_CACHE <- new.env()
get_data <- function(name) {
  if (is.null(DATA_CACHE[[name]])) {
    DATA_CACHE[[name]] <- if (startsWith(name, "Breast")) read_breast(name) else read_gbm(name)
  }
  DATA_CACHE[[name]]
}

# ---- output helpers ---------------------------------------------------------------------------
num <- function(v) ifelse(is.na(v), "NA", sprintf("%.17g", v))

write_named <- function(v, path, value = "value") {
  nm <- names(v)
  if (is.null(nm)) nm <- rep(NA_character_, length(v))
  tab <- data.frame(gene = ifelse(is.na(nm), "NA", nm), value = num(as.numeric(v)), stringsAsFactors = FALSE)
  colnames(tab)[2] <- value
  write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

write_matrix <- function(m, path) {
  tab <- data.frame(row = rownames(m), matrix(num(m), nrow = nrow(m)), stringsAsFactors = FALSE)
  colnames(tab) <- c("row", colnames(m))
  write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

dump_search <- function(env, res, dir, x, light = FALSE) {
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  if (inherits(res, "error")) {
    writeLines(conditionMessage(res), file.path(dir, "error.txt"))
    return(invisible(FALSE))
  }
  for (it in names(res$D)) write_named(res$D[[it]], file.path(dir, sprintf("D_iter_%s.tsv", it)), "D")
  for (it in names(res$network.growth)) {
    writeLines(as.character(res$network.growth[[it]]), file.path(dir, sprintf("added_iter_%s.txt", it)))
  }
  writeLines(as.character(res$removed), file.path(dir, "removed_by_D.txt"))
  cl <- res$correlationLists
  sizes <- vapply(cl, length, integer(1))
  write_named(sizes, file.path(dir, "final_cluster_sizes.tsv"), "size")
  if (!light) {
    for (it in names(res$O)) write_matrix(res$O[[it]], file.path(dir, sprintf("O_iter_%s.tsv", it)))
    if (length(cl) > 0) {
      # O of the final lists, with the package's own function
      mf <- env$get.matrix.format(x = cl, all.features = rownames(x))
      write_matrix(mf[["O"]], file.path(dir, "O_final.tsv"))
      long <- do.call(rbind, lapply(names(cl), function(s) {
        v <- cl[[s]]
        if (length(v) == 0) return(NULL)
        nm <- names(v)
        data.frame(seed = s, gene = ifelse(is.na(nm), "NA", nm), value = num(as.numeric(v)), stringsAsFactors = FALSE)
      }))
      write.table(long, file.path(dir, "final_neighbourhoods.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    }
    if (length(res$correlationMatrix) > 0) {
      deg <- env$estimate.nodes.degree(correlation.matrix.list = res$correlationMatrix)
      for (it in names(deg)) write_named(deg[[it]], file.path(dir, sprintf("degree_iter_%s.tsv", it)), "degree")
    }
  }
  invisible(TRUE)
}

quietly <- function(expr, log_path) {
  con <- file(log_path, open = "wt")
  sink(con)
  sink(con, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(con)
  })
  tryCatch(withCallingHandlers(expr, warning = function(w) invokeRestart("muffleWarning")),
    error = function(e) e
  )
}

search_args <- function(case, x, covariates) {
  c(
    list(
      x = x, Alist = case$seeds, Blist = case$universe, covariates = covariates,
      scale = case$scale, D.lower = case$D.lower, D.higher = case$D.higher,
      O.lower = case$O.lower, keep.Alist = case$keep.Alist, max.iter = case$max.iter,
      sd.threshold = case$sd.threshold, verbose = TRUE,
      nodes.association.function = case$association
    ),
    case$assoc_args
  )
}

# ---- cases ------------------------------------------------------------------------------------
base_case <- function(...) {
  case <- list(
    kind = "growth", dataset = "Breastdata1", seeds = "cluster", association = "nodes.correlation",
    assoc_args = list(correlation.threshold = 0.5, correlation.method = "spearman"),
    scale = TRUE, D.lower = 0.1, D.higher = 0.24, O.lower = 0.01, keep.Alist = TRUE,
    max.iter = 10, sd.threshold = 0, patches = character(), covariates = character(),
    covariate_names = character(), universe = "all", inject_na = FALSE, light = FALSE, note = ""
  )
  given <- list(...)
  for (name in names(given)) case[[name]] <- given[[name]] # no recursive merge of assoc_args
  case
}

glm_args <- function(...) utils::modifyList(list(significance.threshold = 1e-05, glm.family = "gaussian", mc.cores = cores), list(...))

CASES <- list(
  base_case(
    id = "breast_single_pass", max.iter = 1,
    note = "one pass with the thresholds of SEARCH.randomise, no sd prefilter"
  ),
  base_case(
    id = "breast_growth",
    note = "growth to convergence with the thresholds of SEARCH.randomise, no sd prefilter"
  ),
  base_case(
    id = "breast_doc_example", seeds = "first10", D.lower = 0, D.higher = 0, O.lower = 0.1, max.iter = 3,
    patches = c("omax_full"),
    note = "the correlation example of SEARCH.randomise.Rd, run through SEARCH()"
  ),
  base_case(
    id = "breast_sd_prefilter", sd.threshold = 0.1, patches = c("sd_names"),
    note = "all the defaults of SEARCH.randomise, sd prefilter included"
  ),
  base_case(
    id = "breast_seeds_not_immune", keep.Alist = FALSE, patches = c("omax_full", "drop_pruned_lists"),
    note = "keep.seeds = FALSE: every seed can be pruned on Omax"
  ),
  base_case(
    id = "breast_doc_not_immune", seeds = "first10", D.lower = 0, D.higher = 0, O.lower = 0.1,
    max.iter = 3, keep.Alist = FALSE, patches = c("omax_full", "drop_pruned_lists"),
    note = "the documented example with keep.seeds = FALSE"
  ),
  base_case(
    id = "breast_pearson",
    assoc_args = list(correlation.threshold = 0.5, correlation.method = "pearson"),
    note = "Pearson instead of Spearman"
  ),
  base_case(
    id = "breast_unscaled", scale = FALSE, D.lower = 1, D.higher = 2, max.iter = 3,
    note = "scale = FALSE: D is a count"
  ),
  base_case(
    id = "breast_universe_subset", universe = "first600",
    note = "universe restricted to the first 600 probesets (the seeds are added by SEARCH itself)"
  ),
  base_case(
    id = "breast_pruned_members", D.higher = 0.1, patches = c("drop_d_pruned_members"),
    note = "a lower D.higher: more genes become seeds and some are later removed on D"
  ),
  base_case(
    id = "breast_omax_pruning", D.higher = 0.1, O.lower = 0.3,
    patches = c("omax_full", "drop_o_pruned_members", "drop_d_pruned_members"),
    note = "a higher O.lower: new seeds are removed on Omax"
  ),
  base_case(
    id = "breast_missing_values", inject_na = TRUE, D.higher = 0.1,
    patches = c("cor_na_names", "drop_d_pruned_members"),
    note = "as breast_pruned_members with one missing value in 25 probesets that are not seeds"
  ),
  base_case(
    id = "breast_glm_doc", seeds = "first10", association = "nodes.glm", assoc_args = glm_args(),
    D.lower = 0, D.higher = 0, O.lower = 0.1, max.iter = 2, patches = c("omax_full", "glm_single_name"),
    note = "the GLM example of SEARCH.randomise.Rd, run through SEARCH()"
  ),
  base_case(
    id = "breast_fastlm_doc", seeds = "first10", association = "nodes.fastLm", assoc_args = glm_args(),
    D.lower = 0, D.higher = 0, O.lower = 0.1, max.iter = 2, patches = c("omax_full", "fastlm_single_name"),
    note = "the same example with nodes.fastLm"
  ),
  base_case(
    id = "breast_glm_covariates_plain", association = "nodes.glm",
    assoc_args = glm_args(significance.threshold = 0.001, formula.string = "t_rfs + e_rfs"),
    covariates = c("t.rfs", "e.rfs"), covariate_names = c("t_rfs", "e_rfs"), max.iter = 1,
    note = "the same two covariates under names without dots"
  ),
  base_case(
    id = "breast_glm_covariates", association = "nodes.glm",
    assoc_args = glm_args(significance.threshold = 0.001, formula.string = "t.rfs + e.rfs"),
    covariates = c("t.rfs", "e.rfs"), max.iter = 1,
    note = "gaussian GLM with two covariates given through formula.string"
  ),
  base_case(
    id = "breast_glm_covariates_ignored", association = "nodes.glm",
    assoc_args = glm_args(significance.threshold = 0.001), covariates = c("t.rfs", "e.rfs"),
    max.iter = 1,
    note = "covariates passed without formula.string: nodes.glm leaves them out of the model"
  ),
  base_case(
    id = "breast_glm_poisson", association = "nodes.glm",
    assoc_args = glm_args(glm.family = "poisson"), max.iter = 1,
    note = "poisson family at the documented threshold: no association passes on these data"
  ),
  base_case(
    id = "breast_glm_poisson_loose", association = "nodes.glm",
    assoc_args = glm_args(glm.family = "poisson", significance.threshold = 0.9), max.iter = 1,
    note = "poisson family (z statistics) with a threshold loose enough to give neighbourhoods"
  ),
  base_case(
    id = "breast_glm_gamma", association = "nodes.glm",
    assoc_args = glm_args(glm.family = "Gamma"), max.iter = 1,
    note = "Gamma family"
  ),
  base_case(
    id = "breast_glm_bh", association = "nodes.glm",
    assoc_args = glm_args(significance.threshold = 0.01, adjust.p = "BH"), max.iter = 1,
    note = "gaussian GLM with Benjamini-Hochberg adjustment of the p-values of each seed"
  ),
  base_case(
    id = "breast_fastlm", association = "nodes.fastLm", assoc_args = glm_args(), max.iter = 1,
    note = "nodes.fastLm (RcppArmadillo least squares)"
  ),
  base_case(
    id = "it_growth", dataset = "IT", seeds = "gbm", max.iter = 3, light = TRUE,
    note = "GBM IT cohort, ten hypoxia seeds, three growth iterations"
  ),
  base_case(
    id = "it_sd_prefilter", dataset = "IT", seeds = "gbm", max.iter = 1, sd.threshold = 0.1,
    patches = c("sd_names"), light = TRUE,
    note = "GBM IT cohort, one pass with the sd prefilter"
  ),
  base_case(
    id = "breast_bootstrap", kind = "bootstrap", max.iter = 2, n.boot = 5, sampling = "bootstrap",
    sampling.seed = 1, n.batches.null = 2, n.boot.null = 5,
    note = "SEARCH.randomise with bootstrap, null runs with random seeds, summarise.iterations"
  ),
  base_case(
    id = "breast_bootstrap_6iter", kind = "bootstrap", max.iter = 6, n.boot = 8, sampling = "bootstrap",
    sampling.seed = 1, n.batches.null = 1, n.boot.null = 3,
    note = "as breast_bootstrap with up to six growth iterations: replicates stop at different iterations"
  ),
  base_case(
    id = "breast_bootstrap_defaults", kind = "bootstrap", max.iter = 1, sd.threshold = 0.1,
    n.boot = 10, sampling = "bootstrap", sampling.seed = 1, n.batches.null = 2, n.boot.null = 5,
    patches = c("sd_names"),
    note = "SEARCH.randomise with every default (the search2014 preset), bootstrap of ten"
  ),
  base_case(
    id = "breast_jackknife", kind = "bootstrap", max.iter = 1, n.boot = 4, sampling = "jacknife",
    sampling.seed = 7, n.batches.null = 1, n.boot.null = 3,
    note = "SEARCH.randomise with sampling = jacknife"
  )
)

resolve_case <- function(case) {
  x <- get_data(case$dataset)
  case$seeds <- switch(case$seeds,
    cluster = seed_file("examples/breast/seeds_breast_cluster.csv"),
    first10 = rownames(x)[1:10],
    gbm = seed_file("examples/gbm/seeds_buffa2010_setA.csv"),
    stop("unknown seed choice")
  )
  case$universe <- switch(case$universe,
    all = rownames(x),
    first600 = rownames(x)[1:600],
    stop("unknown universe choice")
  )
  case
}

# One missing value (third sample) in every 50th probeset from the 25th, seeds excluded.
NA_ROWS <- function(x, seeds) setdiff(rownames(x)[seq(25, nrow(x), by = 50)], seeds)
case_data <- function(case) {
  x <- get_data(case$dataset)
  if (isTRUE(case$inject_na)) x[NA_ROWS(x, case$seeds), 3] <- NA
  x
}

covariates_of <- function(case, x) {
  if (length(case$covariates) == 0) return(NULL)
  ann <- read_breast_annotation(case$dataset, colnames(x))
  cov <- ann[, case$covariates, drop = FALSE]
  if (length(case$covariate_names) > 0) colnames(cov) <- case$covariate_names
  cov
}

run_growth_case <- function(case, dir) {
  x <- case_data(case)
  covariates <- covariates_of(case, x)
  variants <- list(verbatim = character())
  if (length(case$patches) > 0) variants$patched <- case$patches
  for (v in names(variants)) {
    vdir <- file.path(dir, paste0("r_", v))
    unlink(vdir, recursive = TRUE)
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    env <- load_sources(variants[[v]])
    t0 <- proc.time()[["elapsed"]]
    res <- quietly(do.call(env$SEARCH, search_args(case, x, covariates)), file.path(vdir, "log.txt"))
    writeLines(sprintf("%.2f", proc.time()[["elapsed"]] - t0), file.path(vdir, "seconds.txt"))
    dump_search(env, res, vdir, x, light = case$light)
  }
}

run_bootstrap_case <- function(case, dir) {
  x <- case_data(case)
  n <- ncol(x)
  k <- length(case$seeds)
  # the draws of .SEARCH.randomise.inner, reproduced to hand the same samples to Python
  draw <- function(i) {
    set.seed(i + as.numeric(case$sampling.seed) - 1)
    idx <- 1:n
    if (case$sampling == "bootstrap") sample(idx, replace = TRUE) else sample(idx, size = n - 1, replace = FALSE)
  }
  variants <- list(verbatim = character())
  if (length(case$patches) > 0) variants$patched <- case$patches
  # random seed sets of the null runs (the package left this to the user)
  set.seed(20141)
  null_seeds <- lapply(seq_len(case$n.batches.null), function(b) sample(setdiff(rownames(x), case$seeds), k))
  write_json(null_seeds, file.path(dir, "null_seed_sets.json"), auto_unbox = TRUE)
  draws <- t(vapply(seq_len(max(case$n.boot, case$n.boot.null)), draw, integer(if (case$sampling == "bootstrap") n else n - 1)))
  write.table(draws, file.path(dir, "sample_draws_1based.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

  for (v in names(variants)) {
    vdir <- file.path(dir, paste0("r_", v))
    unlink(vdir, recursive = TRUE)
    dir.create(vdir, recursive = TRUE, showWarnings = FALSE)
    env <- load_sources(variants[[v]])
    randomise <- function(seeds, output.directory, n.boot) {
      do.call(env$SEARCH.randomise, c(
        list(
          data = x, seed.nodes = seeds, universe.nodes = case$universe, covariates = NULL,
          scale = case$scale, D.lower = case$D.lower, D.higher = case$D.higher,
          O.lower = case$O.lower, keep.seeds = case$keep.Alist, max.iter = case$max.iter,
          verbose = FALSE, sd.threshold = case$sd.threshold, n.boot = n.boot,
          sampling = case$sampling, sampling.seed = case$sampling.seed,
          output.directory = output.directory, nodes.association.function = case$association
        ),
        case$assoc_args
      ))
    }
    real_dir <- file.path(vdir, "real")
    null_root <- file.path(vdir, "null")
    runs <- quietly({
      randomise(case$seeds, real_dir, case$n.boot)
      for (b in seq_len(case$n.batches.null)) {
        randomise(null_seeds[[b]], file.path(null_root, sprintf("batch_%d_0_FALSE_0_0", b)), case$n.boot.null)
      }
    }, file.path(vdir, "log.txt"))
    if (inherits(runs, "error")) {
      writeLines(conditionMessage(runs), file.path(vdir, "error.txt"))
      next
    }
    res <- quietly(
      env$summarise.iterations(
        output.directory = paste0(real_dir, "/"), max.iter = case$max.iter, n.boot = case$n.boot,
        output.directory.NULL = paste0(null_root, "/"), n.boot.NULL = case$n.boot.null,
        batches.NULL = case$n.batches.null
      ),
      file.path(vdir, "log_summarise.txt")
    )
    if (inherits(res, "error")) {
      writeLines(conditionMessage(res), file.path(vdir, "summarise_error.txt"))
    } else {
      for (it in names(res$real)) {
        write_matrix(res$real[[it]], file.path(vdir, sprintf("summary_real_iter_%s.tsv", it)))
        if (!is.null(res$`NULL`[[it]])) write_matrix(res$`NULL`[[it]], file.path(vdir, sprintf("summary_null_iter_%s.tsv", it)))
      }
    }
    # the same replicates through SEARCH() directly, to compare D without the five-decimal
    # rounding of the files
    for (i in seq_len(case$n.boot)) {
      xi <- x[, draws[i, ]]
      single <- quietly(do.call(env$SEARCH, search_args(case, xi, NULL)), file.path(vdir, sprintf("log_rep_%d.txt", i)))
      dump_search(env, single, file.path(vdir, sprintf("rep_%d", i)), xi, light = TRUE)
    }
  }
}

# ---- main -------------------------------------------------------------------------------------
for (case in CASES) {
  if (length(only) > 0 && !(case$id %in% only)) next
  if (!startsWith(case$dataset, "Breast") && !nzchar(gbm_dir)) {
    cat(sprintf("%-28s skipped: SEARCH_GBM_DATA not set\n", case$id))
    next
  }
  if (case$association == "nodes.fastLm" && !have_fastlm) {
    cat(sprintf("%-28s skipped: RcppArmadillo not installed\n", case$id))
    next
  }
  t0 <- proc.time()[["elapsed"]]
  case <- resolve_case(case)
  dir <- file.path(out_root, case$id)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  spec <- case
  x_all <- get_data(case$dataset)
  # NULL: the universe is every row of the matrix, in matrix order
  if (length(case$universe) == nrow(x_all)) spec$universe <- NULL
  spec$na_rows <- if (isTRUE(case$inject_na)) NA_ROWS(x_all, case$seeds) else character()
  spec$assoc_args$mc.cores <- NULL
  write_json(spec, file.path(dir, "case.json"), auto_unbox = TRUE, digits = NA, null = "null", pretty = TRUE)
  if (case$kind == "growth") run_growth_case(case, dir) else run_bootstrap_case(case, dir)
  cat(sprintf("%-28s done in %.1f s\n", case$id, proc.time()[["elapsed"]] - t0))
}
