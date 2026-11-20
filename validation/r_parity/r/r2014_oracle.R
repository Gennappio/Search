#!/usr/bin/env Rscript
# Oracle for the verbatim port (search.r2014): runs the SEARCH 0.0.2 package, unmodified, on the
# cases of a JSON file and writes everything it returns, or the error it stops with.
#
# The package sources are read from reference/r_sources/SEARCH/R as they are. Two things are
# replaced, as in search2014_run.R: multicore's mclapply by parallel's, and the plotting function
# by an empty one (it needs the igraph of 2014).
#
# Usage: Rscript r2014_oracle.R <repository root> <cases file> <output directory>
#
# The cases file is JSON (read with jsonlite) or, with the extension .R, an R source file that
# defines the list `cases`. The second form and the JSON written here need no package, so the
# script also runs in the R of 2014, where the package was developed.
#
# Numbers are written with 17 significant digits, which is enough to read back the same double;
# NA and NaN are kept apart; a missing name is written as null.

suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: r2014_oracle.R <repository root> <cases.json> <output directory>")
repo <- normalizePath(args[1])
cases_file <- normalizePath(args[2])
out_dir <- args[3]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_dir <- normalizePath(out_dir)
src_dir <- file.path(repo, "reference", "r_sources", "SEARCH", "R")

have_fastlm <- suppressWarnings(require("RcppArmadillo", quietly = TRUE, character.only = TRUE))

SRC_FILES <- c(
  "SEARCH.R", "nodes.correlation.R", "nodes.glm.R", "nodes.fastLm.R", "getO.R", "getOmax.R",
  "removeNode.R", "removeByName.R", "get.matrix.format.R", "estimate.nodes.degree.R",
  "write.data.files.R", "SEARCH.randomise.R", "SEARCH-internal.R", "summarise.iterations.R"
)
pkg <- new.env(parent = globalenv())
for (f in SRC_FILES) sys.source(file.path(src_dir, f), envir = pkg, keep.source = FALSE)
pkg$make.default.plots <- function(...) invisible(NULL)

# ---- serialisation ------------------------------------------------------------------------------
source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "r2014_json.R"))
named_vector <- function(v) {
  list(has_names = !is.null(names(v)), names = chr(names(v)), values = num(v), is_null = is.null(v))
}
matrix_out <- function(m) {
  list(
    dim = raw_json(paste0("[", nrow(m), ",", ncol(m), "]")), has_rownames = !is.null(rownames(m)), rownames = chr(rownames(m)),
    has_colnames = !is.null(colnames(m)), colnames = chr(colnames(m)), values = num(as.vector(m))
  )
}
per_iteration <- function(l, f) {
  if (length(l) == 0) return(raw_json("{}"))
  lapply(l, f)
}
search_out <- function(res, degrees = NULL) {
  cl <- res$correlationLists
  list(
    removed = chr(res$removed),
    correlationLists = list(has_names = !is.null(names(cl)), names = chr(names(cl)), items = lapply(unname(cl), named_vector)),
    D = per_iteration(res$D, named_vector),
    network.growth = per_iteration(res$network.growth, function(v) list(is_null = is.null(v), values = chr(v))),
    O = per_iteration(res$O, matrix_out),
    correlationMatrix = per_iteration(res$correlationMatrix, matrix_out),
    other.data = per_iteration(res$other.data, matrix_out),
    nodes.degree = per_iteration(degrees, named_vector)
  )
}

# ---- inputs -------------------------------------------------------------------------------------
DATA_CACHE <- new.env()
dataset <- function(name) {
  if (is.null(DATA_CACHE[[name]])) {
    if (substr(name, 1, 6) == "Breast") {
      # exactly as in the examples of SEARCH.randomise.Rd
      DATA_CACHE[[name]] <- as.matrix(read.table(
        file = file.path(repo, "tests", "data", "breast", name, "mRNA_abundance.txt.gz"),
        row.names = 1, header = TRUE, sep = "\t"
      ))
    } else {
      gbm <- Sys.getenv("SEARCH_GBM_DATA", unset = "")
      if (!nzchar(gbm)) stop("SEARCH_GBM_DATA is not set")
      m <- as.matrix(read.csv(file.path(gbm, "DATA", "NORMALIZED", paste0(name, "_GBM_RNA_Deseq2_vst_Combat.csv")), row.names = 1, check.names = FALSE))
      storage.mode(m) <- "double"
      DATA_CACHE[[name]] <- m
    }
  }
  DATA_CACHE[[name]]
}
case_matrix <- function(case) {
  d <- case$data
  if (!is.null(d$dataset)) {
    x <- dataset(d$dataset)
    if (!is.null(d$rows)) x <- x[unlist(d$rows) + 1, , drop = FALSE]
    if (!is.null(d$cols)) x <- x[, unlist(d$cols) + 1, drop = FALSE]
  } else {
    v <- suppressWarnings(as.numeric(vapply(d$values, function(s) if (is.null(s)) NA_character_ else as.character(s), character(1))))
    x <- matrix(v, nrow = d$nrow, ncol = d$ncol, byrow = TRUE)
    rownames(x) <- unlist(d$rownames)
    if (!is.null(d$colnames)) colnames(x) <- unlist(d$colnames)
  }
  for (cell in case$na_cells) x[cell[[1]] + 1, cell[[2]] + 1] <- NA
  x
}
strings <- function(v) {
  if (is.null(v)) return(NULL)
  vapply(v, function(s) if (is.null(s)) NA_character_ else as.character(s), character(1))
}
case_covariates <- function(case, x) {
  cv <- case$covariates
  if (is.null(cv)) return(NULL)
  df <- as.data.frame(lapply(cv$columns, function(col) {
    v <- unlist(lapply(col$values, function(s) if (is.null(s)) NA else s))
    if (identical(col$type, "numeric")) suppressWarnings(as.numeric(v)) else as.character(v)
  }), stringsAsFactors = FALSE, check.names = FALSE)
  colnames(df) <- vapply(cv$columns, function(col) col$name, character(1))
  rownames(df) <- colnames(x)
  df
}

quietly <- function(expr) {
  sink_file <- file("/dev/null", open = "wt")
  sink(sink_file)
  sink(sink_file, type = "message")
  on.exit({
    sink(type = "message")
    sink()
    close(sink_file)
  })
  warnings_seen <- character()
  value <- tryCatch(
    withCallingHandlers(expr, warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    }),
    error = function(e) e
  )
  list(value = value, warnings = unique(warnings_seen))
}

read_gz <- function(path) {
  con <- gzfile(path, open = "rb")
  on.exit(close(con))
  rawToChar(readBin(con, what = "raw", n = 5e8))
}
files_out <- function(dir) {
  files <- sort(list.files(dir, recursive = TRUE))
  setNames(lapply(files, function(f) {
    p <- file.path(dir, f)
    if (grepl("\\.gz$", f)) read_gz(p) else rawToChar(readBin(p, what = "raw", n = file.info(p)$size))
  }), files)
}
summary_out <- function(s) {
  dims <- function(v) raw_json(paste0("[", paste(sprintf("%.17g", as.numeric(v)), collapse = ","), "]"))
  lapply(s, function(one) {
    lapply(one, function(el) if (is.list(el)) (if (length(el) == 0) raw_json("{}") else lapply(el, dims)) else dims(el))
  })
}

run_search <- function(case) {
  x <- case_matrix(case)
  universe <- if (is.null(case$universe)) rownames(x) else strings(case$universe)
  a <- c(
    list(x = x, Alist = strings(case$seeds), Blist = universe, covariates = case_covariates(case, x)),
    case$args, case$assoc_args
  )
  run <- quietly(do.call(pkg$SEARCH, a))
  res <- run$value
  if (inherits(res, "error")) return(list(error = conditionMessage(res), warnings = as.list(run$warnings)))
  deg <- quietly(pkg$estimate.nodes.degree(correlation.matrix.list = res$correlationMatrix))
  degrees <- if (inherits(deg$value, "error")) NULL else deg$value
  c(list(error = NULL, warnings = as.list(run$warnings)), search_out(res, degrees))
}

run_randomise <- function(case, dir) {
  x <- case_matrix(case)
  unlink(dir, recursive = TRUE)
  dir.create(dir, recursive = TRUE)
  # R 3.6.0 changed sample(); before it there is one kind only, the one now called "Rounding"
  if (HAS_SAMPLE_KIND) {
    suppressWarnings(RNGkind(sample.kind = if (is.null(case$sample_kind)) "Rejection" else case$sample_kind))
    on.exit(suppressWarnings(RNGkind(sample.kind = "Rejection")))
  }
  out <- list(error = NULL, sample_kind = if (HAS_SAMPLE_KIND) RNGkind()[3] else "Rounding")
  real_dir <- if (is.null(case$real_dir)) "real/" else case$real_dir
  universe <- if (is.null(case$universe)) rownames(x) else strings(case$universe)
  runs <- c(list(list(dir = real_dir, seeds = case$seeds, n.boot = case$args$n.boot, n.boot.start = case$args$n.boot.start)), case$null_runs)
  summaries <- list()
  for (r in runs) {
    a <- c(
      list(data = x, seed.nodes = strings(r$seeds), universe.nodes = universe, covariates = case_covariates(case, x),
        output.directory = paste0(dir, "/", r$dir)),
      case$args, case$assoc_args
    )
    a$n.boot <- r$n.boot
    if (!is.null(r$n.boot.start)) a$n.boot.start <- r$n.boot.start
    run <- quietly(do.call(pkg$SEARCH.randomise, a))
    if (inherits(run$value, "error")) {
      out$error <- conditionMessage(run$value)
      out$error_run <- r$dir
      break
    }
    summaries[[r$dir]] <- summary_out(run$value)
  }
  out$summaries <- summaries
  if (is.null(out$error) && !is.null(case$summarise)) {
    s <- case$summarise
    run <- quietly(pkg$summarise.iterations(
      output.directory = paste0(dir, "/", real_dir), max.iter = s$max.iter, n.boot = s$n.boot,
      output.directory.NULL = paste0(dir, "/null/"), n.boot.NULL = s$n.boot.NULL, batches.NULL = s$batches.NULL
    ))
    if (inherits(run$value, "error")) {
      out$summarise_error <- conditionMessage(run$value)
    } else {
      out$summarise <- list(
        real = per_iteration(run$value$real, matrix_out),
        null = per_iteration(run$value$`NULL`, matrix_out)
      )
    }
  }
  out$files <- files_out(dir)
  out
}

HAS_SAMPLE_KIND <- "sample.kind" %in% names(formals(RNGkind))
if (grepl("\\.R$", cases_file)) {
  source(cases_file, local = TRUE)
} else {
  cases <- jsonlite::fromJSON(cases_file, simplifyVector = FALSE)
}
results <- vector("list", length(cases))
for (i in seq_along(cases)) {
  case <- cases[[i]]
  if (identical(case$args$nodes.association.function, "nodes.fastLm") && !have_fastlm) {
    results[[i]] <- list(id = case$id, skipped = "RcppArmadillo not installed")
    next
  }
  res <- if (identical(case$kind, "randomise")) run_randomise(case, file.path(out_dir, "work", case$id)) else run_search(case)
  results[[i]] <- c(list(id = case$id), res)
}
info <- list(
  r_version = R.version.string, platform = R.version$platform,
  sample_kind = if (HAS_SAMPLE_KIND) RNGkind()[3] else "Rounding",
  long_double_bytes = .Machine$sizeof.longdouble, locale_collate = Sys.getlocale("LC_COLLATE"),
  date = format(Sys.time(), "%Y-%m-%d")
)
writeLines(to_json(list(info = info, results = results)), file.path(out_dir, "oracle.json"))
cat(sprintf("%d cases written to %s\n", length(cases), file.path(out_dir, "oracle.json")))
