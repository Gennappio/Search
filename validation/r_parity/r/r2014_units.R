#!/usr/bin/env Rscript
# Unit oracle for the verbatim port (search.r2014): one call of one function per case.
#
# A case names a function of the SEARCH 0.0.2 package, read unmodified from
# reference/r_sources/SEARCH/R, or one of the rules of R the port rests on (RULES below), and
# gives its arguments. The script calls it and writes what it returns, or the message of the
# error it stops with, and the files it wrote.
#
# Usage: Rscript r2014_units.R <repository root> <cases.json> <output.json>
#
# Values travel as tagged JSON (decode and encode below). Numbers go in as hexadecimal floating
# point, which R reads exactly, and come out with 17 significant digits, which is enough to read
# back the same double; NA and NaN are kept apart; a missing name is null.
#
# As in r2014_oracle.R two things of the package are replaced: multicore's mclapply by
# parallel's, and the plotting function by an empty one (it needs the igraph of 2014).

suppressPackageStartupMessages(library(parallel))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("usage: r2014_units.R <repository root> <cases.json> <output.json>")
repo <- normalizePath(args[1])
cases_file <- normalizePath(args[2])
out_file <- args[3]
work_root <- tempfile("r2014_units_")
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

source(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))), "r2014_json.R"))

# ---- values in and out ----------------------------------------------------------------------------
strs <- function(v) vapply(v, function(s) if (is.null(s)) NA_character_ else as.character(s), character(1))
nums <- function(v) suppressWarnings(as.numeric(strs(v)))
with_names <- function(out, names) {
  if (!is.null(names)) names(out) <- strs(names)
  out
}
decode <- function(v) {
  kind <- v$t
  if (kind == "null") return(NULL)
  if (kind == "num") return(with_names(nums(v$v), v$names))
  if (kind == "chr") return(with_names(strs(v$v), v$names))
  if (kind == "lgl") return(vapply(v$v, function(s) if (is.null(s)) NA else isTRUE(s), logical(1)))
  if (kind == "mat") {
    return(matrix(nums(v$v), nrow = v$dim[[1]], ncol = v$dim[[2]], dimnames = list(
      if (is.null(v$rownames)) NULL else strs(v$rownames),
      if (is.null(v$colnames)) NULL else strs(v$colnames)
    )))
  }
  if (kind == "list") return(with_names(lapply(v$items, decode), v$names))
  if (kind == "df") {
    df <- as.data.frame(lapply(v$columns, function(col) {
      if (identical(col$type, "numeric")) nums(col$v) else strs(col$v)
    }), stringsAsFactors = FALSE, check.names = FALSE)
    colnames(df) <- vapply(v$columns, function(col) col$name, character(1))
    rownames(df) <- strs(v$rownames)
    return(df)
  }
  stop("decode: unknown kind ", kind)
}
dims <- function(m) raw_json(paste0("[", nrow(m), ",", ncol(m), "]"))
encode <- function(v) {
  if (is.null(v)) return(list(t = "null"))
  if (is.data.frame(v)) return(list(t = "other", class = "data.frame"))
  if (is.matrix(v)) {
    if (is.character(v)) {
      return(list(t = "cmat", dim = dims(v), v = chr(as.vector(v)), rownames = chr(rownames(v)), colnames = chr(colnames(v))))
    }
    return(list(t = "mat", dim = dims(v), v = num(as.vector(v)), rownames = chr(rownames(v)), colnames = chr(colnames(v))))
  }
  if (is.list(v)) return(list(t = "list", items = lapply(unname(v), encode), names = chr(names(v))))
  if (is.character(v)) return(list(t = "chr", v = chr(v), names = chr(names(v))))
  if (is.numeric(v) || is.logical(v)) return(list(t = "num", v = num(v), names = chr(names(v))))
  list(t = "other", class = class(v)[1])
}

# ---- the functions called -------------------------------------------------------------------------
WORK <- NULL
write_gz <- function(path, text) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- gzfile(path, open = "wb")
  on.exit(close(con))
  writeBin(charToRaw(text), con)
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

PACKAGE <- list(
  "getO" = pkg$getO,
  "getOmax" = pkg$getOmax,
  "removeNode" = pkg$removeNode,
  "removeByName" = pkg$removeByName,
  "get.matrix.format" = pkg$get.matrix.format,
  "estimate.nodes.degree" = pkg$estimate.nodes.degree,
  "nodes.correlation" = pkg$nodes.correlation,
  "nodes.glm" = pkg$nodes.glm,
  "nodes.fastLm" = pkg$nodes.fastLm,
  "SEARCH" = pkg$SEARCH,
  "write.data.files" = function(results.obj, iteration) {
    pkg$write.data.files(results.obj = results.obj, data.directory = paste0(WORK, "/"), iteration = iteration)
  },
  ".SEARCH.randomise.inner" = function(...) pkg$.SEARCH.randomise.inner(..., output.directory = WORK),
  "SEARCH.randomise" = function(..., directory = NULL) {
    pkg$SEARCH.randomise(..., output.directory = if (is.null(directory)) WORK else file.path(WORK, directory))
  },
  "summarise.iterations" = function(files, ...) {
    for (name in names(files)) write_gz(file.path(WORK, name), files[[name]])
    pkg$summarise.iterations(output.directory = paste0(WORK, "/real/"), output.directory.NULL = paste0(WORK, "/null/"), ...)
  }
)
# The rules of R the port reproduces, each as the expression the package uses it in.
RULES <- list(
  "v" = function(v) v,
  "v[names] <- values" = function(v, names, values) {
    v[names] <- values
    v
  },
  "v[names]" = function(v, names) v[names],
  "v[v < t]" = function(v, t) v[v < t],
  "v[v > t]" = function(v, t) v[v > t],
  "v[abs(v) > t]" = function(v, t) v[abs(v) > t],
  "names(v) <- names" = function(v, names) {
    names(v) <- names
    v
  },
  "union(a, b)" = function(a, b) union(a, b),
  "intersect(a, b)" = function(a, b) intersect(a, b),
  "setdiff(a, b)" = function(a, b) setdiff(a, b),
  "unique(a)" = function(a) unique(a),
  "a %in% b" = function(a, b) a %in% b,
  "c(a, b)" = function(a, b) c(a, b),
  "sort(v, decreasing = TRUE)" = function(v) sort(v, decreasing = TRUE),
  "as.numeric(sprintf(\"%.3f\", v))" = function(v) as.numeric(sprintf("%.3f", v)),
  "a:b" = function(a, b) a:b,
  "write.table(x)" = function(x) {
    f <- file.path(WORK, "table.txt")
    write.table(x = x, file = f, row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
    rawToChar(readBin(f, what = "raw", n = file.info(f)$size))
  },
  "read.table(file)" = function(text) {
    f <- file.path(WORK, "table.txt.gz")
    write_gz(f, text)
    as.matrix(read.table(file = f, row.names = 1, header = TRUE, sep = "\t"))
  }
)
UNITS <- c(PACKAGE, RULES)
WRITES_FILES <- c("write.data.files", ".SEARCH.randomise.inner", "SEARCH.randomise")

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

run_case <- function(case, index) {
  f <- UNITS[[case$fn]]
  if (is.null(f)) stop("unknown unit: ", case$fn)
  if (identical(case$fn, "nodes.fastLm") && !have_fastlm) return(list(id = case$id, skipped = "RcppArmadillo not installed"))
  WORK <<- file.path(work_root, sprintf("case_%05d", index))
  dir.create(WORK, recursive = TRUE)
  on.exit(unlink(WORK, recursive = TRUE))
  call_args <- lapply(case$args, decode)
  # the package draws its samples with sample(); every case starts from the same state
  set.seed(1)
  run <- quietly(do.call(f, call_args))
  out <- list(id = case$id, error = NULL, value = NULL)
  if (inherits(run$value, "error")) out$error <- conditionMessage(run$value) else out$value <- encode(run$value)
  out$warnings <- as.list(run$warnings)
  if (case$fn %in% WRITES_FILES) out$files <- files_out(WORK)
  out
}

cases <- jsonlite::fromJSON(cases_file, simplifyVector = FALSE)
results <- vector("list", length(cases))
for (i in seq_along(cases)) results[[i]] <- run_case(cases[[i]], i)
unlink(work_root, recursive = TRUE)
info <- list(
  r_version = R.version.string, platform = R.version$platform, sample_kind = RNGkind()[3],
  long_double_bytes = .Machine$sizeof.longdouble, locale_collate = Sys.getlocale("LC_COLLATE"),
  date = format(Sys.time(), "%Y-%m-%d")
)
writeLines(to_json(list(info = info, results = results)), out_file)
cat(sprintf("%d cases written to %s\n", length(cases), out_file))
