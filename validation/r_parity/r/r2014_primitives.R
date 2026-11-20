#!/usr/bin/env Rscript
# The building blocks of R that the SEARCH package leans on, asked directly: the random draws of
# sample(), how a number is read from text and written to it, round(), sd(), cor() and p.adjust().
# The verbatim port (search.r2014) reimplements each of them; py/check_r2014_primitives.py
# compares. Run in a second R, this script also shows which of them changed between versions and
# platforms.
#
# Usage: Rscript r2014_primitives.R <output directory>
# No package is needed. Everything is drawn from fixed seeds.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: r2014_primitives.R <output directory>")
out_dir <- args[1]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
here <- (function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
})()
source(file.path(here, "r2014_json.R"))
HAS_SAMPLE_KIND <- "sample.kind" %in% names(formals(RNGkind))
quiet <- function(expr) suppressWarnings(expr)

# ---- sample() ------------------------------------------------------------------------------------
kinds <- if (HAS_SAMPLE_KIND) c("Rejection", "Rounding") else "Rounding"
rng <- list()
for (kind in kinds) {
  if (HAS_SAMPLE_KIND) quiet(RNGkind("Mersenne-Twister", "Inversion", kind))
  for (seed in c(1, 2, 7, 42, 5784, 1000001, 2147483647, -5, 0)) {
    for (n in c(1, 2, 3, 5, 77, 95, 279, 1000, 65536, 65537)) {
      set.seed(seed); boot <- sample(1:n, replace = TRUE)
      set.seed(seed); jack <- if (n > 1) sample(1:n, size = n - 1, replace = FALSE) else integer(0)
      set.seed(seed); u <- runif(3)
      rng[[length(rng) + 1]] <- list(kind = kind, seed = seed, n = n,
        boot_head = num(head(boot, 200)), boot_sum = sum(as.numeric(boot) * seq_along(boot)),
        jack_head = num(head(jack, 200)), jack_sum = sum(as.numeric(jack) * seq_along(jack)), unif = num(u))
    }
  }
}
if (HAS_SAMPLE_KIND) quiet(RNGkind("Mersenne-Twister", "Inversion", "Rejection"))

# ---- numbers to text and back ----------------------------------------------------------------------
set.seed(5)
v <- c()
for (K in c(1:300, 640, 1000, 1024, 3200, 21254)) v <- c(v, round((0:min(K, 300)) / K, 5), round(-1 / K, 5))
v <- c(v, round(runif(3000), 3), round(rnorm(3000), 3), round(runif(2000) * 1e-3, 5), 10^(-12:12), 123456 * 10^(-10:10),
  9996, 99999, 100000, 999999.5, 1e15, 1e16, 123456789012345678, 0.1 + 0.2, 1 / 3, 2 / 3, 1e-300, 1e300, 100, 1234.5,
  1e-5, 1e-4, 1e-3, 99.99999, 0.999995, 15, 150000, 1500000, 0.5, -0.5, 3e-5, 3.1e-5, 0.00031, 123456.7)
v <- c(v, rnorm(10000) * 10^sample(-20:20, 10000, TRUE), runif(10000), exp(rnorm(3000, 0, 20)))
f <- tempfile()
write.table(data.frame(x = v), f, row.names = FALSE, col.names = FALSE, quote = FALSE)
written <- readLines(f)

texts <- c()
for (digits in c(6, 12, 15, 16, 17, 18, 20)) {
  x <- ifelse(runif(3000) < 0.8, runif(3000, 0, 20), runif(3000, -1e-4, 1e6))
  texts <- c(texts, sprintf(paste0("%.", digits, "g"), x))
}
read_back <- as.numeric(texts)

# ---- round() -----------------------------------------------------------------------------------------
r <- c()
for (K in c(1:400, 640, 1000, 1024, 3200, 6400, 12800, 21254)) r <- c(r, (0:min(K, 400)) / K, -1 / K, (0:min(K, 50)) / (K - 0.999999))
r <- c(r, rnorm(10000), runif(10000) * 1e-4, runif(2000) * 1e6, 0.5 + (0:1000) * 1e-6, 5e-6 * (1:2000), NA, NaN, Inf, -Inf, 0,
  1e-300, 1e300, 2.5, 0.000005, 0.000015, 0.000025, 0.000035, 0.003125, 0.015625)

# ---- sd() and cor() ------------------------------------------------------------------------------------
set.seed(11)
methods <- c("pearson", "spearman", "kendall")
uses <- c("everything", "all.obs", "complete.obs", "na.or.complete", "pairwise.complete.obs")
cors <- list()
for (rep in 1:60) {
  g <- sample(3:14, 1); n <- sample(c(2, 3, 5, 8, 13, 21, 77), 1); k <- sample(1:3, 1)
  x <- matrix(rnorm(g * n, 5, 2), g, n)
  if (rep %% 2 == 0) x <- round(x)                                # ties
  if (rep %% 3 == 0) x[sample(length(x), sample(1:4, 1))] <- NA
  if (rep %% 5 == 0) x[sample(g, 1), ] <- 3                       # a constant row
  if (rep %% 7 == 0) x[sample(g, 1), ] <- NA                      # a row of missing values
  rownames(x) <- paste0("g", 1:g)
  tx <- t(x); seeds <- sample(rownames(x), k)
  for (m in methods) for (u in uses) {
    res <- tryCatch(quiet(cor(tx[, seeds], tx, method = m, use = u)), error = function(e) conditionMessage(e))
    cors[[length(cors) + 1]] <- list(g = g, n = n, x = num(as.vector(x)), seeds = num(match(seeds, rownames(x)) - 1),
      method = m, use = u, error = if (is.character(res)) res else NULL,
      r = if (is.character(res)) NULL else num(as.vector(matrix(res, nrow = k))), sd = num(apply(x, 1, sd)))
  }
}

# ---- p.adjust() ------------------------------------------------------------------------------------------
padj <- list()
for (rep in 1:200) {
  n <- sample(0:12, 1); p <- runif(n)^sample(1:6, 1)
  if (n > 2 && rep %% 3 == 0) p[sample(n, 1)] <- NA
  if (n > 3 && rep %% 4 == 0) p[sample(n, 2)] <- p[1]
  if (n > 3 && rep %% 5 == 0) p[sample(n, 1)] <- 0
  for (m in p.adjust.methods) padj[[length(padj) + 1]] <- list(p = num(p), method = m, out = num(p.adjust(p, m)))
}

# ---- rowMeans(): the accumulator of the platform ----------------------------------------------------------
# D of a gene in ten replicates, as read from the files (five decimals), and their mean. A mean of
# equal values is not always that value, and whether it is depends on the precision of the sum.
rm_in <- matrix(round(sample(0:60, 10 * 4000, TRUE) / sample(c(3, 6, 7, 9, 10, 11, 13), 4000, TRUE), 5), 4000, 10)
rm_in[1:2000, ] <- rm_in[1:2000, 1]                                  # the same value in every replicate
rm_in[sample(length(rm_in), 2000)] <- NA

# ---- a constant seed in a regression (R defect R11) ---------------------------------------------------------
set.seed(1)
n <- 20; g <- 200
x <- matrix(round(rnorm(g * n, 8, 1.5), 2), g, n)
constant <- lapply(c(8.7, 8, 5.48947589862524), function(value) {
  y <- rep(value, n)
  p1 <- sapply(1:g, function(i) coef(summary(glm(y ~ x[i, ], family = "gaussian")))[2, 4])
  p2 <- sapply(1:g, function(i) { xi <- x[i, n:1]; yr <- y[n:1]; coef(summary(glm(yr ~ xi, family = "gaussian")))[2, 4] })
  list(value = value, below_0.01_first_order = sum(p1 < 0.01, na.rm = TRUE), below_0.01_reversed_order = sum(p2 < 0.01, na.rm = TRUE),
    below_0.01_in_both = sum(p1 < 0.01 & p2 < 0.01, na.rm = TRUE), genes = g)
})

info <- list(r_version = R.version.string, platform = R.version$platform, long_double_bytes = .Machine$sizeof.longdouble,
  sample_kinds = chr(kinds), date = format(Sys.time(), "%Y-%m-%d"))
writeLines(to_json(list(
  info = info, rng = rng,
  format = list(value = num(v), text = chr(written)),
  read = list(text = chr(texts), value = num(read_back)),
  round = list(value = num(r), five = num(round(r, 5)), three = num(round(r, 3))),
  cor = cors, padj = padj,
  row_means = list(nrow = nrow(rm_in), ncol = ncol(rm_in), value = num(as.vector(rm_in)), mean = num(rowMeans(rm_in, na.rm = TRUE))),
  constant_seed = constant
)), file.path(out_dir, "primitives.json"))
cat("written", file.path(out_dir, "primitives.json"), "\n")
