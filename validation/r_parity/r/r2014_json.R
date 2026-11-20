# JSON output in base R, for the scripts that must also run in the R of 2014 (no package).
# Lists become objects (with names) or arrays (without); num() and chr() render whole vectors.

raw_json <- function(text) structure(text, class = "json_raw")
json_string <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\n", s, fixed = TRUE)
  s <- gsub("\t", "\\t", s, fixed = TRUE)
  s <- gsub("\r", "\\r", s, fixed = TRUE)
  paste0("\"", s, "\"")
}
to_json <- function(v) {
  if (is.null(v)) return("null")
  if (inherits(v, "json_raw")) return(unclass(v))
  if (is.list(v)) {
    if (length(v) == 0) return(if (is.null(names(v))) "[]" else "{}")
    parts <- vapply(v, to_json, character(1), USE.NAMES = FALSE)
    if (!is.null(names(v))) return(paste0("{", paste(paste0(json_string(names(v)), ":", parts), collapse = ","), "}"))
    return(paste0("[", paste(parts, collapse = ","), "]"))
  }
  if (length(v) != 1) stop("to_json: a vector outside a list")
  if (is.na(v)) return("null")
  if (is.character(v)) return(json_string(v))
  if (is.logical(v)) return(if (v) "true" else "false")
  sprintf("%.17g", v)
}
# numbers as strings of 17 significant digits, enough to read back the same double;
# NA and NaN kept apart
num <- function(v) {
  v <- as.numeric(v)
  if (length(v) == 0) return(raw_json("[]"))
  raw_json(paste0("[", paste0("\"", ifelse(is.na(v) & !is.nan(v), "NA", sprintf("%.17g", v)), "\"", collapse = ","), "]"))
}
# strings, a missing one as null; NULL (no names at all) as null
chr <- function(v) {
  if (is.null(v)) return(NULL)
  v <- as.character(v)
  if (length(v) == 0) return(raw_json("[]"))
  raw_json(paste0("[", paste(ifelse(is.na(v), "null", json_string(v)), collapse = ","), "]"))
}
script_dir <- function() {
  arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  dirname(normalizePath(sub("^--file=", "", arg[1])))
}
