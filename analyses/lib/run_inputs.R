# run_inputs.R — pick the run whose classifications an analysis script reads, and
# say so out loud.
#
# Why this exists. analyses/derived/varviz_classifications.tsv sat six months stale,
# dated 2026-05-02, while the authoritative run was 2026-09-28. Its row count still
# matched exactly at 90,701, so the usual check passed; the contents did not match.
# Six scripts read that path with no way to point them at a run, three of which
# write the published Panel A, Panel B and Panel C reports.
#
# Those scripts did print the path they read, so the failure was not silence. It was
# that nothing compared the file against the runs sitting next to it. A path plus a
# date is only useful if something says the date is wrong.
#
# The default is kept, so every documented invocation still works. What changes is
# that a stale fallback announces itself, and that a run can be selected without
# editing the script.
#
#   VARVIZ_CLASSIFICATIONS  explicit path to a consolidated classifications TSV
#   REVAL_DIR               a run's classifications/ directory of <gene>__dual.tsv
#   VARVIZ_RUN              a run's summary.tsv; its sibling classifications/ is used
#
# In that order of precedence, most specific first.

# A run directory holds one checkpoint per gene, but every caller reads a single
# path, so consolidate to one temp file here rather than in each script. na =
# character() on the way in and na = "" on the way out: read.delim's default turns an
# empty pm1_pathway into NA and writing it back produces the literal "NA", which
# breaks the nchar() > 0 filter every PM1 tally uses.
.consolidate_checkpoints <- function(dir) {
  f <- list.files(dir, pattern = "__dual\\.tsv$", full.names = TRUE)
  if (!length(f)) return(NA_character_)
  parts <- lapply(f, function(p)
    utils::read.delim(p, stringsAsFactors = FALSE, na.strings = character(),
                      colClasses = "character", check.names = FALSE))
  keep <- Reduce(intersect, lapply(parts, names))
  if (!length(keep)) return(NA_character_)
  all <- do.call(rbind, lapply(parts, function(d) d[, keep, drop = FALSE]))
  out <- tempfile(pattern = "varviz_classifications_", fileext = ".tsv")
  utils::write.table(all, out, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  out
}

# Run directories newer than the file actually being read. Named rather than
# counted, because "something newer exists" is not actionable and "ps_am792_v123_final
# is newer" is.
.newer_runs <- function(path, run_glob) {
  if (is.na(path) || !file.exists(path)) return(character(0))
  mt <- file.info(path)$mtime
  cands <- Sys.glob(run_glob)
  cands <- cands[dir.exists(cands)]
  if (!length(cands)) return(character(0))
  stamp <- function(d) {
    inner <- c(file.path(d, "summary.tsv"), file.path(d, "classifications"))
    inner <- inner[file.exists(inner)]
    if (!length(inner)) return(as.POSIXct(NA))
    max(file.info(inner)$mtime, na.rm = TRUE)
  }
  ts <- vapply(cands, function(d) as.numeric(stamp(d)), numeric(1))
  ok <- !is.na(ts) & ts > as.numeric(mt)
  basename(cands[ok])
}

resolve_checkpoint_dir <- function(default = "analyses/classifications") {
  rd <- Sys.getenv("REVAL_DIR", "")
  if (nzchar(rd)) return(rd)
  vr <- Sys.getenv("VARVIZ_RUN", "")
  if (nzchar(vr)) return(file.path(dirname(vr), "classifications"))
  default
}

resolve_classifications <- function(default = "analyses/derived/varviz_classifications.tsv",
                                    tag = "run", quiet = FALSE,
                                    run_glob = "analyses/ps_*") {
  path <- NA_character_; src <- NA_character_
  explicit <- Sys.getenv("VARVIZ_CLASSIFICATIONS", "")
  rd       <- Sys.getenv("REVAL_DIR", "")
  vr       <- Sys.getenv("VARVIZ_RUN", "")
  if (nzchar(explicit)) {
    path <- explicit; src <- "VARVIZ_CLASSIFICATIONS"
  } else if (nzchar(rd)) {
    path <- .consolidate_checkpoints(rd); src <- "REVAL_DIR"
  } else if (nzchar(vr)) {
    path <- .consolidate_checkpoints(file.path(dirname(vr), "classifications"))
    if (is.na(path) && file.exists(vr)) path <- vr
    src <- "VARVIZ_RUN"
  } else {
    path <- default; src <- "default"
  }
  newer <- if (identical(src, "default")) .newer_runs(path, run_glob) else character(0)
  stale <- length(newer) > 0L
  if (!quiet) {
    mt <- if (!is.na(path) && file.exists(path))
      format(file.info(path)$mtime, "%Y-%m-%d %H:%M") else "missing"
    message(sprintf("[%s] classifications: %s  (%s, via %s)", tag, path, mt, src))
    if (stale) {
      message(sprintf("[%s] STALE: %d run director%s newer than this file: %s",
                      tag, length(newer), if (length(newer) == 1L) "y is" else "ies are",
                      paste(newer, collapse = ", ")))
      message(sprintf("[%s]   set REVAL_DIR=analyses/<run>/classifications to use one",
                      tag))
    }
  }
  list(path = path, source = src, stale = stale, newer_runs = newer)
}
