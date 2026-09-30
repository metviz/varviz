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

# Read the dual-pass classifications as one row per variant.
#
# The table is keyed to universe rows, and the universe has one row per
# (variant, source, scoreset): 90,701 rows over 54,454 distinct variants on the
# 14-gene run, with 20,944 variants appearing more than once. A variant's
# classification does not depend on which scoreset it arrived from, and the run
# confirms it: of those 20,944, zero disagree on the Pass-Full class, the Pass-Blind
# class or the points.
#
# So a per-variant metric computed without deduplicating weights a variant by how
# many scoresets happen to contain it. That produced two answers for one quantity:
# Panel C's Spearman rho reads 0.847/0.991 deduplicated over 4,905 variants and
# 0.811/0.883 undeduplicated over 15,303 rows, on the same run.
#
# distinct = FALSE is for a caller whose unit of analysis is not the variant. Panel B
# is the case: one variant genuinely has a different DMS score in each scoreset. Its
# scoreset rows come from the universe rather than from here, so deduplicating this
# table does not collapse them and Panel B needs no exception.
#
# Disagreeing duplicates raise. Today none disagree; if a run produces a variant with
# two classifications that is a defect, and keeping the first silently would hide it
# behind a deduplication step.
read_classifications <- function(path, distinct = TRUE, quiet = FALSE,
                                 key = c("gene", "p_notation")) {
  d <- utils::read.delim(path, stringsAsFactors = FALSE, na.strings = character(),
                         check.names = FALSE)
  n_in <- nrow(d)
  if (!all(key %in% names(d))) {
    if (!quiet) message(sprintf("[classifications] %d rows; no %s column, left as read",
                                n_in, paste(setdiff(key, names(d)), collapse = "/")))
    return(d)
  }
  if (!distinct) {
    if (!quiet) message(sprintf("[classifications] %d rows, not deduplicated", n_in))
    return(d)
  }
  k <- do.call(paste, c(unname(d[key]), list(sep = "\037")))
  # Only the columns a classification actually asserts. Provenance columns such as
  # study legitimately differ between duplicate rows and must not count as conflict.
  assert_cols <- intersect(c("varviz_classification_full", "varviz_classification_blind",
                             "varviz_pts_full", "varviz_pts_blind", "pm1_pathway",
                             "tags_full", "tags_blind"), names(d))
  if (length(assert_cols)) {
    sig <- do.call(paste, c(unname(d[assert_cols]), list(sep = "\037")))
    conflict <- unique(k[duplicated(cbind(k, sig)) == FALSE &
                         k %in% k[duplicated(k)]])
    per_key <- tapply(sig, k, function(x) length(unique(x)))
    bad <- names(per_key)[per_key > 1L]
    if (length(bad)) {
      ex <- head(gsub("\037", " ", bad), 3)
      stop(sprintf(paste0("%d variant(s) carry more than one classification in %s ",
                          "(e.g. %s). Deduplicating would hide a real inconsistency."),
                   length(bad), path, paste(ex, collapse = "; ")), call. = FALSE)
    }
  }
  out <- d[!duplicated(k), , drop = FALSE]
  if (!quiet) {
    message(sprintf("[classifications] %d rows in, %d distinct variants out (%d duplicates collapsed)",
                    n_in, nrow(out), n_in - nrow(out)))
  }
  out
}

# ---------------------------------------------------------------------------
# Report provenance
# ---------------------------------------------------------------------------
# A report states numbers without saying where they came from, so
# analyses/derived/varviz_classifications.tsv could sit six months stale behind 21
# newer runs and nothing in panel_a.html would show it. The resolver prints its path
# to stdout, which is gone the moment the script exits.
#
# These build a stamp a reader can check months later: which run, which universe and
# its fingerprint, which commit, and whether the run was selected deliberately or
# fallen back to. Same content the harness writes per gene in <gene>__manifest.json,
# so a report and a run can be matched up.

.short_sha <- function(path, n = 12L) {
  if (is.na(path) || !file.exists(path)) return(NA_character_)
  out <- tryCatch(system2("sha256sum", shQuote(path), stdout = TRUE, stderr = FALSE),
                  error = function(e) NA_character_, warning = function(w) NA_character_)
  if (length(out) == 0L || is.na(out[1])) return(NA_character_)
  substr(sub(" .*$", "", out[1]), 1L, n)
}

.git_state <- function() {
  sha <- tryCatch(system2("git", c("rev-parse", "--short=8", "HEAD"),
                          stdout = TRUE, stderr = FALSE), error = function(e) NA_character_)
  dirty <- tryCatch(length(system2("git", c("status", "--porcelain"),
                                   stdout = TRUE, stderr = FALSE)) > 0L,
                    error = function(e) NA)
  list(sha = if (length(sha)) sha[1] else NA_character_, dirty = dirty)
}

run_provenance_lines <- function(classifications = NULL,
                                 universe = "analyses/derived/variant_universe_gnomad.tsv") {
  cls <- if (is.null(classifications)) resolve_classifications(quiet = TRUE) else classifications
  g <- .git_state()
  # The run directory, not the temp file the resolver consolidates into: a reader
  # needs the run's name, and /tmp/RtmpXXXX/... says nothing.
  run_label <- switch(cls$source,
    REVAL_DIR              = Sys.getenv("REVAL_DIR", cls$path),
    VARVIZ_RUN             = Sys.getenv("VARVIZ_RUN", cls$path),
    VARVIZ_CLASSIFICATIONS = Sys.getenv("VARVIZ_CLASSIFICATIONS", cls$path),
    cls$path)
  ln <- c(
    sprintf("run: %s (via %s)", run_label, cls$source),
    sprintf("universe: %s  sha256:%s", universe, .short_sha(universe)),
    sprintf("code: commit %s%s, R %s", g$sha,
            if (isTRUE(g$dirty)) " (working tree dirty)" else "",
            paste(R.version$major, R.version$minor, sep = ".")),
    sprintf("generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))
  )
  if (isTRUE(cls$stale)) {
    ln <- append(ln, sprintf("STALE: %d run director%s newer than this input: %s",
                             length(cls$newer_runs),
                             if (length(cls$newer_runs) == 1L) "y is" else "ies are",
                             paste(utils::head(cls$newer_runs, 6), collapse = ", ")),
                 after = 1L)
  }
  ln
}

run_provenance_html <- function(...) {
  ln <- run_provenance_lines(...)
  htmltools::tags$footer(
    style = paste("margin-top:2em;padding-top:0.6em;border-top:1px solid #bbb;",
                  "font:11px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace;color:#555"),
    htmltools::tags$strong("Provenance"),
    htmltools::tags$ul(
      style = "margin:0.3em 0 0 1em;padding:0;list-style:none",
      lapply(ln, function(x) htmltools::tags$li(x))
    )
  )
}

run_provenance_md <- function(...) {
  c("", "---", "", "**Provenance**", "",
    paste0("- ", run_provenance_lines(...)))
}
