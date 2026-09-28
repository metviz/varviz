# Guards the 1000 Genomes and ExAC allele frequencies coming from the local
# dbNSFP row rather than from a network call that never happened.
#
# Why this exists. Pop_1KG and Pop_ExAC were empty for every variant of every run.
# server.R initialised them to NA, logged
#
#   [dbNSFP] Population fields absent from dbnsfp — will fetch via rsID
#
# and then called fetch_pop_by_rsid(), which is defined nowhere in the codebase.
# The call sits inside tryCatch(..., error = function(e) NULL), so it failed
# silently on every variant and the block that would have filled the columns
# never ran. The log line described a fetch that could not occur.
#
# The fix needs no network. dbNSFP 4.9a already carries both frequencies, and the
# harness preloads the full 458-column row per variant; the positional mapper in
# analyses/lib/local_predictors.R simply did not expose them. Verified against the
# real chr1 header: column 213 is 1000Gp3_AF and column 235 is ExAC_AF.
#
# Sourcing them locally also keeps the run reproducible, which a per-variant rsID
# lookup would not: the columns would then depend on whatever a remote service
# returned on the day.
#
# Run from project root: Rscript analyses/tests/test_dbnsfp_population.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
suppressMessages(source(resolve("analyses/lib/local_predictors.R"), local = FALSE))

fails <- 0L
check <- function(ok, msg) {
  # isTRUE: a missing list entry yields NULL and is.na(NULL) is logical(0), which
  # would abort the run instead of reporting the failure.
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# A dbNSFP row in its real shape: 458 columns, "." for everything absent. Same
# synthetic-row approach as analyses/tests/test_dbnsfp_aa_lookup.R, so no dbNSFP
# file and no network are needed.
N_COLS <- 458L
COL_1000G <- 213L
COL_EXAC  <- 235L
mkrow <- function(...) {
  p <- rep(".", N_COLS)
  kv <- list(...)
  for (nm in names(kv)) p[as.integer(nm)] <- kv[[nm]]
  p
}

check(exists(".dbnsfp_row_to_hit"),
      ".dbnsfp_row_to_hit is available from analyses/lib/local_predictors.R")

cat("-- the frequencies are exposed at the verified column positions --\n")
hit <- .dbnsfp_row_to_hit(mkrow("213" = "0.0123", "235" = "0.00456"))
d <- hit$dbnsfp
check(!is.null(d$onekg), "the hit exposes a onekg entry")
check(!is.null(d$exac),  "the hit exposes an exac entry")
check(isTRUE(all.equal(d$onekg$af, 0.0123)),
      "1000Gp3_AF is read from column 213")
check(isTRUE(all.equal(d$exac$af, 0.00456)),
      "ExAC_AF is read from column 235")

cat("-- absent values are NA, not zero --\n")
# The distinction matters: 0 means observed-and-absent, NA means not-measured,
# and a frequency criterion must not read one as the other.
d0 <- .dbnsfp_row_to_hit(mkrow())$dbnsfp
check(is.na(d0$onekg$af), "a '.' 1000Gp3_AF gives NA")
check(is.na(d0$exac$af),  "a '.' ExAC_AF gives NA")
de <- .dbnsfp_row_to_hit(mkrow("213" = "", "235" = ""))$dbnsfp
check(is.na(de$onekg$af), "an empty 1000Gp3_AF gives NA")
check(is.na(de$exac$af),  "an empty ExAC_AF gives NA")
dz <- .dbnsfp_row_to_hit(mkrow("213" = "0", "235" = "0"))$dbnsfp
check(identical(dz$onekg$af, 0), "a genuine zero stays 0 and does not become NA")

cat("-- multi-transcript rows take the first value, as every other field does --\n")
dm <- .dbnsfp_row_to_hit(mkrow("213" = "0.5;0.9", "235" = ".;0.25"))$dbnsfp
check(isTRUE(all.equal(dm$onekg$af, 0.5)), "a ';' list takes the first 1000Gp3_AF value")
check(isTRUE(all.equal(dm$exac$af, 0.25)), "a leading '.' is skipped for ExAC_AF")

cat("-- the undefined-function path is gone from server.R --\n")
src <- readLines(resolve("server.R"), warn = FALSE)
# Strip comments before looking: the removal is documented in a comment that
# names the function, and a historical note is not a call. What must not exist is
# executable code invoking a function that is defined nowhere.
code <- sub("#.*$", "", src)
check(!any(grepl("fetch_pop_by_rsid", code)),
      "no executable line calls fetch_pop_by_rsid, which was never defined")
check(!any(grepl("Population fields absent from dbnsfp", src, fixed = TRUE)),
      "the log line claiming a fetch that never happened is gone")
check(any(grepl('safe_extract_num\\(d, "onekg", "af"\\)', src)),
      "server.R reads the 1000 Genomes frequency from the dbNSFP hit")
check(any(grepl('safe_extract_num\\(d, "exac", "af"\\)', src)),
      "server.R reads the ExAC frequency from the dbNSFP hit")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
