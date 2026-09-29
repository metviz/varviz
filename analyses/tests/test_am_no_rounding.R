# Guards against rounding an AlphaMissense score before comparing it to a
# calibrated cut-point.
#
# Why this exists. The AlphaFold-CSV path did
#
#   am_exact_score <- round(as.numeric(afs_data$am_pathogenicity[idx]), 3)
#
# and every threshold comparison then ran on the rounded value. Rounding to three
# decimals moves anything within half a thousandth of a cut-point up to it, so a
# score genuinely below a tier was scored as meeting it. The tiers affected are
# 0.990, 0.972, 0.906 and 0.792 for PP3 and 0.564 and 0.340 for the class label,
# and the error only ever promotes.
#
# Observed on KCNQ1 p.L342F, whose CSV value is 0.9899: it rounded to 0.990,
# satisfied the >= 0.990 rung and received PP3_strong, a tier the Bergquist
# calibration reserves for 0.990 and above.
#
# Scale at the time of the fix: the CSV path supplied 78,510 of the 90,701
# records in the 14-gene universe (86.6%), and 1,019 of them sat exactly on a
# cut-point after rounding, 622 at 0.990 alone.
#
# The fix keeps the value unrounded for every comparison and rounds only where it
# is displayed. Rounding for presentation is fine; rounding before a threshold
# decision is not.
#
# Run from project root: Rscript analyses/tests/test_am_no_rounding.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src <- readLines(resolve("server.R"), warn = FALSE)

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

cat("-- the defect arithmetic, so the guard states what it is guarding --\n")
# Not a test of server.R: a demonstration that the two orderings disagree, which
# is the whole reason the extraction must not round.
check(isTRUE(round(0.9899, 3) >= 0.990),
      "round(0.9899, 3) >= 0.990 is TRUE  (the promoting behaviour)")
check(isFALSE(0.9899 >= 0.990),
      "0.9899 >= 0.990 is FALSE           (the correct behaviour)")

cat("-- the extraction no longer rounds --\n")
code <- sub("#.*$", "", src)
ext <- grep("am_exact_score <- ", code)
check(length(ext) >= 1L, "am_exact_score is assigned in server.R")
assigns <- code[ext]
check(!any(grepl("round\\(", assigns)),
      "no am_exact_score assignment wraps the value in round()")
check(any(grepl("as\\.numeric\\(afs_data\\$am_pathogenicity", assigns)),
      "the CSV value is still read from afs_data$am_pathogenicity")

cat("-- the class boundaries compare the unrounded value --\n")
# am_exact_class uses 0.564 and 0.340; those are cut-points too and were subject
# to the same promotion.
cls <- grep("am_exact_class <- if", code)
check(length(cls) == 1L, "am_exact_class is assigned once")
if (length(cls) == 1L) {
  blk <- paste(code[cls:(cls + 2L)], collapse = " ")
  check(grepl("am_exact_score >= 0.564", blk, fixed = TRUE),
        "the pathogenic boundary compares am_exact_score directly")
  check(!grepl("round(am_exact_score", blk, fixed = TRUE),
        "the class boundary does not round before comparing")
}

cat("-- display still rounds --\n")
# The export column is where rounding belongs: a 16-digit float in a TSV is
# noise, and nothing downstream makes a decision on the displayed value.
exp_line <- grep("AM_Score\\s*=", code)
check(length(exp_line) >= 1L, "an AM_Score export column exists")
check(any(grepl("round\\(", code[exp_line])),
      "the AM_Score export rounds for display")

cat("-- the scoring path takes the unrounded value --\n")
raw_line <- grep("if \\(!is\\.na\\(raw\\)\\) raw else am_exact_score", code)
check(length(raw_line) == 1L,
      "am_sc falls back to am_exact_score, which is now unrounded")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
