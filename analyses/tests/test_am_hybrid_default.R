# Guards a latent configuration bug: the two places that read
# getOption("varviz.am_hybrid") disagreed on the default.
#
#   server.R:5881  .am_hybrid <- isTRUE(getOption("varviz.am_hybrid", TRUE))
#   server.R:5972  ... !isTRUE(getOption("varviz.am_hybrid", FALSE)) ...
#
# With varviz.am_calibrated = TRUE and varviz.am_hybrid left unset, the PP3 side
# read the default as TRUE and ran the hybrid ladder (supporting at the developer
# threshold 0.564, ignoring am_calibrated entirely), while the BP4 side read the
# same unset option as FALSE and moved the benign cut to the calibrated 0.169.
# That is a run scoring PP3 off one calibration and BP4 off another.
#
# No shipped run hit it, because analyses/ps_final_harness.R sets both flags
# explicitly. It becomes reachable the moment a caller sets only am_calibrated.
#
# This asserts on the shipped source rather than on a re-expression of the rule:
# every getOption("varviz.am_hybrid", <default>) in server.R must carry the same
# default literal. Re-typing the two expressions here and comparing them would
# only prove that a private copy is self-consistent (cf. commit 4c93a92, where a
# test re-implemented the Pathogenic cap and passed while the shipped rule was
# dead code).
#
# Run from project root: Rscript analyses/tests/test_am_hybrid_default.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src_path <- resolve("server.R")

src <- readLines(src_path, warn = FALSE)

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# Every read of the option, with its default argument captured.
pat  <- 'getOption\\(\\s*"varviz\\.am_hybrid"\\s*,\\s*([A-Za-z]+)\\s*\\)'
hits <- regmatches(src, gregexpr(pat, src))
lines_with <- which(lengths(hits) > 0)

defaults <- unlist(lapply(src[lines_with], function(l) {
  m <- regmatches(l, gregexpr(pat, l))[[1]]
  sub(paste0("^.*", "varviz\\.am_hybrid\"\\s*,\\s*"), "", sub("\\s*\\)$", "", m))
}))

cat(sprintf("varviz.am_hybrid read at %d line(s): %s\n",
            length(lines_with), paste(lines_with, collapse = ", ")))
cat(sprintf("defaults found: %s\n", paste(defaults, collapse = ", ")))

check(length(lines_with) >= 2,
      "the option is read in at least two places (otherwise this guard is vacuous)")

check(length(unique(defaults)) == 1L,
      sprintf("all reads of varviz.am_hybrid use one default (found: %s)",
              paste(unique(defaults), collapse = " vs ")))

check(all(defaults == "TRUE"),
      "the shared default is TRUE, matching the documented post-1.2.1 behaviour")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
