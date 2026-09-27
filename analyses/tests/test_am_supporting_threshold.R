# Guards the AlphaMissense supporting rung and the marker the boundary cap needs.
#
# The supporting rung used to sit at AlphaMissense's developer threshold of
# 0.564. The Bergquist 2025 ClinGen SVI calibration places 0.564 inside the
# indeterminate interval and begins supporting evidence at 0.792, so the engine
# was awarding a point at a score its own cited calibration treats as
# uninformative. "No reliable evidence" is an acceptable result, so the rung
# moved to the calibrated 0.792.
#
# Why this is a code change and not a flag change. Setting
# varviz.am_calibrated = TRUE with varviz.am_hybrid = FALSE also moves the
# supporting rung to 0.792, but it does two further things that were not wanted:
#
#   1. It lets AlphaMissense award PP3_strong at 0.990 from a standing start.
#      Under the hybrid ladder a variant with no other PP3 can only reach
#      supporting, so the calibrated branch is MORE permissive at the top of the
#      scale, not less.
#   2. It silently disables the Pathogenic-boundary cap. `.am_raised_from` is
#      assigned in exactly one place, inside the hybrid branch, and
#      cap_am_pathogenic() returns its input untouched when that marker is NA.
#      Routing through the calibrated branch therefore leaves the cap inert
#      while every log line and option still claims it is enabled.
#
# Assertion 4 below is the regression guard for (2). It is the one that matters:
# the marker assignment must stay inside the hybrid branch, because nothing else
# in the file sets it.
#
# Asserts on the shipped source rather than on a re-expression of the rule.
# Re-typing the comparisons here would only prove a private copy is
# self-consistent, which is the failure mode of commit 4c93a92.
#
# Run from project root: Rscript analyses/tests/test_am_supporting_threshold.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src <- readLines(resolve("server.R"), warn = FALSE)

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# --- the module-level constant ----------------------------------------------
const <- grep("^AM_PP3_SUPPORTING\\s*<-", src)
check(length(const) == 1L,
      sprintf("AM_PP3_SUPPORTING is defined exactly once at module level (found %d)",
              length(const)))
if (length(const) == 1L) {
  val <- suppressWarnings(as.numeric(
    sub("^AM_PP3_SUPPORTING\\s*<-\\s*([0-9.]+).*$", "\\1", src[const])))
  check(isTRUE(all.equal(val, 0.792)),
        sprintf("AM_PP3_SUPPORTING is the calibrated 0.792 (found %s)", format(val)))
}

# --- the hybrid branch uses it, and no longer hardcodes 0.564 ---------------
beg <- grep("^\\s*if \\(\\.am_hybrid\\) \\{", src)
check(length(beg) == 1L, "the hybrid branch appears exactly once")

if (length(beg) == 1L) {
  # The branch ends where the calibrated branch begins.
  end <- beg - 1L + grep("\\}\\s*else if \\(\\.am_cal\\)", src[beg:length(src)])[1]
  branch <- src[beg:end]

  check(any(grepl("AM_PP3_SUPPORTING", branch)),
        "the hybrid supporting rung reads AM_PP3_SUPPORTING")

  check(!any(grepl("am_sc >= 0\\.564", branch)),
        "the hybrid branch no longer compares am_sc against a hardcoded 0.564")

  # THE REGRESSION GUARD: the boundary cap is inert without this assignment,
  # and this branch is the only place in the file that performs it.
  check(any(grepl("\\.am_raised_from <- \\.lvl0", branch)),
        "the hybrid branch still records .am_raised_from, so the Pathogenic cap can fire")

  # The raise rungs are unchanged; only the supporting rung moved.
  for (thr in c("0\\.990", "0\\.972", "0\\.906")) {
    check(any(grepl(paste0("am_sc >= ", thr), branch)),
          sprintf("the raise rung at %s is unchanged", gsub("\\\\", "", thr)))
  }
}

# --- the marker is still assigned in exactly one place ----------------------
assigns <- grep("\\.am_raised_from <- \\.lvl0", src)
check(length(assigns) == 1L,
      sprintf(".am_raised_from is assigned in exactly one place (found %d)", length(assigns)))

# --- the last-resort branch keeps the developer threshold -------------------
# The `else` branch is the documented pre-1.2.1 behaviour and is not part of
# this change, so it should still carry 0.564.
check(any(grepl("am_sc >= 0\\.564 && pp3_level\\(acmg_tags\\) < 1L", src)),
      "the pre-1.2.1 last-resort branch still uses the developer threshold 0.564")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
