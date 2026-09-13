# Self-check for the 3-point PP3 rung.
#
# The 2015 guidelines have no strength between Moderate (2) and Strong (4), but
# Pejaver 2022 and Bergquist 2025 both report calibrated 3-point intervals
# because the point-based system is expected to adopt one. REVEL's interval is
# 0.879-0.931 and AlphaMissense's is 0.972-0.989; both previously collapsed into
# the Moderate tag below them.
#
# The rung is gated on varviz.pp3_3pt (default FALSE) so it can be measured
# before it is adopted. With it off, a score inside a 3-point interval must award
# exactly what the engine awarded before: the moderate tag.
#
# Run from project root: Rscript analyses/tests/test_pp3_3pt.R

suppressMessages(source("server.R"))

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

check(identical(unname(ACMG_TAG_PTS[["PP3_moderate_plus"]]), 3),
      "PP3_moderate_plus is worth 3 points")
check(ACMG_TAG_PTS[["PP3_moderate"]] == 2 && ACMG_TAG_PTS[["PP3_strong"]] == 4,
      "it sits between Moderate (2) and Strong (4)")

# The tag must score, not just exist: a tag absent from the points table would
# silently contribute nothing, which is how an invented tier hides.
r <- classify_acmg(c("PP3_moderate_plus", "PM2", "PP2"))
check(r$pts == 5L, "PP3_moderate_plus + PM2 + PP2 = 5 pts, so the tag scores")
check(identical(r$classification, "VUS-High"), "5 pts lands in VUS-High")

# One point below Strong, which is enough to decide a class at the boundary:
# PM1_strong + PM2 + PP2 lands on 6, so the PP3 tier chosen sets 9 vs 10.
check(classify_acmg(c("PM1","PM2","PP2","PP3_moderate_plus"))$pts == 7L,
      "PM1 + PM2 + PP2 + 3-point = 7 pts")
check(classify_acmg(c("PM1","PM2","PP2","PP3_strong"))$pts == 8L,
      "the same tags with PP3_strong = 8 pts, one more")
check(identical(classify_acmg(c("PM1_strong","PM2","PP2","PP3_moderate_plus"))$classification,
                "Likely Pathogenic"),
      "at the boundary the 3-point tier gives 9 pts = Likely Pathogenic")
check(identical(classify_acmg(c("PM1_strong","PM2","PP2","PP3_strong"))$classification,
                "Pathogenic"),
      "the same tags with PP3_strong give 10 pts = Pathogenic")

# The blinding helper strips ClinVar criteria only; PP3 is predictor evidence
# and must survive under either tier.
source("analyses/lib/clinvar_blind.R")
b <- acmg_blind(c("PP3_moderate_plus", "PM2", "PP2"), "ccrs")
check("PP3_moderate_plus" %in% b$tags, "the 3-point tag survives Pass-Blind")

# Default OFF: the rung must not change any call until it is switched on.
check(!isTRUE(getOption("varviz.pp3_3pt", FALSE)),
      "varviz.pp3_3pt defaults to FALSE")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
