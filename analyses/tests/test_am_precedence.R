# Self-check for AlphaMissense PP3 precedence.
#
# BRCA1 p.C27W is an override case. REVEL sets PP3_moderate; the
# calibrated AlphaMissense ladder raises it to PP3_strong, which is the +2 that
# carries the variant from Likely Pathogenic (8) to Pathogenic (10).
#
# Run from project root: Rscript analyses/tests/test_am_precedence.R
#
# Deferential mode must leave REVEL's tag alone and reproduce the release call,
# while still using calibrated thresholds. If it does not, the 2h run measures
# nothing.
suppressMessages(source("server.R"))
cat("option plumbing:\n")
for (cal in c(FALSE, TRUE)) for (def in c(FALSE, TRUE)) {
  options(varviz.am_calibrated = cal, varviz.am_no_override = def)
  cat(sprintf("  am_calibrated=%-5s am_no_override=%-5s -> reads %s / %s\n", cal, def,
      getOption("varviz.am_calibrated"), getOption("varviz.am_no_override")))
}
# Replicate the ladder exactly as server.R runs it, over the three modes.
ladder <- function(am_sc, revel_lvl, cal, def, hyb = FALSE) {
  lvl <- revel_lvl
  pl  <- function() lvl
  if (hyb) {
    if (pl() == 0L) {
      if (am_sc >= 0.564) lvl <- 1L
    } else {
      if      (am_sc >= 0.990 && pl() < 4L) lvl <- 4L
      else if (am_sc >= 0.972 && pl() < 3L) lvl <- 3L
      else if (am_sc >= 0.906 && pl() < 2L) lvl <- 2L
    }
  } else if (cal) {
    if (!def || pl() == 0L) {
      if      (am_sc >= 0.990 && pl() < 4L) lvl <- 4L
      else if (am_sc >= 0.972 && pl() < 3L) lvl <- 3L
      else if (am_sc >= 0.906 && pl() < 2L) lvl <- 2L
      else if (am_sc >= 0.792 && pl() < 1L) lvl <- 1L
    }
  } else if (am_sc >= 0.564 && pl() < 1L) lvl <- 1L
  lvl
}
fails <- 0L
chk <- function(ok,msg){cat(sprintf("  [%s] %s\n", if(ok)"ok" else "FAIL", msg)); if(!ok) fails<<-fails+1L}
cat("\nladder, AM=0.995 (strong band):\n")
chk(ladder(0.995, 2L, TRUE,  FALSE) == 4L, "override ON : REVEL moderate(2) is raised to strong(4)")
chk(ladder(0.995, 2L, TRUE,  TRUE)  == 2L, "override OFF: REVEL moderate(2) is left alone")
chk(ladder(0.995, 0L, TRUE,  TRUE)  == 4L, "override OFF: with no other PP3, AM still speaks -> strong(4)")
chk(ladder(0.995, 2L, FALSE, FALSE) == 2L, "uncalibrated: AM never overrides, stays moderate(2)")
cat("\nladder, AM=0.70 (below calibrated supporting, above developer 0.564):\n")
chk(ladder(0.70, 0L, TRUE,  TRUE)  == 0L, "calibrated: 0.70 is indeterminate, no tag")
chk(ladder(0.70, 0L, FALSE, FALSE) == 1L, "uncalibrated: 0.70 clears 0.564, supporting")
cat("\nhybrid: developer threshold at the supporting rung, calibrated rungs above:\n")
chk(ladder(0.995, 2L, FALSE, FALSE, TRUE) == 4L,
    "AM 0.995 raises REVEL moderate(2) to strong(4) - the productive case")
chk(ladder(0.70,  0L, FALSE, FALSE, TRUE) == 1L,
    "AM 0.70 with no other PP3 still gives supporting, as the release does")
chk(ladder(0.50,  0L, FALSE, FALSE, TRUE) == 0L,
    "AM 0.50 is below the developer threshold, no tag")
chk(ladder(0.995, 0L, FALSE, FALSE, TRUE) == 1L,
    "AM 0.995 alone gives only supporting - the unproductive case is NOT taken")
chk(ladder(0.92,  1L, FALSE, FALSE, TRUE) == 2L,
    "AM 0.92 raises a supporting tag to moderate")
chk(ladder(0.85,  2L, FALSE, FALSE, TRUE) == 2L,
    "AM 0.85 is below the moderate rung, leaves REVEL moderate alone")

cat(sprintf("\n%s: %d failure(s)\n", if (fails==0L) "PASS" else "FAIL", fails))
quit(status = if (fails==0L) 0L else 1L)
