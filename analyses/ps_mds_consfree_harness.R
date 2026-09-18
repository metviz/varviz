# ps_mds_consfree — MDS corroborate + conservation-release build of the harness.
#
# RETIRED 2026-09-06 (engine v2.0.0). varviz.mds_frees_cons is a no-op now that
# MDS corroboration takes max(pathway PM1, MDS) instead of the sum, so this
# configuration is identical to ps_mds_corroborate. Both are superseded by
# ps_final (see analyses/REPRODUCIBILITY.md "Run configurations"). Kept only so
# the historical analyses/ps_mds_consfree/ directory has its provenance.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# Adds to the corroborate build: where MDS corroborates a pathway that had spent
# cons_strong to reach PM1_strong, the conservation flag is released so the
# conservation tier of PP3 is no longer suppressed. PM1 stays Strong on
# pathway(+2) + MDS tier(>=+2); conservation returns to PP3 where it belongs.
#
# Run from project root: Rscript analyses/ps_mds_consfree_harness.R
# Output: analyses/ps_mds_consfree/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL   <- "ps_mds_consfree"
HARNESS_OUTDIR  <- "analyses/ps_mds_consfree"
HARNESS_OPTIONS <- list(
  varviz.mds_tiered     = TRUE,  # LR tiers: +2 / +3 / +4 at MDS <= -4 / -8 / -12
  varviz.mds_pm1        = TRUE,  # MDS evaluated for every variant (originate + corroborate)
  varviz.mds_frees_cons = TRUE   # MDS carries the Strong upgrade; conservation released to PP3
)

source("analyses/05_classify_harness.R")
