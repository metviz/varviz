# ps_mds_corroborate — MDS corroborate build of the dual-pass harness.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# server.R Path 4 now scores MDS for every variant, originating PM1 where no
# other path fired and adding its LR tier to the base strength (capped at
# Strong) where one did. The ClinVar-hotspot upgrade records its pre-upgrade
# tag so strip_clinvar_tags() can demote it under Pass-Blind.
#
# Run from project root: Rscript analyses/ps_mds_corroborate_harness.R
# Output: analyses/ps_mds_corroborate/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL   <- "ps_mds_corroborate"
HARNESS_OUTDIR  <- "analyses/ps_mds_corroborate"
HARNESS_OPTIONS <- list(
  varviz.mds_tiered = TRUE,   # LR tiers: +2 / +3 / +4 at MDS <= -4 / -8 / -12
  varviz.mds_pm1    = TRUE    # MDS evaluated for every variant (originate + corroborate)
)

source("analyses/05_classify_harness.R")
