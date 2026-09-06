# ps_tiered — tiered-MDS build of the dual-pass classification harness.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# Turns on the exploratory MDS Strong tier: PM1_strong (+4) at MDS <= -12
# (LR+ 39.0 genome-wide, 32.1 on 1-star held-out). Everything from -4 to -12
# stays Moderate; the former "Moderate-plus" (+3) tier was retired 2026-09-06.
#
# Run from project root: Rscript analyses/ps_tiered_harness.R
# Output: analyses/ps_tiered/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL   <- "ps_tiered"
HARNESS_OUTDIR  <- "analyses/ps_tiered"
HARNESS_OPTIONS <- list(varviz.mds_tiered = TRUE)

source("analyses/05_classify_harness.R")
