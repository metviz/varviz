# ps_tiered — tiered-MDS build of the dual-pass classification harness.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# Turns on the second MDS tier: PM1_moderate_plus (+3) at MDS <= -8, where LR+
# is ~11 on the ClinGen/Pejaver scale — between Moderate (4.3) and Strong (18.7).
#
# Run from project root: Rscript analyses/ps_tiered_harness.R
# Output: analyses/ps_tiered/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL   <- "ps_tiered"
HARNESS_OUTDIR  <- "analyses/ps_tiered"
HARNESS_OPTIONS <- list(varviz.mds_tiered = TRUE)

source("analyses/05_classify_harness.R")
