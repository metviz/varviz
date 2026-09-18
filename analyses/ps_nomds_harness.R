# ps_nomds — MDS-off ablation of the dual-pass classification harness.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# server.R's PM1 Path 4 (MDS) is disabled, so PM1 can fire only from the
# domain, functional-site, PTM and ClinVar-hotspot pathways. Read against
# ps_baseline this isolates what MDS contributes.
#
# Run from project root: Rscript analyses/ps_nomds_harness.R
# Output: analyses/ps_nomds/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL   <- "ps_nomds"
HARNESS_OUTDIR  <- "analyses/ps_nomds"
HARNESS_OPTIONS <- list(varviz.mds_pm1 = FALSE)

source("analyses/05_classify_harness.R")
