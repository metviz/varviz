# ps_baseline — dual-pass classification harness against the pre-PS-hotspot engine.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# This run drives the server.R snapshot taken before the PS/CCRS hotspot work
# landed, so its numbers are the baseline the later ablations are read against.
# The pinned snapshot is a working-tree backup, not a tracked file: if it is
# missing, restore it from git history rather than substituting current server.R.
#
# Run from project root: Rscript analyses/ps_baseline_harness.R
# Output: analyses/ps_baseline/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL  <- "ps_baseline"
HARNESS_SERVER <- "server.R.bak.20260729_214027.pre-ps-hotspot"
HARNESS_OUTDIR <- "analyses/ps_baseline"

source("analyses/05_classify_harness.R")
