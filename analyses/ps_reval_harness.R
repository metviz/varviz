# ps_reval — re-evaluation run of the dual-pass classification harness.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# Engine settings are the defaults: this re-runs the current server.R over the
# same universe into its own directory, so a re-classified set can be diffed
# against an earlier one without overwriting it.
#
# Run from project root: Rscript analyses/ps_reval_harness.R
# Output: analyses/ps_reval/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL  <- "ps_reval"
HARNESS_OUTDIR <- "analyses/ps_reval"

source("analyses/05_classify_harness.R")
