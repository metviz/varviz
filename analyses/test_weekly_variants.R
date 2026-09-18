# weekly_test — dual-pass harness over the weekly regression variant set.
#
# Configuration only; the harness itself is analyses/05_classify_harness.R.
# Engine settings are the defaults; only the universe and the output directory
# differ, so week-to-week drift in the live APIs shows up as a diff of
# analyses/weekly_test/summary.tsv rather than of a published run.
#
# Run from project root: Rscript analyses/test_weekly_variants.R
# Output: analyses/weekly_test/classifications/<gene>__dual.tsv + summary.tsv

HARNESS_LABEL    <- "weekly_test"
HARNESS_UNIVERSE <- "analyses/weekly_test/universe.tsv"
HARNESS_OUTDIR   <- "analyses/weekly_test"

source("analyses/05_classify_harness.R")
