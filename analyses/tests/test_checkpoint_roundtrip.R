# Self-check for the per-gene checkpoint round-trip used by every harness.
#
# A concatenation pass re-reads each cached gene's TSV and writes them back out
# as summary.tsv. That path must be idempotent: a gene computed once and then
# concatenated N times must produce identical bytes every time.
#
# It was not. Two readr defaults corrupted the pm1_pathway column, which every
# PM1 tally filters with nchar(x) > 0:
#   1. read_tsv's default na = "NA" turned the empty pm1_pathway of a no-PM1
#      variant into a real NA, which write_tsv then wrote as the string "NA".
#      5,714 cells per benchmark run.
#   2. na = character() alone is not enough: a column blank for EVERY row of a
#      gene (DDR2, which has no PM1 evidence for either of its two variants) is
#      type-guessed as logical, so write_tsv serialised it as "NA" again.
#
# This test writes a checkpoint holding both shapes, reads it with the harness
# idiom, writes it back, and asserts the bytes are unchanged.
#
# Run from project root: Rscript analyses/tests/test_checkpoint_roundtrip.R

suppressPackageStartupMessages(library(readr))

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# The exact read used by run_one() in ps_final_harness.R,
# analyses/05_classify_harness.R and analyses/repro/08_casestudy_harness.R.
read_checkpoint <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, na = character())
  lg <- vapply(df, is.logical, logical(1))
  if (any(lg)) df[lg] <- lapply(df[lg], function(x) rep("", length(x)))
  df
}

tmp  <- tempfile(fileext = ".tsv")
tmp2 <- tempfile(fileext = ".tsv")
on.exit(unlink(c(tmp, tmp2)), add = TRUE)

# A checkpoint with the two failure shapes: pm1_pathway blank on SOME rows,
# and all_blank blank on EVERY row.
ckpt <- data.frame(
  gene         = rep("DDR2", 3),
  p_notation   = c("p.G505S", "p.S461L", "p.R752C"),
  pm1_pathway  = c("", "", "uniprot_domain"),
  all_blank    = c("", "", ""),
  tags_full    = c("PP2, PP3", "PP2, PP3", "PM1, PP2"),
  pts_full     = c(2L, 2L, 3L),
  stringsAsFactors = FALSE
)
write_tsv(ckpt, tmp)
on_disk <- readLines(tmp)

df <- read_checkpoint(tmp)

check(identical(df$pm1_pathway, c("", "", "uniprot_domain")),
      "a partly-blank character column keeps its empty strings")
check(is.character(df$all_blank) && identical(df$all_blank, c("", "", "")),
      "an all-blank column stays character, not logical NA")
check(is.integer(df$pts_full) || is.numeric(df$pts_full),
      "a numeric column keeps its numeric type")
check(sum(nchar(df$pm1_pathway) > 0) == 1L,
      "the nchar() > 0 filter every PM1 tally uses counts 1, not 3")

write_tsv(df, tmp2)
check(identical(readLines(tmp2), on_disk),
      "write(read(checkpoint)) is byte-identical to the checkpoint")
check(!any(grepl("\tNA\t|\tNA$", readLines(tmp2))),
      "no literal NA is introduced anywhere in the round-trip")

# Two round-trips must be as stable as one -- a concat pass can run repeatedly.
write_tsv(read_checkpoint(tmp2), tmp2)
check(identical(readLines(tmp2), on_disk), "a second round-trip is also stable")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
