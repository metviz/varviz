#!/usr/bin/env Rscript
# Which run on disk reproduces which manuscript claim?
#
# The manuscript reports benchmark statistics from more than one pipeline run.
# This script scores every candidate run against every claim in 00_claims.tsv so
# the authoritative run can be chosen deliberately rather than inferred.
#
# Run from the repository root:  Rscript analyses/repro/03_which_run.R

RUNS <- c(
  derived_mds   = "analyses/derived/varviz_classifications_mds.tsv",
  derived_base  = "analyses/derived/varviz_classifications.tsv",
  ps_reval      = "analyses/ps_reval/summary.tsv",
  ps_baseline   = "analyses/ps_baseline/summary.tsv",
  ps_nomds      = "analyses/ps_nomds/summary.tsv"
)
CLAIMS <- "analyses/repro/00_claims.tsv"
OUT    <- "analyses/repro/out/which_run.tsv"

cl <- read.delim(CLAIMS, sep = "\t", quote = "", stringsAsFactors = FALSE)
tmp <- tempfile()
mat <- matrix("", nrow = nrow(cl), ncol = length(RUNS),
              dimnames = list(cl$id, names(RUNS)))

for (rn in names(RUNS)) {
  p <- RUNS[[rn]]
  if (!file.exists(p)) { mat[, rn] <- "-"; next }
  system2("Rscript", c("analyses/repro/01_recompute.R", p, tmp), stdout = NULL)
  cp <- read.delim(tmp, sep = "\t", stringsAsFactors = FALSE)
  v  <- setNames(cp$value, cp$id)
  for (i in seq_len(nrow(cl))) {
    id <- cl$id[i]; rep <- cl$reported[i]
    if (!id %in% names(v)) { mat[i, rn] <- "."; next }
    tol <- if (rep == round(rep) && abs(rep) >= 10) 0.5 else 0.1
    mat[i, rn] <- if (abs(v[[id]] - rep) <= tol) "MATCH" else "no"
  }
}

w <- max(nchar(rownames(mat)))
cat("\nWhich run reproduces which claim   (MATCH / no / . = not computable for that run)\n\n")
cat(sprintf(paste0("%-", w, "s"), ""), sprintf("%13s", colnames(mat)), "\n")
for (i in seq_len(nrow(mat)))
  cat(sprintf(paste0("%-", w, "s"), rownames(mat)[i]),
      sprintf("%13s", mat[i, ]), "\n")

cat("\nclaims matched per run:\n")
for (rn in colnames(mat))
  cat(sprintf("  %-14s %d / %d\n", rn, sum(mat[, rn] == "MATCH"), nrow(mat)))

write.table(data.frame(id = rownames(mat), mat, check.names = FALSE),
            OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("\n[repro] wrote %s\n", OUT))
