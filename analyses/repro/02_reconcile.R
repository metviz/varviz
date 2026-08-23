#!/usr/bin/env Rscript
# Join the manuscript's reported values (00_claims.tsv) to the recomputed ones
# (out/computed.tsv) and flag every disagreement.
#
# Tolerance: percentages must agree to 0.1 (the precision the manuscript reports);
# counts must agree exactly; ratios to 0.1.
#
# Run from the repository root:  Rscript analyses/repro/02_reconcile.R

CLAIMS <- "analyses/repro/00_claims.tsv"
# regenerated.tsv, not computed.tsv: 07_regenerate.R rebuilds every claim from
# the authoritative run in one pass, and its ids match 00_claims.tsv exactly.
# out/computed.tsv was written 2026-08-21 from ps_mds_corroborate, which lost
# KCNQ1 and TP53 to silent UniProt failures -- reconciling against it compares
# the manuscript to contaminated numbers.
COMP   <- "analyses/repro/out/regenerated.tsv"
OUT    <- "analyses/repro/out/reconciliation.tsv"

cl <- read.delim(CLAIMS, sep = "\t", quote = "", stringsAsFactors = FALSE)
cp <- read.delim(COMP,   sep = "\t", quote = "", stringsAsFactors = FALSE)

m <- merge(cl, cp, by = "id", all.x = TRUE, sort = FALSE)
m <- m[match(cl$id, m$id), ]

is_count <- m$reported == round(m$reported) & abs(m$reported) >= 10
tol <- ifelse(is_count, 0.5, 0.1)
m$delta  <- m$value - m$reported
m$status <- ifelse(is.na(m$value), "NOT COMPUTED",
             ifelse(abs(m$delta) <= tol, "match", "MISMATCH"))

fmt <- function(x) ifelse(is.na(x), "-", formatC(x, format = "f", digits = 2))
cat("\n", strrep("=", 96), "\n", sep = "")
cat(sprintf("%-22s %-26s %10s %10s %9s  %s\n",
            "id", "location", "reported", "computed", "delta", "status"))
cat(strrep("-", 96), "\n", sep = "")
for (i in seq_len(nrow(m))) {
  cat(sprintf("%-22s %-26s %10s %10s %9s  %s\n",
              m$id[i], substr(m$location[i], 1, 26),
              fmt(m$reported[i]), fmt(m$value[i]), fmt(m$delta[i]), m$status[i]))
}
cat(strrep("=", 96), "\n", sep = "")
nbad <- sum(m$status != "match")
cat(sprintf("%d of %d claims reconcile; %d need attention\n",
            sum(m$status == "match"), nrow(m), nbad))

write.table(m[, c("id","location","claim","reported","value","delta","status")],
            OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("[repro] wrote %s\n", OUT))
quit(status = if (nbad > 0) 1 else 0)
