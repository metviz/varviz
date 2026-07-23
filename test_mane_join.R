# Self-check for the MANE join + numbering-consistency guard.
# Run: Rscript test_mane_join.R   (after analyses/16_build_mane_join.R)
#
# Reads the produced artifact (offline) and checks the join landed and the
# numbering_ok flag means what the ACMG matching will trust it to mean.
TSV <- "analyses/derived/gene_mane.tsv"
if (!file.exists(TSV)) {
  cat("skip: run analyses/16_build_mane_join.R first (", TSV, " missing)\n", sep = "")
  quit(status = 0)
}
g <- read.delim(TSV, stringsAsFactors = FALSE)
row <- function(sym) g[g$gene_name == sym, ][1, ]

# ── Known-good genes: UniProt canonical == MANE protein ───────────────────
for (sym in c("GCK", "CASR", "SNCA", "TP53")) {
  r <- row(sym)
  stopifnot(
    !is.na(r$numbering_ok), isTRUE(r$numbering_ok),
    r$uniprot_len == r$mane_prot_len,
    grepl("^NM_", r$mane_refseq_nuc), grepl("^NP_", r$mane_refseq_prot)
  )
}
# GCK anchors to the exact transcript that carried the isoform-collision bug.
stopifnot(row("GCK")$mane_refseq_nuc == "NM_000162.5",
          row("GCK")$uniprot_len == 465)

# ── Known divergences must be flagged, not silently trusted ───────────────
# TTN: UniProt canonical is a shorter titin isoform than MANE.
# NACA: gene_data holds the short NACA accession (215) vs MANE 2078.
# KRAS: 189 vs 188 — the 4A/4B C-terminal difference (conservative flag).
for (sym in c("TTN", "NACA", "KRAS")) {
  r <- row(sym)
  stopifnot(!is.na(r$numbering_ok), isFALSE(as.logical(r$numbering_ok)),
            r$uniprot_len != r$mane_prot_len)
}

# ── Flag semantics hold across the whole table ────────────────────────────
both_known <- !is.na(g$uniprot_len) & !is.na(g$mane_prot_len)
# TRUE/FALSE exactly when both lengths known; NA exactly when one is missing.
stopifnot(
  all(!is.na(g$numbering_ok[both_known])),
  all(is.na(g$numbering_ok[!both_known])),
  # and the flag is the length equality, nothing else
  all(g$numbering_ok[both_known] == (g$uniprot_len[both_known] == g$mane_prot_len[both_known]))
)

# ── Coverage sanity: the vast majority join to a MANE transcript ──────────
# read.delim restores written NA (na="") as "", so test emptiness, not is.na.
cov <- mean(nzchar(g$mane_refseq_nuc) & !is.na(g$mane_refseq_nuc))
stopifnot(cov > 0.95)

cat(sprintf("mane join: all checks pass (%d genes, %.1f%% MANE-matched, %d flagged)\n",
            nrow(g), 100 * cov, sum(g$numbering_ok == FALSE, na.rm = TRUE)))
