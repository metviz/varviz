# Self-check for the MANE join + accession anchoring + numbering guard.
# Run: Rscript test_mane_join.R   (after analyses/16_build_mane_join.R)
#
# Reads the produced artifact (offline) and checks that (a) the join landed,
# (b) MANE anchoring picks the numbering-correct UniProt accession, and
# (c) numbering_ok means what the ACMG matching will trust it to mean.
TSV <- "analyses/derived/gene_mane.tsv"
if (!file.exists(TSV)) {
  cat("skip: run analyses/16_build_mane_join.R first (", TSV, " missing)\n", sep = "")
  quit(status = 0)
}
g <- read.delim(TSV, stringsAsFactors = FALSE)
row <- function(sym) g[g$gene_name == sym, ][1, ]

# ── Known-good genes: gene_data accession already IS the MANE entry ────────
for (sym in c("GCK", "CASR", "SNCA", "TP53")) {
  r <- row(sym)
  stopifnot(
    isTRUE(r$accession_agrees), r$uniprot_id == r$mane_uniprot_id,
    isTRUE(r$numbering_ok), r$mane_uniprot_len == r$mane_prot_len,
    grepl("^NM_", r$mane_refseq_nuc), grepl("^NP_", r$mane_refseq_prot)
  )
}
stopifnot(row("GCK")$mane_refseq_nuc == "NM_000162.5",
          row("GCK")$mane_uniprot_len == 465)

# ── Anchoring CORRECTS a wrong stored accession ───────────────────────────
# gene_data held the short NACA isoform (Q13765/215); the MANE transcript ties
# to E9PAV3/2078, which matches the MANE protein — so anchoring makes NACA safe.
naca <- row("NACA")
stopifnot(
  isFALSE(as.logical(naca$accession_agrees)),
  naca$uniprot_id == "Q13765", naca$mane_uniprot_id == "E9PAV3",
  naca$mane_uniprot_len == 2078, naca$mane_prot_len == 2078,
  isTRUE(naca$numbering_ok)               # fixed once the right accession is used
)
# Every correction row must genuinely disagree and name a different accession.
corr <- g[which(g$accession_agrees == FALSE), ]
stopifnot(nrow(corr) >= 10,
          all(corr$uniprot_id != corr$mane_uniprot_id | is.na(corr$uniprot_id)))

# ── Residual mismatch is TRUE isoform divergence, not a wrong accession ────
# TTN: gene_data accession already agrees, but UniProt canonical (34350) is a
# shorter titin isoform than MANE (35991) — only alignment fixes this.
ttn <- row("TTN")
stopifnot(isTRUE(ttn$accession_agrees),
          isFALSE(as.logical(ttn$numbering_ok)),
          ttn$mane_uniprot_len != ttn$mane_prot_len)
# KRAS: 189 vs 188 (4A/4B). Flagged, conservative — the G12/G13/Q61 hotspots
# are N-terminal and unaffected, but the length guard still fires.
kras <- row("KRAS")
stopifnot(isFALSE(as.logical(kras$numbering_ok)), kras$mane_uniprot_len == 189)

# ── Flag semantics hold across the whole table ────────────────────────────
both_known <- !is.na(g$mane_uniprot_len) & !is.na(g$mane_prot_len)
stopifnot(
  all(!is.na(g$numbering_ok[both_known])),
  all(is.na(g$numbering_ok[!both_known])),
  all(g$numbering_ok[both_known] ==
        (g$mane_uniprot_len[both_known] == g$mane_prot_len[both_known]))
)

# ── Coverage: read.delim restores written NA ("") as "", so test emptiness ──
cov_mane   <- mean(nzchar(g$mane_refseq_nuc) & !is.na(g$mane_refseq_nuc))
cov_anchor <- mean(nzchar(g$mane_uniprot_id) & !is.na(g$mane_uniprot_id))
stopifnot(cov_mane > 0.95, cov_anchor > 0.90)

cat(sprintf("mane join: all checks pass (%d genes, %.1f%% anchored, %d corrected, %d flagged)\n",
            nrow(g), 100 * cov_anchor,
            sum(g$accession_agrees == FALSE, na.rm = TRUE),
            sum(g$numbering_ok == FALSE, na.rm = TRUE)))
