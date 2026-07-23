#!/usr/bin/env Rscript
# 17_merge_mane_into_rdata.R — fold the MANE columns into data/VarViz.RData.
#
# Adds columns only; the app's fetch key stays uniprot_id, so this is inert to
# the running app until server.R is switched to prefer mane_uniprot_id. That
# switch is deliberately NOT done here — the 15 accession corrections are meant
# to be reviewed by hand first (see the TSV written below).
#
# Backs up VarViz.RData before writing, and preserves ccrs_data untouched.
#
# Usage:  Rscript analyses/17_merge_mane_into_rdata.R

RDATA <- "data/VarViz.RData"
AUG   <- "analyses/derived/gene_data_with_mane.rds"
stopifnot(file.exists(RDATA), file.exists(AUG))

# ── Back up first ─────────────────────────────────────────────────────────
stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
bak <- paste0(RDATA, ".bak.", stamp, ".pre-mane")
file.copy(RDATA, bak, overwrite = FALSE)
message("[merge] backed up ", RDATA, " -> ", bak)

# ── Load, swap gene_data, keep everything else ────────────────────────────
load(RDATA)                                     # gene_data, ccrs_data
orig <- gene_data
aug  <- readRDS(AUG)

# Safety: same genes, same order, original columns unchanged.
stopifnot(
  nrow(aug) == nrow(orig),
  identical(aug$gene_name, orig$gene_name),
  identical(aug$uniprot_id, orig$uniprot_id),
  all(colnames(orig) %in% colnames(aug))
)
new_cols <- setdiff(colnames(aug), colnames(orig))
message("[merge] adding columns: ", paste(new_cols, collapse = ", "))

gene_data <- aug
save(gene_data, ccrs_data, file = RDATA)        # ccrs_data preserved verbatim
message("[merge] wrote ", RDATA, " (", round(file.size(RDATA)/1e6, 1), " MB)")

# ── Verify the round-trip ─────────────────────────────────────────────────
chk <- new.env(); load(RDATA, envir = chk)
stopifnot(
  identical(nrow(chk$ccrs_data), nrow(ccrs_data)),   # untouched
  all(new_cols %in% colnames(chk$gene_data)),
  nrow(chk$gene_data) == nrow(orig)
)
message("[merge] verified: ccrs_data intact, gene_data now ",
        ncol(chk$gene_data), " columns")

# ── Review TSVs at repo root (non-ignored, easy to open) ──────────────────
cols <- c("gene_name", "uniprot_id", "uniprot_len",
          "mane_uniprot_id", "mane_uniprot_len", "accession_agrees",
          "mane_refseq_nuc", "mane_refseq_prot", "mane_ensembl",
          "mane_prot_len", "numbering_ok")
corr <- gene_data[which(gene_data$accession_agrees == FALSE), cols]
corr <- corr[order(abs(corr$uniprot_len - corr$mane_uniprot_len), decreasing = TRUE), ]
f1 <- paste0("VarViz_MANE_accession_corrections_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write.table(corr, f1, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

flag <- gene_data[which(gene_data$numbering_ok == FALSE), cols]
flag <- flag[order(abs(flag$mane_uniprot_len - flag$mane_prot_len), decreasing = TRUE), ]
f2 <- paste0("VarViz_MANE_numbering_flags_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write.table(flag, f2, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

message(sprintf("[merge] review tables: %s (%d rows), %s (%d rows)",
                f1, nrow(corr), f2, nrow(flag)))
