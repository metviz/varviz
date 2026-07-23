#!/usr/bin/env Rscript
# 16_build_mane_join.R — attach MANE Select transcripts to gene_data and flag
# genes whose UniProt canonical numbering may not match the transcript numbering
# ClinVar / dbNSFP report against.
#
# Why. VarViz matches variants by integer residue number. UniProt-based layers
# (AlphaFold pLDDT, AlphaMissense, UniProt features, ConSurf) count residues on
# the UniProt CANONICAL isoform; ClinVar and dbNSFP count them on a RefSeq
# transcript (MANE Select for most clinically-used genes). When those two
# sequences differ in length, residue N in one is not residue N in the other,
# and PM1 / PS1 / PM5 / PP3 silently line up against the wrong residue with no
# error raised. Storing the MANE transcript alone does NOT fix that — remapping
# needs sequence alignment — but comparing the two lengths ENUMERATES the genes
# at risk, which is the whole point of this guard.
#
# Inputs (staged under analyses/tmp/mane/ by the caller, see header of run):
#   mane_summary.txt.gz : NCBI MANE summary (symbol -> NM_/NP_/ENST)
#   np_len.tsv          : NP_ accession -> protein length, parsed from the MANE
#                         RefSeq protein FASTA
# UniProt canonical lengths are streamed once from the reviewed human proteome.
#
# Outputs (analyses/derived/):
#   gene_mane.tsv            : per-gene join + numbering_ok flag (all genes)
#   gene_data_with_mane.rds  : gene_data + the new columns, for a later merge
#                              into data/VarViz.RData (kept separate on purpose —
#                              this script never clobbers the app's RData)
#   mane_numbering_mismatch.tsv : the at-risk subset (length disagreement)
#
# Usage:  Rscript analyses/16_build_mane_join.R

suppressMessages({library(httr2); library(jsonlite)})

TMP <- "analyses/tmp/mane"
DER <- "analyses/derived"
dir.create(DER, showWarnings = FALSE, recursive = TRUE)
SUMMARY <- file.path(TMP, "mane_summary.txt.gz")
NPLEN   <- file.path(TMP, "np_len.tsv")
stopifnot(file.exists(SUMMARY), file.exists(NPLEN))

message("[mane] loading gene_data")
load("data/VarViz.RData")                       # gene_data
gene_data$gene_name  <- as.character(gene_data$gene_name)
gene_data$uniprot_id <- as.character(gene_data$uniprot_id)

# ── MANE summary: one MANE Select row per gene symbol ─────────────────────
message("[mane] reading MANE summary")
mane <- read.delim(gzfile(SUMMARY), stringsAsFactors = FALSE, check.names = FALSE)
colnames(mane)[colnames(mane) == "#NCBI_GeneID"] <- "NCBI_GeneID"
sel <- mane[mane$MANE_status == "MANE Select",
            c("symbol", "RefSeq_nuc", "RefSeq_prot", "Ensembl_nuc", "Ensembl_prot")]
sel <- sel[!duplicated(sel$symbol), ]           # symbol is unique for Select

# ── MANE protein length from the RefSeq FASTA (NP_ -> length) ──────────────
nplen <- read.delim(NPLEN, header = FALSE, col.names = c("np", "len"),
                     stringsAsFactors = FALSE)
np_len <- setNames(nplen$len, nplen$np)
sel$mane_prot_len <- unname(np_len[sel$RefSeq_prot])

# ── UniProt reviewed human: accession, canonical length, MANE-Select ENST ──
# One stream, cached so re-runs are offline. The MANE-Select cross-reference is
# what lets us ANCHOR gene->UniProt on the transcript instead of trusting the
# accession already in gene_data (which is sometimes a wrong/short isoform, e.g.
# NACA held Q13765/215aa where the MANE entry is E9PAV3/2078aa).
UP_CACHE <- file.path(TMP, "uniprot_mane.tsv")
if (file.exists(UP_CACHE)) {
  message("[mane] uniprot stream from cache")
  upl <- read.delim(UP_CACHE, stringsAsFactors = FALSE)
} else {
  message("[mane] streaming UniProt reviewed human (accession,length,MANE) (one-time)")
  rows <- list(); cursor <- NULL; page <- 0L
  repeat {
    req <- request("https://rest.uniprot.org/uniprotkb/search") |>
      req_url_query(query = "reviewed:true AND organism_id:9606",
                    fields = "accession,length,xref_mane-select",
                    format = "tsv", size = 500)
    if (!is.null(cursor)) req <- req_url_query(req, cursor = cursor)
    resp <- req_perform(req)
    df  <- read.delim(text = resp_body_string(resp), stringsAsFactors = FALSE)
    if (nrow(df) == 0) break
    rows[[length(rows) + 1L]] <- df
    page <- page + 1L
    link <- resp_headers(resp)[["link"]]   # Link: rel="next" carries the cursor
    if (is.null(link) || !grepl("rel=\"next\"", link)) break
    cursor <- sub(".*[?&]cursor=([^&>]+).*", "\\1", link)
    if (page %% 5L == 0L) message("  ", page * 500L, " streamed")
  }
  upl <- do.call(rbind, rows)
  colnames(upl) <- c("uniprot_id", "uniprot_len", "mane_enst")
  write.table(upl, UP_CACHE, sep = "\t", row.names = FALSE, quote = FALSE)
}
up_len <- setNames(as.integer(upl$uniprot_len), upl$uniprot_id)

# ENST(MANE) -> UniProt accession, matched on the versionless base so a version
# skew between UniProt's record and the MANE summary can't break the join.
base_enst <- function(x) sub("\\..*$", "", trimws(sub(";.*$", "", x)))
have_mane <- nzchar(trimws(upl$mane_enst))
enst2acc <- setNames(upl$uniprot_id[have_mane], base_enst(upl$mane_enst[have_mane]))
enst2len <- setNames(as.integer(upl$uniprot_len[have_mane]), base_enst(upl$mane_enst[have_mane]))

# ── Join everything onto gene_data ────────────────────────────────────────
m <- match(gene_data$gene_name, sel$symbol)
gene_data$mane_refseq_nuc  <- sel$RefSeq_nuc[m]
gene_data$mane_refseq_prot <- sel$RefSeq_prot[m]
gene_data$mane_ensembl     <- sel$Ensembl_nuc[m]
gene_data$mane_prot_len    <- sel$mane_prot_len[m]

# The UniProt accession UniProt itself ties to this gene's MANE transcript, and
# that entry's canonical length. This is the numbering-correct accession.
me <- base_enst(gene_data$mane_ensembl)
gene_data$mane_uniprot_id  <- unname(enst2acc[me])
gene_data$mane_uniprot_len <- unname(enst2len[me])

# Existing accession's own canonical length, kept for comparison.
gene_data$uniprot_len <- unname(up_len[gene_data$uniprot_id])

# Does the accession already in gene_data agree with the MANE-anchored one?
gene_data$accession_agrees <- with(gene_data,
  ifelse(is.na(mane_uniprot_id) | is.na(uniprot_id), NA,
         uniprot_id == mane_uniprot_id))

# numbering_ok compares the MANE-ANCHORED UniProt entry to the MANE protein.
# TRUE  = both lengths known and equal (residue numbers align by construction)
# FALSE = both known and DIFFERENT (UniProt canonical is a genuinely different
#         isoform than MANE — needs alignment, not just a better accession)
# NA    = a length is missing (no MANE gene, or MANE maps to no reviewed entry)
gene_data$numbering_ok <- with(gene_data,
  ifelse(is.na(mane_uniprot_len) | is.na(mane_prot_len), NA,
         mane_uniprot_len == mane_prot_len))

# ── Write outputs ─────────────────────────────────────────────────────────
out_cols <- c("gene_name", "uniprot_id", "uniprot_len",
              "mane_uniprot_id", "mane_uniprot_len", "accession_agrees",
              "mane_refseq_nuc", "mane_refseq_prot", "mane_ensembl",
              "mane_prot_len", "numbering_ok")
write.table(gene_data[, out_cols], file.path(DER, "gene_mane.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
saveRDS(gene_data, file.path(DER, "gene_data_with_mane.rds"))

mismatch <- gene_data[which(gene_data$numbering_ok == FALSE), out_cols]
mismatch <- mismatch[order(abs(mismatch$mane_uniprot_len - mismatch$mane_prot_len),
                           decreasing = TRUE), ]
write.table(mismatch, file.path(DER, "mane_numbering_mismatch.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# Genes whose stored accession disagrees with the MANE-anchored one — the map
# a later merge would correct.
corrected <- gene_data[which(gene_data$accession_agrees == FALSE), out_cols]
write.table(corrected, file.path(DER, "mane_accession_corrections.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# ── Report ────────────────────────────────────────────────────────────────
n <- nrow(gene_data)
has_mane   <- sum(!is.na(gene_data$mane_refseq_nuc))
has_anchor <- sum(!is.na(gene_data$mane_uniprot_id))
has_both   <- sum(!is.na(gene_data$numbering_ok))
n_ok       <- sum(gene_data$numbering_ok, na.rm = TRUE)
n_bad      <- sum(gene_data$numbering_ok == FALSE, na.rm = TRUE)
n_diff     <- sum(gene_data$accession_agrees == FALSE, na.rm = TRUE)
message("\n[mane] summary")
message(sprintf("  genes in gene_data              : %s", format(n, big.mark = ",")))
message(sprintf("  matched to a MANE Select        : %s (%.1f%%)", format(has_mane, big.mark=","), 100*has_mane/n))
message(sprintf("  MANE-anchored UniProt found     : %s (%.1f%%)", format(has_anchor, big.mark=","), 100*has_anchor/n))
message(sprintf("  stored accession disagrees      : %s  <- corrected by anchoring", format(n_diff, big.mark=",")))
message(sprintf("  both lengths known (anchored)   : %s", format(has_both, big.mark=",")))
message(sprintf("  numbering OK after anchoring    : %s (%.1f%% of comparable)", format(n_ok, big.mark=","), 100*n_ok/max(has_both,1)))
message(sprintf("  numbering MISMATCH remaining    : %s  <- true UniProt/MANE isoform divergence", format(n_bad, big.mark=",")))
message(sprintf("  written: %s/{gene_mane.tsv, gene_data_with_mane.rds,", DER))
message("           mane_numbering_mismatch.tsv, mane_accession_corrections.tsv}")
if (n_diff > 0) {
  message("\n  sample accession corrections (gene_data -> MANE-anchored):")
  cs <- head(corrected[order(abs(corrected$uniprot_len - corrected$mane_uniprot_len),
                             decreasing = TRUE), ], 12)
  for (i in seq_len(nrow(cs)))
    message(sprintf("    %-12s %-10s(%s) -> %-10s(%s)  MANE prot=%s",
                    cs$gene_name[i], cs$uniprot_id[i], cs$uniprot_len[i],
                    cs$mane_uniprot_id[i], cs$mane_uniprot_len[i], cs$mane_prot_len[i]))
}
if (n_bad > 0) {
  message("\n  residual mismatches (MANE-anchored UniProt canonical vs MANE protein):")
  top <- head(mismatch, 12)
  for (i in seq_len(nrow(top)))
    message(sprintf("    %-12s uniprot=%-5s mane=%-5s  d=%+d",
                    top$gene_name[i], top$mane_uniprot_len[i], top$mane_prot_len[i],
                    top$mane_uniprot_len[i] - top$mane_prot_len[i]))
}
