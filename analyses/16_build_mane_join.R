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

# ── UniProt canonical length: one stream of the reviewed human proteome ────
# Cached so re-runs are offline. ~20k Swiss-Prot accessions, 500 per page.
UP_CACHE <- file.path(TMP, "uniprot_len.tsv")
if (file.exists(UP_CACHE)) {
  message("[mane] uniprot lengths from cache")
  upl <- read.delim(UP_CACHE, stringsAsFactors = FALSE)
} else {
  message("[mane] streaming UniProt reviewed human lengths (one-time)")
  rows <- list(); cursor <- NULL; page <- 0L
  repeat {
    req <- request("https://rest.uniprot.org/uniprotkb/search") |>
      req_url_query(query = "reviewed:true AND organism_id:9606",
                    fields = "accession,length", format = "tsv", size = 500)
    if (!is.null(cursor)) req <- req_url_query(req, cursor = cursor)
    resp <- req_perform(req)
    txt <- resp_body_string(resp)
    df  <- read.delim(text = txt, stringsAsFactors = FALSE)
    if (nrow(df) == 0) break
    rows[[length(rows) + 1L]] <- df
    page <- page + 1L
    # UniProt paginates with a Link: rel="next" header carrying the cursor.
    link <- resp_headers(resp)[["link"]]
    if (is.null(link) || !grepl("rel=\"next\"", link)) break
    cursor <- sub(".*[?&]cursor=([^&>]+).*", "\\1", link)
    if (page %% 5L == 0L) message("  ", page * 500L, " streamed")
  }
  upl <- do.call(rbind, rows)
  colnames(upl) <- c("uniprot_id", "uniprot_len")
  write.table(upl, UP_CACHE, sep = "\t", row.names = FALSE, quote = FALSE)
}
up_len <- setNames(as.integer(upl$uniprot_len), upl$uniprot_id)

# ── Join everything onto gene_data ────────────────────────────────────────
m <- match(gene_data$gene_name, sel$symbol)
gene_data$mane_refseq_nuc  <- sel$RefSeq_nuc[m]
gene_data$mane_refseq_prot <- sel$RefSeq_prot[m]
gene_data$mane_ensembl     <- sel$Ensembl_nuc[m]
gene_data$mane_prot_len    <- sel$mane_prot_len[m]
gene_data$uniprot_len      <- unname(up_len[gene_data$uniprot_id])

# numbering_ok: TRUE  = both lengths known and equal (safe to match by residue)
#               FALSE = both known and DIFFERENT (at-risk, needs alignment)
#               NA    = a length is missing (no MANE gene, or TrEMBL accession)
gene_data$numbering_ok <- with(gene_data,
  ifelse(is.na(uniprot_len) | is.na(mane_prot_len), NA,
         uniprot_len == mane_prot_len))

# ── Write outputs ─────────────────────────────────────────────────────────
out_cols <- c("gene_name", "uniprot_id", "uniprot_len",
              "mane_refseq_nuc", "mane_refseq_prot", "mane_ensembl",
              "mane_prot_len", "numbering_ok")
write.table(gene_data[, out_cols], file.path(DER, "gene_mane.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")
saveRDS(gene_data, file.path(DER, "gene_data_with_mane.rds"))

mismatch <- gene_data[which(gene_data$numbering_ok == FALSE), out_cols]
mismatch <- mismatch[order(abs(mismatch$uniprot_len - mismatch$mane_prot_len),
                           decreasing = TRUE), ]
write.table(mismatch, file.path(DER, "mane_numbering_mismatch.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# ── Report ────────────────────────────────────────────────────────────────
n <- nrow(gene_data)
has_mane <- sum(!is.na(gene_data$mane_refseq_nuc))
has_both <- sum(!is.na(gene_data$numbering_ok))
n_ok     <- sum(gene_data$numbering_ok, na.rm = TRUE)
n_bad    <- sum(gene_data$numbering_ok == FALSE, na.rm = TRUE)
message("\n[mane] summary")
message(sprintf("  genes in gene_data          : %s", format(n, big.mark = ",")))
message(sprintf("  matched to a MANE Select    : %s (%.1f%%)", format(has_mane, big.mark=","), 100*has_mane/n))
message(sprintf("  both lengths known          : %s", format(has_both, big.mark=",")))
message(sprintf("  numbering OK (len match)    : %s (%.1f%% of comparable)", format(n_ok, big.mark=","), 100*n_ok/max(has_both,1)))
message(sprintf("  numbering MISMATCH          : %s", format(n_bad, big.mark=",")))
message(sprintf("  written: %s/gene_mane.tsv, gene_data_with_mane.rds, mane_numbering_mismatch.tsv", DER))
if (n_bad > 0) {
  message("\n  top mismatches (UniProt canonical vs MANE protein length):")
  top <- head(mismatch, 15)
  for (i in seq_len(nrow(top)))
    message(sprintf("    %-12s uniprot=%-5s mane=%-5s  d=%+d",
                    top$gene_name[i], top$uniprot_len[i], top$mane_prot_len[i],
                    top$uniprot_len[i] - top$mane_prot_len[i]))
}
