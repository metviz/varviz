#!/usr/bin/env Rscript
# Per-gene MANE Select residue maps: position -> reference amino acid.
#
# Why this exists. analyses/check_universe_numbering.py first used the
# per-accession AlphaMissense CSV as the reference sequence, because it enumerates
# every substitution of every residue and so states the reference residue at each
# position. That is the UniProt canonical sequence, and for KRAS it is the wrong
# one: the universe's 18,979 KRAS variants match MANE ENST00000311936 at 100.0%
# and P01116 at only 90.2%, differing at 22 positions between 151 and 189 where
# the two KRAS splice isoforms diverge. Checked against the CSV, KRAS looked a
# tenth off-sequence; checked against MANE it is exact.
#
# The map is built from dbNSFP without a network call. Column 5 (aaref) is the
# row's reference residue and column 12 (aapos) gives that residue's position in
# each transcript listed in column 15, so selecting the MANE transcript's element
# yields position -> residue on the MANE protein.
#
# Output is regenerable and lives under analyses/derived/, which is gitignored:
# it is derived data, not code.
#
# Run from project root:
#   Rscript analyses/build_mane_residue_maps.R
suppressMessages({ library(readr) })
source("analyses/lib/local_predictors.R")

COL_AAREF <- 5L; COL_AAPOS <- 12L; COL_GENE <- 13L; COL_TX <- 15L
OUT  <- "analyses/derived/mane_residue_maps.tsv"
UNIV <- Sys.getenv("VARVIZ_UNIVERSE", "analyses/derived/variant_universe_gnomad.tsv")

genes <- sort(unique(read_tsv(UNIV, show_col_types = FALSE)$gene))
out <- list()
for (g in genes) {
  mt <- .dbnsfp_mane_tx(g)
  if (is.na(mt)) { message("skip ", g, ": no MANE in gene_mane.tsv"); next }
  ex <- tryCatch(readRDS(file.path("analyses/cache/ensembl", paste0(g, "_exons.rds"))),
                 error = function(e) NULL)
  if (is.null(ex)) { message("skip ", g, ": no exon cache"); next }
  ch <- as.character(ex$chr[1])
  rows <- tryCatch(system(sprintf("tabix %s %s:%d-%d", shQuote(DBNSFP_DEFAULT_PATH),
                                  ch, min(ex$genomic_start), max(ex$genomic_end)),
                          intern = TRUE, ignore.stderr = TRUE),
                   error = function(e) character(0))
  if (!length(rows)) { message("skip ", g, ": no dbNSFP rows"); next }
  res <- character(0)
  for (l in rows) {
    f <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(f) < COL_TX) next
    aaref <- f[COL_AAREF]
    if (is.na(aaref) || nchar(aaref) != 1L || aaref == ".") next
    tx <- strsplit(f[COL_TX], ";", fixed = TRUE)[[1]]
    ap <- strsplit(f[COL_AAPOS], ";", fixed = TRUE)[[1]]
    gn <- strsplit(f[COL_GENE], ";", fixed = TRUE)[[1]]
    i  <- which(tx == mt)
    if (!length(i) || i[1] > length(ap) || i[1] > length(gn)) next
    # Guard against an overlapping gene sharing the row.
    if (!identical(gn[i[1]], g)) next
    p <- suppressWarnings(as.integer(ap[i[1]]))
    if (is.na(p)) next
    k <- as.character(p)
    # First write wins: every genomic encoding of a residue reports the same
    # reference amino acid, so a later row cannot disagree.
    if (is.na(res[k])) res[k] <- aaref
  }
  if (!length(res)) { message("skip ", g, ": MANE not listed in any row"); next }
  pos <- as.integer(names(res))
  ord <- order(pos)
  out[[g]] <- data.frame(gene = g, mane_tx = mt, pos = pos[ord],
                         aa = unname(res)[ord], stringsAsFactors = FALSE)
  message(sprintf("  %-9s %-18s positions=%-5d max=%d", g, mt, length(res), max(pos)))
}

df <- do.call(rbind, out)
dir.create("analyses/derived", showWarnings = FALSE, recursive = TRUE)
write_tsv(df, OUT)
cat(sprintf("\nwrote %s (%d rows, %d genes)\n", OUT, nrow(df), length(out)))
