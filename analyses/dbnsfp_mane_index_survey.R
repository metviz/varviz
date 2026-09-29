#!/usr/bin/env Rscript
# Where does the MANE transcript sit in dbNSFP's per-transcript lists, and does
# the per-UniProt AlphaMissense CSV agree with the MANE element?
#
# dbNSFP stores several predictors as ';'-separated per-transcript lists, and
# .dbnsfp_row_to_hit() reads them with num_first()/chr_first(), i.e. the first
# non-'.' element. That is only correct when MANE happens to be listed first.
# Observed on KCNH2 p.A1017D, where MANE is the second of two transcripts and the
# first element is a different isoform's score.
#
# This surveys all 14 benchmark genes and reports, per gene:
#   - whether the MANE transcript appears in dbNSFP's list at all
#   - at which index, and how often that index is 1 (where the current read is
#     accidentally correct)
#   - whether the AlphaFold/EBI CSV value equals the MANE element's value
#
# The last column is the one that decides the fix. Where the CSV equals the MANE
# element, the CSV is already a MANE-correct source for AlphaMissense and no
# transcript indexing is needed for it. Where it does not, the CSV is on another
# sequence, which analyses/derived/gene_mane.tsv records as numbering_ok = FALSE.
#
# Run from project root:
#   Rscript analyses/dbnsfp_mane_index_survey.R
suppressMessages({ library(readr); library(dplyr) })
source("analyses/lib/local_predictors.R")

COL_AAPOS <- 12L; COL_AAALT <- 6L; COL_TX <- 15L; COL_AM <- 138L
N_SAMPLE  <- as.integer(Sys.getenv("N_SAMPLE", "400"))

UID <- c(BAP1="Q92560", BRCA1="P38398", CASR="P41180", DDR2="Q16832",
         GCK="P35557", KCNH2="Q12809", KCNQ1="P51787", KRAS="P01116",
         LDLR="P01130", NUDT15="Q9NV35", PTEN="P60484", SLC13A5="Q86YT5",
         TP53="P04637", TSHR="P16473")

mane <- read_tsv("analyses/derived/gene_mane.tsv", show_col_types = FALSE)
univ <- read_tsv(Sys.getenv("VARVIZ_UNIVERSE", "analyses/derived/variant_universe_gnomad.tsv"),
                 show_col_types = FALSE)

out <- list()
for (g in names(UID)) {
  mrow <- mane[match(g, mane$gene_name), ]
  mt   <- sub("\\..*$", "", mrow$mane_ensembl)
  nok  <- identical(toupper(trimws(as.character(mrow$numbering_ok))), "TRUE")
  ex <- tryCatch(readRDS(file.path("analyses/cache/ensembl", paste0(g, "_exons.rds"))),
                 error = function(e) NULL)
  if (is.null(ex)) { message("skip ", g); next }
  chrom <- as.character(ex$chr[1])
  rows <- tryCatch(system(sprintf("tabix %s %s:%d-%d", shQuote(DBNSFP_DEFAULT_PATH),
                                  chrom, min(ex$genomic_start), max(ex$genomic_end)),
                          intern = TRUE, ignore.stderr = TRUE),
                   error = function(e) character(0))
  am_csv <- tryCatch(load_alphamissense_csv(UID[[g]]), error = function(e) NULL)
  if (!length(rows)) { message("skip ", g, ": no dbNSFP rows"); next }

  v <- univ$p_notation[univ$gene == g]
  if (length(v) > N_SAMPLE) { set.seed(1); v <- sample(v, N_SAMPLE) }
  key <- function(p) list(pos = sub("^p\\.[A-Z*]([0-9]+).*$", "\\1", p),
                          alt = sub("^p\\.[A-Z*][0-9]+([A-Z*])$", "\\1", p))

  n_seen <- 0L; n_listed <- 0L; n_first <- 0L; n_csv_eq_mane <- 0L; n_csv_cmp <- 0L
  idxs <- integer(0); n_multi <- 0L
  for (p_not in v) {
    k <- key(p_not)
    for (l in rows) {
      p <- strsplit(l, "\t", fixed = TRUE)[[1]]
      if (length(p) < COL_AM) next
      if (!grepl(paste0("(^|;)", k$pos, "($|;)"), p[COL_AAPOS])) next
      if (p[COL_AAALT] != k$alt) next
      tx <- strsplit(p[COL_TX], ";", fixed = TRUE)[[1]]
      am <- strsplit(p[COL_AM], ";", fixed = TRUE)[[1]]
      n_seen <- n_seen + 1L
      if (length(tx) > 1L) n_multi <- n_multi + 1L
      i <- which(tx == mt)
      if (length(i)) {
        n_listed <- n_listed + 1L; idxs <- c(idxs, i[1])
        if (i[1] == 1L) n_first <- n_first + 1L
        mv <- suppressWarnings(as.numeric(am[min(i[1], length(am))]))
        cv <- tryCatch(load_alphamissense_csv(UID[[g]]), error = function(e) NULL)
        cs <- if (is.null(am_csv)) NULL else lookup_alphamissense(am_csv, p_not)
        if (!is.null(cs) && !is.null(cs$score) && !is.na(mv)) {
          n_csv_cmp <- n_csv_cmp + 1L
          if (isTRUE(all.equal(as.numeric(cs$score), mv, tolerance = 1e-6))) n_csv_eq_mane <- n_csv_eq_mane + 1L
        }
      }
      break
    }
  }
  out[[g]] <- tibble(gene = g, mane_tx = mt, numbering_ok = nok,
                     sampled = length(v), matched = n_seen,
                     multi_tx = n_multi, mane_listed = n_listed,
                     mane_is_first = n_first,
                     median_mane_idx = if (length(idxs)) median(idxs) else NA_real_,
                     max_mane_idx = if (length(idxs)) max(idxs) else NA_integer_,
                     csv_compared = n_csv_cmp, csv_eq_mane = n_csv_eq_mane,
                     pct_csv_eq_mane = if (n_csv_cmp) round(100*n_csv_eq_mane/n_csv_cmp, 1) else NA_real_)
  message(sprintf("  %-9s matched=%-5d mane_listed=%-5d first=%-5d csv=mane %s%%",
                  g, n_seen, n_listed, n_first,
                  ifelse(n_csv_cmp, round(100*n_csv_eq_mane/n_csv_cmp,1), "n/a")))
}

res <- bind_rows(out)
cat("\n=============== MANE index survey (dbNSFP per-transcript lists) ===============\n")
as.data.frame(res) |> print(row.names = FALSE)
dir.create("analyses/reports", showWarnings = FALSE, recursive = TRUE)
write_tsv(res, "analyses/reports/dbnsfp_mane_index_survey.tsv")
cat("\nwrote analyses/reports/dbnsfp_mane_index_survey.tsv\n")
