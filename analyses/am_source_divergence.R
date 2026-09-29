#!/usr/bin/env Rscript
# How much does the choice of AlphaMissense source change the PP3 tier?
#
# The engine resolves AlphaMissense from two places: the dbNSFP 4.9a column and
# the per-UniProt AlphaFold/EBI substitution CSV. They are not the same numbers.
# KCNQ1 p.F270L is 0.9868 in dbNSFP and 0.9949 in the CSV, which straddles the
# 0.990 Strong cut-point, so the same variant scores PP3_moderate or PP3_strong
# depending only on which source answered.
#
# This quantifies the exposure over the whole benchmark universe: for every
# variant, bucket the value each source would give into the tier its cut-point
# implies, and count the disagreements and their direction.
#
# Buckets follow the shipped ladder: >= 0.990 strong, >= 0.972 three-point (which
# is demoted to moderate unless varviz.pp3_3pt), >= 0.906 moderate, >= 0.792
# supporting, below that nothing. Bucketing the AM value in isolation is not the
# same as the final PP3 tag, which also depends on what other predictors
# contributed; it measures the exposure attributable to the source choice, which
# is the question.
#
# Run from project root:
#   Rscript analyses/am_source_divergence.R
suppressMessages({ library(readr); library(dplyr) })
source("analyses/lib/local_predictors.R")

UNIV <- Sys.getenv("VARVIZ_UNIVERSE", "analyses/derived/variant_universe_gnomad.tsv")
AMDIR <- "analyses/raw/alphamissense"

UID <- c(BAP1="Q92560", BRCA1="P38398", CASR="P41180", DDR2="Q16832",
         GCK="P35557", KCNH2="Q12809", KCNQ1="P51787", KRAS="P01116",
         LDLR="P01130", NUDT15="Q9NV35", PTEN="P60484", SLC13A5="Q86YT5",
         TP53="P04637", TSHR="P16473")

THR <- c(strong = 0.990, three_pt = 0.972, moderate = 0.906, supporting = 0.792)
bucket <- function(x) {
  ifelse(is.na(x), "absent",
  ifelse(x >= THR[["strong"]],     "strong",
  ifelse(x >= THR[["three_pt"]],   "three_pt",
  ifelse(x >= THR[["moderate"]],   "moderate",
  ifelse(x >= THR[["supporting"]], "supporting", "none")))))
}
RANK <- c(absent = -1L, none = 0L, supporting = 1L, moderate = 2L,
          three_pt = 3L, strong = 4L)

univ <- read_tsv(UNIV, show_col_types = FALSE)
out <- list()

for (g in sort(unique(univ$gene))) {
  uid <- UID[[g]]
  if (is.null(uid) || is.na(uid)) { message("skip ", g, ": no accession"); next }
  ex <- tryCatch(readRDS(file.path("analyses/cache/ensembl", paste0(g, "_exons.rds"))),
                 error = function(e) NULL)
  if (is.null(ex)) { message("skip ", g, ": no exon cache"); next }
  chrom <- as.character(ex$chr[1])
  env <- tryCatch(load_dbnsfp_for_region(chrom, min(ex$genomic_start), max(ex$genomic_end)),
                  error = function(e) NULL)
  am_csv <- tryCatch(load_alphamissense_csv(uid, am_dir = AMDIR), error = function(e) NULL)
  if (is.null(env) || is.null(am_csv)) { message("skip ", g, ": dbnsfp/csv unavailable"); next }

  v <- univ$p_notation[univ$gene == g]
  pos <- suppressWarnings(as.integer(sub("^p\\.[A-Z*]([0-9]+).*$", "\\1", v)))
  ref <- sub("^p\\.([A-Z*])[0-9]+[A-Z*]$", "\\1", v)
  alt <- sub("^p\\.[A-Z*][0-9]+([A-Z*])$", "\\1", v)

  d_am <- vapply(seq_along(v), function(i) {
    if (is.na(pos[i])) return(NA_real_)
    h <- tryCatch(lookup_dbnsfp_by_aa(env, chrom, pos[i], alt[i], aa_ref = ref[i]),
                  error = function(e) NULL)
    if (is.null(h)) NA_real_ else {
      x <- h$dbnsfp$alphamissense$am_pathogenicity
      if (is.null(x) || length(x) == 0) NA_real_ else as.numeric(x[1])
    }
  }, numeric(1))

  c_am <- vapply(v, function(p) {
    r <- tryCatch(lookup_alphamissense(am_csv, p), error = function(e) NULL)
    if (is.null(r) || is.null(r$score)) NA_real_ else as.numeric(r$score)
  }, numeric(1))

  out[[g]] <- tibble(gene = g, p_notation = v, dbnsfp = d_am, csv = c_am)
  message(sprintf("  %-9s n=%-6d dbNSFP=%-6d csv=%-6d",
                  g, length(v), sum(!is.na(d_am)), sum(!is.na(c_am))))
}

df <- bind_rows(out) |>
  mutate(b_db = bucket(dbnsfp), b_csv = bucket(csv),
         r_db = RANK[b_db], r_csv = RANK[b_csv])

both <- df |> filter(!is.na(dbnsfp), !is.na(csv))
cat("\n================ AlphaMissense source divergence ================\n")
cat(sprintf("universe variants            : %d\n", nrow(df)))
cat(sprintf("present in dbNSFP            : %d\n", sum(!is.na(df$dbnsfp))))
cat(sprintf("present in AlphaFold CSV     : %d\n", sum(!is.na(df$csv))))
cat(sprintf("present in BOTH (comparable) : %d\n", nrow(both)))

cat(sprintf("\nvalue differs at all         : %d (%.1f%% of comparable)\n",
            sum(both$dbnsfp != both$csv), 100*mean(both$dbnsfp != both$csv)))
cat(sprintf("median absolute difference   : %.4f\n", median(abs(both$dbnsfp - both$csv))))
cat(sprintf("max absolute difference      : %.4f\n", max(abs(both$dbnsfp - both$csv))))

disagree <- both |> filter(b_db != b_csv)
cat(sprintf("\nTIER DISAGREEMENT            : %d (%.2f%% of comparable)\n",
            nrow(disagree), 100*nrow(disagree)/nrow(both)))
cat(sprintf("  CSV gives the HIGHER tier  : %d\n", sum(disagree$r_csv > disagree$r_db)))
cat(sprintf("  dbNSFP gives HIGHER tier   : %d\n", sum(disagree$r_db > disagree$r_csv)))

cat("\ntransitions (dbNSFP tier -> CSV tier), most common first:\n")
disagree |> count(b_db, b_csv, sort = TRUE) |> as.data.frame() |> print(row.names = FALSE)

cat("\nby gene:\n")
both |> group_by(gene) |>
  summarise(n = n(), disagree = sum(b_db != b_csv),
            pct = round(100*disagree/n, 2),
            csv_higher = sum(RANK[b_csv] > RANK[b_db]),
            db_higher  = sum(RANK[b_db] > RANK[b_csv]), .groups = "drop") |>
  arrange(desc(disagree)) |> as.data.frame() |> print(row.names = FALSE)

cat("\nvariants present in only ONE source (source choice decides whether AM exists):\n")
cat(sprintf("  dbNSFP only : %d\n", sum(!is.na(df$dbnsfp) & is.na(df$csv))))
cat(sprintf("  CSV only    : %d\n", sum(is.na(df$dbnsfp) & !is.na(df$csv))))

dir.create("analyses/reports", showWarnings = FALSE, recursive = TRUE)
write_tsv(df, "analyses/reports/am_source_divergence.tsv")
cat("\nwrote analyses/reports/am_source_divergence.tsv\n")
