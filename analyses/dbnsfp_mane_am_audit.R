#!/usr/bin/env Rscript
# How much of the published AlphaMissense input rests on the wrong transcript?
#
# dbNSFP 4.9a stores AlphaMissense, SIFT, PolyPhen2 and PROVEAN as ';'-separated
# per-transcript lists, and .dbnsfp_row_to_hit() reads them with num_first(), i.e.
# the first non-'.' element. That is only the MANE transcript's value when MANE
# happens to be listed first, which analyses/reports/dbnsfp_mane_index_survey.tsv
# shows is the minority for most genes (GCK 0 of 88, TP53 0 of 91, KCNH2 3 of 55,
# BRCA1 16 of 110).
#
# This audit compares, for every variant in the benchmark universe, the value the
# first element gives against the value the MANE transcript's own element gives,
# and buckets both into the PP3 tier their cut-point implies. It then intersects
# that with what the published run actually consumed: AM_Source in the master table
# records whether the score came from the dbNSFP/MyVariant path (exposed to the
# defect) or from the per-UniProt AlphaFold CSV (already MANE-correct, because the
# CSV is keyed by UniProt accession and has no transcript dimension).
#
# Indexing by transcript is safe here: the transcript list and each predictor list
# were verified to have equal length in every row examined, and the MANE
# transcript's own aapos element is used for the position match, because aapos is
# per-transcript too and other transcripts carry a constant offset (CASR: +10 and
# -77 relative to MANE).
#
# Run from project root:
#   Rscript analyses/dbnsfp_mane_am_audit.R
suppressMessages({ library(readr); library(dplyr) })
source("analyses/lib/local_predictors.R")

COL_AAREF <- 5L; COL_AAALT <- 6L; COL_AAPOS <- 12L; COL_GENE <- 13L
COL_TX <- 15L; COL_AM <- 138L

MASTER <- Sys.getenv("VARVIZ_MASTER",
                     "analyses/ps_am792_v123_final/master_all_variants.tsv")
UNIV   <- Sys.getenv("VARVIZ_UNIVERSE",
                     "analyses/derived/variant_universe_gnomad.tsv")

# The shipped ladder. Bucketing the AlphaMissense value in isolation is not the
# final PP3 tag, which also depends on the other predictors; it measures the
# exposure attributable to transcript selection, which is the question here.
THR  <- c(strong = 0.990, three_pt = 0.972, moderate = 0.906, supporting = 0.792)
RANK <- c(absent = -1L, none = 0L, supporting = 1L, moderate = 2L,
          three_pt = 3L, strong = 4L)
bucket <- function(x) {
  ifelse(is.na(x), "absent",
  ifelse(x >= THR[["strong"]],     "strong",
  ifelse(x >= THR[["three_pt"]],   "three_pt",
  ifelse(x >= THR[["moderate"]],   "moderate",
  ifelse(x >= THR[["supporting"]], "supporting", "none")))))
}

mane <- read_tsv("analyses/derived/gene_mane.tsv", show_col_types = FALSE)
univ <- read_tsv(UNIV, show_col_types = FALSE)

out <- list()
for (g in sort(unique(univ$gene))) {
  mrow <- mane[match(g, mane$gene_name), ]
  if (!nrow(mrow) || is.na(mrow$mane_ensembl)) { message("skip ", g, ": no MANE"); next }
  mt <- sub("\\..*$", "", mrow$mane_ensembl)
  ex <- tryCatch(readRDS(file.path("analyses/cache/ensembl", paste0(g, "_exons.rds"))),
                 error = function(e) NULL)
  if (is.null(ex)) { message("skip ", g, ": no exon cache"); next }
  ch <- as.character(ex$chr[1])
  rows <- tryCatch(system(sprintf("tabix %s %s:%d-%d", shQuote(DBNSFP_DEFAULT_PATH),
                                  ch, min(ex$genomic_start), max(ex$genomic_end)),
                          intern = TRUE, ignore.stderr = TRUE),
                   error = function(e) character(0))
  if (!length(rows)) { message("skip ", g, ": no dbNSFP rows"); next }

  # key -> c(first_element, mane_element). Keyed on the protein substitution, so
  # the several nucleotide rows encoding one substitution collapse to one entry.
  env <- new.env(hash = TRUE, parent = emptyenv())
  for (l in rows) {
    f <- strsplit(l, "\t", fixed = TRUE)[[1]]
    if (length(f) < COL_AM) next
    if (!grepl(g, f[COL_GENE], fixed = TRUE)) next
    am <- strsplit(f[COL_AM], ";", fixed = TRUE)[[1]]
    if (!length(am)) next
    tx <- strsplit(f[COL_TX], ";", fixed = TRUE)[[1]]
    ap <- strsplit(f[COL_AAPOS], ";", fixed = TRUE)[[1]]
    i  <- which(tx == mt)
    if (!length(i) || i[1] > length(ap)) next
    pos <- ap[i[1]]
    if (is.na(pos) || pos == ".") next
    k <- paste0(f[COL_AAREF], pos, f[COL_AAALT])
    if (!is.null(env[[k]])) next
    nz <- suppressWarnings(as.numeric(am[am != "." & am != ""]))
    env[[k]] <- c(if (length(nz)) nz[1] else NA_real_,
                  suppressWarnings(as.numeric(am[min(i[1], length(am))])))
  }

  v   <- univ$p_notation[univ$gene == g]
  ref <- sub("^p\\.([A-Z*])[0-9]+[A-Z*]$", "\\1", v)
  pos <- sub("^p\\.[A-Z*]([0-9]+)[A-Z*]$", "\\1", v)
  alt <- sub("^p\\.[A-Z*][0-9]+([A-Z*])$", "\\1", v)
  keys <- paste0(ref, pos, alt)
  got <- lapply(keys, function(k) { x <- env[[k]]; if (is.null(x)) c(NA_real_, NA_real_) else x })

  out[[g]] <- tibble(gene = g, Variant = v,
                     am_first = vapply(got, `[`, numeric(1), 1L),
                     am_mane  = vapply(got, `[`, numeric(1), 2L))
  message(sprintf("  %-9s n=%-6d dbNSFP_hit=%-6d", g, length(v), sum(!is.na(out[[g]]$am_mane))))
}

df <- bind_rows(out) |>
  mutate(b_first = bucket(am_first), b_mane = bucket(am_mane),
         r_first = RANK[b_first], r_mane = RANK[b_mane])

cat("\n============ first-element vs MANE-element AlphaMissense ============\n")
both <- df |> filter(!is.na(am_first), !is.na(am_mane))
cat(sprintf("universe variants resolved in dbNSFP : %d\n", nrow(both)))
cat(sprintf("value differs                        : %d (%.1f%%)\n",
            sum(both$am_first != both$am_mane), 100*mean(both$am_first != both$am_mane)))
cat(sprintf("TIER differs                         : %d (%.2f%%)\n",
            sum(both$b_first != both$b_mane), 100*mean(both$b_first != both$b_mane)))
d <- both |> filter(b_first != b_mane)
cat(sprintf("  first-element gives HIGHER tier    : %d\n", sum(d$r_first > d$r_mane)))
cat(sprintf("  MANE gives HIGHER tier             : %d\n", sum(d$r_mane > d$r_first)))
cat(sprintf("max |difference|                     : %.4f\n",
            if (nrow(both)) max(abs(both$am_first - both$am_mane)) else NA))

cat("\nby gene:\n")
both |> group_by(gene) |>
  summarise(n = n(), val_diff = sum(am_first != am_mane),
            tier_diff = sum(b_first != b_mane),
            pct = round(100*tier_diff/n, 2), .groups = "drop") |>
  arrange(desc(tier_diff)) |> as.data.frame() |> print(row.names = FALSE)

# The part that decides whether published numbers move: only rows whose score the
# run actually took from the dbNSFP/MyVariant path were exposed.
if (file.exists(MASTER)) {
  m <- read_tsv(MASTER, show_col_types = FALSE,
                col_select = c("gene", "Variant", "AM_Score", "AM_Source")) |>
       distinct(gene, Variant, .keep_all = TRUE)
  j <- df |> inner_join(m, by = c("gene", "Variant")) |>
       mutate(b_used = bucket(AM_Score), r_used = RANK[b_used])
  cat(sprintf("\n=========== exposure in the published run (%s) ===========\n", basename(MASTER)))
  cat(sprintf("joined rows                          : %d\n", nrow(j)))
  print(as.data.frame(j |> count(AM_Source)), row.names = FALSE)
  ex <- j |> filter(AM_Source != "AlphaFold", !is.na(am_mane))
  cat(sprintf("\nrows sourced from dbNSFP/MyVariant    : %d\n", nrow(ex)))
  cat(sprintf("  used value != MANE value           : %d\n",
              sum(abs(ex$AM_Score - ex$am_mane) > 5e-5, na.rm = TRUE)))
  cat(sprintf("  TIER would change under MANE       : %d\n", sum(ex$b_used != ex$b_mane)))
  chg <- ex |> filter(b_used != b_mane)
  if (nrow(chg)) {
    cat("\n  transitions (used tier -> MANE tier):\n")
    chg |> count(b_used, b_mane, sort = TRUE) |> as.data.frame() |> print(row.names = FALSE)
    cat("\n  by gene:\n")
    chg |> count(gene, sort = TRUE) |> as.data.frame() |> print(row.names = FALSE)
  }
  dir.create("analyses/reports", showWarnings = FALSE, recursive = TRUE)
  write_tsv(j, "analyses/reports/dbnsfp_mane_am_audit.tsv")
  cat("\nwrote analyses/reports/dbnsfp_mane_am_audit.tsv\n")
} else {
  message("master table not found at ", MASTER, "; skipped exposure section")
}
