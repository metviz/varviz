#!/usr/bin/env Rscript
# How many variants move on ANY transcript-aligned predictor under MANE-aware
# selection, not just on AlphaMissense?
#
# analyses/dbnsfp_mane_am_audit.R measured AlphaMissense alone. The fix in
# analyses/lib/local_predictors.R also moves SIFT, PolyPhen2 HDIV and HVAR,
# FATHMM, PROVEAN, MetaRNN and REVEL, every one verified aligned with
# Ensembl_transcriptid over 5,256 multi-transcript rows. REVEL and MetaRNN carry
# the most weight because both feed PP3 directly: REVEL through its own
# calibrated intervals and MetaRNN through the meta-predictor consensus rule.
#
# The comparison drives the real lookup path twice per variant, with
# varviz.dbnsfp_mane TRUE and FALSE, so it measures shipped behaviour rather than
# a re-expression of the rule.
#
# Run from project root:
#   Rscript analyses/dbnsfp_mane_predictor_delta.R
suppressMessages({ library(readr); library(dplyr) })
source("analyses/lib/local_predictors.R")

univ <- read_tsv(Sys.getenv("VARVIZ_UNIVERSE",
                 "analyses/derived/variant_universe_gnomad.tsv"), show_col_types = FALSE)

PRED <- list(AM = c("alphamissense","am_pathogenicity"), REVEL = c("revel","score"),
             MetaRNN = c("metarnn","score"), SIFT = c("sift","score"),
             PP2HDIV = c("polyphen2","hdiv","score"), PP2HVAR = c("polyphen2","hvar","score"),
             FATHMM = c("fathmm","score"), PROVEAN = c("provean","score"))
dig <- function(h, path) {
  if (is.null(h)) return(NA_real_)
  r <- h$dbnsfp
  for (k in path) { if (is.null(r)) return(NA_real_); r <- r[[k]] }
  if (is.null(r) || length(r) == 0L) NA_real_ else as.numeric(r[1])
}
moved <- function(x, y) {
  if (is.na(x) && is.na(y)) FALSE
  else if (is.na(x) || is.na(y)) TRUE
  else !isTRUE(all.equal(x, y))
}

out <- list()
for (g in sort(unique(univ$gene))) {
  ex <- tryCatch(readRDS(file.path("analyses/cache/ensembl", paste0(g, "_exons.rds"))),
                 error = function(e) NULL)
  if (is.null(ex)) { message("skip ", g); next }
  ch  <- as.character(ex$chr[1])
  env <- tryCatch(load_dbnsfp_for_region(ch, min(ex$genomic_start), max(ex$genomic_end)),
                  error = function(e) NULL)
  if (is.null(env)) { message("skip ", g, ": no dbNSFP"); next }
  v   <- univ$p_notation[univ$gene == g]
  pos <- suppressWarnings(as.integer(sub("^p\\.[A-Z*]([0-9]+)[A-Z*]$", "\\1", v)))
  ref <- sub("^p\\.([A-Z*])[0-9]+[A-Z*]$", "\\1", v)
  alt <- sub("^p\\.[A-Z*][0-9]+([A-Z*])$", "\\1", v)
  acc <- vector("list", length(v))
  for (k in seq_along(v)) {
    if (is.na(pos[k])) next
    options(varviz.dbnsfp_mane = TRUE)
    a <- tryCatch(lookup_dbnsfp_by_aa(env, ch, pos[k], alt[k], aa_ref = ref[k]),
                  error = function(e) NULL)
    if (is.null(a)) next
    options(varviz.dbnsfp_mane = FALSE)
    b <- tryCatch(lookup_dbnsfp_by_aa(env, ch, pos[k], alt[k], aa_ref = ref[k]),
                  error = function(e) NULL)
    acc[[k]] <- c(list(gene = g, Variant = v[k]),
                  setNames(lapply(PRED, function(pa) moved(dig(a, pa), dig(b, pa))),
                           names(PRED)))
  }
  options(varviz.dbnsfp_mane = TRUE)
  acc <- acc[!vapply(acc, is.null, logical(1))]
  if (length(acc)) out[[g]] <- bind_rows(acc)
  message(sprintf("  %-9s resolved=%d", g, length(acc)))
}

df <- bind_rows(out)
cat("\n====== predictors moved by MANE-aware selection ======\n")
cat(sprintf("variants resolved in dbNSFP : %d\n\n", nrow(df)))
cat(sprintf("%-9s %9s %8s\n", "predictor", "moved", "pct"))
for (nm in names(PRED)) {
  n <- sum(df[[nm]], na.rm = TRUE)
  cat(sprintf("%-9s %9d %7.2f%%\n", nm, n, 100*n/nrow(df)))
}
any_ch <- Reduce(`|`, df[names(PRED)])
cat(sprintf("\nANY predictor moved         : %d (%.2f%%)\n", sum(any_ch), 100*mean(any_ch)))
cat("\nby gene:\n")
df |> mutate(any_ch = any_ch) |> group_by(gene) |>
  summarise(n = n(), moved = sum(any_ch), pct = round(100*moved/n, 2), .groups = "drop") |>
  arrange(desc(moved)) |> as.data.frame() |> print(row.names = FALSE)

dir.create("analyses/reports", showWarnings = FALSE, recursive = TRUE)
write_tsv(df, "analyses/reports/dbnsfp_mane_predictor_delta.tsv")
cat("\nwrote analyses/reports/dbnsfp_mane_predictor_delta.tsv\n")
