# Panel A (VariBench clinical) recomputed for an arbitrary classification dir.
# Regenerates Table 2 Panel A AUROC/MCC (Full + Blind) after the split-roles
# PM1 engine change. Compares old-engine (analyses/classifications) vs
# new-engine (analyses/ps_reval/classifications) on the identical VariBench set.
# ISOLATION: read-only analysis; sources only the metric library, no app edits.
suppressMessages({ library(dplyr); library(readr) })
source("analyses/lib/metric_suite.R")

UNIVERSE <- "analyses/derived/variant_universe.tsv"
univ <- read_tsv(UNIVERSE, show_col_types = FALSE) |> filter(source == "VariBench")

panel_a <- function(dir, tag) {
  fs <- list.files(dir, pattern = "__dual\\.tsv$", full.names = TRUE)
  cls <- bind_rows(lapply(fs, read_tsv, show_col_types = FALSE))
  joined <- cls |>
    inner_join(univ, by = c("gene", "p_notation")) |>
    mutate(truth = case_when(
      grepl("[Pp]athogenic", label) ~ "Pathogenic",
      grepl("[Bb]enign",     label) ~ "Benign",
      TRUE                          ~ NA_character_))
  mf <- compute_metric_suite(joined$truth, joined$varviz_classification_full,  joined$varviz_pts_full)
  mb <- compute_metric_suite(joined$truth, joined$varviz_classification_blind, joined$varviz_pts_blind)
  get <- function(m, k) { v <- m$estimate[m$metric == k]; if (length(v)) v[1] else NA_real_ }
  cat(sprintf("\n=== Panel A [%s]  (%s)\n", tag, dir))
  cat(sprintf("  VariBench variants matched: %d  (P=%d  B=%d)\n",
              nrow(joined), sum(joined$truth=="Pathogenic",na.rm=TRUE), sum(joined$truth=="Benign",na.rm=TRUE)))
  for (k in c("AUROC","MCC","Sensitivity","Specificity")) {
    f <- get(mf,k); b <- get(mb,k)
    cat(sprintf("  %-12s Full=%.3f  Blind=%.3f  Δ(Blind-Full)=%+.3f\n", k, f, b, b-f))
  }
  invisible(list(full=mf, blind=mb, n=nrow(joined)))
}

panel_a("analyses/classifications",          "OLD engine (pre split-roles)")
panel_a("analyses/ps_reval/classifications", "NEW engine (split-roles)")
