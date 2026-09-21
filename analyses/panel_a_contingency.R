#!/usr/bin/env Rscript
# Panel A as a full 3x2 contingency table, with both the conditional metrics
# and the unconditional recall.
#
# compute_metric_suite() drops VUS calls before computing sensitivity,
# specificity and MCC (metric_suite.R: `collapsed != "VUS"`). That is a defensible
# convention -- a VUS is not a wrong answer, it is a refusal to answer -- but it
# means "sensitivity 1.000" is conditional on the engine having returned a
# non-VUS call. Under blinding the refusals are numerous, so the conditional
# and unconditional figures diverge sharply.
#
# This prints both, plus the 3x2 table they are derived from, so the
# denominator is explicit rather than implied.
#
# Usage:
#   Rscript analyses/panel_a_contingency.R [classification-dir]
suppressMessages({ library(dplyr); library(readr) })

DIR  <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(DIR)) DIR <- "analyses/ps_final_v231/classifications"
UNIV <- "analyses/derived/variant_universe.tsv"

univ <- read_tsv(UNIV, show_col_types = FALSE) |> filter(source == "VariBench")
fs   <- list.files(DIR, pattern = "__dual\\.tsv$", full.names = TRUE)
cls  <- bind_rows(lapply(fs, read_tsv, show_col_types = FALSE)) |>
  distinct(gene, p_notation, .keep_all = TRUE)

joined <- cls |>
  inner_join(univ |> distinct(gene, p_notation, .keep_all = TRUE),
             by = c("gene", "p_notation")) |>
  mutate(truth = case_when(grepl("[Pp]athogenic", label) ~ "Pathogenic",
                           grepl("[Bb]enign",     label) ~ "Benign",
                           TRUE                          ~ NA_character_)) |>
  filter(!is.na(truth))

collapse3 <- function(x) ifelse(x %in% c("Pathogenic", "Likely Pathogenic"), "P/LP",
                         ifelse(x %in% c("Benign", "Likely Benign"),        "B/LB", "VUS"))

cat(sprintf("Panel A contingency — %s\n", DIR))
cat(sprintf("VariBench variants matched: %d  (Pathogenic %d, Benign %d)\n\n",
            nrow(joined), sum(joined$truth == "Pathogenic"), sum(joined$truth == "Benign")))

for (pas in c("full", "blind")) {
  col  <- paste0("varviz_classification_", pas)
  pred <- collapse3(joined[[col]])
  tab  <- table(factor(joined$truth, levels = c("Pathogenic", "Benign")),
                factor(pred, levels = c("P/LP", "VUS", "B/LB")))
  cat(sprintf("--- Pass-%s ---\n", ifelse(pas == "full", "Full", "Blind")))
  print(tab)
  P <- tab["Pathogenic", ]; B <- tab["Benign", ]
  # conditional: VUS dropped, as compute_metric_suite() does
  cond_sens <- P["P/LP"] / (P["P/LP"] + P["B/LB"])
  cond_spec <- B["B/LB"] / (B["B/LB"] + B["P/LP"])
  # unconditional: VUS counted as a miss against the full reference class
  unc_sens  <- P["P/LP"] / sum(P)
  unc_spec  <- B["B/LB"] / sum(B)
  cat(sprintf("  conditional  (VUS excluded): sensitivity %.3f  specificity %.3f\n",
              cond_sens, cond_spec))
  cat(sprintf("  unconditional (VUS counted): sensitivity %.3f  specificity %.3f\n",
              unc_sens, unc_spec))
  cat(sprintf("  reference-pathogenic called VUS: %d of %d (%.1f%%)\n\n",
              P["VUS"], sum(P), 100 * P["VUS"] / sum(P)))
}
