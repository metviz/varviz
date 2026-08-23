#!/usr/bin/env Rscript
# Panel A (VariBench clinical concordance) recomputed from the authoritative run.
#
# Uses analyses/lib/metric_suite.R unchanged, so Sens/Spec/MCC/AUROC/VUS-rate are
# computed by exactly the same code as analyses/06_panel_a_clinical.R. The only
# difference is the input: the authoritative derived_mds run, scored for three
# arms (Pass-Full, Pass-Blind, Pass-Blind + MDS) instead of two.
#
# Note on Sens/Spec: metric_suite.R excludes VUS calls before computing them, so
# both are conditional on the engine making a non-VUS call. The VUS rate is the
# other half of that picture and must be read alongside them.
#
# Run from the repository root:  Rscript analyses/repro/04_panel_a.R

suppressMessages({ library(dplyr); library(readr) })
source("analyses/lib/metric_suite.R")

CLS  <- "analyses/derived/varviz_classifications_mds.tsv"
UNIV <- "analyses/derived/variant_universe.tsv"
OUT  <- "analyses/repro/out/panel_a.tsv"

cls  <- read_tsv(CLS,  show_col_types = FALSE)
univ <- read_tsv(UNIV, show_col_types = FALSE)

joined <- cls |>
  inner_join(univ |> filter(source == "VariBench"), by = c("gene", "p_notation")) |>
  mutate(truth = case_when(
    grepl("[Pp]athogenic", label) ~ "Pathogenic",
    grepl("[Bb]enign",     label) ~ "Benign",
    TRUE                          ~ NA_character_))

cat(sprintf("[panel_a] VariBench-labelled variants: %d  (P %d / B %d)\n",
            nrow(joined), sum(joined$truth == "Pathogenic", na.rm = TRUE),
            sum(joined$truth == "Benign", na.rm = TRUE)))

arms <- list(
  Full        = list(cls = joined$varviz_classification_full,      pts = joined$varviz_pts_full),
  Blind       = list(cls = joined$varviz_classification_blind,     pts = joined$varviz_pts_blind),
  `Blind+MDS` = list(cls = joined$varviz_classification_blind_mds, pts = joined$varviz_pts_blind_mds)
)

res <- bind_rows(lapply(names(arms), function(a) {
  compute_metric_suite(truth = joined$truth,
                       pred_class = arms[[a]]$cls,
                       pred_score = arms[[a]]$pts) |> mutate(arm = a)
}))

wide <- res |>
  select(metric, arm, estimate) |>
  tidyr::pivot_wider(names_from = arm, values_from = estimate)

cat("\n")
print(as.data.frame(wide), digits = 4, row.names = FALSE)

au <- function(a) res$estimate[res$metric == "AUROC" & res$arm == a]
cat(sprintf("\nDelta AUROC (arm - Pass-Full):  Blind %+.3f   Blind+MDS %+.3f\n",
            au("Blind") - au("Full"), au("Blind+MDS") - au("Full")))
cat(sprintf("MDS recovers %+.3f AUROC over baseline Pass-Blind\n",
            au("Blind+MDS") - au("Blind")))

write.table(res, OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("[repro] wrote %s\n", OUT))
