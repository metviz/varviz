#!/usr/bin/env Rscript
# Panel A across harness runs, using analyses/lib/metric_suite.R unchanged.
#
# With MDS ungated, the pathway is part of the engine rather than a bolt-on arm,
# so the ablation comparator is the MDS-off run rather than a separate column of
# the same file.
#
# Run from the repository root:  Rscript analyses/repro/06_panel_a_runs.R

suppressMessages({ library(dplyr); library(readr) })
source("analyses/lib/metric_suite.R")

RUNS <- c(`MDS corroborate` = "analyses/ps_final/summary.tsv",
          `MDS off`         = "analyses/ps_nomds_v2/summary.tsv")
univ <- read_tsv("analyses/derived/variant_universe.tsv", show_col_types = FALSE) |>
  filter(source == "VariBench")

res <- bind_rows(lapply(names(RUNS), function(rn) {
  cls <- read_tsv(RUNS[[rn]], show_col_types = FALSE)
  # de-duplicated basis: one row per distinct variant (duplicates are
  # variant x MaveDB-scoreset records and classify identically)
  ded <- cls |> distinct(gene, p_notation, .keep_all = TRUE)
  bind_rows(lapply(c("rows", "variants"), function(basis) {
    tbl <- if (basis == "rows") cls else ded
    j <- tbl |> inner_join(univ, by = c("gene", "p_notation")) |>
      mutate(truth = case_when(grepl("[Pp]athogenic", label) ~ "Pathogenic",
                               grepl("[Bb]enign",     label) ~ "Benign",
                               TRUE ~ NA_character_))
    bind_rows(lapply(c("full", "blind"), function(a) {
      compute_metric_suite(truth = j$truth,
                           pred_class = j[[paste0("varviz_classification_", a)]],
                           pred_score = j[[paste0("varviz_pts_", a)]]) |>
        mutate(run = rn, basis = basis, arm = a, n_joined = nrow(j))
    }))
  }))
}))

show <- res |>
  filter(metric %in% c("Sensitivity","Specificity","MCC","AUROC","VUS_rate")) |>
  select(basis, run, arm, metric, estimate) |>
  tidyr::pivot_wider(names_from = metric, values_from = estimate) |>
  arrange(desc(basis), run, arm)
print(as.data.frame(show), digits = 4, row.names = FALSE)

cat("\njoined N by basis: ",
    paste(unique(paste0(res$basis, "=", res$n_joined)), collapse = "  "), "\n")
write.table(res, "analyses/repro/out/panel_a_runs.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("[repro] wrote analyses/repro/out/panel_a_runs.tsv\n")
