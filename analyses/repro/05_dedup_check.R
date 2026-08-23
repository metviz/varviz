#!/usr/bin/env Rscript
# Is the benchmark counting variants, or variant x scoreset rows?
#
# variant_universe.tsv is built from MaveDB (89,957 rows) + VariBench (744).
# MaveDB contributes the same protein variant once per scoreset, so the 90,701
# "variants" are rows, not distinct variants. This script quantifies the effect
# and recomputes Panel A on distinct variants.
#
# Run from the repository root:  Rscript analyses/repro/05_dedup_check.R

suppressMessages({ library(dplyr); library(readr) })
source("analyses/lib/metric_suite.R")

cls  <- read_tsv("analyses/derived/varviz_classifications_mds.tsv", show_col_types = FALSE)
univ <- read_tsv("analyses/derived/variant_universe.tsv",           show_col_types = FALSE)

cat(sprintf("classification rows           : %d\n", nrow(cls)))
cat(sprintf("distinct gene + p_notation    : %d\n",
            n_distinct(paste(cls$gene, cls$p_notation))))

# Are duplicate rows for one variant classified identically? If so, dedup is lossless.
incons <- cls |>
  group_by(gene, p_notation) |>
  summarise(k = n(),
            n_full  = n_distinct(varviz_classification_full),
            n_blind = n_distinct(varviz_classification_blind_mds), .groups = "drop") |>
  filter(k > 1, n_full > 1 | n_blind > 1)
cat(sprintf("duplicate keys classified inconsistently: %d\n", nrow(incons)))

dedup <- cls |> distinct(gene, p_notation, .keep_all = TRUE)

panel <- function(cl_tbl, tag) {
  j <- cl_tbl |>
    inner_join(univ |> filter(source == "VariBench"), by = c("gene", "p_notation")) |>
    mutate(truth = case_when(grepl("[Pp]athogenic", label) ~ "Pathogenic",
                             grepl("[Bb]enign",     label) ~ "Benign",
                             TRUE ~ NA_character_))
  cat(sprintf("\n%s: %d joined rows (P %d / B %d)\n", tag, nrow(j),
              sum(j$truth == "Pathogenic", na.rm = TRUE),
              sum(j$truth == "Benign",     na.rm = TRUE)))
  bind_rows(lapply(c("full", "blind", "blind_mds"), function(a) {
    compute_metric_suite(truth = j$truth,
                         pred_class = j[[paste0("varviz_classification_", a)]],
                         pred_score = j[[paste0("varviz_pts_",           a)]]) |>
      mutate(arm = a, basis = tag)
  }))
}

res <- bind_rows(panel(cls, "as-published (rows)"), panel(dedup, "de-duplicated (variants)"))

cmp <- res |>
  filter(metric %in% c("Sensitivity","Specificity","MCC","AUROC","VUS_rate")) |>
  select(metric, arm, basis, estimate) |>
  tidyr::pivot_wider(names_from = basis, values_from = estimate) |>
  mutate(delta = `de-duplicated (variants)` - `as-published (rows)`)
cat("\n")
print(as.data.frame(cmp), digits = 4, row.names = FALSE)

write.table(res, "analyses/repro/out/panel_a_dedup.tsv", sep = "\t",
            row.names = FALSE, quote = FALSE)
cat("\n[repro] wrote analyses/repro/out/panel_a_dedup.tsv\n")
