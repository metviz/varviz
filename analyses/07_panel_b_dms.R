# Phase 3.2 — Panel B: MaveDB DMS concordance across BENCHMARK_GENES.
#
# DMS truth is functional and orthogonal to ClinVar circularity, so this is a
# single-pass analysis (Pass-Full only). For each MaveDB study URN, score_raw
# is collapsed into a binary functional truth (LOF -> "Pathogenic",
# WT-like -> "Benign"). Variants in the middle 50% per study are excluded
# from the binary metric calc but counted in the descriptive bin plot.
#
# Threshold scheme (deviation from plan):
# The plan listed paper-specific thresholds keyed by paper-name strings
# (giacomelli_2018, findlay_2018, ...), but our universe carries MaveDB
# scoreset URN suffixes (e.g. 00000003-a-2) — there is no URN -> paper
# mapping in the ingest yet. Until that mapping is added, we use a
# data-driven per-study quartile split: bottom quartile = LOF, top quartile
# = WT-like, middle 50% = Partial. Convention: lower score = more damaging,
# the dominant convention across MaveDB function/fitness assays. Inverted
# assays (e.g. abundance-as-stability where higher = more stable) would
# require per-URN sign correction.
#
# Run from project root: `Rscript analyses/07_panel_b_dms.R`

suppressMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
})
source("analyses/lib/metric_suite.R")

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

CLASS_SUMMARY  <- "analyses/derived/varviz_classifications.tsv"
CHECKPOINT_DIR <- "analyses/classifications"
UNIVERSE       <- "analyses/derived/variant_universe.tsv"
FIG_OUT        <- "analyses/figures/panel_b_dms.pdf"
HTML_OUT       <- "analyses/reports/panel_b.html"
dir.create(dirname(FIG_OUT),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(HTML_OUT), recursive = TRUE, showWarnings = FALSE)

# Load classifications: prefer summary, fall back to per-gene checkpoints.
cls <- if (file.exists(CLASS_SUMMARY)) {
  cat(sprintf("[panel_b] Reading summary %s\n", CLASS_SUMMARY))
  read_tsv(CLASS_SUMMARY, show_col_types = FALSE)
} else {
  ckpt_files <- list.files(CHECKPOINT_DIR, pattern = "__dual\\.tsv$", full.names = TRUE)
  if (length(ckpt_files) == 0L) stop("No classification output yet — run 05_classify_harness.R first.")
  cat(sprintf("[panel_b] Reading %d per-gene checkpoints\n", length(ckpt_files)))
  bind_rows(lapply(ckpt_files, read_tsv, show_col_types = FALSE))
}
univ <- read_tsv(UNIVERSE, show_col_types = FALSE)

# Restrict to MaveDB rows (have score_raw populated) and join with classifications.
mavedb <- univ |> filter(source == "MaveDB", !is.na(score_raw))
joined <- mavedb |> inner_join(cls, by = c("gene", "p_notation"))
cat(sprintf("[panel_b] MaveDB variants with classifications: %d (over %d studies)\n",
            nrow(joined), n_distinct(joined$study)))

# Per-study quartile-based binning. score column is named score_raw.x in the
# universe-side of the join (universe.tsv has score_raw and the cls has none).
score_col <- if ("score_raw" %in% names(joined)) {
  "score_raw"
} else if ("score_raw.x" %in% names(joined)) {
  "score_raw.x"
} else {
  stop("No score_raw column found after join")
}
joined$score_dms <- joined[[score_col]]

joined <- joined |>
  group_by(study) |>
  mutate(
    q25      = quantile(score_dms, 0.25, na.rm = TRUE),
    q75      = quantile(score_dms, 0.75, na.rm = TRUE),
    dms_bin  = case_when(
      score_dms <= q25                    ~ "LOF",
      score_dms >= q75                    ~ "WT-like",
      score_dms >  q25 & score_dms <  q75 ~ "Partial",
      TRUE                                ~ NA_character_
    ),
    truth_dms = case_when(
      dms_bin == "LOF"     ~ "Pathogenic",
      dms_bin == "WT-like" ~ "Benign",
      TRUE                 ~ NA_character_
    )
  ) |> ungroup() |>
  mutate(dms_bin = factor(dms_bin, levels = c("LOF", "Partial", "WT-like")))

cat("\n[panel_b] DMS bin counts (per-study quartile split):\n")
print(joined |> count(dms_bin, useNA = "ifany"))

# --- Metric suite: per-study + pooled ---
metrics_per_study <- joined |>
  group_by(gene, study) |>
  group_modify(~ compute_metric_suite(
    truth      = .x$truth_dms,
    pred_class = .x$varviz_classification_full,
    pred_score = .x$varviz_pts_full
  )) |> ungroup()

metrics_pooled <- compute_metric_suite(
  truth      = joined$truth_dms,
  pred_class = joined$varviz_classification_full,
  pred_score = joined$varviz_pts_full
)

cat("\n--- Panel B per-study metrics ---\n"); print(metrics_per_study, n = Inf)
cat("\n--- Panel B pooled metrics ---\n");    print(metrics_pooled)

# --- Descriptive bar plot: per-bin VarViz Pathogenic-call rate ---
joined <- joined |> mutate(
  varviz_path = varviz_classification_full %in% c("Pathogenic", "Likely Pathogenic")
)

summary_df <- joined |>
  filter(!is.na(dms_bin)) |>
  group_by(gene, study, dms_bin) |>
  summarize(n      = n(),
            n_path = sum(varviz_path),
            pct_path = if (n > 0) n_path / n else NA_real_,
            .groups = "drop")

cat("\n[panel_b] Per-study, per-bin VarViz P-call rate:\n")
print(summary_df, n = Inf)

# Headline pooled metrics for the figure subtitle
auc_row  <- metrics_pooled[metrics_pooled$metric == "AUROC", ]
sens_row <- metrics_pooled[metrics_pooled$metric == "Sensitivity", ]
spec_row <- metrics_pooled[metrics_pooled$metric == "Specificity", ]

subtitle <- if (nrow(metrics_pooled) > 0 && !is.na(auc_row$estimate %||% NA)) {
  sprintf("Pooled AUROC=%.3f (95%% CI %.3f-%.3f) | Sens=%.3f | Spec=%.3f",
          auc_row$estimate, auc_row$lower, auc_row$upper,
          sens_row$estimate %||% NA, spec_row$estimate %||% NA)
} else {
  "Awaiting more harness output (need >=1 LOF + >=1 WT-like with non-VUS classification)"
}

if (nrow(summary_df) > 0) {
  p <- ggplot(summary_df, aes(dms_bin, pct_path, fill = dms_bin)) +
    geom_col() +
    geom_text(aes(label = sprintf("%d/%d", n_path, n)), vjust = -0.3, size = 3) +
    facet_wrap(~ gene + study, ncol = 4, labeller = label_both) +
    scale_fill_manual(values = c("LOF" = "#ef4444", "Partial" = "#f59e0b", "WT-like" = "#10b981")) +
    scale_y_continuous(labels = scales::percent_format(), limits = c(0, 1.10)) +
    labs(title    = "Panel B - VarViz Pathogenic call rate by MaveDB DMS bin",
         subtitle = subtitle,
         x        = "DMS bin (per-study quartile split)",
         y        = "VarViz P / LP (%)") +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none",
          strip.text = element_text(size = 8))

  n_studies <- n_distinct(summary_df$study)
  fig_h <- min(14, 2 + 1.6 * ceiling(n_studies / 4))
  ggsave(FIG_OUT, p, width = 14, height = fig_h, units = "in", limitsize = FALSE)
} else {
  cat("[panel_b] No MaveDB rows in classifications yet — skipping figure.\n")
}

# HTML companion report
htmltools::save_html(
  htmltools::tagList(
    htmltools::h2("Panel B - MaveDB DMS concordance (multi-study)"),
    htmltools::p(sprintf("Studies: %d   Variants joined to classifications: %d",
                         n_distinct(joined$study), nrow(joined))),
    htmltools::h3("Pooled metric suite"),
    htmltools::tags$pre(paste(capture.output(print(metrics_pooled)),    collapse = "\n")),
    htmltools::h3("Per-study metric suite"),
    htmltools::tags$pre(paste(capture.output(print(metrics_per_study, n = Inf)), collapse = "\n")),
    htmltools::h3("Per-study, per-bin VarViz P-call rate"),
    htmltools::tags$pre(paste(capture.output(print(summary_df, n = Inf)), collapse = "\n"))
  ),
  HTML_OUT
)

cat(sprintf("\n[panel_b] Wrote %s and %s\n", FIG_OUT, HTML_OUT))
