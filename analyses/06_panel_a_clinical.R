# Phase 3.1 — Panel A: VariBench clinical concordance, dual-pass (MF1+MF3).
#
# A1 = Pass-Full vs VariBench (acknowledged circular)
# A2 = Pass-Blind vs VariBench (rigorous primary metric)
#
# Reads dual-pass classifications. Prefers analyses/derived/varviz_classifications.tsv
# (final summary) but falls back to concatenating analyses/classifications/*__dual.tsv
# when the harness is partial — useful for incremental sanity checks.
#
# Run from project root: `Rscript analyses/06_panel_a_clinical.R`

suppressMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
})
source("analyses/lib/metric_suite.R")

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

CLASS_SUMMARY <- "analyses/derived/varviz_classifications.tsv"
CHECKPOINT_DIR <- "analyses/classifications"
UNIVERSE       <- "analyses/derived/variant_universe.tsv"
FIG_OUT        <- "analyses/figures/panel_a_clinical.pdf"
HTML_OUT       <- "analyses/reports/panel_a.html"
dir.create(dirname(FIG_OUT),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(HTML_OUT), recursive = TRUE, showWarnings = FALSE)

# Load classifications: prefer summary, fall back to per-gene checkpoints.
cls <- if (file.exists(CLASS_SUMMARY)) {
  cat(sprintf("[panel_a] Reading summary %s\n", CLASS_SUMMARY))
  read_tsv(CLASS_SUMMARY, show_col_types = FALSE)
} else {
  ckpt_files <- list.files(CHECKPOINT_DIR, pattern = "__dual\\.tsv$", full.names = TRUE)
  if (length(ckpt_files) == 0L) stop("No classification output yet — run 05_classify_harness.R first.")
  cat(sprintf("[panel_a] Reading %d per-gene checkpoints\n", length(ckpt_files)))
  bind_rows(lapply(ckpt_files, read_tsv, show_col_types = FALSE))
}
univ <- read_tsv(UNIVERSE, show_col_types = FALSE)

# Restrict to VariBench rows (have a clinical label) and join with classifications.
joined <- cls |>
  inner_join(univ |> filter(source == "VariBench"),
             by = c("gene", "p_notation")) |>
  mutate(truth = case_when(
    grepl("[Pp]athogenic", label) ~ "Pathogenic",
    grepl("[Bb]enign",     label) ~ "Benign",
    TRUE                          ~ NA_character_
  ))

cat(sprintf("[panel_a] VariBench-labeled variants with classifications: %d\n", nrow(joined)))
cat(sprintf("  Pathogenic: %d   Benign: %d\n",
            sum(joined$truth == "Pathogenic", na.rm = TRUE),
            sum(joined$truth == "Benign",     na.rm = TRUE)))

# A1 — Pass-Full
metrics_full <- compute_metric_suite(
  truth      = joined$truth,
  pred_class = joined$varviz_classification_full,
  pred_score = joined$varviz_pts_full
) |> mutate(pass = "Full")

# A2 — Pass-Blind
metrics_blind <- compute_metric_suite(
  truth      = joined$truth,
  pred_class = joined$varviz_classification_blind,
  pred_score = joined$varviz_pts_blind
) |> mutate(pass = "Blind")

# When the harness output is partial and the available variants are all one
# truth class (or all VUS-classified), compute_metric_suite returns 0 rows.
# Defer pivot/delta in that case and surface a clear note.
have_metrics <- nrow(metrics_full) > 0 && nrow(metrics_blind) > 0
delta <- NULL
if (have_metrics) {
  both <- bind_rows(metrics_full, metrics_blind) |>
    pivot_wider(names_from = pass, values_from = c(estimate, lower, upper))
  delta <- both |>
    filter(metric %in% c("Sensitivity", "Specificity", "AUROC", "MCC")) |>
    mutate(delta = estimate_Full - estimate_Blind)
}

cat("\n--- Panel A1: Pass-Full (acknowledged circular) ---\n");      print(metrics_full)
cat("\n--- Panel A2: Pass-Blind (rigorous primary metric) ---\n");   print(metrics_blind)
if (have_metrics) {
  cat("\n--- Δ (Full − Blind): drop attributable to ClinVar criteria ---\n"); print(delta)
} else {
  cat("\n[panel_a] No evaluable rows yet (need >=1 P + >=1 B with non-VUS classification).\n")
  cat("  This is expected for partial harness output — re-run after more genes complete.\n")
}

cm_full  <- build_confusion_3x2(joined$truth, joined$varviz_classification_full)
cm_blind <- build_confusion_3x2(joined$truth, joined$varviz_classification_blind)
cat("\nA1 confusion (Pass-Full):\n");  print(cm_full)
cat("\nA2 confusion (Pass-Blind):\n"); print(cm_blind)

# Figure: side-by-side confusion heatmaps + Δ headline
make_cm_df <- function(cm, label) {
  as_tibble(as.data.frame(cm)) |> mutate(pass = label)
}
cm_long <- bind_rows(
  make_cm_df(cm_full,  "A1: Pass-Full"),
  make_cm_df(cm_blind, "A2: Pass-Blind")
)

if (!is.null(delta) && nrow(delta) > 0) {
  dsens  <- delta$delta[delta$metric == "Sensitivity"]
  dauroc <- delta$delta[delta$metric == "AUROC"]
  dmcc   <- delta$delta[delta$metric == "MCC"]
  subtitle <- sprintf("ΔSens=%+.3f  ΔAUROC=%+.3f  ΔMCC=%+.3f  (Full − Blind; positive = circular gain)",
                      dsens %||% NA, dauroc %||% NA, dmcc %||% NA)
} else {
  subtitle <- "Awaiting more harness output (need >=1 P + >=1 B with non-VUS classification)"
}

p <- ggplot(cm_long, aes(Truth, VarViz, fill = Freq)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = Freq), size = 5) +
  scale_fill_gradient(low = "#f1f5f9", high = "#003f5c") +
  facet_wrap(~ pass, ncol = 2) +
  labs(title = NULL, subtitle = NULL, x = "VariBench truth", y = "VarViz call") +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none",
        panel.grid = element_blank(),
        strip.text = element_text(face = "bold"))

ggsave(FIG_OUT, p, width = 10, height = 5, units = "in")
ggsave(sub("\\.pdf$", ".png", FIG_OUT), p, width = 10, height = 5, units = "in", dpi = 300)

# HTML companion report
htmltools::save_html(
  htmltools::tagList(
    htmltools::h2("Panel A — VariBench clinical concordance (dual-pass)"),
    htmltools::p(sprintf("VariBench-labeled variants with classifications: %d (P=%d, B=%d)",
                         nrow(joined),
                         sum(joined$truth == "Pathogenic", na.rm = TRUE),
                         sum(joined$truth == "Benign",     na.rm = TRUE))),
    htmltools::h3("A1 — Pass-Full (acknowledged circular)"),
    htmltools::tags$pre(paste(capture.output(print(metrics_full)),  collapse = "\n")),
    htmltools::tags$pre(paste(capture.output(print(cm_full)),       collapse = "\n")),
    htmltools::h3("A2 — Pass-Blind (rigorous primary metric)"),
    htmltools::tags$pre(paste(capture.output(print(metrics_blind)), collapse = "\n")),
    htmltools::tags$pre(paste(capture.output(print(cm_blind)),      collapse = "\n")),
    htmltools::h3("Δ (Full − Blind)"),
    htmltools::tags$pre(paste(
      if (!is.null(delta) && nrow(delta) > 0) capture.output(print(delta))
      else "Awaiting more harness output (need >=1 P + >=1 B with non-VUS classification).",
      collapse = "\n"))
  ),
  HTML_OUT
)

cat(sprintf("\n[panel_a] Wrote %s and %s\n", FIG_OUT, HTML_OUT))
