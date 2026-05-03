# Phase 3.3b — Panel C: CAPS analysis (orthogonal population-genetic track).
#
# Restrict to gnomAD-present variants from the universe; bin each by VarViz
# classification (both Pass-Full and Pass-Blind for symmetry); compute
# observed singleton proportion + Wilson 95% CI per bin; CAPS =
# (observed - expected) / expected. Spearman rho between ordinal bin rank
# and CAPS tests whether the classifier is monotonic with negative
# selection signal.
#
# Deviation from plan: plan called for Gudkov 2025 published per-consequence
# expected proportions. Their GitHub repo (mgudVCCRI/PopGenVariantFiltering)
# provides the Snakemake pipeline that DERIVES these from gnomAD VCFs but
# does not ship pre-computed values. Running the pipeline requires gnomAD
# VCFs (~TB) — out of scope. Pragmatic substitute: use the pooled baseline
# singleton fraction across all gnomAD-present universe variants as the
# "expected". Same interpretation: per-bin enrichment vs. an
# annotation-agnostic background. The CAPS sign and Spearman monotonicity
# remain meaningful; only the absolute CAPS magnitude shifts under a
# different baseline.
#
# Run from project root: `Rscript analyses/08_panel_c_caps.R`

suppressMessages({
  library(dplyr); library(readr); library(ggplot2); library(tidyr)
})
source("analyses/lib/caps_compute.R")

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

CLASS_SUMMARY  <- "analyses/derived/varviz_classifications.tsv"
CHECKPOINT_DIR <- "analyses/classifications"
GNOMAD_IN      <- "analyses/derived/variant_universe_gnomad.tsv"
FIG_OUT        <- "analyses/figures/panel_c_caps.pdf"
HTML_OUT       <- "analyses/reports/panel_c.html"
dir.create(dirname(FIG_OUT),  recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(HTML_OUT), recursive = TRUE, showWarnings = FALSE)

# Load classifications (prefer summary, fall back to checkpoints).
cls <- if (file.exists(CLASS_SUMMARY)) {
  cat(sprintf("[panel_c] Reading summary %s\n", CLASS_SUMMARY))
  read_tsv(CLASS_SUMMARY, show_col_types = FALSE)
} else {
  ckpt_files <- list.files(CHECKPOINT_DIR, pattern = "__dual\\.tsv$", full.names = TRUE)
  if (length(ckpt_files) == 0L) stop("No classification output yet — run 05_classify_harness.R first.")
  cat(sprintf("[panel_c] Reading %d per-gene checkpoints\n", length(ckpt_files)))
  bind_rows(lapply(ckpt_files, read_tsv, show_col_types = FALSE))
}
gn <- read_tsv(GNOMAD_IN, show_col_types = FALSE)

# Restrict to gnomAD-present variants joined with both passes.
joined <- cls |>
  inner_join(
    gn |> select(gene, p_notation, gnomad_present, gnomad_singleton),
    by = c("gene", "p_notation")
  ) |>
  filter(gnomad_present)

cat(sprintf("[panel_c] gnomAD-present variants with classifications: %d\n", nrow(joined)))
cat(sprintf("  Singletons in pool: %d (%.1f%%)\n",
            sum(joined$gnomad_singleton, na.rm = TRUE),
            100 * mean(joined$gnomad_singleton, na.rm = TRUE)))

bin_order <- c("Benign", "Likely Benign", "VUS-Low", "VUS-Mid", "VUS-High",
               "Likely Pathogenic", "Pathogenic")

# Pooled baseline = overall singleton fraction across all gnomAD-present
# variants in the analysis (annotation-agnostic). Documented deviation
# from Gudkov 2025 per-consequence baselines.
baseline <- if (nrow(joined) > 0) {
  mean(joined$gnomad_singleton, na.rm = TRUE)
} else NA_real_
cat(sprintf("[panel_c] Pooled-baseline expected singleton fraction: %.4f\n", baseline %||% NA))

# Helper: compute_caps_table for one classification pass.
caps_for_pass <- function(class_col_name, pass_label) {
  if (is.na(baseline)) return(NULL)
  df <- joined |>
    mutate(bin = .data[[class_col_name]]) |>
    filter(!is.na(bin), bin %in% bin_order)
  if (nrow(df) == 0) return(NULL)
  expected_per_bin <- setNames(rep(baseline, length(bin_order)), bin_order)
  tab <- compute_caps_table(
    df               = df,
    bin_col          = "bin",
    singleton_col    = "gnomad_singleton",
    expected_per_bin = expected_per_bin
  )
  tab$bin <- factor(tab$bin, levels = bin_order)
  tab <- tab[!is.na(tab$bin), , drop = FALSE]
  tab <- tab[order(tab$bin), , drop = FALSE]
  tab$pass <- pass_label
  tab
}

caps_full  <- caps_for_pass("varviz_classification_full",  "Full")
caps_blind <- caps_for_pass("varviz_classification_blind", "Blind")

if (is.null(caps_full) && is.null(caps_blind)) {
  cat("[panel_c] No gnomAD-present classifications available yet — skipping figure.\n")
  htmltools::save_html(
    htmltools::tagList(
      htmltools::h2("Panel C — CAPS analysis"),
      htmltools::p("Awaiting harness output. No gnomAD-present rows in classifications yet.")
    ),
    HTML_OUT
  )
  cat(sprintf("[panel_c] Wrote %s (placeholder)\n", HTML_OUT))
  quit(save = "no")
}

cat("\n--- Panel C CAPS table (Pass-Full) ---\n");  if (!is.null(caps_full))  print(caps_full)
cat("\n--- Panel C CAPS table (Pass-Blind) ---\n"); if (!is.null(caps_blind)) print(caps_blind)

# Spearman rho per pass: ordinal bin rank vs CAPS
spearman_rho <- function(tab) {
  if (is.null(tab) || nrow(tab) < 3) return(NA_real_)
  ranks <- as.integer(tab$bin)
  caps_vals <- tab$caps
  if (sum(!is.na(caps_vals)) < 3) return(NA_real_)
  suppressWarnings(cor(ranks, caps_vals, method = "spearman", use = "complete.obs"))
}
rho_full  <- spearman_rho(caps_full)
rho_blind <- spearman_rho(caps_blind)
cat(sprintf("\n[panel_c] Spearman rho (bin rank vs CAPS): Full=%.3f  Blind=%.3f\n",
            rho_full %||% NA, rho_blind %||% NA))

# Combined plot: side-by-side dual-pass
both <- bind_rows(caps_full, caps_blind)
both$pass <- factor(both$pass, levels = c("Full", "Blind"),
                    labels = c("Pass-Full", "Pass-Blind"))

n_total_full  <- if (!is.null(caps_full))  sum(caps_full$n)  else 0
n_total_blind <- if (!is.null(caps_blind)) sum(caps_blind$n) else 0

p <- ggplot(both, aes(bin, caps, fill = bin)) +
  geom_col() +
  geom_errorbar(aes(ymin = (ci_lower - p_expected) / p_expected,
                    ymax = (ci_upper - p_expected) / p_expected),
                width = 0.3, linewidth = 0.4, color = "gray30") +
  geom_text(aes(label = sprintf("n=%d", n), y = pmax(caps, 0, na.rm = TRUE) + 0.05),
            size = 2.8, na.rm = TRUE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  facet_wrap(~ pass, ncol = 2) +
  scale_fill_manual(values = c(
    "Benign"            = "#10b981",
    "Likely Benign"     = "#34d399",
    "VUS-Low"           = "#a3e635",
    "VUS-Mid"           = "#facc15",
    "VUS-High"          = "#f97316",
    "Likely Pathogenic" = "#f43f5e",
    "Pathogenic"        = "#dc2626"
  )) +
  labs(title = NULL, subtitle = NULL,
       x = "VarViz classification (Benign -> Pathogenic)",
       y = "CAPS = (observed - baseline) / baseline") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 30, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"))

ggsave(FIG_OUT, p, width = 11, height = 5, units = "in")
ggsave(sub("\\.pdf$", ".png", FIG_OUT), p, width = 11, height = 5, units = "in", dpi = 300)

htmltools::save_html(
  htmltools::tagList(
    htmltools::h2("Panel C - CAPS by VarViz classification bin (gnomAD-present)"),
    htmltools::p(sprintf("Pooled baseline singleton fraction: %.4f", baseline %||% NA)),
    htmltools::p(sprintf("Pass-Full Spearman rho: %.3f  |  Pass-Blind Spearman rho: %.3f",
                         rho_full %||% NA, rho_blind %||% NA)),
    htmltools::h3("Pass-Full"),
    htmltools::tags$pre(paste(capture.output(print(caps_full %||% "")),  collapse = "\n")),
    htmltools::h3("Pass-Blind"),
    htmltools::tags$pre(paste(capture.output(print(caps_blind %||% "")), collapse = "\n"))
  ),
  HTML_OUT
)

cat(sprintf("\n[panel_c] Wrote %s and %s\n", FIG_OUT, HTML_OUT))
