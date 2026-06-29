#!/usr/bin/env Rscript
# 12_dolphin_delta_report.R — Pass-Blind vs Pass-Blind+DOLPHIN delta on the
# Panel A (VariBench clinical) and Panel C (CAPS population-genetic) metrics.
#
# Panel B (MaveDB DMS) is unchanged because it scores Pass-FULL output, which
# DOLPHIN does not alter — DOLPHIN augments only Pass-Blind PM1 where the
# CCRS / UniProt-domain / ClinVar-hotspot pathways did not fire.

suppressPackageStartupMessages({
  library(data.table)
  source("analyses/lib/metric_suite.R", local = FALSE)
  source("analyses/lib/caps_compute.R", local = FALSE)
})

CLASS_TSV       <- "analyses/derived/varviz_classifications_dolphin.tsv"
UNIVERSE_TSV    <- "analyses/derived/variant_universe.tsv"        # for Panel A (VariBench-labelled subset)
UNIVERSE_GN_TSV <- "analyses/derived/variant_universe_gnomad.tsv" # for Panel C (gnomAD-present subset)

cat("============================================================\n")
cat(" DOLPHIN delta report: Pass-Blind vs Pass-Blind + DOLPHIN\n")
cat("============================================================\n")

cls <- fread(CLASS_TSV, sep = "\t", header = TRUE)

# Sanity check: how many rows were re-binned by DOLPHIN augmentation?
rebinned <- sum(cls$varviz_classification_blind != cls$varviz_classification_blind_dolphin)
cat("\n[delta] universe rows:", nrow(cls), "\n")
cat("[delta] rebinned by DOLPHIN:", rebinned,
    sprintf("(%.2f%%)\n", 100 * rebinned / nrow(cls)))

# ---------------------------------------------------------------------------
# Panel A: VariBench clinical concordance
# ---------------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Panel A: VariBench clinical concordance\n")
cat(strrep("=", 60), "\n", sep = "")

uni_full <- fread(UNIVERSE_TSV, sep = "\t", header = TRUE)
vb <- uni_full[source == "VariBench" & label %in% c("Pathogenic", "Benign"),
               .(gene, p_notation, label)]
vb <- unique(vb, by = c("gene", "p_notation"))

setkey(cls, gene, p_notation)
setkey(vb,  gene, p_notation)
joined <- vb[cls, nomatch = 0]

cat("[A] VariBench-labelled rows in universe:", nrow(vb), "\n")
cat("[A] joined to engine:", nrow(joined), "\n")
cat("[A] label split: ",
    sum(joined$label == "Pathogenic"), " P / ",
    sum(joined$label == "Benign"),     " B\n", sep = "")

cat("\n--- Pass-Blind (original) ---\n")
m_blind <- compute_metric_suite(
  truth      = joined$label,
  pred_class = joined$varviz_classification_blind,
  pred_score = joined$varviz_pts_blind
)
print(m_blind)

cat("\n--- Pass-Blind + DOLPHIN ---\n")
m_dolphin <- compute_metric_suite(
  truth      = joined$label,
  pred_class = joined$varviz_classification_blind_dolphin,
  pred_score = joined$varviz_pts_blind_dolphin
)
print(m_dolphin)

# ---------------------------------------------------------------------------
# Panel C: CAPS population-genetic concordance
# ---------------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(" Panel C: CAPS singleton-enrichment concordance\n")
cat(strrep("=", 60), "\n", sep = "")

uni <- fread(UNIVERSE_GN_TSV, sep = "\t", header = TRUE)

# Distinct rows; first occurrence wins (some variants appear in multiple sets)
uni <- unique(uni, by = c("gene", "p_notation"))
setkey(uni, gene, p_notation)

joined_c <- uni[cls, nomatch = 0]
joined_c <- joined_c[gnomad_present == TRUE | gnomad_present == "TRUE" | gnomad_ac >= 1]

cat("[C] gnomAD-present rows:", nrow(joined_c), "\n")

baseline <- mean(joined_c$gnomad_singleton, na.rm = TRUE)
cat("[C] pooled-baseline singleton fraction:",
    sprintf("%.4f", baseline), "\n")

bin_order <- c("Benign", "Likely Benign", "VUS-Low", "VUS-Mid", "VUS-High",
               "Likely Pathogenic", "Pathogenic")

caps_for_pass <- function(col, label) {
  expected_per_bin <- setNames(rep(baseline, length(bin_order)), bin_order)
  tab <- compute_caps_table(
    df               = data.frame(
      bin              = joined_c[[col]],
      gnomad_singleton = joined_c$gnomad_singleton
    ),
    bin_col          = "bin",
    singleton_col    = "gnomad_singleton",
    expected_per_bin = expected_per_bin
  )
  tab <- tab[match(bin_order, tab$bin), ]
  tab <- tab[!is.na(tab$bin), ]
  cat("[C/", label, "] per-bin table:\n", sep = "")
  print(tab, row.names = FALSE, digits = 3)
  ranks <- match(tab$bin, bin_order)
  rho <- suppressWarnings(cor(ranks, tab$caps, method = "spearman",
                              use = "complete.obs"))
  cat("[C/", label, "] Spearman ρ (bin rank vs CAPS) = ",
      sprintf("%.3f", rho), "\n", sep = "")
  list(table = tab, rho = rho)
}

cat("\n--- Pass-Blind (original) ---\n")
c_blind <- caps_for_pass("varviz_classification_blind", "Blind")

cat("\n--- Pass-Blind + DOLPHIN ---\n")
c_dolphin <- caps_for_pass("varviz_classification_blind_dolphin", "Dolphin")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
cat("\n", strrep("=", 60), "\n", sep = "")
cat(" SUMMARY\n")
cat(strrep("=", 60), "\n", sep = "")

get_metric <- function(m, k) {
  # compute_metric_suite returns a tibble with columns (metric, estimate, lower, upper).
  # Lookup the row whose metric value matches k and return its estimate.
  if (is.null(m) || !"metric" %in% colnames(m)) return(NA_real_)
  v <- m$estimate[m$metric == k]
  if (length(v) == 0) return(NA_real_)
  suppressWarnings(as.numeric(v[1]))
}

summary_tbl <- data.frame(
  panel = c("A: VariBench AUROC",
            "A: VariBench MCC",
            "A: VariBench Sensitivity",
            "A: VariBench Specificity",
            "A: VariBench VUS rate",
            "C: CAPS Spearman ρ"),
  blind   = c(get_metric(m_blind, "AUROC"),
              get_metric(m_blind, "MCC"),
              get_metric(m_blind, "Sensitivity"),
              get_metric(m_blind, "Specificity"),
              get_metric(m_blind, "VUS_rate"),
              c_blind$rho),
  dolphin = c(get_metric(m_dolphin, "AUROC"),
              get_metric(m_dolphin, "MCC"),
              get_metric(m_dolphin, "Sensitivity"),
              get_metric(m_dolphin, "Specificity"),
              get_metric(m_dolphin, "VUS_rate"),
              c_dolphin$rho),
  stringsAsFactors = FALSE
)
summary_tbl$delta <- summary_tbl$dolphin - summary_tbl$blind
print(summary_tbl, row.names = FALSE, digits = 4)

cat("\nDONE.\n")
