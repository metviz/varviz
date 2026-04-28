# Phase 3 driver — run Panels A / B / C in sequence and emit a unified summary.
#
# Reads dual-pass classifications (analyses/derived/varviz_classifications.tsv
# or per-gene checkpoints if summary not yet present), then sources each panel
# script in turn. Each panel writes its own PDF + HTML; this script just
# orchestrates and surfaces a cross-panel summary at the end.
#
# Run from project root: `Rscript analyses/09_run_all_panels.R`

cat("============================================================\n")
cat("  Pass 2 panel runner — A (clinical) + B (DMS) + C (CAPS)\n")
cat("============================================================\n\n")

panels <- list(
  list(name = "Panel A — VariBench clinical (dual-pass)", script = "analyses/06_panel_a_clinical.R"),
  list(name = "Panel B — MaveDB DMS concordance",         script = "analyses/07_panel_b_dms.R"),
  list(name = "Panel C — CAPS analysis (dual-pass)",      script = "analyses/08_panel_c_caps.R")
)

results <- list()
for (p in panels) {
  cat(sprintf("\n>>> %s\n>>> %s\n\n", p$name, strrep("-", nchar(p$name))))
  t0 <- Sys.time()
  ok <- tryCatch({
    source(p$script, local = TRUE)
    TRUE
  }, error = function(e) {
    cat(sprintf("    PANEL FAILED: %s\n", conditionMessage(e)))
    FALSE
  })
  results[[p$script]] <- list(name = p$name, ok = ok,
                              elapsed = as.numeric(Sys.time() - t0, units = "secs"))
}

cat("\n============================================================\n")
cat("  Run summary\n")
cat("============================================================\n")
for (r in results) {
  cat(sprintf("  %-50s  %s  %.1fs\n",
              r$name,
              if (r$ok) "OK     " else "FAILED ",
              r$elapsed))
}
cat("\nFigures: analyses/figures/{panel_a_clinical,panel_b_dms,panel_c_caps}.pdf\n")
cat("Reports: analyses/reports/{panel_a,panel_b,panel_c}.html\n")
