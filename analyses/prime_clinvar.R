#!/usr/bin/env Rscript
# Prime analyses/cache/clinvar/<GENE>_<clinsig>.rds for the benchmark genes.
#
# NCBI eutils truncates under concurrent load instead of erroring: a gene comes
# back with a fraction of its records and the run continues. Six concurrent runs
# produced five different ClinVar record counts for PTEN, which silently thinned
# PM5 and PS1 and collapsed the ClinVar density behind the PM1 hotspot pathway.
#
# Fetching once, serially and politely, then reading from disk removes the
# failure and pins every later run to one ClinVar snapshot, without which two
# runs are not comparable.
#
# Truncation only ever returns FEWER rows than a healthy response, never more,
# and the log counts ("Total records") are raw ClinVar records rather than the
# filtered missense rows this function returns, so no absolute floor is
# available. Instead each gene is fetched several times and the largest result
# kept: agreement across attempts is evidence the response was complete, and a
# disagreement is reported so a short prime is visible rather than silent.
suppressMessages(source("server.R"))

DIR <- CLINVAR_OVERRIDE_DIR
dir.create(DIR, recursive = TRUE, showWarnings = FALSE)

# VARVIZ_PRIME_GENES overrides the default benchmark list, so the RASopathy and
# external-163 cohorts can be primed with the same guarantees. Already-primed
# genes are skipped, so overlapping lists cost nothing.
.env_genes <- Sys.getenv("VARVIZ_PRIME_GENES", "")
genes <- if (nzchar(.env_genes)) {
  trimws(strsplit(.env_genes, "[,[:space:]]+")[[1]])
} else c("BAP1","BRCA1","CASR","DDR2","GCK","KCNH2","KCNQ1","KRAS",
         "LDLR","NUDT15","PTEN","SLC13A5","TP53","TSHR")
genes <- genes[nzchar(genes)]
ATTEMPTS <- as.integer(Sys.getenv("VARVIZ_PRIME_ATTEMPTS", "3"))

disagreed <- character(0)
for (g in genes) for (sig in c("path", "benign")) {
  out <- file.path(DIR, paste0(g, "_", sig, ".rds"))
  if (file.exists(out)) { cat(sprintf("%-9s %-7s cached\n", g, sig)); next }
  best <- NULL; seen <- integer(0)
  for (attempt in seq_len(ATTEMPTS)) {
    # Drop this session's entry so a retry is a real fetch, not a replay of the
    # truncated response we are retrying because of.
    api_cache$clinvar[[paste0(g, ":", sig)]] <- NULL
    df <- tryCatch(extract_clinvar(g, sig), error = function(e) NULL)
    n  <- if (is.null(df)) 0L else nrow(df)
    seen <- c(seen, n)
    if (!is.null(df) && (is.null(best) || n > nrow(best))) best <- df
    Sys.sleep(4)
  }
  if (is.null(best) || nrow(best) == 0L) {
    cat(sprintf("%-9s %-7s EMPTY after %d attempts, not written\n", g, sig, ATTEMPTS))
    disagreed <- c(disagreed, paste0(g, "_", sig, " (empty)"))
    next
  }
  saveRDS(best, out)
  agree <- length(unique(seen)) == 1L
  cat(sprintf("%-9s %-7s %s %5d rows   attempts: %s\n", g, sig,
              if (agree) "OK  " else "KEPT", nrow(best), paste(seen, collapse = "/")))
  if (!agree) disagreed <- c(disagreed, paste0(g, "_", sig, " (", paste(seen, collapse = "/"), ")"))
  Sys.sleep(3)
}
if (length(disagreed)) {
  cat("\nATTEMPTS DISAGREED (largest kept; check before trusting these genes):\n  ",
      paste(disagreed, collapse = ", "), "\n")
} else cat("\nall genes primed, every attempt agreed\n")
cat("PRIME_DONE\n")
