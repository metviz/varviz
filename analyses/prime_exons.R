#!/usr/bin/env Rscript
# Prime analyses/cache/ensembl/<GENE>_exons.rds.
#
# rest.ensembl.org returns sporadic timeouts under sustained load. A failed exon
# lookup leaves chrom = NA, which disables the dbNSFP and REVEL preloads at once,
# so conservation reads as GERP=NA and PM1/PP3 collapse. The harness refuses to
# checkpoint such a gene, which is the right behaviour but aborts the whole run:
# three external-163 genes (KREMEN1, OTOF, PCDH15) failed a 44-gene run this way
# after the other 41 had completed.
#
# server.R prefers a primed copy over the network (fetch_ensembl_exons ->
# ENSEMBL_OVERRIDE_DIR), so priming takes the flaky remote off the critical path
# for every later run and pins the gene's coordinates across runs.
#
# Retries with backoff, one gene at a time. Writes only a non-empty table, so a
# partial fetch never becomes a poisoned cache.
#
#   VARVIZ_PRIME_GENES="GENE1 GENE2" Rscript analyses/prime_exons.R
suppressMessages(source("server.R"))

DIR <- ENSEMBL_OVERRIDE_DIR
dir.create(DIR, recursive = TRUE, showWarnings = FALSE)

.env_genes <- Sys.getenv("VARVIZ_PRIME_GENES", "")
genes <- if (nzchar(.env_genes)) {
  trimws(strsplit(.env_genes, "[,[:space:]]+")[[1]])
} else c("BAP1","BRCA1","CASR","DDR2","GCK","KCNH2","KCNQ1","KRAS",
         "LDLR","NUDT15","PTEN","SLC13A5","TP53","TSHR","SNCA","SLC16A2")
genes <- unique(genes[nzchar(genes)])

failed <- character(0)
for (g in genes) {
  out <- file.path(DIR, paste0(g, "_exons.rds"))
  if (file.exists(out)) { cat(sprintf("%-10s cached\n", g)); next }
  ok <- FALSE
  for (attempt in 1:5) {
    df <- tryCatch(fetch_ensembl_exons(g), error = function(e) NULL)
    if (!is.null(df) && nrow(df) > 0) {
      saveRDS(df, out)
      cat(sprintf("%-10s OK   %3d exons, chr%s (attempt %d)\n",
                  g, nrow(df), as.character(df$chr[1]), attempt))
      ok <- TRUE; break
    }
    Sys.sleep(attempt * 8)
  }
  if (!ok) { cat(sprintf("%-10s FAILED after 5 attempts\n", g)); failed <- c(failed, g) }
  Sys.sleep(2)
}
if (length(failed)) cat("\nNOT PRIMED:", paste(failed, collapse = ", "), "\n")
cat("PRIME_EXONS_DONE\n")
