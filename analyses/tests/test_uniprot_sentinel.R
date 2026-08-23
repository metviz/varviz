# Regression guard for the silent UniProt corruption that spoiled three benchmark
# runs (ps_mds_corroborate, ps_mds_consfree, and the first ps_final).
#
# Both load-bearing rest.uniprot.org fetches return NULL on failure so the Shiny
# app degrades gracefully. That NULL strips the domain/site PM1 pathways and
# leaves MDS unable to score, while the gene's ROW COUNT stays correct -- so
# row-count validation cannot see it. server.R therefore records the failed
# accession in options(varviz.uniprot_failed), and ps_final_harness.R refuses to
# write a checkpoint for any gene that grew that list.
#
# This asserts the recording actually happens, for BOTH fetches, by pointing them
# at a dead host to force the connection-level failure ("Failed to perform HTTP
# request") that these runs actually hit -- the class that req_retry() ignores
# unless retry_on_failure = TRUE.
#
# Run from project root: Rscript analyses/tests/test_uniprot_sentinel.R
# Slow by design (~5 min): the retry backoff it exercises is real.

# Both UniProt fetches must record into varviz.uniprot_failed on failure, so the
# harness can refuse a checkpoint. Force failure by pointing them at a dead host.
suppressMessages(library(httr2)); suppressMessages(library(jsonlite))
suppressMessages(source("server.R"))

stopifnot(is.null(getOption("varviz.uniprot_failed")))

# Unreachable host -> connection failure, the exact class that broke these runs.
trace_fail <- function(fn) {
  body(fn) <- parse(text = gsub("https://rest\\.uniprot\\.org/uniprotkb/",
                                "https://127.0.0.1:9/", paste(deparse(body(fn)), collapse="\n")))
  fn
}
fetch_pfam_bad     <- trace_fail(extract_pfam)
fetch_geneinfo_bad <- trace_fail(extract_gene_info_uniprot)

options(varviz.uniprot_failed = NULL)
invisible(suppressMessages(fetch_pfam_bad("TEST_PFAM")))
r1 <- getOption("varviz.uniprot_failed", character(0))
cat("pfam fetch     ->", if ("TEST_PFAM" %in% r1) "RECORDED ok" else "NOT RECORDED **FAIL**", "\n")

invisible(suppressMessages(fetch_geneinfo_bad("TEST_GENEINFO", "TESTGENE")))
r2 <- getOption("varviz.uniprot_failed", character(0))
cat("geneinfo fetch ->", if ("TEST_GENEINFO" %in% r2) "RECORDED ok" else "NOT RECORDED **FAIL**", "\n")
cat("accumulates    ->", if (length(r2) == 2) "ok (both retained)" else paste("**FAIL**", length(r2)), "\n")
stopifnot("TEST_PFAM" %in% r2, "TEST_GENEINFO" %in% r2, length(r2) == 2)
cat("\nALL SENTINEL CHECKS PASSED\n")
