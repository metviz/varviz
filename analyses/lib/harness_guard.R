# Failure sentinels for the classification harnesses.
#
# Every harness fetches gene-level sources (Pfam, UniProt, gnomAD, ClinVar,
# CCRS, AlphaFold, mean pathogenicity, gene info) once per gene inside
# tryCatch(..., error = function(e) NULL). A source that FAILED used to be
# indistinguishable from a source that is genuinely ABSENT for that gene, so a
# transient outage dropped the domain/site PM1 pathways while the row count
# stayed correct and the gene was checkpointed as a real result. That is the
# silent contamination documented for ps_mds_corroborate / ps_mds_consfree
# (analyses/repro/README.md "OPEN ISSUE").
#
# harness_fetch() keeps the tryCatch-to-NULL behaviour but records the failure
# in options(varviz.harness_fetch_failed) as "<gene>:<source>". server.R
# maintains its own sentinels for the fetches it owns (varviz.uniprot_failed,
# varviz.clinvar_batch_failed) and ps_final_harness.R adds
# varviz.localpred_failed. sentinel_snapshot() reads all four families so
# run_one() can diff before/after a gene and refuse the checkpoint if anything
# new failed.
#
# Sourced by analyses/05_classify_harness.R and analyses/repro/08_casestudy_harness.R.

HARNESS_SENTINEL_OPTIONS <- c("varviz.harness_fetch_failed",
                              "varviz.uniprot_failed",
                              "varviz.localpred_failed",
                              "varviz.clinvar_batch_failed")

# Evaluate `expr`; on error print the cause, record "<gene>:<source>" and
# return NULL. A NULL *result* without an error is passed through unflagged --
# that is a legitimate "nothing for this gene".
harness_fetch <- function(source, gene, expr) {
  tryCatch(expr, error = function(e) {
    cat(sprintf("    [%s] FETCH FAILED %s: %s\n", gene, source, conditionMessage(e)))
    key <- paste0(gene, ":", source)
    options(varviz.harness_fetch_failed =
              union(getOption("varviz.harness_fetch_failed", character(0)), key))
    NULL
  })
}

# Flat character vector of every recorded sentinel across all families.
sentinel_snapshot <- function() {
  as.character(unlist(lapply(HARNESS_SENTINEL_OPTIONS,
                             function(o) getOption(o, character(0))),
                      use.names = FALSE))
}

# Sentinels present in `after` that were not in `before`.
new_sentinels <- function(before, after) setdiff(after, before)
