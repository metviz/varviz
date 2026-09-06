# Self-check for analyses/lib/harness_guard.R.
#
# The harnesses fetch eight gene-level sources per gene inside tryCatch(...,
# error = function(e) NULL). A failed source used to look identical to an
# absent one, so a UniProt or Ensembl outage was checkpointed as a real
# classification (the ps_mds_corroborate / ps_mds_consfree contamination).
# harness_fetch() records every failure in an option() sentinel and
# sentinel_snapshot() / new_sentinels() let run_one() refuse the checkpoint.
#
# Run from project root: Rscript analyses/tests/test_harness_guard.R

source("analyses/lib/harness_guard.R")

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

reset <- function() options(varviz.harness_fetch_failed = NULL,
                            varviz.uniprot_failed = NULL,
                            varviz.localpred_failed = NULL,
                            varviz.clinvar_batch_failed = NULL)

# 1. success passes the value through and records nothing
reset()
v <- harness_fetch("pfam", "GENE1", 42)
check(identical(v, 42), "harness_fetch returns the value on success")
check(is.null(getOption("varviz.harness_fetch_failed")), "no sentinel on success")

# 2. failure returns NULL and records gene:source
reset()
out <- capture.output(v <- harness_fetch("uniprot", "GENE1", stop("boom")))
check(is.null(v), "harness_fetch returns NULL on error")
check(identical(getOption("varviz.harness_fetch_failed"), "GENE1:uniprot"),
      "failure recorded as gene:source")
check(any(grepl("GENE1.*uniprot.*boom", out)), "failure message names gene, source and error")

# 3. repeated failures accumulate without duplicates
out <- capture.output({
  harness_fetch("gnomad", "GENE1", stop("x"))
  harness_fetch("gnomad", "GENE1", stop("x"))
})
check(identical(sort(getOption("varviz.harness_fetch_failed")),
                c("GENE1:gnomad", "GENE1:uniprot")), "sentinels accumulate, deduplicated")

# 4. snapshot / diff sees every sentinel family, including server.R's own
reset()
before <- sentinel_snapshot()
options(varviz.uniprot_failed = "P12345")
out <- capture.output(harness_fetch("ccrs", "GENE2", stop("y")))
after <- sentinel_snapshot()
d <- new_sentinels(before, after)
check(setequal(d, c("P12345", "GENE2:ccrs")), "new_sentinels reports server.R and harness sentinels")
check(length(new_sentinels(after, after)) == 0L, "no diff when nothing new failed")

# 5. the classic silent NULL: value NULL but no error is NOT a failure
reset()
v <- harness_fetch("af", "GENE3", NULL)
check(is.null(v) && is.null(getOption("varviz.harness_fetch_failed")),
      "a legitimate NULL result is not flagged")

reset()
cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
