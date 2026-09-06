# Offline self-check for analyses/lib/dolphin_bulk.R failure handling.
#
# fetch_dolphin() returns NULL when the API call fails (timeout, 5xx, 429
# budget exhausted). fetch_dolphin_gene() used to turn that into pm1=FALSE and
# write it to the resumable per-gene TSV, so a network outage was cached as
# "DOLPHIN says no PM1" and never retried. Now a failed fetch is written as
# pm1=NA and NA rows are redone (and replaced) on the next run.
#
# Plain R, no testthat (analyses/tests/test_dolphin_bulk.R needs testthat and
# the network). Run from project root:
#   Rscript analyses/tests/test_dolphin_bulk_failure.R

source("analyses/lib/dolphin_bulk.R")

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok" else "FAIL", msg))
  if (!isTRUE(ok)) fails <<- fails + 1L
}

tmp <- tempfile("dolphin_fail_")
dir.create(tmp, recursive = TRUE)
real_fetch <- get("fetch_dolphin", envir = globalenv())
on.exit({
  unlink(tmp, recursive = TRUE)
  assign("fetch_dolphin", real_fetch, envir = globalenv())
}, add = TRUE)
variants <- c("p.A53T", "p.E46K")
enst <- "ENST00000394986"

# Run 1: API down -> NULL from fetch_dolphin for every variant.
assign("fetch_dolphin", function(...) NULL, envir = globalenv())
df1 <- suppressMessages(fetch_dolphin_gene("SNCA", variants, ensembl = enst,
                                           cache_dir = tmp, rate_sleep_s = 0,
                                           verbose = FALSE))
check(nrow(df1) == 2, "run 1: one row per variant")
check(all(is.na(df1$pm1)), "run 1: failed fetch -> pm1 is NA, not FALSE")
on_disk <- read.table(file.path(tmp, "SNCA.tsv"), sep = "\t", header = TRUE,
                      stringsAsFactors = FALSE, quote = "", comment.char = "")
check(all(is.na(on_disk$pm1)), "run 1: NA survives the TSV round-trip")

# Run 2: API back; A53T is in-domain with PM1, E46K is outside any domain.
calls <- character(0)
assign("fetch_dolphin", function(ensembl, p_notation, ...) {
  calls <<- c(calls, p_notation)
  if (p_notation == "p.A53T") list(results = list(list(acmg = "PM1", pfid = "PF01387")))
  else list(results = "mutation not within domain")
}, envir = globalenv())
df2 <- suppressMessages(fetch_dolphin_gene("SNCA", variants, ensembl = enst,
                                           cache_dir = tmp, rate_sleep_s = 0,
                                           verbose = FALSE))
check(setequal(calls, variants), "run 2: both NA rows are retried")
check(nrow(df2) == 2, "run 2: retried rows replace the NA rows, no duplicates")
check(isTRUE(df2$pm1[df2$p_notation == "p.A53T"]), "run 2: in-domain PM1 -> TRUE")
check(identical(df2$pm1[df2$p_notation == "p.E46K"], FALSE), "run 2: outside domain -> FALSE")

# Run 3: fully resolved cache -> no network calls.
calls <- character(0)
df3 <- suppressMessages(fetch_dolphin_gene("SNCA", variants, ensembl = enst,
                                           cache_dir = tmp, verbose = FALSE))
check(length(calls) == 0, "run 3: resolved rows are not re-fetched")
check(nrow(df3) == 2 && !anyNA(df3$pm1), "run 3: cache is complete and NA-free")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
