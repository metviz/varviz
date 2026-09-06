# Offline self-check for analyses/lib/dolphin.R cache behaviour on failure.
#
# fetch_dolphin() and dolphin_canonical_enst() used to write their failure
# value (NA) into the session cache exactly like a real response, so one
# timeout made every later lookup of that key in the same harness run return
# "no answer" without touching the network. A failure must not be cached;
# only real answers are.
#
# Deterministic without a stubbed network: a sub-millisecond timeout makes
# httr::GET fail before any handshake can complete. Plain R, no testthat.
# Run from project root: Rscript analyses/tests/test_dolphin_failure_cache.R

source("analyses/lib/dolphin.R")   # standalone -> .dolphin_local_cache is used

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok" else "FAIL", msg))
  if (!isTRUE(ok)) fails <<- fails + 1L
}
cache_keys <- function() ls(.dolphin_local_cache, all.names = TRUE)

# 1. failed variant fetch -> NULL, nothing cached
r <- suppressMessages(fetch_dolphin(ensembl = "ENST00000394986", p_notation = "p.A53T",
                                    timeout_s = 0.001))
check(is.null(r), "fetch_dolphin returns NULL on network failure")
check(length(cache_keys()) == 0, "failed variant fetch is NOT cached")

# 2. failed ENST lookup -> NA_character_, nothing cached
e <- suppressMessages(dolphin_canonical_enst("SNCA", timeout_s = 0.001))
check(identical(e, NA_character_), "dolphin_canonical_enst returns NA on failure")
check(length(cache_keys()) == 0, "failed ENST lookup is NOT cached")

# 3. real answers still cache and short-circuit the network
.dolphin_cache_set("ENST00000394986:A53T", list(results = list(list(acmg = "PM1"))))
.dolphin_cache_set("ENST:SNCA", "ENST00000394986")
r2 <- fetch_dolphin(ensembl = "ENST00000394986", p_notation = "p.A53T", timeout_s = 0.001)
check(isTRUE(dolphin_fires_pm1(r2)), "cached success is served without the network")
e2 <- dolphin_canonical_enst("SNCA", timeout_s = 0.001)
check(identical(e2, "ENST00000394986"), "cached ENST is served without the network")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
