# Self-check for the harness's degraded-predictor-cache guard.
# Run from project root: Rscript analyses/tests/test_harness_cache_guard.R
#
# Why this exists. The guard decided whether a run was degraded by asking
# whether the bulk per-chromosome preload object was NULL. That is not the same
# question as whether the predictor values reached the result. REVEL and
# AlphaMissense also arrive through dbNSFP, so a gene can have a NULL preload
# and a complete set of values.
#
# On 2026-10-11 that gap made the guard refuse two reproducibility runs whose
# data was complete: 4 of 16 RASopathy genes and most of the 44 external ones
# were flagged, and in every one REVEL and AlphaMissense were fully populated
# and identical to the original run. Only eleven REVEL slices exist under
# analyses/raw/revel/by_chrom, so any cohort reaching another chromosome trips
# the old check.
#
# A guard that fires on complete data teaches people to pass
# VARVIZ_ALLOW_DEGRADED=1 by reflex, and that flag cannot tell "preload absent,
# values complete" from "values genuinely missing" -- the case the guard exists
# for, where SNCA scored PP3 Supporting instead of Moderate because the
# meta-predictors really were absent.
#
# So: the preload state stays recorded as a note, and the GATE reads the result
# table. Offline; no network, no fixtures.

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
source(resolve("analyses/lib/harness_guard.R"))

fails <- 0L; checks <- 0L
check <- function(label, ok) {
  ok <- isTRUE(ok); checks <<- checks + 1L
  if (!ok) fails <<- fails + 1L
  cat(if (ok) "  ok   " else "  FAIL ", label, "\n", sep = "")
}
reset <- function() options(varviz.localpred_degraded = character(0),
                            varviz.localpred_preload_note = character(0))

# A result table in the real column shape: the harness reads these names off
# build_variant_table()'s output.
vt <- function(revel, am, n = 3)
  data.frame(Variant = paste0("p.A", seq_len(n), "G"),
             REVEL_Score = revel, AM_Score = am,
             stringsAsFactors = FALSE)

# ── The gate reads values, not preloads ───────────────────────────────────
reset()
check("complete values with a NULL preload do NOT gate",
      isFALSE(record_value_gaps("LZTR1", vt(c(0.91, 0.42, 0.77), c(0.99, 0.31, 0.88)))))
check("and nothing is recorded as a gap",
      length(cache_degradations()) == 0L)

reset()
check("a predictor absent from every row DOES gate",
      isTRUE(record_value_gaps("FAM20A", vt(c(NA, NA, NA), c(0.5, 0.6, 0.7)))))
check("the gap names the gene and the predictor",
      any(grepl("FAM20A", cache_degradations())) &&
      any(grepl("REVEL",  cache_degradations())))
check("a predictor that did arrive is not named",
      !any(grepl("AlphaMissense", cache_degradations())))

reset()
check("empty strings count as missing, not as values",
      isTRUE(record_value_gaps("X", vt(c("", "", ""), c("0.5", "0.6", "0.7")))))
reset()
check("the literal 'NA' string counts as missing",
      isTRUE(record_value_gaps("X", vt(c("NA", "NA", "NA"), c("0.5", "0.6", "0.7")))))

reset()
check("a PARTIAL gap does not gate the run",
      isFALSE(record_value_gaps("Y", vt(c(0.9, NA, 0.7), c(0.9, 0.8, 0.7)))))
reset()
check("both predictors missing gates and names both",
      isTRUE(record_value_gaps("Z", vt(c(NA, NA, NA), c(NA, NA, NA)))) &&
      any(grepl("REVEL", cache_degradations())) &&
      any(grepl("AlphaMissense", cache_degradations())))

reset()
check("a missing column is treated as missing values",
      isTRUE(record_value_gaps("W", data.frame(Variant = "p.A1G", AM_Score = 0.9))))

# ── The preload check survives, demoted to a note ─────────────────────────
reset()
record_cache_degradation("LZTR1", new.env(), NULL, list(1), 100, 100)
check("a NULL preload records a note",
      length(preload_notes()) >= 1L && any(grepl("LZTR1", preload_notes())))
check("and that note is NOT a gate",
      length(cache_degradations()) == 0L)

# ── The real failure the guard exists for still trips it ──────────────────
# SNCA: REVEL (chr4 absent) and AlphaMissense never loaded, the values were
# genuinely absent, and the meta-predictors that earn PP3 Moderate were gone.
reset()
record_cache_degradation("SNCA", new.env(), NULL, NULL, 140, 140)
gated <- record_value_gaps("SNCA", vt(c(NA, NA, NA), c(NA, NA, NA)))
check("the SNCA case -- preload AND values absent -- still gates", isTRUE(gated))

# ── Wiring: the harness must gate on the value check ──────────────────────
h <- readLines(resolve("analyses/ps_final_harness.R"), warn = FALSE)
code <- h[!grepl("^\\s*#", h)]
check("the harness calls record_value_gaps()", any(grepl("record_value_gaps\\(", code)))
# min() of an empty match is Inf, which would pass this vacuously.
.rv <- grep("record_value_gaps\\(", code); .bv <- grep("build_variant_table\\(", code)
check("it is called after build_variant_table()",
      length(.rv) > 0 && length(.bv) > 0 && min(.rv) > min(.bv))
check("cache_degradations() still drives the refusal",
      any(grepl("cache_degradations\\(\\)", code)))
check("the preload note is surfaced too", any(grepl("preload_notes\\(\\)", code)))

# Every harness that gates on cache_degradations() must also populate it.
# Demoting record_cache_degradation() to a note silently disarms any harness
# that was relying on it, which is a worse outcome than the false positives
# this change set out to fix.
for (hp in c("analyses/ps_final_harness.R", "analyses/05_classify_harness.R")) {
  src <- readLines(resolve(hp), warn = FALSE)
  src <- src[!grepl("^\\s*#", src)]
  if (!any(grepl("cache_degradations\\(\\)", src))) next
  check(paste(basename(hp), "gates and therefore calls record_value_gaps()"),
        any(grepl("record_value_gaps\\(", src)))
}

cat(sprintf("\nharness cache guard: %d check(s), %d failure(s)\n", checks, fails))
quit(status = if (fails) 1L else 0L)
