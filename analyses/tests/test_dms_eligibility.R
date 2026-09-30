# Guards which scoresets are eligible to serve as pathogenicity truth at all.
#
# Why this exists. Polarity and numbering ask whether a scoreset is read correctly.
# Eligibility asks a prior question: does its axis measure pathogenicity?
#
# The six KRAS scoresets of urn:mavedb:00000115-a report binding free energy changes
# against six partners (RAF1, SOS1, PIK3CG, RALGDS and the DARPins K27 and K55), from
# PMID 38109937. A larger ddG is weaker binding. That is a real measurement, but it is
# not a pathogenicity axis: KRAS pathogenic variants are gain-of-function acting through
# impaired GTP hydrolysis, and several of them retain or increase effector binding.
#
# Measured on the corrected universe two independent ways:
#
#   Established oncogenic variants (G12D, G12V, G12C, G12A, G13D, Q61H, Q61L) sit at
#   percentiles 6-90 across the six scoresets, median 37, with only 3 of 42 anchor
#   observations in either extreme decile. A pathogenicity axis concentrates them in one
#   tail.
#
#   Using VariBench clinical labels, which are independent of the DMS score and of every
#   predictor, labelled-Pathogenic variants land at median percentile 42-53 in the six
#   KRAS scoresets. Their median distance from random placement is 3, against 26 for the
#   twelve other scoresets, and they occupy the bottom six of eighteen.
#
# They also dominate by weight: six of twenty scoresets from a single experiment, which
# is 30% of an unweighted median.
#
# This is a scope decision, not a defect, so it is recorded per scoreset with its reason
# rather than hardcoded in a filter.
#
# GCK 00000096-a-1 scores weakly on the clinical measure (distance 11) and is NOT
# excluded. The reason first written here, a weak labelled subset, was wrong. The measure
# is invalid for GCK: the gene has both loss-of-function disease (GCK-MODY) and
# gain-of-function disease (hyperinsulinemic hypoglycemia), so its pathogenic variants
# belong at opposite tails. PMID 37101203 names W99R as an HH variant that grew faster
# than wild type, and the data places it at the 99th percentile. The axis is sound, and
# distance-from-random is therefore supporting evidence for an exclusion rather than the
# test for one.
#
# Run from project root: Rscript analyses/tests/test_dms_eligibility.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
suppressMessages(source(resolve("analyses/lib/dms_polarity.R"), local = FALSE))

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

cat("-- the helpers exist --\n")
check(exists("load_dms_eligibility"),    "load_dms_eligibility is available")
check(exists("dms_ineligible_reason"),   "dms_ineligible_reason is available")

cat("-- the recorded exclusions are exactly the six KRAS ddG scoresets --\n")
POL <- resolve("analyses/derived/mavedb_scoreset_polarity.tsv")
el <- load_dms_eligibility(POL)
KRAS6 <- paste0("00000115-a-", 1:6)
# A present table that lacks the column must fail rather than skip: an empty map reads
# as "nothing excluded", which would pass every assertion below for the wrong reason.
if (file.exists(POL)) {
  check(length(el) > 0L,
        "the polarity table carries an ineligible_reason column (regenerate with analyses/determine_mavedb_polarity.py)")
}
if (length(el)) {
  excluded <- names(el)[nzchar(el)]
  check(setequal(excluded, KRAS6),
        sprintf("exactly the six KRAS scoresets are excluded (found %d: %s)",
                length(excluded), paste(excluded, collapse = ", ")))
  for (s in KRAS6) {
    r <- dms_ineligible_reason(el, s)
    check(nzchar(r), sprintf("%s carries an exclusion reason", s))
    check(grepl("binding|ddG|free energy", r, ignore.case = TRUE),
          sprintf("%s names the measured quantity in its reason", s))
  }
  # The ones that must NOT be swept in.
  for (s in c("00000059-a-1", "00000093-a-1", "00000096-a-1", "00001216-a-1")) {
    check(!nzchar(dms_ineligible_reason(el, s)),
          sprintf("%s stays eligible", s))
  }
}

cat("-- an unlisted scoreset is eligible, so a new one is not silently dropped --\n")
check(!nzchar(dms_ineligible_reason(el, "99999999-z-9")),
      "an unknown scoreset is eligible by default")
check(!nzchar(dms_ineligible_reason(setNames(character(0), character(0)), "00000115-a-1")),
      "with no table loaded nothing is excluded")

cat("-- a missing table degrades rather than erroring --\n")
e0 <- load_dms_eligibility(file.path(tempdir(), "definitely_absent.tsv"))
check(length(e0) == 0L, "an absent table yields an empty map")

cat("-- the reason travels with the exclusion, not only the fact of it --\n")
tmp <- file.path(tempdir(), paste0("elig_", as.integer(runif(1, 1e6, 9e6)), ".tsv"))
on.exit(unlink(tmp), add = TRUE)
writeLines(c("gene\tstudy\tineligible_reason",
             "KRAS\t00000115-a-1\tmeasures binding free energy, not pathogenicity",
             "GCK\t00000096-a-1\t"), tmp)
e <- load_dms_eligibility(tmp)
check(identical(dms_ineligible_reason(e, "00000115-a-1"),
                "measures binding free energy, not pathogenicity"),
      "the reason round-trips verbatim")
check(!nzchar(dms_ineligible_reason(e, "00000096-a-1")),
      "an empty reason means eligible")

cat("-- every scoreset consumer honours eligibility --\n")
for (f in c("analyses/07_panel_b_dms.R", "analyses/panel_b_dms_reval.R",
            "analyses/panel_b_auroc_summary.R", "analyses/10_compose_supp_s2.R",
            "analyses/repro/07_panel_b_dms_run.R")) {
  p <- resolve(f)
  if (!file.exists(p)) { cat(sprintf("  [skip] %s absent\n", basename(f))); next }
  code <- sub("#.*$", "", readLines(p, warn = FALSE))
  check(any(grepl("dms_ineligible_reason\\(|load_dms_eligibility\\(", code)),
        sprintf("%s honours recorded eligibility", basename(f)))
}

cat("-- the curation records the reason and its evidence --\n")
d <- resolve("analyses/determine_mavedb_polarity.py")
if (file.exists(d)) {
  src <- paste(readLines(d, warn = FALSE), collapse = "\n")
  check(grepl("ineligible_reason", src), "the determination script emits ineligible_reason")
  check(grepl("38109937", src), "the KRAS publication is cited beside the exclusion")
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
