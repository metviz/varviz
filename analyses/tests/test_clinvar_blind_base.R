#!/usr/bin/env Rscript
# Base-R checks for strip_clinvar_tags(). No testthat dependency, so this runs
# anywhere. Covers the pre-existing contract plus the composite-pathway and
# hotspot-upgrade demotion added with the ungated MDS pathway.
#
# Run from the repository root:  Rscript analyses/tests/test_clinvar_blind_base.R

source("analyses/lib/clinvar_blind.R")

ok <- TRUE
chk <- function(label, got, want) {
  pass <- setequal(got, want)
  if (!pass) ok <<- FALSE
  cat(sprintf("  [%s] %s\n", if (pass) "PASS" else "FAIL", label))
  if (!pass) cat(sprintf("        got:  %s\n        want: %s\n",
                         paste(sort(got), collapse = ","),
                         paste(sort(want), collapse = ",")))
}

cat("strip_clinvar_tags()\n")

chk("direct ClinVar tags removed",
    strip_clinvar_tags(c("PM2", "PP2", "PS1", "PM5", "PP5", "BP6")),
    c("PM2", "PP2"))

chk("CCRS-derived PM1 survives blinding",
    strip_clinvar_tags(c("PM1", "PM2"), "ccrs"), c("PM1", "PM2"))

chk("UniProt-domain PM1 survives blinding",
    strip_clinvar_tags(c("PM1", "PM2"), "uniprot_domain"), c("PM1", "PM2"))

chk("hotspot-only PM1 is dropped",
    strip_clinvar_tags(c("PM1", "PM2"), "clinvar_hotspot"), "PM2")

chk("hotspot-only PM1_strong is dropped too",
    strip_clinvar_tags(c("PM1_strong", "PM2"), "clinvar_hotspot"), "PM2")

# --- composite pathways introduced by the ungated MDS route ------------------
chk("MDS-corroborated CCRS Strong survives intact",
    strip_clinvar_tags(c("PM1_strong", "PM2"), "ccrs+mds"),
    c("PM1_strong", "PM2"))

chk("MDS-originated PM1 survives",
    strip_clinvar_tags(c("PM1_moderate_plus", "PM2"), "mds"),
    c("PM1_moderate_plus", "PM2"))

chk("hotspot upgrade demoted to its pre-upgrade strength",
    strip_clinvar_tags(c("PM1_strong", "PM2"), "ccrs+hotspot_upgrade(PM1)"),
    c("PM1", "PM2"))

chk("hotspot upgrade over a tiered MDS base demotes to that tier",
    strip_clinvar_tags(c("PM1_strong", "PM2"),
                       "ccrs+mds+hotspot_upgrade(PM1_moderate_plus)"),
    c("PM1_moderate_plus", "PM2"))

chk("no pathway supplied leaves tags alone",
    strip_clinvar_tags(c("PM1", "PM2")), c("PM1", "PM2"))

chk("NA pathway is tolerated",
    strip_clinvar_tags(c("PM1", "PM2"), NA_character_), c("PM1", "PM2"))

cat(if (ok) "\nall checks pass\n" else "\nFAILURES ABOVE\n")
quit(status = if (ok) 0 else 1)
