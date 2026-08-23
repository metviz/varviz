# strip_clinvar_tags() — withhold ClinVar-derived ACMG tags from a tags vector
# before passing it to classify_acmg(). Used by the dual-pass classification
# harness (analyses/05_classify_harness.R) to produce the Pass-Blind label.
#
# Reference: docs/plans/2026-04-26-sk-revisions-pass2.md "Methodological
# Foundations" §MF1. The deeper circularity in benchmarking VarViz against
# ClinVar truth is that PS1 / PM5 / PM1-hotspot / PP5 / BP6 are all derived
# from ClinVar lookups — so the classifier reads ClinVar as both input AND
# truth label. Pass-Blind defuses this by withholding those tags.
#
# PP2/BP1 (ClinGen+GenCC validity) are gene-level priors, NOT per-variant
# ClinVar lookups, and are intentionally NOT stripped: removing them would
# defeat the purpose of having any classifier rather than just isolate the
# circularity.

CLINVAR_DIRECT_TAGS <- c(
  "PS1", "PS1_moderate", "PS1_supporting",
  "PM5", "PM5_supporting",
  "PP5", "PP5_moderate", "PP5_supporting",
  "BP6", "BP6_moderate", "BP6_supporting"
)

PM1_TAGS_ALL <- c("PM1", "PM1_moderate_plus", "PM1_strong")

strip_clinvar_tags <- function(tags_vec, pm1_pathway = character(0)) {
  out <- setdiff(tags_vec, CLINVAR_DIRECT_TAGS)
  if (length(pm1_pathway) == 0L) return(out)
  pw <- pm1_pathway[1]
  if (is.na(pw)) return(out)

  # PM1 fires from several independent pathways in server.R: UniProt functional
  # site, CCRS percentile, UniProt domain, MDS (Pfam alignment), and the ClinVar
  # 15-residue hotspot. Only the hotspot is ClinVar-circular; the rest stay.
  # server.R records the firing pathway, joining co-firing routes with "+"
  # (e.g. "ccrs+mds"), so match on substrings rather than equality.

  # Case 1 — PM1 originated from the hotspot alone: drop it entirely.
  if (identical(pw, "clinvar_hotspot")) {
    out <- setdiff(out, PM1_TAGS_ALL)
    return(out)
  }

  # Case 2 — a non-circular pathway fired and the hotspot then upgraded it to
  # PM1_strong. The base strength is legitimate under blinding, the upgrade is
  # not, so restore the pre-upgrade tag that server.R recorded in the pathway
  # string. Leaving PM1_strong in place would carry ClinVar-derived strength
  # into the arm defined as ClinVar-blind.
  m <- regmatches(pw, regexec("\\+hotspot_upgrade\\(([^)]+)\\)", pw))[[1]]
  if (length(m) == 2L && "PM1_strong" %in% out) {
    out <- setdiff(out, PM1_TAGS_ALL)
    out <- c(out, m[2])
  }
  out
}
