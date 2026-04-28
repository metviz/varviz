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
  "PM5",
  "PP5", "PP5_moderate", "PP5_supporting",
  "BP6", "BP6_moderate", "BP6_supporting"
)

strip_clinvar_tags <- function(tags_vec, pm1_pathway = character(0)) {
  out <- setdiff(tags_vec, CLINVAR_DIRECT_TAGS)
  # PM1 fires from three independent pathways in server.R: ClinVar 15-residue
  # hotspot, CCRS percentile, UniProt domain. Only the hotspot path is
  # ClinVar-circular; CCRS and UniProt domain pathways stay.
  if ("PM1" %in% out && length(pm1_pathway) > 0L &&
      "clinvar_hotspot" %in% pm1_pathway) {
    out <- setdiff(out, "PM1")
  }
  out
}
