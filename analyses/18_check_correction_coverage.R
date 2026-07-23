#!/usr/bin/env Rscript
# 18_check_correction_coverage.R — is each MANE accession correction safe to
# apply without losing an upstream data layer?
#
# Switching gene_data's accession changes which AlphaFold / AlphaMissense files
# VarViz fetches (both are served per UniProt accession on its CANONICAL isoform,
# AF-<acc>-F1-*). A correction that points at the biologically right protein but
# one that AlphaFold/AlphaMissense never modelled would trade a real coverage
# gain for a coverage loss. This checks AF pLDDT + AlphaMissense presence for the
# OLD and NEW accession of every correction and emits a safe-to-apply allowlist.
#
# Coordinate note: this does NOT switch VarViz to MANE numbering. Each accession
# is used on its own canonical, which is what all the UniProt-anchored layers
# already agree on. MANE only told us which accession is the right one.
#
# Usage:  Rscript analyses/18_check_correction_coverage.R
suppressMessages(library(httr2))

CORR <- Sys.glob("VarViz_MANE_accession_corrections_*.tsv")
CORR <- CORR[order(CORR, decreasing = TRUE)][1]
stopifnot(!is.na(CORR), file.exists(CORR))
d <- read.delim(CORR, stringsAsFactors = FALSE)

# HEAD a URL; TRUE only on 200. AlphaFold confidence files are versioned
# (v6 current, older entries lag), so accept any version. AlphaMissense is either
# aa-substitutions.csv or the -v1 variant.
head_ok <- function(url) {
  r <- try(request(url) |> req_method("HEAD") |> req_timeout(20) |>
             req_error(is_error = function(x) FALSE) |> req_perform(), silent = TRUE)
  !inherits(r, "try-error") && resp_status(r) == 200
}
af_ok <- function(acc) {
  if (is.na(acc) || !nzchar(acc)) return(FALSE)
  for (v in c("v6","v5","v4","v3","v2"))
    if (head_ok(sprintf("https://alphafold.ebi.ac.uk/files/AF-%s-F1-confidence_%s.json", acc, v)))
      return(TRUE)
  FALSE
}
am_ok <- function(acc) {
  if (is.na(acc) || !nzchar(acc)) return(FALSE)
  for (f in c("aa-substitutions.csv","aa-substitutions-v1.csv"))
    if (head_ok(sprintf("https://alphafold.ebi.ac.uk/files/AF-%s-F1-%s", acc, f)))
      return(TRUE)
  FALSE
}

res <- data.frame()
for (i in seq_len(nrow(d))) {
  old <- d$uniprot_id[i]; new <- d$mane_uniprot_id[i]
  oAF <- af_ok(old); oAM <- am_ok(old); nAF <- af_ok(new); nAM <- am_ok(new)
  aligns <- isTRUE(as.logical(d$numbering_ok[i]))   # new canonical == MANE protein
  # Verdict:
  #  APPLY  — new has AF and (AM present, or AM absent on BOTH so nothing is lost)
  #  HOLD   — new drops AlphaMissense that the old accession served, or lacks AF
  verdict <- if (!nAF) "HOLD" else
             if (!nAM && oAM) "HOLD" else "APPLY"
  reason <- if (!nAF) "new accession has no AlphaFold model" else
            if (!nAM && oAM) "new accession drops AlphaMissense the old one served" else
            if (!nAM && !oAM) "no AlphaMissense either way; AF + everything else gained" else
            if (aligns) "full coverage; canonical == MANE (numbering aligns)" else
            "full coverage; canonical != MANE length (still the correct protein)"
  res <- rbind(res, data.frame(
    gene = d$gene_name[i], old_acc = old, new_acc = new,
    old_len = d$uniprot_len[i], new_len = d$mane_uniprot_len[i],
    mane_prot_len = d$mane_prot_len[i], numbering_aligns = aligns,
    old_AF = oAF, old_AM = oAM, new_AF = nAF, new_AM = nAM,
    verdict = verdict, reason = reason, stringsAsFactors = FALSE))
}
res <- res[order(res$verdict, res$gene), ]

OUT <- paste0("VarViz_MANE_correction_coverage_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write.table(res, OUT, sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# Machine-readable allowlist of the accessions that are safe to swap.
allow <- res[res$verdict == "APPLY", c("gene", "old_acc", "new_acc")]
AOUT <- paste0("VarViz_MANE_correction_allowlist_", format(Sys.Date(), "%Y%m%d"), ".tsv")
write.table(allow, AOUT, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("\n%d corrections: %d APPLY, %d HOLD\n",
            nrow(res), sum(res$verdict == "APPLY"), sum(res$verdict == "HOLD")))
print(res[, c("gene","old_acc","new_acc","new_AF","new_AM","verdict","reason")], row.names = FALSE)
cat(sprintf("\nwritten: %s (full), %s (allowlist)\n", OUT, AOUT))
