#!/usr/bin/env Rscript
# 14_validate_pfam_pssm.R — do our offline Pfam PSSM deltas reproduce DOLPHIN?
#
# Ground truth is the five score_delta values DOLPHIN returned for the variants
# cited in the manuscript, captured from its API on 2026-07-22 (the service was
# serving correctly; only its TLS certificate had lapsed, see
# analyses/lib/dolphin.R). If this script reproduces them from alignments we
# fetch ourselves, the PM1 domain evidence no longer depends on that API.
#
# Usage:  Rscript analyses/14_validate_pfam_pssm.R
# Alignments are cached under analyses/tmp/pfam/ after the first run.

source("analyses/lib/pfam_pssm.R")

CACHE <- "analyses/tmp/pfam"
dir.create(CACHE, showWarnings = FALSE, recursive = TRUE)

# DOLPHIN ground truth: pfam, human row prefix, residue, wt, mut, delta
TRUTH <- list(
  list(pf = "PF01387", row = "SYUA_HUMAN", label = "SNCA G14R",  pos = 14,   wt = "G", mut = "R", dolphin = -9.161),
  list(pf = "PF01387", row = "SYUA_HUMAN", label = "SNCA E46K",  pos = 46,   wt = "E", mut = "K", dolphin = -8.873),
  list(pf = "PF01387", row = "SYUA_HUMAN", label = "SNCA G51D",  pos = 51,   wt = "G", mut = "D", dolphin = -5.415),
  list(pf = "PF01387", row = "SYUA_HUMAN", label = "SNCA A53T",  pos = 53,   wt = "A", mut = "T", dolphin =  1.165),
  list(pf = "PF01094", row = "CASR_HUMAN", label = "CASR G143E", pos = 143,  wt = "G", mut = "E", dolphin = -5.959)
)

fetch_alignment <- function(pf) {
  dest <- file.path(CACHE, paste0(pf, "_full.sto.gz"))
  if (!file.exists(dest)) {
    url <- sprintf(
      "https://www.ebi.ac.uk/interpro/wwwapi/entry/pfam/%s/?annotation=alignment:full", pf)
    message("  downloading ", pf, " full alignment ...")
    utils::download.file(url, dest, mode = "wb", quiet = TRUE)
  }
  dest
}

pssm_cache <- new.env(parent = emptyenv())
get_pssm <- function(pf) {
  if (!is.null(pssm_cache[[pf]])) return(pssm_cache[[pf]])
  aln <- read_stockholm(fetch_alignment(pf))
  message(sprintf("  %s: %s sequences, %s columns", pf,
                  format(length(aln), big.mark = ","), nchar(aln[[1]])))
  obj <- list(aln = aln, M = build_pssm(aln))
  pssm_cache[[pf]] <- obj
  obj
}

message("[pfam-pssm] validating against DOLPHIN ground truth\n")
ours <- numeric(0); theirs <- numeric(0); labs <- character(0)
cat(sprintf("%-12s %-9s %-9s %-10s %-9s %s\n",
            "variant", "pfam", "column", "ours", "DOLPHIN", "diff"))
for (t in TRUTH) {
  p   <- get_pssm(t$pf)
  hid <- grep(t$row, names(p$aln), value = TRUE)[1]
  col <- alignment_column_for_residue(p$aln, hid, t$pos)
  d   <- pssm_delta(p$M, col, t$wt, t$mut)
  ours <- c(ours, d); theirs <- c(theirs, t$dolphin); labs <- c(labs, t$label)
  cat(sprintf("%-12s %-9s %-9s %-10s %-9s %s\n",
              t$label, t$pf, col, round(d, 3), t$dolphin, round(d - t$dolphin, 3)))
}

ok_sign <- sum(sign(ours) == sign(theirs))
cat(sprintf("\n  sign agreement : %d/%d\n", ok_sign, length(ours)))
cat(sprintf("  Pearson  r     : %.3f\n", cor(ours, theirs)))
cat(sprintf("  Spearman rho   : %.3f\n", cor(ours, theirs, method = "spearman")))
cat(sprintf("  mean |diff|    : %.3f\n", mean(abs(ours - theirs))))
cat(sprintf("  max  |diff|    : %.3f\n", max(abs(ours - theirs))))

# A constant offset is expected and benign: DOLPHIN built its matrices from
# Pfam 33.1 (May 2020) filtered to eukaryotes, whereas we use the current full
# alignment unfiltered. What matters for a PM1 call is the ordering and the
# sign, not the absolute value.
cat(sprintf("\n  mean signed offset (ours - DOLPHIN): %+.3f\n", mean(ours - theirs)))
cat(sprintf("  sd of offset                        : %.3f\n", sd(ours - theirs)))
