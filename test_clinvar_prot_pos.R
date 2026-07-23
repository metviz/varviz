# Self-check for ClinVar protein-position parsing.
# Run: Rscript test_clinvar_prot_pos.R
#
# Regression guard for the isoform-collision bug: ClinVar's protein_change field
# lists the same change across every overlapping isoform ("D77G, D78G, D79G"),
# so parsing its first number filed p.Asp78Gly at residue 77 and mis-fired
# PS1/PM5 against a real codon-77 query. The canonical position lives in the
# title's (p.Xxx###Yyy). This test slices the real parse block out of server.R
# (a trusted in-repo build artifact, never user input) and drives it with the
# exact ClinVar strings, so it fails if the source ordering regresses.
src <- readLines("server.R", warn = FALSE)
beg <- grep("^         aa_source <- NULL$", src)
end <- grep("^         # Also check molecular_consequence_list for variant type$", src) - 1L
stopifnot(length(beg) == 1, length(end) == 1, end > beg)
block <- paste(src[beg:end], collapse = "\n")

# Reproduce the two enclosing statements the block relies on.
parse_pos <- function(var_title, prot_change) {
  e <- new.env(parent = globalenv())
  e$var_title   <- var_title
  e$prot_change <- prot_change
  e$rec         <- list()          # variation_set_name absent
  e$prot_pos <- NA; e$var_type <- "other"
  eval(parse(text = block), envir = e)
  # position extraction, verbatim from the lines just after the block
  if (!is.null(e$aa_source) && nchar(e$aa_source) > 0) {
    pm <- regmatches(e$aa_source, regexpr("\\d+", e$aa_source))
    if (length(pm) > 0) e$prot_pos <- as.numeric(pm)
  }
  e$prot_pos
}

# ── Real ClinVar records (fetched from eutils 2026-07-23) ─────────────────
# The user's query residue and the collider that used to land on it.
stopifnot(
  # p.Leu77Pro — canonical codon 77, protein_change "L76P, L77P, L78P"
  parse_pos("NM_000162.5(GCK):c.230T>C (p.Leu77Pro)", "L76P, L77P, L78P") == 77,
  # p.Asp78Gly — canonical codon 78, protein_change "D77G, D78G, D79G".
  # THE BUG: first-number-of-protein_change filed this at 77.
  parse_pos("NM_000162.5(GCK):c.233A>G (p.Asp78Gly)", "D77G, D78G, D79G") == 78,
  parse_pos("NM_000162.5(GCK):c.232G>A (p.Asp78Asn)", "D77N, D78N, D79N") == 78,
  # Larger position, same multi-isoform shape
  parse_pos("NM_000162.5(GCK):c.1133C>T (p.Ala378Val)",
            "A378V, A379V, A41V, A56V, A377V") == 378,
  # Frameshift: position still from the canonical title
  parse_pos("NM_000162.5(GCK):c.755_779dup (p.Phe260fs)", "F259fs, F260fs, F261fs") == 260
)

# ── Fallbacks still work when the title carries no p. change ───────────────
stopifnot(
  # No title p. → first token of protein_change, NOT the first number of the
  # whole comma list.
  parse_pos("NM_000162.5(GCK):c.1000A>G", "D77G, D78G, D79G") == 77,
  # No p. anywhere (CNV / splice) → NA, so the row is dropped upstream.
  is.na(parse_pos("GRCh38/hg38 7p14.1-13(chr7:38177999-45304100)x1", ""))
)

# ── Ordering guard: title parsed before protein_change ────────────────────
# The whole fix is that the canonical title is consulted first; if a later edit
# put protein_change back on top, the collision returns.
title_line <- grep("aa_source <- sub\\(\".\\*\\\\\\\\\\(p", src)[1]
pc_line    <- grep("aa_source <- trimws\\(strsplit\\(prot_change", src)[1]
stopifnot(!is.na(title_line), !is.na(pc_line), title_line < pc_line)

cat("clinvar prot_pos parsing: all checks pass\n")
