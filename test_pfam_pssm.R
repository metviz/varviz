# Self-check for the Pfam-alignment PSSM (Corcuff et al. 2023 construction).
# Run: Rscript test_pfam_pssm.R
# ponytail: synthetic alignments only — no network, no Pfam download — so this
# stays runnable anywhere. The real-data agreement with DOLPHIN is recorded in
# analyses/14_validate_pfam_pssm.R.
source("analyses/lib/pfam_pssm.R")

# ── Stockholm parsing ─────────────────────────────────────────────────────
tmp <- tempfile(fileext = ".sto")
writeLines(c(
  "# STOCKHOLM 1.0",
  "#=GF ID   Test",
  "SEQ_A/1-4    ACDE",
  "SEQ_B/1-4    ACDE",
  "#=GC RF      xxxx",
  "",
  "SEQ_A/1-4    FGHI",          # second block: Stockholm wraps, must concatenate
  "SEQ_B/1-4    FGHI",
  "//"), tmp)
aln <- read_stockholm(tmp)
stopifnot(
  length(aln) == 2L,
  aln[["SEQ_A/1-4"]] == "ACDEFGHI",   # blocks joined in order
  !any(grepl("^#", names(aln)))       # #=GC annotation rows excluded
)

# ── residue -> alignment column, honouring gaps and the id's start offset ──
tmp2 <- tempfile(fileext = ".sto")
writeLines(c("# STOCKHOLM 1.0",
             "P1/1-4       A-CD-E",
             "P2/10-13     AACD-E",
             "//"), tmp2)
a2 <- read_stockholm(tmp2)
stopifnot(
  alignment_column_for_residue(a2, "P1/1-4", 1) == 1L,   # A at col 1
  alignment_column_for_residue(a2, "P1/1-4", 2) == 3L,   # gap at col 2 skipped
  alignment_column_for_residue(a2, "P1/1-4", 3) == 4L,
  # start offset respected: residue 10 is the FIRST residue of this row
  alignment_column_for_residue(a2, "P2/10-13", 10) == 1L,
  alignment_column_for_residue(a2, "P2/10-13", 12) == 3L,
  is.na(alignment_column_for_residue(a2, "P1/1-4", 99)),   # past the end
  is.na(alignment_column_for_residue(a2, "NOPE/1-4", 1))   # unknown id
)

# ── PSSM construction ─────────────────────────────────────────────────────
# A column where every sequence carries the same residue must score that
# residue positively (observed >> background) and any substitution negatively.
tmp3 <- tempfile(fileext = ".sto")
writeLines(c("# STOCKHOLM 1.0",
             paste0("S", 1:40, "/1-2", "   ", "WA"),
             "//"), tmp3)
a3 <- read_stockholm(tmp3)
M3 <- build_pssm(a3)
stopifnot(
  nrow(M3) == 2L, ncol(M3) == 20L,
  identical(colnames(M3), AA20),
  M3[1, "W"] > 0,                       # W is enriched vs its 1.08% background
  M3[1, "A"] < 0,                       # A is absent at this column
  # delta direction: substituting away from the conserved residue is negative,
  # which is the sign convention DOLPHIN reports.
  pssm_delta(M3, 1, "W", "A") < 0,
  pssm_delta(M3, 1, "A", "W") > 0,
  # the identity substitution is exactly zero
  abs(pssm_delta(M3, 1, "W", "W")) < 1e-12
)

# ── pseudocount keeps unobserved residues finite ──────────────────────────
# Without the +c*f_l term an unobserved amino acid gives log(0) = -Inf and the
# delta becomes unusable. This is the whole reason the paper adds it.
stopifnot(
  is.finite(M3[1, "C"]),
  all(is.finite(M3)),
  # a rarely-seen residue still scores below a conserved one
  M3[1, "C"] < M3[1, "W"]
)

# ── guards ────────────────────────────────────────────────────────────────
stopifnot(
  is.na(pssm_delta(M3, NA_integer_, "W", "A")),
  is.na(pssm_delta(M3, 99L, "W", "A")),        # column out of range
  is.na(pssm_delta(M3, 1, "X", "A")),          # non-standard residue
  is.na(pssm_delta(M3, 1, "W", "*"))
)

cat("pfam PSSM: all checks pass\n")
