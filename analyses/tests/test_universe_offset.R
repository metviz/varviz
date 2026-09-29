# Guards the numbering correction applied to MaveDB scoresets in the universe build.
#
# Why this exists. analyses/check_universe_numbering.py found 6,186 of 87,974
# universe variants (7.03%) not numbered on their gene's reference sequence, all
# from MaveDB, and two scoresets account for every one of them:
#
#   TP53  00000059-a-1  7.2% on reference numbering; +101 fixes 2,754 of 2,754.
#                       Its positions span 1-191, so +101 maps them to 102-292,
#                       a domain construct numbered from its own residue 1.
#   BRCA1 00000003-a-2  7.6% on reference numbering; +1 fixes 3,420 of 3,429.
#
# Left uncorrected, those variants were CLASSIFIED at the wrong residues, which is
# why both scoresets had to be dropped from Panel B entirely rather than merely
# re-oriented. The offsets are recorded in
# analyses/derived/mavedb_scoreset_offsets.tsv.
#
# The correction belongs in the universe build, not in the MaveDB pull: the raw
# canonical table should keep what MaveDB published, and the derived universe should
# be on reference numbering. p_notation_raw preserves the construct numbering so a
# variant can still be traced back to its scoreset.
#
# Run from project root: Rscript analyses/tests/test_universe_offset.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
suppressMessages(source(resolve("analyses/lib/dms_polarity.R"), local = FALSE))

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

cat("-- the helper exists --\n")
check(exists("shift_p_notation"), "shift_p_notation is available")

cat("-- a position shifts while the residues are untouched --\n")
check(identical(shift_p_notation("p.A101H", 1L), "p.A102H"),
      "p.A101H at +1 becomes p.A102H")
check(identical(shift_p_notation("p.R74H", 101L), "p.R175H"),
      "TP53 construct p.R74H at +101 becomes p.R175H, the known hotspot")
check(identical(shift_p_notation("p.G144S", 101L), "p.G245S"),
      "construct p.G144S at +101 becomes p.G245S")
check(identical(shift_p_notation("p.M1V", -0L), "p.M1V"), "an offset of zero changes nothing")
check(identical(shift_p_notation("p.R175H", -101L), "p.R74H"), "a negative offset shifts down")

cat("-- stop codons and vectorisation --\n")
check(identical(shift_p_notation("p.R213*", 1L), "p.R214*"),
      "a stop-gain alt survives the shift")
check(identical(shift_p_notation(c("p.A1G", "p.A2G", "p.A3G"), 10L),
                c("p.A11G", "p.A12G", "p.A13G")),
      "a scalar offset recycles over a vector of notations")
check(identical(shift_p_notation(c("p.A1G", "p.A1G"), c(0L, 5L)),
                c("p.A1G", "p.A6G")),
      "a per-row offset vector is applied element-wise")

cat("-- missing and malformed input fail closed --\n")
# Passing an unparseable value through unshifted would silently reintroduce the
# defect for exactly the rows least likely to be noticed.
check(is.na(shift_p_notation(NA_character_, 1L)), "NA in gives NA out")
check(is.na(shift_p_notation("p.Arg175His", 1L)),
      "a three-letter notation is refused rather than passed through unshifted")
check(is.na(shift_p_notation("not a variant", 1L)), "junk gives NA")
check(is.na(shift_p_notation("p.A101H", NA_integer_)), "an NA offset gives NA")

cat("-- the shift never produces a position below 1 --\n")
# A construct offset larger than the position means the recorded offset is wrong for
# that row; producing p.A-3G would corrupt a join key silently.
check(is.na(shift_p_notation("p.A2G", -5L)),
      "an offset that would push the position below 1 gives NA")

cat("-- the recorded offsets are the two known ones --\n")
off <- load_dms_offsets(resolve("analyses/derived/mavedb_scoreset_offsets.tsv"))
if (length(off)) {
  check(identical(dms_offset_for(off, "00000059-a-1"), 101L), "TP53 00000059-a-1 is +101")
  check(identical(dms_offset_for(off, "00000003-a-2"), 1L),   "BRCA1 00000003-a-2 is +1")
  check(identical(dms_offset_for(off, "00000096-a-1"), 0L),   "a clean scoreset is 0")
  nz <- names(off)[off != 0L]
  check(setequal(nz, c("00000059-a-1", "00000003-a-2")),
        sprintf("exactly two scoresets carry an offset (found %d)", length(nz)))
}

cat("-- corrected positions land on the reference residue --\n")
# The check that matters: after shifting, the reference amino acid in the notation
# must equal the residue at that position in the gene's MANE protein.
mp <- resolve("analyses/derived/mane_residue_maps.tsv")
if (file.exists(mp)) {
  m <- utils::read.delim(mp, stringsAsFactors = FALSE)
  ref_at <- function(g, p) {
    r <- m[m$gene == g & m$pos == p, "aa"]
    if (length(r)) r[1] else NA_character_
  }
  for (spec in list(c("TP53", "p.R74H", "101"), c("TP53", "p.G144S", "101"),
                    c("TP53", "p.R147W", "101"))) {
    g <- spec[1]; before <- spec[2]; k <- as.integer(spec[3])
    after <- shift_p_notation(before, k)
    pos <- as.integer(sub("^p\\.[A-Z*]([0-9]+)[A-Z*]$", "\\1", after))
    r   <- sub("^p\\.([A-Z*])[0-9]+[A-Z*]$", "\\1", after)
    check(identical(ref_at(g, pos), r),
          sprintf("%s %s -> %s sits on reference residue %s", g, before, after, r))
  }
} else {
  cat("  [skip] mane_residue_maps.tsv absent; build with analyses/build_mane_residue_maps.R\n")
}

cat("-- the universe build applies the offsets and keeps the construct numbering --\n")
ub <- resolve("analyses/03_build_universe.R")
if (file.exists(ub)) {
  code <- sub("#.*$", "", readLines(ub, warn = FALSE))
  check(any(grepl("dms_polarity\\.R", code)), "03_build_universe.R sources the offsets loader")
  check(any(grepl("shift_p_notation\\(", code)), "03_build_universe.R shifts MaveDB positions")
  check(any(grepl("dms_offset_for\\(", code)), "03_build_universe.R reads the per-scoreset offset")
  # Traceability: without the original notation a corrected row cannot be matched
  # back to the scoreset it came from.
  check(any(grepl("p_notation_raw", code)), "03_build_universe.R preserves p_notation_raw")
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
