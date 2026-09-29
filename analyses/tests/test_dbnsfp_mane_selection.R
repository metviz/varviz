# Guards MANE-aware selection of dbNSFP's per-transcript predictor columns.
#
# Why this exists. dbNSFP 4.9a stores several predictors as ';'-separated
# per-transcript lists, and .dbnsfp_row_to_hit() read them with num_first(), the
# first non-'.' element. That is the MANE transcript's value only when MANE happens
# to be listed first, and analyses/reports/dbnsfp_mane_index_survey.tsv shows that
# is the minority for half the benchmark genes: GCK 0 of 88, TP53 0 of 91, DDR2
# 0 of 2, TSHR 0 of 4, KCNH2 3 of 55, BRCA1 16 of 110, KCNQ1 68 of 119. The MANE
# index also varies row to row (TP53 rows place it 3rd of 12, 4th of 16, 7th of 20,
# 10th of 21), so it has to be resolved per row from Ensembl_transcriptid.
#
# Measured cost of the defect over the 14-gene universe
# (analyses/dbnsfp_mane_am_audit.R): of 26,468 variants resolved in dbNSFP, 6,491
# (24.5%) get a different AlphaMissense value from the first element than from
# MANE's, and 1,285 (4.85%) land in a different PP3 tier, 996 of them a tier
# lower than MANE would give. Max absolute difference 0.6252.
#
# Which columns may be indexed was established empirically over 5,256
# multi-transcript rows across eight genes, comparing each column's list length
# against the Ensembl_transcriptid list length:
#
#   transcript-aligned, 0 mismatches : aapos, Ensembl_proteinid, HGVSp_VEP, SIFT,
#                                      Polyphen2 HDIV, Polyphen2 HVAR, FATHMM,
#                                      PROVEAN, MetaRNN, REVEL, AlphaMissense
#   single-valued, never indexed     : LRT, MetaSVM, MetaLR, CADD, DANN, PhyloP,
#                                      PhastCons, GERP, SiPhy, 1000Gp3_AF, ExAC_AF
#   MIXED, must NOT be indexed       : MutationTaster (357 aligned, 331 single,
#                                      4,568 neither -- its values are per model,
#                                      not per transcript)
#
# Indexing a column that is not transcript-aligned would silently pair a
# transcript's position with another axis's value, which is the same class of
# error this change exists to remove. MutationTaster is the one column that looks
# indexable and is not.
#
# Run from project root: Rscript analyses/tests/test_dbnsfp_mane_selection.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
suppressMessages(source(resolve("analyses/lib/local_predictors.R"), local = FALSE))

fails <- 0L
check <- function(ok, msg) {
  # isTRUE: a missing list entry yields NULL and is.na(NULL) is logical(0), which
  # would abort the run instead of reporting the failure.
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

N_COLS <- 458L
mkrow <- function(...) {
  p <- rep(".", N_COLS)
  kv <- list(...)
  for (nm in names(kv)) p[as.integer(nm)] <- kv[[nm]]
  p
}

# Three transcripts, MANE second. Every transcript-aligned column carries three
# values whose middle entry is the answer; the non-aligned columns carry one.
MANE <- "ENST00000000002"
row3 <- mkrow(
  "5"   = "E",                       # aaref, not per-transcript
  "6"   = "K",                       # aaalt, not per-transcript
  "12"  = "600;604;517",             # aapos
  "13"  = "GENEX;GENEX;GENEX",
  "15"  = "ENST00000000001;ENST00000000002;ENST00000000003",
  "16"  = "ENSP01;ENSP02;ENSP03",
  "24"  = "p.Glu600Lys;p.Glu604Lys;p.Glu517Lys",
  "38"  = "0.11;0.22;0.33",          # SIFT score
  "40"  = "T;D;T",                   # SIFT pred
  "44"  = "0.41;0.42;0.43",          # PP2 HDIV score
  "46"  = "B;D;B",                   # PP2 HDIV pred
  "47"  = "0.51;0.52;0.53",          # PP2 HVAR score
  "49"  = "B;D;B",                   # PP2 HVAR pred
  "50"  = "0.9",                     # LRT, single-valued
  "52"  = "D",
  "54"  = "0.7;0.8;0.9;0.95",        # MutationTaster, NOT transcript-aligned
  "56"  = "A;D;N;P",
  "62"  = "-1.1;-2.2;-3.3",          # FATHMM score
  "64"  = "T;D;T",
  "65"  = "-1.5;-2.5;-3.5",          # PROVEAN score
  "67"  = "N;D;N",
  "70"  = "0.66",                    # MetaSVM, single-valued
  "73"  = "0.77",                    # MetaLR, single-valued
  "77"  = "0.81;0.82;0.83",          # MetaRNN score
  "83"  = "0.71;0.72;0.73",          # REVEL score
  "138" = "0.8652;0.4761;0.1234",    # AlphaMissense score
  "140" = "P;A;B",                   # AlphaMissense pred
  "155" = "24.5",                    # CADD phred, single-valued
  "159" = "0.98",                    # DANN, single-valued
  "213" = "0.0123",                  # 1000Gp3_AF, single-valued
  "235" = "0.00456"                  # ExAC_AF, single-valued
)

cat("-- the MANE element is selected for every transcript-aligned column --\n")
d <- .dbnsfp_row_to_hit(row3, mane_tx = MANE)$dbnsfp
check(isTRUE(all.equal(d$alphamissense$am_pathogenicity, 0.4761)),
      "AlphaMissense takes the MANE element, not the first")
check(identical(d$alphamissense$am_class, "A"),
      "AlphaMissense class takes the MANE element")
check(isTRUE(all.equal(d$sift$score, 0.22)),        "SIFT score takes the MANE element")
check(identical(d$sift$pred, "D"),                  "SIFT pred takes the MANE element")
check(isTRUE(all.equal(d$polyphen2$hdiv$score, 0.42)), "PolyPhen2 HDIV takes the MANE element")
check(identical(d$polyphen2$hdiv$pred, "D"),        "PolyPhen2 HDIV pred takes the MANE element")
check(isTRUE(all.equal(d$polyphen2$hvar$score, 0.52)), "PolyPhen2 HVAR takes the MANE element")
check(isTRUE(all.equal(d$fathmm$score, -2.2)),      "FATHMM takes the MANE element")
check(isTRUE(all.equal(d$provean$score, -2.5)),     "PROVEAN takes the MANE element")
check(isTRUE(all.equal(d$metarnn$score, 0.82)),     "MetaRNN takes the MANE element")
check(isTRUE(all.equal(d$revel$score, 0.72)),       "REVEL takes the MANE element")

cat("-- columns that are not transcript-aligned are left alone --\n")
# Indexing these would pair a transcript position with a different axis. The
# guard is that they still read the first non-'.' entry.
check(isTRUE(all.equal(d$lrt$score, 0.9)),          "LRT is single-valued and unchanged")
check(isTRUE(all.equal(d$metasvm$score, 0.66)),     "MetaSVM is single-valued and unchanged")
check(isTRUE(all.equal(d$metalr$score, 0.77)),      "MetaLR is single-valued and unchanged")
check(isTRUE(all.equal(d$cadd$phred, 24.5)),        "CADD is single-valued and unchanged")
check(isTRUE(all.equal(d$dann$score, 0.98)),        "DANN is single-valued and unchanged")
check(isTRUE(all.equal(d$onekg$af, 0.0123)),        "1000Gp3_AF is single-valued and unchanged")
check(isTRUE(all.equal(d$exac$af, 0.00456)),        "ExAC_AF is single-valued and unchanged")
check(isTRUE(all.equal(d$mutationtaster$score, 0.7)),
      "MutationTaster keeps the FIRST element: its list is per model, not per transcript")

cat("-- MANE absent from the row falls back to the first element --\n")
# A transcript set that does not contain MANE at all. Returning nothing would
# lose the variant; the documented pre-change behaviour is the safe fallback.
d2 <- .dbnsfp_row_to_hit(row3, mane_tx = "ENST99999999999")$dbnsfp
check(isTRUE(all.equal(d2$alphamissense$am_pathogenicity, 0.8652)),
      "AlphaMissense falls back to the first element when MANE is not listed")
check(isTRUE(all.equal(d2$revel$score, 0.71)),
      "REVEL falls back to the first element when MANE is not listed")

cat("-- no MANE supplied reproduces the pre-change behaviour exactly --\n")
# The regression guard. Every existing caller that does not pass mane_tx, and
# every existing test, must see byte-identical output.
d3 <- .dbnsfp_row_to_hit(row3)$dbnsfp
check(isTRUE(all.equal(d3$alphamissense$am_pathogenicity, 0.8652)),
      "without mane_tx AlphaMissense is the first element")
check(isTRUE(all.equal(d3$sift$score, 0.11)), "without mane_tx SIFT is the first element")
check(isTRUE(all.equal(d3$revel$score, 0.71)), "without mane_tx REVEL is the first element")

cat("-- a '.' at the MANE position falls back rather than returning NA --\n")
# dbNSFP leaves a per-transcript entry as '.' when that tool could not score that
# transcript. Dropping to NA would lose evidence the row does carry.
row_dot <- row3
row_dot[38]  <- "0.11;.;0.33"
row_dot[138] <- ".;.;0.1234"
d4 <- .dbnsfp_row_to_hit(row_dot, mane_tx = MANE)$dbnsfp
check(isTRUE(all.equal(d4$sift$score, 0.11)),
      "a '.' at the MANE position falls back to the first scored element")
check(isTRUE(all.equal(d4$alphamissense$am_pathogenicity, 0.1234)),
      "a leading '.' pair falls through to the only scored element")

cat("-- a single-transcript row is unaffected --\n")
row1 <- mkrow("12" = "100", "13" = "GENEX", "15" = MANE,
              "38" = "0.05", "138" = "0.9", "213" = "0.5")
d5 <- .dbnsfp_row_to_hit(row1, mane_tx = MANE)$dbnsfp
check(isTRUE(all.equal(d5$sift$score, 0.05)), "single-transcript SIFT is unchanged")
check(isTRUE(all.equal(d5$alphamissense$am_pathogenicity, 0.9)),
      "single-transcript AlphaMissense is unchanged")
check(isTRUE(all.equal(d5$onekg$af, 0.5)), "single-transcript 1000Gp3_AF is unchanged")

cat("-- the gene to MANE transcript map resolves from gene_mane.tsv --\n")
# Resolved from the row's own genename column so no caller signature changes.
check(exists(".dbnsfp_mane_tx"), ".dbnsfp_mane_tx is available")
if (exists(".dbnsfp_mane_tx")) {
  # Spot values read from analyses/derived/gene_mane.tsv this session.
  check(identical(.dbnsfp_mane_tx("KCNH2"), "ENST00000262186"),
        "KCNH2 resolves to ENST00000262186, which MyVariant lists second")
  check(identical(.dbnsfp_mane_tx("TP53"), "ENST00000269305"), "TP53 resolves to its MANE")
  check(identical(.dbnsfp_mane_tx("CASR"), "ENST00000639785"), "CASR resolves to its MANE")
  check(is.na(.dbnsfp_mane_tx("NOT_A_GENE")), "an unknown gene gives NA, not an error")
  check(is.na(.dbnsfp_mane_tx(NA_character_)), "NA in gives NA out")
}

cat("-- the row mapper resolves MANE from the row when none is supplied --\n")
# .dbnsfp_row_to_hit(p, mane_tx = NA) means 'work it out from column 13'. This is
# what makes every existing call site MANE-aware without touching it.
kc <- row3
kc[13] <- "KCNH2;KCNH2;KCNH2"
kc[15] <- "ENST00000330883;ENST00000262186;ENST00000000003"
d6 <- .dbnsfp_row_to_hit(kc, mane_tx = NA)$dbnsfp
check(isTRUE(all.equal(d6$alphamissense$am_pathogenicity, 0.4761)),
      "a KCNH2 row auto-resolves to the MANE element 0.4761")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
