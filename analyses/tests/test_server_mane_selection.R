# Guards MANE-aware selection on the LIVE application path, and the two defects
# that sat beside it.
#
# Why this exists. analyses/lib/local_predictors.R now reads dbNSFP's
# per-transcript columns at the MANE index, but the live app does not go through
# that mapper. It reads a MyVariant.info response, and MyVariant IS dbNSFP: it
# returns the same per-transcript arrays. Confirmed live against
# myvariant.info/v1/query for KCNH2:
#
#   dbnsfp.ensembl.transcriptid : ["ENST00000330883", "ENST00000262186"]
#   dbnsfp.alphamissense.score  : [0.8652, 0.4761]
#   dbnsfp.sift.score           : [0.009, 0.0]
#   dbnsfp.polyphen2.hvar.score : [0.971, 0.952]
#
# KCNH2's MANE Select is ENST00000262186, the SECOND element, so MANE's
# AlphaMissense is 0.4761 and the first element 0.8652 belongs to another
# isoform. safe_extract() ends with result[1], so every predictor on the live
# path took the wrong transcript for the genes where MANE is not listed first:
# GCK, TP53, DDR2, TSHR, KCNH2, BRCA1 and KCNQ1 of the fourteen benchmark genes.
#
# Two further defects found while tracing that path:
#
# 1. The AlphaMissense key. parse_dbnsfp_scores read
#    safe_extract_num(d, "alphamissense", "am_pathogenicity"). MyVariant names
#    that field "score"; only the local positional mapper names it
#    am_pathogenicity. So the live per-variant path got NA for AlphaMissense and
#    silently fell back to the AlphaFold CSV. Both spellings must resolve.
#
# 2. AM_Source. The export sets AM_Source to "MyVariant" when
#    !is.na(am$s), but the sv() helper returns s = "" (an empty string, not NA)
#    for an absent score, and is.na("") is FALSE. 10,570 of the 90,701 rows in
#    analyses/ps_am792_v123_final/master_all_variants.tsv carry
#    AM_Source = "MyVariant" with an empty AM_Score as a result: a provenance
#    column naming a source that supplied nothing.
#
# Run from project root: Rscript analyses/tests/test_server_mane_selection.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src <- readLines(resolve("server.R"), warn = FALSE)

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# Eval the shipped definitions rather than a copy, the idiom of
# analyses/tests/test_run_manifest.R and test_acmg_posterior.R. Commit 4c93a92 is
# the cautionary case: a test that re-implemented the rule passed while the
# shipped rule was dead code.
load_fn <- function(name) {
  beg <- grep(sprintf("^%s <- function", name), src)
  if (length(beg) != 1L) return(FALSE)
  rest <- src[beg:length(src)]
  end <- beg - 1L + which(rest == "}")[1]
  if (is.na(end)) return(FALSE)
  eval(parse(text = paste(src[beg:end], collapse = "\n")), envir = globalenv())
  TRUE
}

cat("-- the helpers exist and are defined once --\n")
for (fn in c("mane_tx_for_gene", "dbnsfp_mane_index", "safe_extract_tx",
             "safe_extract_num_tx")) {
  n <- length(grep(sprintf("^%s <- function", fn), src))
  check(n == 1L, sprintf("%s is defined exactly once (found %d)", fn, n))
}
ok_all <- all(vapply(c("safe_extract", "safe_extract_tx", "safe_extract_num_tx",
                       "mane_tx_for_gene", "dbnsfp_mane_index"), load_fn, logical(1)))
if (!ok_all) {
  cat(sprintf("\nFAIL: %d failure(s) (helpers could not be loaded)\n", fails + 1L))
  quit(status = 1L)
}

# A stand-in for the app's gene_data, so the test needs no 21 MB RData load.
gd <- data.frame(
  gene_name    = c("KCNH2", "TP53", "GCK"),
  mane_ensembl = c("ENST00000262186.10", "ENST00000269305.9", "ENST00000403799.8"),
  stringsAsFactors = FALSE)

cat("-- the gene to MANE transcript map --\n")
check(identical(mane_tx_for_gene("KCNH2", gd), "ENST00000262186"),
      "KCNH2 resolves to ENST00000262186 with the version suffix dropped")
check(is.na(mane_tx_for_gene("NOT_A_GENE", gd)), "an unknown gene gives NA")
check(is.na(mane_tx_for_gene(NA_character_, gd)), "NA in gives NA out")
check(is.na(mane_tx_for_gene("KCNH2", data.frame(gene_name = "KCNH2"))),
      "a gene_data without mane_ensembl degrades to NA, it does not error")

# The live KCNH2 response, verbatim field names and order from myvariant.info.
kcnh2 <- list(dbnsfp = list(
  genename = "KCNH2",
  ensembl  = list(transcriptid = c("ENST00000330883", "ENST00000262186"),
                  proteinid    = c("ENSP00000328531", "ENSP00000262186")),
  alphamissense = list(score = c(0.8652, 0.4761), pred = c("P", "A")),
  sift          = list(score = c(0.009, 0.0),     pred = c("D", "D")),
  polyphen2     = list(hvar  = list(score = c(0.971, 0.952), pred = c("D", "D"))),
  revel         = list(score = c(0.71, 0.72)),
  metasvm       = list(score = 0.66)))

cat("-- the MANE index is found in a live MyVariant hit --\n")
i <- dbnsfp_mane_index(kcnh2$dbnsfp, gd)
check(identical(as.integer(i), 2L), "KCNH2's MANE is the second transcript listed")
check(is.na(dbnsfp_mane_index(NULL, gd)), "a NULL hit gives NA, not an error")
check(is.na(dbnsfp_mane_index(list(genename = "KCNH2"), gd)),
      "a hit with no transcript list gives NA")
check(is.na(dbnsfp_mane_index(list(genename = "ZZZ",
       ensembl = list(transcriptid = c("ENST1", "ENST2"))), gd)),
      "a gene absent from gene_data gives NA")

cat("-- per-transcript arrays are read at the MANE index --\n")
d <- kcnh2$dbnsfp
check(isTRUE(all.equal(safe_extract_num_tx(d, i, "alphamissense", "score"), 0.4761)),
      "AlphaMissense reads 0.4761 (MANE), not 0.8652 (first)")
check(isTRUE(all.equal(safe_extract_num_tx(d, i, "sift", "score"), 0.0)),
      "SIFT reads the MANE element")
check(isTRUE(all.equal(safe_extract_num_tx(d, i, "polyphen2", "hvar", "score"), 0.952)),
      "PolyPhen2 HVAR reads the MANE element")
check(isTRUE(all.equal(safe_extract_num_tx(d, i, "revel", "score"), 0.72)),
      "REVEL reads the MANE element")
check(identical(safe_extract_tx(d, i, "alphamissense", "pred"), "A"),
      "the AlphaMissense class reads the MANE element")

cat("-- a single-valued column is unaffected by an index past its end --\n")
# MetaSVM, MetaLR, LRT, CADD, DANN and the conservation and frequency columns
# carry one value for the whole row; an index of 2 must not blank them.
check(isTRUE(all.equal(safe_extract_num_tx(d, i, "metasvm", "score"), 0.66)),
      "MetaSVM still reads its single value when the MANE index is 2")

cat("-- NA index reproduces safe_extract exactly --\n")
# The regression guard: with no MANE resolved, behaviour is the old first-element
# read, so a gene outside gene_data is never made worse.
check(isTRUE(all.equal(safe_extract_num_tx(d, NA_integer_, "alphamissense", "score"), 0.8652)),
      "an NA index falls back to the first element")
check(identical(safe_extract_tx(d, NA_integer_, "sift", "pred"),
                safe_extract(d, "sift", "pred")),
      "an NA index matches safe_extract on the same path")

cat("-- an unscored MANE element falls back rather than going NA --\n")
dot <- d
dot$alphamissense$score <- c(0.8652, ".")
check(isTRUE(all.equal(safe_extract_num_tx(dot, i, "alphamissense", "score"), 0.8652)),
      "a '.' at the MANE position falls back to the first scored element")

cat("-- both AlphaMissense spellings resolve --\n")
# MyVariant says "score"; the local positional mapper says "am_pathogenicity".
code <- sub("#.*$", "", src)
am_lines <- grep("alphamissense", code)
check(any(grepl('"score"', code[am_lines])),
      'some AlphaMissense read uses the MyVariant key "score"')
check(any(grepl('"am_pathogenicity"', code[am_lines])),
      'some AlphaMissense read still uses the local key "am_pathogenicity"')
parse_blk <- grep("am_score *<- ", code)
check(length(parse_blk) >= 1L, "am_score is assigned in server.R")
if (length(parse_blk)) {
  win <- paste(code[parse_blk[1]:min(parse_blk[1] + 3L, length(code))], collapse = " ")
  check(grepl("am_pathogenicity", win, fixed = TRUE) && grepl('"score"', win, fixed = TRUE),
        "the am_score assignment tries both spellings")
}

cat("-- the parse path reads predictors at the MANE index --\n")
check(any(grepl("dbnsfp_mane_index\\(", code)),
      "parse_dbnsfp_scores computes a MANE index")
for (p in c("sift", "revel", "alphamissense")) {
  hits <- grep(sprintf('safe_extract_num_tx\\([^)]*"%s"', p), code)
  check(length(hits) >= 1L, sprintf("%s is read through the MANE-aware extractor", p))
}

cat("-- AM_Source no longer names a source that supplied nothing --\n")
# is.na("") is FALSE, so the old guard passed for every absent score.
srcline <- grep("AM_Source *=", code)
check(length(srcline) >= 1L, "AM_Source is assigned")
if (length(srcline)) {
  blk <- paste(code[srcline], collapse = " ")
  check(!grepl("!is.na(am$s)) \"MyVariant\"", blk, fixed = TRUE),
        "the bare !is.na(am$s) guard is gone")
  check(grepl("nzchar", blk, fixed = TRUE),
        "AM_Source tests the score with nzchar, which rejects an empty string")
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
