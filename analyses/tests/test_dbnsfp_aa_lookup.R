# Self-check for lookup_dbnsfp_by_aa() reference-residue matching.
#
# dbNSFP rows carry a ";"-separated multi-transcript aapos list, so one row can
# be indexed at several protein positions. Keying only on (aapos, aaalt) lets a
# query for G51D match a row that is really N51D in an alternate transcript --
# a different substitution whose predictor scores then drive PP3.
#
# Real collision this was found on (SNCA, chr4, canonical P37840):
#   query      G51D  -- REVEL 0.72, MetaSVM 0.456 (D), MetaLR 0.679 (D)
#   row "51_D" N51D  -- REVEL ".",  MetaSVM -0.1978 (T), MetaLR 0.4472 (T)
# The wrong row calls the variant tolerated, so the meta-predictor rule never
# upgrades PP3 to Moderate and the variant loses a point.
#
# Plain R, no network, no dbNSFP file: the index is built from synthetic rows in
# the real 458-column shape. Run from project root:
#   Rscript analyses/tests/test_dbnsfp_aa_lookup.R

source("analyses/lib/local_predictors.R")

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok" else "FAIL", msg))
  if (!isTRUE(ok)) fails <<- fails + 1L
}

# One dbNSFP row: 458 tab columns. Only the fields the index and the converter
# read are populated; the rest are the "." missing marker.
mkrow <- function(pos, ref, alt, aapos, aaref, aaalt, revel, metasvm, metasvm_pred) {
  p <- rep(".", 458)
  p[1] <- "4"; p[2] <- as.character(pos); p[3] <- ref; p[4] <- alt
  p[5] <- aaref; p[6] <- aaalt; p[12] <- aapos
  p[83] <- as.character(revel)
  p[70] <- as.character(metasvm); p[72] <- metasvm_pred
  p
}

# Two rows that collide on the (51, D) protein key:
#   the canonical G51D, and an N51D from an alternate transcript numbering.
row_n51d <- mkrow(89822359, "T", "C", "51;65;65;51;65", "N", "D", ".",   -0.1978, "T")
row_g51d <- mkrow(89828149, "G", "A", "51;51;65",       "G", "D", 0.72,   0.456,  "D")

build_env <- function(rows) {
  env   <- new.env(hash = TRUE, parent = emptyenv())
  env_p <- new.env(hash = TRUE, parent = emptyenv())
  for (p in rows) {
    assign(paste0(p[2], "_", p[3], "_", p[4]), p, envir = env)
    for (ap in unique(strsplit(p[12], ";", fixed = TRUE)[[1]])) {
      # Mirror the real loader: write both the (pos, ref, alt) key and the
      # legacy (pos, alt) key, first-write-wins on each.
      for (k in c(paste0(ap, "_", p[5], "_", p[6]), paste0(ap, "_", p[6]))) {
        if (!exists(k, envir = env_p, inherits = FALSE)) assign(k, p, envir = env_p)
      }
    }
  }
  attr(env, "dbnsfp_chrom") <- "4"
  attr(env, "protein_env")  <- env_p
  env
}
# N51D is indexed FIRST, so a key that ignores aaref would return it for G51D.
e <- build_env(list(row_n51d, row_g51d))

# 1. the queried substitution wins, not the position-only collision
h <- lookup_dbnsfp_by_aa(e, "4", 51, "D", aa_ref = "G")
check(!is.null(h), "G51D resolves to a row")
if (!is.null(h)) {
  check(isTRUE(all.equal(h$dbnsfp$revel$score, 0.72)),
        "G51D gets its own REVEL 0.72, not the N51D row's missing value")
  check(identical(h$dbnsfp$metasvm$pred, "D"),
        "G51D gets its own MetaSVM damaging call, not the N51D tolerated call")
}

# 2. the other substitution at the same position still resolves to itself
h2 <- lookup_dbnsfp_by_aa(e, "4", 51, "D", aa_ref = "N")
check(!is.null(h2) && identical(h2$dbnsfp$metasvm$pred, "T"),
      "N51D still resolves to the N51D row")

# 3. a reference residue that matches nothing returns no hit rather than a
#    neighbouring row -- silence is recoverable, a wrong score is not
check(is.null(lookup_dbnsfp_by_aa(e, "4", 51, "D", aa_ref = "W")),
      "unmatched reference residue returns NULL, never a different variant")

# 4. omitting aa_ref keeps the old position-only behaviour for callers that
#    genuinely cannot supply it
check(!is.null(lookup_dbnsfp_by_aa(e, "4", 51, "D")),
      "aa_ref omitted still resolves (back-compatible)")

# 5. wrong chromosome is refused
check(is.null(lookup_dbnsfp_by_aa(e, "7", 51, "D", aa_ref = "G")),
      "chromosome mismatch returns NULL")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
