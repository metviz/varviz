# dms_polarity.R — per-scoreset DMS score polarity for the Panel B binarisation.
#
# A deep mutational scan produces one number per variant, and nothing in the
# number says which direction is damaging. It depends on what was measured. A
# growth assay on a tumour suppressor enriches variants that LOST function, so
# higher is worse. An abundance assay depletes them, so lower is worse. A binding
# free-energy assay reports a larger value for weaker binding, so higher is worse.
# A scoreset published in ACMG points is higher for more evidence of pathogenicity.
#
# analyses/07_panel_b_dms.R applied one convention to all 21 scoresets and said so:
# "Convention: lower score = more damaging ... Inverted assays ... would require
# per-URN sign correction", alongside "there is no URN -> paper mapping in the
# ingest yet". Nine of the 21 are inverted, so each reported 1 - its true AUROC,
# and those nine are exactly the nine below 0.5 in the published table.
#
# The directions come from analyses/determine_mavedb_polarity.py, which decides
# them from assay semantics and verified publications, and cross-checks them
# against clinical labels that are independent of both the DMS measurement and
# every predictor. Polarity is deliberately NOT inferred from agreement with a
# predictor: orienting a truth scale by whether it agrees with the predictor under
# evaluation, then scoring that predictor against it, manufactures the result.
#
# This rule lives in one file because it previously lived in three, which is how a
# fix applied to one script drifts away from the others. Same shape as the problem
# commit fbf663d addressed for the classification harness: configurable, not copied.

# The documented convention, used for any scoreset the table does not list. It is
# the majority case and it is what the panel did before per-scoreset polarity
# existed, so an unlisted scoreset behaves exactly as it used to.
DMS_POLARITY_DEFAULT <- "lower"
DMS_POLARITY_VALID   <- c("lower", "higher")

# study -> "lower" | "higher". A missing or malformed table yields an empty map
# rather than an error, so every caller degrades to the previous behaviour instead
# of losing the panel.
load_dms_polarity <- function(path = "analyses/derived/mavedb_scoreset_polarity.tsv") {
  empty <- setNames(character(0), character(0))
  p <- if (file.exists(path)) path else file.path("..", path)
  if (!file.exists(p)) return(empty)
  tb <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(tb) || !all(c("study", "damaging_direction") %in% names(tb))) return(empty)
  s <- as.character(tb$study)
  d <- as.character(tb$damaging_direction)
  ok <- !is.na(s) & nzchar(s) & !is.na(d) & d %in% DMS_POLARITY_VALID
  if (!any(ok)) return(empty)
  setNames(d[ok], s[ok])
}

dms_direction_for <- function(pol, study) {
  if (length(pol) == 0L) return(DMS_POLARITY_DEFAULT)
  i <- match(as.character(study)[1], names(pol))
  if (is.na(i)) DMS_POLARITY_DEFAULT else unname(pol[[i]])
}

# Binary functional truth: 1 = loss of function, 0 = WT-like, NA = the middle 50%
# the panel excludes. score, q25, q75 and direction are recycled to a common
# length so this drops into a grouped mutate().
#
# Tail precedence matches the case_when() it replaces: under "lower" the
# bottom-quartile branch is evaluated first, so a degenerate scoreset with
# q25 == q75 keeps calling that variant loss-of-function rather than WT-like.
# Under "higher" the mirror branch takes precedence, for the same reason.
dms_truth <- function(score, q25, q75, direction) {
  n <- max(length(score), length(q25), length(q75), length(direction))
  if (n == 0L) return(numeric(0))
  score     <- rep_len(suppressWarnings(as.numeric(score)), n)
  q25       <- rep_len(suppressWarnings(as.numeric(q25)), n)
  q75       <- rep_len(suppressWarnings(as.numeric(q75)), n)
  direction <- rep_len(as.character(direction), n)
  bad <- is.na(direction) | !direction %in% DMS_POLARITY_VALID
  if (any(bad)) {
    # Fail closed. Treating a typo as the default is precisely how an inverted
    # scoreset would be scored the wrong way round again, silently.
    stop(sprintf("unrecognised DMS polarity: %s (expected one of %s)",
                 paste(unique(direction[bad]), collapse = ", "),
                 paste(DMS_POLARITY_VALID, collapse = ", ")), call. = FALSE)
  }
  # A degenerate split is not a split. Where q25 == q75 the two tails overlap and
  # every variant satisfies both comparisons, so the result depends entirely on
  # which branch is evaluated first rather than on the data. BRCA1 00000081-a-1 is
  # the case that exposed this: 5 distinct scores over 1,059 variants with
  # q25 = q75 = 0, so under the original branch order "LOF" silently meant a score
  # of exactly 0 (876 variants) and "WT-like" meant anything above it (183), which
  # is not the documented quartile method and produced a number in either
  # direction. Refuse it instead, and let the caller drop the scoreset.
  usable <- !is.na(q25) & !is.na(q75) & q75 > q25
  lo <- usable & !is.na(score) & score <= q25
  hi <- usable & !is.na(score) & score >= q75
  up <- direction == "higher"
  out <- rep(NA_real_, n)
  # "lower": assign the WT-like tail first so the LOF branch overwrites it.
  out[hi & !up] <- 0; out[lo & !up] <- 1
  # "higher": mirrored precedence.
  out[lo &  up] <- 0; out[hi &  up] <- 1
  out
}

# Convenience for the two panel scripts that label bins rather than code them.
dms_bin_label <- function(score, q25, q75, direction,
                          labels = c(lof = "LOF", wt = "WT-like", mid = "Partial")) {
  t <- dms_truth(score, q25, q75, direction)
  ifelse(is.na(t), labels[["mid"]], ifelse(t == 1, labels[["lof"]], labels[["wt"]]))
}

# A scoreset numbered on a different sequence is a separate problem from polarity,
# and a worse one: its variants were classified at the wrong residues, so no sign
# convention makes the result interpretable. TP53 00000059-a-1 is offset by +101,
# a domain construct numbered from its own residue 1, and BRCA1 00000003-a-2 by +1.
# Offsets come from analyses/check_mavedb_offsets.py. Until the universe is rebuilt
# on reference numbering and re-classified, such a scoreset has to be dropped
# rather than reported.
load_dms_offsets <- function(path = "analyses/derived/mavedb_scoreset_offsets.tsv") {
  empty <- setNames(integer(0), character(0))
  p <- if (file.exists(path)) path else file.path("..", path)
  if (!file.exists(p)) return(empty)
  tb <- tryCatch(utils::read.delim(p, stringsAsFactors = FALSE, check.names = FALSE),
                 error = function(e) NULL)
  if (is.null(tb) || !all(c("study", "offset") %in% names(tb))) return(empty)
  s <- as.character(tb$study)
  k <- suppressWarnings(as.integer(tb$offset))
  ok <- !is.na(s) & nzchar(s) & !is.na(k)
  if (!any(ok)) return(empty)
  setNames(k[ok], s[ok])
}

dms_offset_for <- function(off, study) {
  if (length(off) == 0L) return(0L)
  i <- match(as.character(study)[1], names(off))
  if (is.na(i)) 0L else as.integer(off[[i]])
}

# Move a one-letter protein notation by a recorded scoreset offset, keeping the
# reference and alternate residues. Used by the universe build to put MaveDB
# scoresets on their gene's reference numbering: TP53 00000059-a-1 is a domain
# construct numbered from its own residue 1 and needs +101, BRCA1 00000003-a-2
# needs +1. Left uncorrected those variants are classified at the wrong residues.
#
# Fails closed. A value that is not a well-formed one-letter notation returns NA
# rather than passing through unshifted, because an unshifted row is exactly the
# defect this corrects and it would be the hardest kind to notice. A shift that
# would put the position below 1 also returns NA: that means the recorded offset is
# wrong for the row, and a negative position would silently corrupt a join key.
shift_p_notation <- function(p_notation, offset) {
  n <- max(length(p_notation), length(offset))
  if (n == 0L) return(character(0))
  p <- rep_len(as.character(p_notation), n)
  k <- rep_len(suppressWarnings(as.integer(offset)), n)
  out <- rep(NA_character_, n)
  ok <- !is.na(p) & !is.na(k) & grepl("^p\\.[A-Z*][0-9]+[A-Z*]$", p)
  if (!any(ok)) return(out)
  ref <- sub("^p\\.([A-Z*])[0-9]+[A-Z*]$", "\\1", p[ok])
  pos <- as.integer(sub("^p\\.[A-Z*]([0-9]+)[A-Z*]$", "\\1", p[ok]))
  alt <- sub("^p\\.[A-Z*][0-9]+([A-Z*])$", "\\1", p[ok])
  new <- pos + k[ok]
  res <- ifelse(is.na(new) | new < 1L, NA_character_, paste0("p.", ref, new, alt))
  out[ok] <- res
  out
}
