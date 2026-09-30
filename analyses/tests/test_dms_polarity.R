# Guards per-scoreset DMS score polarity in the Panel B binarisation.
#
# Why this exists. Three scripts collapse a MaveDB score into a binary functional
# truth, and all three carried the same inline rule: bottom quartile = loss of
# function, top quartile = WT-like. analyses/07_panel_b_dms.R:14-19 states the
# assumption and its limit outright: "Convention: lower score = more damaging, the
# dominant convention across MaveDB function/fitness assays. Inverted assays (e.g.
# abundance-as-stability where higher = more stable) would require per-URN sign
# correction."
#
# Nine of the 21 scoresets are inverted, established in
# analyses/determine_mavedb_polarity.py from each scoreset's assay semantics and,
# where enough of them exist, from clinical labels that are independent of both the
# DMS measurement and every predictor:
#
#   TP53  00000059-a-1  cell-growth fitness over the DNA-binding domain. p53 is a
#                       tumour suppressor, so a variant that LOST function grows
#                       faster and scores higher.
#   BRCA1 00000081-a-1  homology-directed repair DEPLETION score.
#   BRCA1 00000093-a-1  values reported in the BASE_ACMG POINTS system, where
#                       positive points are evidence FOR pathogenicity.
#   KRAS  00000115-a-1  binding free energy change against six partners; a larger
#         ... a-6       ddG is weaker binding.
#
# An inverted scoreset reports 1 - its true AUROC. The published table has nine
# scoresets below 0.5 and they are exactly these nine.
#
# The rule lives in one place now. Three copies of it is how a half-fix survives:
# the same shape of defect is what commit fbf663d addressed for the classification
# harness, "configurable, not copied".
#
# Run from project root: Rscript analyses/tests/test_dms_polarity.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
suppressMessages(source(resolve("analyses/lib/dms_polarity.R"), local = FALSE))

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

cat("-- the helpers exist --\n")
check(exists("load_dms_polarity"), "load_dms_polarity is available")
check(exists("dms_truth"),         "dms_truth is available")

cat("-- lower-is-damaging reproduces the previous inline rule exactly --\n")
# truth = 1 for loss of function, 0 for WT-like, NA for the middle 50%.
s <- c(0, 1, 2, 3, 4, 5, 6, 7, 8)
q25 <- 2; q75 <- 6
old <- ifelse(s <= q25, 1, ifelse(s >= q75, 0, NA_real_))
new <- dms_truth(s, q25, q75, "lower")
check(identical(as.numeric(new), as.numeric(old)),
      "lower: identical to bottom-quartile-is-LOF for a full range of scores")

cat("-- higher-is-damaging inverts the two tails and nothing else --\n")
inv <- dms_truth(s, q25, q75, "higher")
check(identical(as.numeric(inv), as.numeric(ifelse(s >= q75, 1, ifelse(s <= q25, 0, NA_real_)))),
      "higher: the top quartile becomes LOF and the bottom becomes WT-like")
check(identical(which(is.na(new)), which(is.na(inv))),
      "the excluded middle 50% is the same set under either direction")
check(sum(new == 1, na.rm = TRUE) == sum(inv == 0, na.rm = TRUE),
      "the class counts swap rather than change size")

cat("-- the quartile boundaries are inclusive on both sides --\n")
# The inline rule used <= q25 and >= q75; a variant sitting exactly on a boundary
# must keep its class, or the corrected panel would silently lose variants.
check(identical(as.numeric(dms_truth(c(2, 6), 2, 6, "lower")), c(1, 0)),
      "lower: a score exactly on q25 is LOF and exactly on q75 is WT-like")
check(identical(as.numeric(dms_truth(c(2, 6), 2, 6, "higher")), c(0, 1)),
      "higher: the same two boundary scores swap class")

cat("-- missing values stay missing --\n")
check(is.na(dms_truth(NA_real_, 2, 6, "lower")),  "an NA score gives NA under lower")
check(is.na(dms_truth(NA_real_, 2, 6, "higher")), "an NA score gives NA under higher")
check(all(is.na(dms_truth(c(1, 5), NA_real_, NA_real_, "lower"))),
      "NA quartiles give NA rather than an error")

cat("-- a degenerate quartile split is refused in either direction --\n")
# q25 == q75 makes both tail comparisons true for every variant, so the answer
# comes from branch order rather than from the data. BRCA1 00000081-a-1 is real:
# 1,059 variants over 5 distinct scores with q25 = q75 = 0. Under the original
# branch order that silently reported "LOF" for a score of exactly 0 and
# "WT-like" for anything above, which is not a quartile split, and it produced a
# published AUROC of 0.379 that was never valid in either direction.
deg <- c(0, 0, 0, 1, 2, 4)
check(all(is.na(dms_truth(deg, 0, 0, "lower"))),
      "q25 == q75 yields all NA under lower rather than a branch-order artefact")
check(all(is.na(dms_truth(deg, 0, 0, "higher"))),
      "q25 == q75 yields all NA under higher too, so the refusal is symmetric")
check(!any(is.na(dms_truth(c(0, 5), 1, 4, "lower"))),
      "a genuine split with q75 > q25 still classifies both tails")

cat("-- an unknown direction fails closed --\n")
# Silently treating a typo as lower-is-damaging is how an inverted scoreset would
# slip back in.
bad <- tryCatch({ dms_truth(s, q25, q75, "sideways"); "no error" },
                error = function(e) "error")
check(identical(bad, "error"), "an unrecognised direction raises rather than defaulting")

cat("-- the polarity table loads, and is keyed by scoreset --\n")
tmp <- file.path(tempdir(), paste0("pol_", as.integer(runif(1, 1e6, 9e6)), ".tsv"))
on.exit(unlink(tmp), add = TRUE)
writeLines(c("gene\tstudy\tdamaging_direction",
             "TP53\t00000059-a-1\thigher",
             "GCK\t00000096-a-1\tlower"), tmp)
p <- load_dms_polarity(tmp)
check(identical(unname(p[["00000059-a-1"]]), "higher"), "an inverted scoreset reads as higher")
check(identical(unname(p[["00000096-a-1"]]), "lower"),  "a normal scoreset reads as lower")

cat("-- a scoreset absent from the table defaults to lower, the documented convention --\n")
check(identical(dms_direction_for(p, "99999999-z-9"), "lower"),
      "an unlisted scoreset falls back to lower-is-damaging")
check(identical(dms_direction_for(p, "00000059-a-1"), "higher"),
      "a listed scoreset uses its recorded direction")

cat("-- a missing table degrades to the old behaviour rather than erroring --\n")
p0 <- load_dms_polarity(file.path(tempdir(), "definitely_absent_polarity.tsv"))
check(length(p0) == 0L, "an absent table yields an empty map")
check(identical(dms_direction_for(p0, "00000059-a-1"), "lower"),
      "with no table every scoreset reads as lower, i.e. the pre-change panel")

cat("-- all three Panel B scripts use the shared rule, not their own copy --\n")
# The reason the defect survived: the same quartile rule was written out three
# times. Any script still computing it inline can be corrected in isolation and
# drift again.
# 10_compose_supp_s2.R carried a fourth copy and is checked when present. It is
# untracked in this repository, so its absence is not a failure.
for (f in c("analyses/07_panel_b_dms.R", "analyses/panel_b_dms_reval.R",
            "analyses/panel_b_auroc_summary.R", "analyses/10_compose_supp_s2.R")) {
  if (!file.exists(resolve(f))) {
    cat(sprintf("  [skip] %s is not present\n", basename(f)))
    next
  }
  src <- readLines(resolve(f), warn = FALSE)
  code <- sub("#.*$", "", src)
  check(any(grepl("dms_polarity\\.R", code)), sprintf("%s sources the shared helper", basename(f)))
  # Either entry point is the shared rule: dms_bin_label() labels the bins and
  # calls dms_truth() itself. What must not survive is a script computing the
  # quartile comparison inline, which the next assertion is the real guard on.
  check(any(grepl("dms_truth\\(|dms_bin_label\\(", code)),
        sprintf("%s goes through the shared rule", basename(f)))
  check(!any(grepl("score_dms *<= *q25", code)),
        sprintf("%s no longer hardcodes score_dms <= q25", basename(f)))
  # A scoreset on different numbering was classified at the wrong residues. If one
  # script drops it and another plots it, the figure and the table disagree about
  # what the panel contains, which is how the three copies diverged before.
  # dms_pending_offset() is the current rule: exclude on the offset still
  # OUTSTANDING, so a scoreset the universe build has already corrected is kept.
  # dms_offset_for() is still valid where the universe is uncorrected.
  check(any(grepl("dms_offset_for\\(|dms_pending_offset\\(", code)),
        sprintf("%s excludes scoresets whose numbering is still uncorrected", basename(f)))
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
