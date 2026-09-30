# Guards deduplication of the dual-pass classifications to one row per variant.
#
# Why this exists. The classifications table is keyed to universe rows, and the
# universe has one row per (variant, source, scoreset). On the 14-gene run that is
# 90,701 rows over 54,454 distinct variants: 20,944 variants appear more than once.
# A variant's classification does not depend on which scoreset it arrived from, and
# checking the run confirms it: of those 20,944, exactly ZERO disagree on the
# Pass-Full class, the Pass-Blind class or the points. The duplicates are pure
# repetition.
#
# So any per-variant metric computed without deduplicating weights a variant by how
# many scoresets happen to contain it. That is what produced two different answers
# for one quantity: on the same run, Panel C's Spearman rho reads 0.847/0.991 from
# panel_c_caps_reval.R, which deduplicates to 4,905 distinct variants, and
# 0.811/0.883 from 08_panel_c_caps.R and 10_compose_supp_s2.R, which do not and use
# 15,303 rows. panel_c.html quoted one and the supplement the other.
#
# Panel B is the exception and is deliberately left alone: its unit of analysis is
# (variant, scoreset), because one variant genuinely has a different DMS score in
# each scoreset. The scoreset rows come from the universe, not from the
# classifications, so deduplicating the classifications does not collapse them.
#
# Run from project root: Rscript analyses/tests/test_classifications_dedup.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
suppressMessages(source(resolve("analyses/lib/run_inputs.R"), local = FALSE))

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

cat("-- the reader exists --\n")
check(exists("read_classifications"), "read_classifications is available")

tmp <- file.path(tempdir(), paste0("dedup_", as.integer(runif(1, 1e6, 9e6))))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

wr <- function(df, name) {
  p <- file.path(tmp, name)
  utils::write.table(df, p, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
  p
}

# Two genes, one variant repeated three times with identical values, as the real
# table has it.
dup <- data.frame(
  gene = c("GENEA", "GENEA", "GENEA", "GENEB", "GENEB"),
  p_notation = c("p.A1G", "p.A1G", "p.A1G", "p.C2T", "p.D3E"),
  varviz_classification_full  = c("VUS-Low", "VUS-Low", "VUS-Low", "Pathogenic", "Benign"),
  varviz_classification_blind = c("VUS-Low", "VUS-Low", "VUS-Low", "VUS-High", "Benign"),
  varviz_pts_full = c(1, 1, 1, 10, -7),
  stringsAsFactors = FALSE)
p_dup <- wr(dup, "dup.tsv")

cat("-- identical duplicates collapse to one row per variant --\n")
d <- read_classifications(p_dup, quiet = TRUE)
check(nrow(d) == 3L, sprintf("5 rows over 3 variants collapse to 3 (got %d)", nrow(d)))
check(nrow(dplyr::distinct(d, gene, p_notation)) == nrow(d),
      "the result has exactly one row per (gene, p_notation)")
check(identical(sort(d$p_notation), c("p.A1G", "p.C2T", "p.D3E")),
      "every distinct variant survives; none is dropped")
# as.numeric: read.delim types a bare 10 as integer, and identical(10L, 10) is FALSE.
check(isTRUE(all.equal(as.numeric(d$varviz_pts_full[d$p_notation == "p.C2T"]), 10)),
      "the retained row keeps its values intact")

cat("-- the collapse is announced, not silent --\n")
# The whole failure mode was two scripts disagreeing without either saying why.
msg <- paste(capture.output(read_classifications(p_dup, quiet = FALSE), type = "message"),
             collapse = " ")
check(grepl("5", msg) && grepl("3", msg),
      "the message reports rows in and variants out")

cat("-- deduplication can be declined where the unit is not the variant --\n")
d2 <- read_classifications(p_dup, distinct = FALSE, quiet = TRUE)
check(nrow(d2) == 5L, "distinct = FALSE preserves every row")

cat("-- duplicates that DISAGREE fail closed --\n")
# Today zero of the 20,944 repeated variants disagree. If a future run produces a
# variant with two different classifications, keeping the first silently would hide a
# real defect behind a deduplication step.
bad <- dup
bad$varviz_classification_full[2] <- "Pathogenic"
p_bad <- wr(bad, "bad.tsv")
res <- tryCatch({ read_classifications(p_bad, quiet = TRUE); "no error" },
                error = function(e) conditionMessage(e))
check(!identical(res, "no error"), "disagreeing duplicates raise rather than collapse")
check(grepl("GENEA", res) || grepl("p.A1G", res),
      "the error names the offending variant, not just the count")

cat("-- a table with no duplicates is unchanged --\n")
uniq <- dup[!duplicated(dup[, c("gene", "p_notation")]), ]
d3 <- read_classifications(wr(uniq, "uniq.tsv"), quiet = TRUE)
check(nrow(d3) == 3L, "an already-distinct table passes through at the same size")

cat("-- every panel script deduplicates the classifications side --\n")
# Both members of each pair. The pairs are the hazard: Panel C diverged because
# 08_panel_c_caps.R did not deduplicate while panel_c_caps_reval.R did, and Panel B
# diverged the same way between 07_panel_b_dms.R and panel_b_dms_reval.R. Panel B
# keeps its (variant, scoreset) unit through the universe side of the join, so
# deduplicating the classifications side is right there too.
for (f in c("analyses/06_panel_a_clinical.R", "analyses/panel_a_fig_reval.R",
            "analyses/07_panel_b_dms.R", "analyses/panel_b_dms_reval.R",
            "analyses/08_panel_c_caps.R", "analyses/panel_c_caps_reval.R",
            "analyses/10_compose_supp_s2.R")) {
  p <- resolve(f)
  if (!file.exists(p)) { cat(sprintf("  [skip] %s absent\n", basename(f))); next }
  code <- sub("#.*$", "", readLines(p, warn = FALSE))
  dd <- grep(paste0("read_classifications\\(|distinct\\(cls|cls \\|> distinct\\(|",
                    "distinct\\(gene, p_notation"), code)
  check(length(dd) > 0L, sprintf("%s deduplicates to one row per variant", basename(f)))
  # Presence is not enough. A dedupe placed after the join has no effect on the
  # result, which is exactly the mistake made while writing this.
  jn <- grep("inner_join\\(|left_join\\(", code)
  if (length(dd) && length(jn)) {
    check(min(dd) < min(jn),
          sprintf("%s deduplicates BEFORE its first join, not after", basename(f)))
  }
}
cat("-- both sides of a per-variant join are deduplicated, not just one --\n")
# Deduplicating one side and joining it to a duplicated other side re-inflates the
# result. S2 deduplicated cls, joined it to an undeduplicated gnomAD universe, and
# still disagreed with Panel C on the Pass-Blind rho.
for (f in c("analyses/08_panel_c_caps.R", "analyses/10_compose_supp_s2.R")) {
  p <- resolve(f)
  if (!file.exists(p)) next
  code <- sub("#.*$", "", readLines(p, warn = FALSE))
  check(any(grepl("distinct\\(gn", code)),
        sprintf("%s deduplicates the gnomAD side too", basename(f)))
}

pr <- resolve("analyses/panel_c_caps_reval.R")
if (file.exists(pr)) {
  code <- sub("#.*$", "", readLines(pr, warn = FALSE))
  check(any(grepl("distinct\\(gene, p_notation", code)),
        "panel_c_caps_reval.R already deduplicated and still does")
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
