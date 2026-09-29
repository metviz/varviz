# Repo-wide invariants for the analysis scripts, so a rule cannot be fixed in one
# place and left wrong in another.
#
# Why this exists. Every defect in the DMS and dbNSFP work was the same mistake:
# the site being looked at was fixed and the others were not enumerated first.
#
#   the MANE per-transcript read was fixed in analyses/lib/local_predictors.R and
#     server.R was left on the first element
#   the polarity rule was applied to three scripts, then a fourth was found in
#     10_compose_supp_s2.R, then a fifth in repro/07_panel_b_dms_run.R
#   the offset exclusion was applied to the AUROC table and not to the two figure
#     scripts, so the table and the figure disagreed about what the panel contained
#
# Each was found by a later grep rather than by a check. These are those greps,
# written down. A new copy of a rule now fails the suite instead of surviving until
# someone notices a number that does not add up.
#
# Run from project root: Rscript analyses/tests/test_analysis_consistency.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
root <- resolve("analyses")
if (!dir.exists(root)) { cat("analyses/ not found\n"); quit(status = 1L) }

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}
warn <- function(msg) cat(sprintf("  [warn] %s\n", msg))

# Live analysis scripts only. Archived runs are frozen history, .bak files are
# backups, and the tests themselves quote the patterns they guard.
scripts <- list.files(root, pattern = "\\.(R|py)$", recursive = TRUE, full.names = TRUE)
scripts <- scripts[!grepl("runs_archive|/tests/|\\.bak|zenodo_deposit", scripts)]
code_of <- function(f) sub("#.*$", "", readLines(f, warn = FALSE))

cat(sprintf("-- %d live analysis scripts --\n", length(scripts)))

cat("-- no script computes the DMS quartile rule inline --\n")
# The rule that was copied five times. Any new copy lands here.
inline <- Filter(function(f) {
  cd <- code_of(f)
  any(grepl("score_dms[[:space:]]*<=[[:space:]]*q25", cd)) ||
  any(grepl("score_dms[[:space:]]*>=[[:space:]]*q75", cd))
}, scripts)
check(length(inline) == 0L,
      sprintf("zero inline quartile comparisons (found %d%s)", length(inline),
              if (length(inline)) paste0(": ", paste(basename(inline), collapse = ", ")) else ""))

cat("-- every script that bins DMS scores uses the shared rule and drops shifted scoresets --\n")
binners <- Filter(function(f) any(grepl("quantile\\([[:space:]]*score_dms", code_of(f))), scripts)
check(length(binners) > 0L, sprintf("found the DMS binning scripts (%d)", length(binners)))
for (f in binners) {
  cd <- code_of(f); b <- basename(f)
  check(any(grepl("dms_polarity\\.R", cd)), sprintf("%s sources dms_polarity.R", b))
  check(any(grepl("dms_bin_label\\(|dms_truth\\(", cd)), sprintf("%s goes through the shared rule", b))
  # A scoreset on different numbering was classified at the wrong residues. One
  # script dropping it while another plots it is how the figure and the table
  # came to disagree.
  check(any(grepl("dms_offset_for\\(", cd)), sprintf("%s excludes shifted scoresets", b))
}

cat("-- the local dbNSFP mapper indexes only transcript-aligned columns --\n")
lp <- resolve("analyses/lib/local_predictors.R")
if (file.exists(lp)) {
  cd <- paste(code_of(lp), collapse = "\n")
  # Verified aligned with Ensembl_transcriptid over 5,256 multi-transcript rows.
  for (col in c("38", "40", "44", "46", "47", "49", "62", "64", "65", "67", "77", "83", "138", "140")) {
    check(grepl(sprintf("(num_tx|chr_tx)\\(%s\\)", col), cd),
          sprintf("column %s is read MANE-aware", col))
  }
  # MutationTaster looks indexable and is not: its list is per prediction model,
  # 357 rows aligned against 4,568 that matched neither count.
  for (col in c("50", "52", "54", "56", "70", "72", "73", "75", "155", "159", "213", "235")) {
    check(grepl(sprintf("(num_first|chr_first)\\(%s\\)", col), cd),
          sprintf("column %s is NOT transcript-indexed", col))
  }
}

cat("-- the live app reads the same columns MANE-aware --\n")
sv <- resolve("server.R")
if (file.exists(sv)) {
  cd <- paste(code_of(sv), collapse = "\n")
  for (p in c("sift", "fathmm", "provean", "metarnn", "revel", "alphamissense")) {
    check(grepl(sprintf('safe_extract(_num)?_tx\\([^)]*"%s"', p), cd),
          sprintf("%s is read through the MANE-aware extractor", p))
  }
  check(grepl('safe_extract_num_tx\\(d, mane_i, "polyphen2"', cd),
        "polyphen2 is read through the MANE-aware extractor")
  # These carry one value per row; indexing them would pair a transcript position
  # with a different axis.
  for (p in c("metasvm", "metalr", "lrt")) {
    check(!grepl(sprintf('safe_extract(_num)?_tx\\([^)]*"%s"', p), cd),
          sprintf("%s is NOT transcript-indexed", p))
  }
}

cat("-- polarity and offset are recorded for every scoreset in the universe --\n")
pol_p <- resolve("analyses/derived/mavedb_scoreset_polarity.tsv")
off_p <- resolve("analyses/derived/mavedb_scoreset_offsets.tsv")
uni_p <- resolve("analyses/derived/variant_universe_gnomad.tsv")
if (all(file.exists(c(pol_p, off_p, uni_p)))) {
  pol <- utils::read.delim(pol_p, stringsAsFactors = FALSE)
  off <- utils::read.delim(off_p, stringsAsFactors = FALSE)
  uni <- utils::read.delim(uni_p, stringsAsFactors = FALSE,
                           colClasses = "character")[, c("source", "study")]
  studies <- sort(unique(uni$study[uni$source == "MaveDB" & nzchar(uni$study)]))
  check(length(studies) > 0L, sprintf("universe lists %d MaveDB scoresets", length(studies)))
  miss_p <- setdiff(studies, pol$study)
  miss_o <- setdiff(studies, off$study)
  check(length(miss_p) == 0L,
        sprintf("every scoreset has a recorded polarity (missing %d%s)", length(miss_p),
                if (length(miss_p)) paste0(": ", paste(miss_p, collapse = ", ")) else ""))
  check(length(miss_o) == 0L,
        sprintf("every scoreset has a recorded offset (missing %d%s)", length(miss_o),
                if (length(miss_o)) paste0(": ", paste(miss_o, collapse = ", ")) else ""))
  check(all(pol$damaging_direction %in% c("lower", "higher")),
        "every recorded direction is one of lower or higher")
  check("direction_evidence" %in% names(pol) && all(nzchar(pol$direction_evidence)),
        "every direction records how it was established")
} else {
  warn("polarity/offset/universe tables absent; regenerate with analyses/check_mavedb_offsets.py and analyses/determine_mavedb_polarity.py")
}

cat("-- run-specific inputs go through the resolver, not a hardcoded path --\n")
# analyses/derived/varviz_classifications.tsv was six months stale while its row
# count still matched the current run, which is what made it dangerous. On the real
# repository 21 run directories are newer than it. Five scripts read it with no way
# to select a run; they now resolve it, and the fallback announces its own staleness.
#
# analyses/repro/03_which_run.R is exempt by design. It scores several candidate runs
# against the manuscript's claims, so the consolidated file being one hardcoded
# candidate among five is the question it answers, not a defect.
CLS <- "analyses/derived/varviz_classifications.tsv"
# 03_which_run.R: its candidate list IS the question it answers.
# run_inputs.R: it defines the resolver, so it necessarily names the default.
EXEMPT <- c("03_which_run.R", "run_inputs.R")
unresolved <- character(0)
for (f in scripts) {
  cd <- code_of(f); b <- basename(f)
  hit <- grep(CLS, cd, fixed = TRUE)
  if (!length(hit)) next
  if (b %in% EXEMPT) next
  if (any(grepl("SUMMARY_OUT|OUT[[:space:]]*<-", cd[hit]))) next   # a writer
  if (any(grepl("resolve_classifications\\(", cd))) next          # resolved
  win <- unique(unlist(lapply(hit, function(i) seq(max(1L, i - 2L), min(length(cd), i + 1L)))))
  if (any(grepl("Sys.getenv", cd[win]))) next
  unresolved <- c(unresolved, b)
}
check(length(unresolved) == 0L,
      sprintf("no script reads the consolidated classifications unresolved (found %d%s)",
              length(unresolved),
              if (length(unresolved)) paste0(": ", paste(unique(unresolved), collapse = ", ")) else ""))

cat("-- the derived classification chain is still a single global path --\n")
# 11_merge_dolphin.R now resolves its input, but the files it and 22_mds_benchmark.R
# write are fixed paths, so a stale _dolphin or _mds file carries the same hazard one
# link down the chain. Reported, not enforced: giving these a run directory changes
# their interfaces and the chain's shape, which is a design decision.
chain <- c("analyses/derived/varviz_classifications_dolphin.tsv",
           "analyses/derived/varviz_classifications_mds.tsv")
downstream <- character(0)
for (f in scripts) {
  cd <- code_of(f)
  if (any(vapply(chain, function(p) any(grepl(p, cd, fixed = TRUE)), logical(1)))) {
    if (!any(grepl("resolve_classifications\\(|Sys.getenv|args\\[1\\]", cd)))
      downstream <- c(downstream, basename(f))
  }
}
if (length(downstream)) {
  warn(sprintf("%d script(s) read the derived chain at a fixed path: %s",
               length(downstream), paste(unique(downstream), collapse = ", ")))
  warn("  regenerate the chain from the intended run before quoting their numbers")
} else {
  check(TRUE, "every derived-chain reader can be pointed at a run")
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
