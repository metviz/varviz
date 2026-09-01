# Self-check for the harness wrappers around analyses/05_classify_harness.R.
#
# Each analyses/ps_*_harness.R and analyses/test_weekly_variants.R is a few
# lines of configuration over one shared harness. This asserts that every
# wrapper still resolves to the run it is named for: the same output paths,
# the same server.R, the same universe and the same varviz.* options the
# hand-copied versions used to hard-code.
#
# --dry-run prints the resolved configuration and exits before server.R is
# sourced, so this needs no network, no data files and about a second.
#
# Run from project root: Rscript analyses/tests/test_harness_shims.R

EXPECTED <- list(
  ps_baseline_harness.R = list(
    server = "server.R.bak.20260729_214027.pre-ps-hotspot",
    universe = "analyses/derived/variant_universe_gnomad.tsv",
    ckpt = "analyses/ps_baseline/classifications",
    summary = "analyses/ps_baseline/summary.tsv",
    options = character(0)),
  ps_nomds_harness.R = list(
    server = "server.R",
    universe = "analyses/derived/variant_universe_gnomad.tsv",
    ckpt = "analyses/ps_nomds/classifications",
    summary = "analyses/ps_nomds/summary.tsv",
    options = c("varviz.mds_pm1 = FALSE")),
  ps_reval_harness.R = list(
    server = "server.R",
    universe = "analyses/derived/variant_universe_gnomad.tsv",
    ckpt = "analyses/ps_reval/classifications",
    summary = "analyses/ps_reval/summary.tsv",
    options = character(0)),
  ps_tiered_harness.R = list(
    server = "server.R",
    universe = "analyses/derived/variant_universe_gnomad.tsv",
    ckpt = "analyses/ps_tiered/classifications",
    summary = "analyses/ps_tiered/summary.tsv",
    options = c("varviz.mds_tiered = TRUE")),
  ps_mds_corroborate_harness.R = list(
    server = "server.R",
    universe = "analyses/derived/variant_universe_gnomad.tsv",
    ckpt = "analyses/ps_mds_corroborate/classifications",
    summary = "analyses/ps_mds_corroborate/summary.tsv",
    options = c("varviz.mds_tiered = TRUE", "varviz.mds_pm1 = TRUE")),
  ps_mds_consfree_harness.R = list(
    server = "server.R",
    universe = "analyses/derived/variant_universe_gnomad.tsv",
    ckpt = "analyses/ps_mds_consfree/classifications",
    summary = "analyses/ps_mds_consfree/summary.tsv",
    options = c("varviz.mds_tiered = TRUE", "varviz.mds_pm1 = TRUE",
                "varviz.mds_frees_cons = TRUE")),
  test_weekly_variants.R = list(
    server = "server.R",
    universe = "analyses/weekly_test/universe.tsv",
    ckpt = "analyses/weekly_test/classifications",
    summary = "analyses/weekly_test/summary.tsv",
    options = character(0))
)

if (!file.exists("analyses/05_classify_harness.R"))
  stop("run from the project root")

grab <- function(out, pattern) {
  hit <- grep(pattern, out, value = TRUE)
  if (length(hit) != 1) stop("expected exactly one line matching ", pattern)
  sub(pattern, "\\1", hit)
}

failures <- 0L
for (script in names(EXPECTED)) {
  want <- EXPECTED[[script]]
  out <- suppressWarnings(system2("Rscript", c(file.path("analyses", script), "--dry-run"),
                                  stdout = TRUE, stderr = TRUE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    cat(sprintf("FAIL %s — exited %d\n%s\n", script, status, paste(out, collapse = "\n")))
    failures <- failures + 1L
    next
  }

  got <- list(
    server   = grab(out, "^\\[harness\\] run=.* server=(.*)$"),
    universe = grab(out, "^\\[harness\\] universe=(.*)$"),
    ckpt     = grab(out, "^\\[harness\\] checkpoints=(.*)$"),
    summary  = grab(out, "^\\[harness\\] summary=(.*)$"))
  opt_line <- grep("^\\[harness\\] options: ", out, value = TRUE)
  got$options <- if (length(opt_line) == 0) character(0)
                 else trimws(strsplit(sub("^\\[harness\\] options: ", "", opt_line), " ; ")[[1]])

  bad <- names(want)[!vapply(names(want),
                             function(k) identical(sort(want[[k]]), sort(got[[k]])),
                             logical(1))]
  if (length(bad)) {
    failures <- failures + 1L
    cat(sprintf("FAIL %s\n", script))
    for (k in bad)
      cat(sprintf("       %-9s want %s\n       %-9s got  %s\n",
                  k, paste(want[[k]], collapse = " ; "),
                  "", paste(got[[k]], collapse = " ; ")))
  } else {
    cat(sprintf("  ok  %-30s -> %s\n", script, got$ckpt))
  }
}

if (failures) stop(sprintf("%d harness wrapper(s) no longer resolve to their run", failures))
cat("harness wrappers: all checks pass\n")
