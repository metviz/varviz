# Guards a provenance stamp on every generated report.
#
# Why this exists. analyses/derived/varviz_classifications.tsv sat six months stale
# while 21 run directories were newer, and nothing noticed because a report states
# numbers without stating where they came from. Reading panel_a.html could not tell
# you which run produced it. The resolver prints its path to stdout, which is gone the
# moment the script exits.
#
# A report that carries its own run, universe fingerprint and code commit can be
# checked by a reader months later. This is the same content the harness already writes
# per gene in <gene>__manifest.json, so a report and a run can be matched up.
#
# Run from project root: Rscript analyses/tests/test_run_provenance.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
suppressMessages(source(resolve("analyses/lib/run_inputs.R"), local = FALSE))

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

VARS <- c("VARVIZ_CLASSIFICATIONS", "REVAL_DIR", "VARVIZ_RUN")
saved <- Sys.getenv(VARS, unset = NA_character_, names = TRUE)
on.exit({
  for (v in VARS) Sys.unsetenv(v)
  for (v in VARS) if (!is.na(saved[[v]])) do.call(Sys.setenv, setNames(list(saved[[v]]), v))
}, add = TRUE)

cat("-- the helpers exist --\n")
check(exists("run_provenance_lines"), "run_provenance_lines is available")
check(exists("run_provenance_html"),  "run_provenance_html is available")
check(exists("run_provenance_md"),    "run_provenance_md is available")

cat("-- the stamp names the run, the inputs and the code --\n")
for (v in VARS) Sys.unsetenv(v)
Sys.setenv(REVAL_DIR = "analyses/ps_am792_v124_refnum/classifications")
ln <- run_provenance_lines()
txt <- paste(ln, collapse = " | ")
check(length(ln) >= 4L, sprintf("at least four provenance lines (got %d)", length(ln)))
check(grepl("ps_am792_v124_refnum", txt),
      "the run directory appears, so a reader can tell which run this is")
check(grepl("REVAL_DIR", txt), "the mechanism that selected the run is named")
check(grepl("universe", txt, ignore.case = TRUE), "the universe is named")
check(grepl("[0-9a-f]{8}", txt), "a fingerprint or commit hash is present")
check(grepl("generated|written", txt, ignore.case = TRUE), "a generation timestamp is present")

cat("-- a stale fallback says so in the report, not only on stdout --\n")
# The whole point: the warning has to survive into the artefact.
for (v in VARS) Sys.unsetenv(v)
st <- paste(run_provenance_lines(), collapse = " | ")
check(grepl("default", st), "falling back to the default is stated")
check(grepl("STALE", st) || grepl("newer", st),
      "a fallback older than a run on disk is flagged in the stamp itself")

cat("-- the HTML form is a real tag, not pasted markup --\n")
Sys.setenv(REVAL_DIR = "analyses/ps_am792_v124_refnum/classifications")
h <- run_provenance_html()
check(inherits(h, "shiny.tag") || inherits(h, "shiny.tag.list"),
      "run_provenance_html returns an htmltools tag")
hs <- paste(as.character(h), collapse = " ")
check(grepl("ps_am792_v124_refnum", hs), "the rendered HTML carries the run")

cat("-- the markdown form is plain lines --\n")
m <- run_provenance_md()
check(is.character(m) && length(m) >= 4L, "run_provenance_md returns character lines")
check(any(grepl("ps_am792_v124_refnum", m)), "the markdown carries the run")

cat("-- the stamp can carry extra lines, and excluded scoresets are among them --\n")
# panel_b.html carried the post-exclusion numbers while mentioning KRAS zero times, so a
# reader could not tell that six of twenty scoresets were absent. The exclusion was
# announced on stdout, which does not survive the script exiting: the same gap the run
# stamp exists to close.
suppressMessages(source(resolve("analyses/lib/dms_polarity.R"), local = FALSE))
check(exists("dms_exclusion_lines"), "dms_exclusion_lines is available")
if (exists("dms_exclusion_lines")) {
  el <- load_dms_eligibility(resolve("analyses/derived/mavedb_scoreset_polarity.tsv"))
  xl <- dms_exclusion_lines(el)
  check(length(xl) >= 1L, "excluded scoresets produce at least one line")
  jx <- paste(xl, collapse = " ")
  check(grepl("00000115-a-1", jx), "the excluded scoresets are named")
  check(grepl("binding|ddG|free energy", jx, ignore.case = TRUE),
        "the reason travels with them, not just the count")
  check(length(dms_exclusion_lines(setNames(character(0), character(0)))) == 0L,
        "no exclusions produces no lines")
}
x <- run_provenance_lines(extra = c("excluded: something"))
check(any(grepl("excluded: something", x)), "run_provenance_lines accepts extra lines")
hx <- paste(as.character(run_provenance_html(extra = "excluded: something")), collapse = " ")
check(grepl("excluded: something", hx), "the HTML footer renders extra lines")
mx <- run_provenance_md(extra = "excluded: something")
check(any(grepl("excluded: something", mx)), "the markdown footer renders extra lines")

cat("-- every scoreset consumer discloses its exclusions in the report --\n")
for (f in c("analyses/07_panel_b_dms.R", "analyses/panel_b_dms_reval.R",
            "analyses/10_compose_supp_s2.R", "analyses/repro/07_panel_b_dms_run.R")) {
  p <- resolve(f)
  if (!file.exists(p)) { cat(sprintf("  [skip] %s absent\n", basename(f))); next }
  code <- sub("#.*$", "", readLines(p, warn = FALSE))
  check(any(grepl("dms_exclusion_lines\\(", code)),
        sprintf("%s puts its exclusions in the report, not only on stdout", basename(f)))
}

cat("-- every report writer stamps its output --\n")
WRITERS <- c("analyses/06_panel_a_clinical.R", "analyses/panel_a_fig_reval.R",
             "analyses/07_panel_b_dms.R", "analyses/panel_b_dms_reval.R",
             "analyses/08_panel_c_caps.R", "analyses/panel_c_caps_reval.R",
             "analyses/10_compose_supp_s2.R", "analyses/panel_b_auroc_summary.R",
             "analyses/repro/07_panel_b_dms_run.R")
for (f in WRITERS) {
  p <- resolve(f)
  if (!file.exists(p)) { cat(sprintf("  [skip] %s absent\n", basename(f))); next }
  code <- sub("#.*$", "", readLines(p, warn = FALSE))
  check(any(grepl("run_provenance_(html|md|lines)\\(", code)),
        sprintf("%s stamps its report", basename(f)))
  # Calling it is not enough: a script that stamps without sourcing the library
  # aborts on its last line, after the figures are already written.
  check(any(grepl("run_inputs\\.R", code)),
        sprintf("%s sources the library the stamp comes from", basename(f)))
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
