# Guards how an analysis script picks the run whose classifications it reads.
#
# Why this exists. analyses/derived/varviz_classifications.tsv was six months
# stale, dated 2026-05-02, while the authoritative run was from 2026-09-28. Its row
# count still matched exactly, 90,701, which is what made it dangerous: the usual
# sanity check passes. Content differed substantially.
#
# Six scripts read that path with no way to point them at a run. They do print what
# they read, so the failure was not silence, it was that nothing said the file was
# older than every run on disk. Three of them write the published Panel A, Panel B
# and Panel C reports.
#
# The resolver honours an explicit override, accepts a run directory, and warns when
# the file it fell back to predates a run that exists. It keeps the default, so every
# documented invocation still works; what changes is that staleness announces itself.
#
# analyses/repro/03_which_run.R is deliberately NOT a caller. Its whole purpose is to
# score several candidate runs against the manuscript's claims, so the consolidated
# file being one hardcoded candidate among five is the point, not a defect.
#
# Run from project root: Rscript analyses/tests/test_run_inputs.R

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
clear <- function() for (v in VARS) Sys.unsetenv(v)
on.exit({
  clear()
  for (v in VARS) if (!is.na(saved[[v]])) do.call(Sys.setenv, setNames(list(saved[[v]]), v))
}, add = TRUE)

tmp <- file.path(tempdir(), paste0("runin_", as.integer(runif(1, 1e6, 9e6))))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

# A run directory in the real shape: <gene>__dual.tsv per gene.
ck <- file.path(tmp, "run_a", "classifications")
dir.create(ck, recursive = TRUE, showWarnings = FALSE)
mk <- function(gene, n) {
  d <- data.frame(gene = rep(gene, n), p_notation = paste0("p.A", seq_len(n), "G"),
                  varviz_classification_full = rep("VUS-Low", n),
                  stringsAsFactors = FALSE)
  write.table(d, file.path(ck, paste0(gene, "__dual.tsv")), sep = "\t",
              quote = FALSE, row.names = FALSE)
}
mk("GENEA", 3); mk("GENEB", 4)
explicit <- file.path(tmp, "explicit.tsv")
write.table(data.frame(gene = "GENEZ", p_notation = "p.M1V"), explicit,
            sep = "\t", quote = FALSE, row.names = FALSE)
dflt <- file.path(tmp, "default.tsv")
write.table(data.frame(gene = "GENED", p_notation = "p.M1L"), dflt,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("-- the helpers exist --\n")
check(exists("resolve_classifications"), "resolve_classifications is available")
check(exists("resolve_checkpoint_dir"),  "resolve_checkpoint_dir is available")

cat("-- an explicit file wins --\n")
clear(); Sys.setenv(VARVIZ_CLASSIFICATIONS = explicit)
r <- resolve_classifications(default = dflt, quiet = TRUE)
check(identical(normalizePath(r$path), normalizePath(explicit)),
      "VARVIZ_CLASSIFICATIONS is used when set")
check(identical(r$source, "VARVIZ_CLASSIFICATIONS"), "the source is reported as the override")

cat("-- a run directory is consolidated to a readable path --\n")
# Callers read a single path, so a run directory has to collapse to one. Doing it
# here rather than in five scripts is the point.
clear(); Sys.setenv(REVAL_DIR = ck)
r <- resolve_classifications(default = dflt, quiet = TRUE)
check(file.exists(r$path), "REVAL_DIR yields an existing path")
check(identical(r$source, "REVAL_DIR"), "the source is reported as REVAL_DIR")
d <- utils::read.delim(r$path, stringsAsFactors = FALSE)
check(nrow(d) == 7L, sprintf("all 7 checkpoint rows are consolidated (got %d)", nrow(d)))
check(setequal(unique(d$gene), c("GENEA", "GENEB")), "both genes survive the consolidation")

cat("-- a run's summary.tsv locates its checkpoints --\n")
clear(); Sys.setenv(VARVIZ_RUN = file.path(tmp, "run_a", "summary.tsv"))
check(identical(normalizePath(resolve_checkpoint_dir(default = tmp)), normalizePath(ck)),
      "VARVIZ_RUN resolves the sibling classifications directory")
r <- resolve_classifications(default = dflt, quiet = TRUE)
check(identical(r$source, "VARVIZ_RUN"), "VARVIZ_RUN is used when no stronger override is set")

cat("-- precedence is explicit file, then run directory, then run summary --\n")
clear()
Sys.setenv(VARVIZ_CLASSIFICATIONS = explicit, REVAL_DIR = ck,
           VARVIZ_RUN = file.path(tmp, "run_a", "summary.tsv"))
check(identical(resolve_classifications(default = dflt, quiet = TRUE)$source,
                "VARVIZ_CLASSIFICATIONS"), "the explicit file outranks both others")
clear(); Sys.setenv(REVAL_DIR = ck, VARVIZ_RUN = file.path(tmp, "run_a", "summary.tsv"))
check(identical(resolve_classifications(default = dflt, quiet = TRUE)$source, "REVAL_DIR"),
      "a run directory outranks a run summary")

cat("-- with nothing set the default is kept, so documented invocations still work --\n")
clear()
r <- resolve_classifications(default = dflt, quiet = TRUE)
check(identical(normalizePath(r$path), normalizePath(dflt)), "the default path is returned")
check(identical(r$source, "default"), "the source is reported as the default")

cat("-- falling back to a file older than a run on disk is reported as stale --\n")
# The whole point. A newer run existing while the default is used is the exact
# situation that produced supplementary numbers from a May run.
newer <- file.path(tmp, "ps_newer")
dir.create(file.path(newer, "classifications"), recursive = TRUE, showWarnings = FALSE)
writeLines("x", file.path(newer, "summary.tsv"))
Sys.setmtime <- function(p, t) Sys.setFileTime(p, t)
Sys.setFileTime(dflt, Sys.time() - 86400 * 30)
Sys.setFileTime(file.path(newer, "summary.tsv"), Sys.time())
clear()
r <- resolve_classifications(default = dflt, quiet = TRUE, run_glob = file.path(tmp, "ps_*"))
check(isTRUE(r$stale), "the fallback is flagged stale when a newer run exists")
check(length(r$newer_runs) >= 1L, "the newer run is named, not just counted")
msg <- paste(capture.output(
  resolve_classifications(default = dflt, quiet = FALSE,
                          run_glob = file.path(tmp, "ps_*")), type = "message"),
  collapse = " ")
check(grepl("STALE", msg), "a loud STALE warning is emitted when not quiet")

cat("-- a fallback newer than every run is not flagged --\n")
Sys.setFileTime(dflt, Sys.time())
r <- resolve_classifications(default = dflt, quiet = TRUE, run_glob = file.path(tmp, "ps_*"))
check(!isTRUE(r$stale), "a current default is not called stale")

cat("-- every panel and dolphin script goes through the resolver --\n")
for (f in c("analyses/06_panel_a_clinical.R", "analyses/07_panel_b_dms.R",
            "analyses/08_panel_c_caps.R", "analyses/11_merge_dolphin.R",
            "analyses/run_dolphin_benchmark.R")) {
  p <- resolve(f)
  if (!file.exists(p)) { cat(sprintf("  [skip] %s absent\n", basename(f))); next }
  code <- sub("#.*$", "", readLines(p, warn = FALSE))
  check(any(grepl("run_inputs\\.R", code)), sprintf("%s sources run_inputs.R", basename(f)))
  check(any(grepl("resolve_classifications\\(", code)),
        sprintf("%s resolves its classifications input", basename(f)))
}

cat("-- the run-selection script stays hardcoded on purpose --\n")
p <- resolve("analyses/repro/03_which_run.R")
if (file.exists(p)) {
  code <- sub("#.*$", "", readLines(p, warn = FALSE))
  check(!any(grepl("resolve_classifications\\(", code)),
        "03_which_run.R does not resolve: its candidate list IS the question it answers")
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
