# Guards the per-gene run manifest.
#
# Why this exists. A run recorded its inputs and its analysis parameters but not
# the code that consumed them: no version, no commit, and none of the engine
# flags. That is the gap that made ps_final_v231 impossible to account for. Its
# caches and its code were provably unchanged, yet 16 KCNQ1 variants scored
# differently, and nothing stored with the run could say which configuration had
# produced it. Identical inputs plus a different flag is a different call, so
# inputs alone are not provenance.
#
# The manifest is deliberately written per gene rather than once per run. Shards
# run concurrently against one output directory, so a single shared file would be
# a write race; per-gene files are disjoint, which is the same reason the
# classification checkpoints are per gene.
#
# Run from project root: Rscript analyses/tests/test_run_manifest.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src <- readLines(resolve("analyses/ps_final_harness.R"), warn = FALSE)

fails <- 0L
check <- function(ok, msg) {
  ok <- isTRUE(ok)
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

beg <- grep("^write_run_manifest <- function", src)
check(length(beg) == 1L,
      sprintf("write_run_manifest is defined exactly once (found %d)", length(beg)))
if (length(beg) != 1L) {
  cat(sprintf("\nFAIL: %d failure(s)\n", fails + 1L)); quit(status = 1L)
}
end <- beg - 1L + which(src[beg:length(src)] == "}")[1]
# readLines/parse/eval of a fixed line range of this repository's own tracked
# source, the same idiom as analyses/tests/test_acmg_posterior.R. Reading it from
# the harness is the point: the assertions must land on the shipped function, not
# on a copy of it in this file.
eval(parse(text = paste(src[beg:end], collapse = "\n")), envir = globalenv())

tmp <- file.path(tempdir(), paste0("manifest_", as.integer(runif(1, 1e6, 9e6))))
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

# A flag whose value must appear, to prove the manifest reads live state rather
# than restating a hardcoded list.
options(varviz.am_pp3_supporting = 0.792, varviz.am_hybrid = TRUE,
        varviz.am_calibrated = FALSE, varviz.mds_pm1 = TRUE)

p <- write_run_manifest("CASR", tmp, universe = resolve("analyses/derived/variant_universe_gnomad.tsv"))
check(file.exists(p), "a manifest file is written")
check(basename(p) == "CASR__manifest.json", "it is named <gene>__manifest.json")

txt <- paste(readLines(p, warn = FALSE), collapse = "\n")
cat("-- code identity --\n")
for (k in c("varviz_version", "git_commit", "git_dirty", "r_version", "written_at")) {
  check(grepl(paste0('"', k, '"'), txt, fixed = TRUE), sprintf("records %s", k))
}
cat("-- engine flags, read from live state --\n")
for (k in c("varviz.am_pp3_supporting", "varviz.am_hybrid", "varviz.am_calibrated",
            "varviz.mds_pm1", "varviz.am_cap_pathogenic", "varviz.pm2_points")) {
  check(grepl(k, txt, fixed = TRUE), sprintf("records %s", k))
}
check(grepl("0.792", txt, fixed = TRUE),
      "the recorded supporting threshold is the live value, not a default")

cat("-- input identity --\n")
for (k in c("universe", "universe_sha256", "gene")) {
  check(grepl(paste0('"', k, '"'), txt, fixed = TRUE), sprintf("records %s", k))
}
# A manifest that names an input without fingerprinting it cannot show the input
# was the same one.
sha <- sub('^.*"universe_sha256"\\s*:\\s*"([0-9a-f]*)".*$', "\\1", gsub("\n", " ", txt))
check(nchar(sha) == 64L, "universe_sha256 is a full sha256, not a placeholder")

cat("-- it is valid JSON --\n")
ok <- tryCatch({ jsonlite::fromJSON(p); TRUE }, error = function(e) FALSE)
check(ok, "the manifest parses as JSON")

cat("-- the harness writes one per gene beside the master table --\n")
body <- paste(src, collapse = "\n")
check(grepl("write_run_manifest(", body, fixed = TRUE),
      "the harness calls write_run_manifest")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
