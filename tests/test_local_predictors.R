library(testthat)
library(data.table)

# testthat::test_file() sets cwd to tests/, plain Rscript leaves it at project
# root. Resolve every path relative to whichever cwd we're under.
resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src_path  <- resolve("analyses/lib/local_predictors.R")
am_dir    <- resolve("analyses/raw/alphamissense")
revel_dir <- resolve("analyses/raw/revel/by_chrom")

source(src_path, local = TRUE)

# ---------------------------------------------------------------------------
# AlphaMissense (small per-UID CSVs, ~8 KB — always run)
# ---------------------------------------------------------------------------

test_that("load_alphamissense_csv returns keyed data.table for PTEN (P60484)", {
  am <- load_alphamissense_csv("P60484", am_dir = am_dir)
  expect_s3_class(am, "data.table")
  expect_setequal(names(am), c("protein_variant", "am_pathogenicity", "am_class"))
  expect_equal(key(am), "protein_variant")
  expect_gt(nrow(am), 1000L)
})

test_that("load_alphamissense_csv returns NULL when the UID has no bulk file", {
  expect_null(load_alphamissense_csv("P16473", am_dir = am_dir))  # TSHR — no AM file
})

test_that("lookup_alphamissense finds a known PTEN variant in 'p.X#Y' form", {
  am <- load_alphamissense_csv("P60484", am_dir = am_dir)
  out <- lookup_alphamissense(am, "p.M1A")
  expect_equal(out$score, 0.7696)
  expect_equal(out$class, "LPath")
})

test_that("lookup_alphamissense accepts 'X#Y' (no 'p.' prefix)", {
  am <- load_alphamissense_csv("P60484", am_dir = am_dir)
  out <- lookup_alphamissense(am, "M1A")
  expect_equal(out$score, 0.7696)
})

test_that("lookup_alphamissense returns NA list for an absent variant", {
  am <- load_alphamissense_csv("P60484", am_dir = am_dir)
  out <- lookup_alphamissense(am, "p.Z999Z")
  expect_true(is.na(out$score))
  expect_true(is.na(out$class))
})

test_that("lookup_alphamissense gracefully handles a NULL am_dt (missing UID)", {
  out <- lookup_alphamissense(NULL, "p.M1A")
  expect_true(is.na(out$score))
  expect_true(is.na(out$class))
})

# ---------------------------------------------------------------------------
# REVEL (per-chrom CSV ~150 MB — skip if bulk data not present, e.g. on CI)
# ---------------------------------------------------------------------------

revel_chr17_path <- file.path(revel_dir, "chr_17.csv")

test_that("load_revel_for_chrom returns keyed data.table for chr 17", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  expect_s3_class(rev, "data.table")
  expect_true(all(c("grch38_pos", "ref", "alt", "REVEL") %in% names(rev)))
  expect_equal(key(rev), c("grch38_pos", "ref", "alt"))
  expect_gt(nrow(rev), 1e6L)
})

test_that("load_revel_for_chrom collapses (pos, ref, alt) duplicates across transcripts", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  dups <- rev[, .N, by = .(grch38_pos, ref, alt)][N > 1L]
  expect_equal(nrow(dups), 0L)
})

test_that("lookup_revel finds a known chr17 variant", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  out <- lookup_revel(rev, "17", 156220L, "G", "A")
  expect_equal(out, 0.073)
  expect_type(out, "double")
})

test_that("lookup_revel returns NA_real_ for an absent (pos, ref, alt)", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  out <- lookup_revel(rev, "17", 1L, "A", "T")
  expect_true(is.na(out))
  expect_type(out, "double")
})

test_that("lookup_revel errors when called with the wrong chromosome", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  expect_error(lookup_revel(rev, "10", 156220L, "G", "A"), "chromosome")
})
