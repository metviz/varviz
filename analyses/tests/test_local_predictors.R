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

# ---------------------------------------------------------------------------
# lookup_revel_by_aa — codon-positions + alt_aa filter for missense lookups
# ---------------------------------------------------------------------------

# Fixture: chr17:156220 G>A produces aaref=T, aaalt=M with REVEL=0.073.
# The codon spans positions 156220..156222 on the plus strand at that locus.

test_that("lookup_revel_by_aa finds REVEL for codon positions + alt_aa", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  out <- lookup_revel_by_aa(rev, "17", c(156220L, 156221L, 156222L), "M")
  expect_equal(out, 0.073)
  expect_type(out, "double")
})

test_that("lookup_revel_by_aa returns NA when no row matches alt_aa", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  out <- lookup_revel_by_aa(rev, "17", c(156220L, 156221L, 156222L), "Z")
  expect_true(is.na(out))
  expect_type(out, "double")
})

test_that("lookup_revel_by_aa returns NA for empty / all-NA codon_positions", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  expect_true(is.na(lookup_revel_by_aa(rev, "17", integer(0), "M")))
  expect_true(is.na(lookup_revel_by_aa(rev, "17", c(NA_integer_, NA_integer_), "M")))
})

test_that("lookup_revel_by_aa errors when called with the wrong chromosome", {
  skip_if_not(file.exists(revel_chr17_path), "REVEL bulk data not on disk")
  rev <- load_revel_for_chrom("17", revel_dir = revel_dir)
  expect_error(
    lookup_revel_by_aa(rev, "10", c(156220L, 156221L, 156222L), "M"),
    "chromosome"
  )
})

# ---------------------------------------------------------------------------
# dbNSFP 4.9a single-source lookup (~40 GB bgzipped — skip if not built)
#
# Fixture: BAP1 p.Q280K — chr3:52405858 G>T (NM_004656.4 / NC_000003.12).
# Expected per-tool values cross-checked against the live VarViz app's
# TSV download for the same variant (rs746045584). Because dbNSFP and
# MyVariant.info both source from dbNSFP, the values should match to
# the precision dbNSFP stores them.
# ---------------------------------------------------------------------------

dbnsfp_path <- resolve("analyses/tmp/dbNSFPv4.9a_hg38_custombuild.bgz")

test_that("load_dbnsfp_for_region returns NULL for nonexistent file", {
  expect_null(load_dbnsfp_for_region("3", 1L, 100L, dbnsfp_path = "/nonexistent.bgz"))
})

test_that("load_dbnsfp_for_region pulls BAP1 region and hashes by (pos,ref,alt)", {
  skip_if_not(file.exists(dbnsfp_path), "dbNSFP build not on disk yet")
  env <- load_dbnsfp_for_region("3", 52401000L, 52410000L, dbnsfp_path = dbnsfp_path)
  expect_true(is.environment(env))
  expect_gt(length(ls(env)), 0L)
  expect_equal(attr(env, "dbnsfp_chrom"), "3")
})

test_that("lookup_dbnsfp_local returns parse_dbnsfp_scores-shaped list for BAP1 p.Q280K", {
  skip_if_not(file.exists(dbnsfp_path), "dbNSFP build not on disk yet")
  env <- load_dbnsfp_for_region("3", 52401000L, 52410000L, dbnsfp_path = dbnsfp_path)
  hit <- lookup_dbnsfp_local(env, "3", 52405858L, "G", "T")
  expect_type(hit, "list")
  expect_named(hit, "dbnsfp")
  d <- hit$dbnsfp
  # Values from dbNSFP 4.9a (Aug 2024). The live app's TSV uses MyVariant.info
  # which caches an older dbNSFP; CADD shifted from 17.77 (old) -> 21.9 (4.9a).
  # REVEL/AM/PhyloP/PhastCons match the live app to dbNSFP precision.
  expect_equal(d$revel$score,                    0.190, tolerance = 0.005)
  expect_equal(d$alphamissense$am_pathogenicity, 0.1204, tolerance = 0.005)
  expect_equal(d$alphamissense$am_class,         "B")
  expect_equal(d$cadd$phred,                     21.9,  tolerance = 0.05)
  expect_equal(d$dann$score,                     0.978, tolerance = 0.005)
  # Four predictors that were NA under the REVEL+AM+UCSC monkey-patch are
  # now populated — this is the apples-to-apples-with-live-app win.
  expect_equal(d$sift$score,                     0.086, tolerance = 0.005)
  expect_equal(d$polyphen2$hdiv$score,           0.008, tolerance = 0.005)
  expect_equal(d$metasvm$score,                  -0.894, tolerance = 0.005)
  expect_equal(d$`gerp++`$rs,                    5.35,  tolerance = 0.05)
  expect_equal(d$phylop$`100way_vertebrate`$score,    5.454, tolerance = 0.005)
  expect_equal(d$phylop$`470way_mammalian`$score,    11.827, tolerance = 0.005)
  expect_equal(d$phastcons$`100way_vertebrate`$score, 1.000, tolerance = 0.005)
})

test_that("lookup_dbnsfp_local returns NULL for absent (pos,ref,alt)", {
  skip_if_not(file.exists(dbnsfp_path), "dbNSFP build not on disk yet")
  env <- load_dbnsfp_for_region("3", 52401000L, 52410000L, dbnsfp_path = dbnsfp_path)
  expect_null(lookup_dbnsfp_local(env, "3", 52405858L, "G", "Z"))  # bogus alt
  expect_null(lookup_dbnsfp_local(env, "3", 1L, "A", "T"))         # outside region
})

test_that("lookup_dbnsfp_local rejects cross-chrom queries silently", {
  skip_if_not(file.exists(dbnsfp_path), "dbNSFP build not on disk yet")
  env <- load_dbnsfp_for_region("3", 52401000L, 52410000L, dbnsfp_path = dbnsfp_path)
  expect_null(lookup_dbnsfp_local(env, "17", 52405858L, "G", "T"))
})

test_that("lookup_dbnsfp_by_aa resolves BAP1 p.Q280K via protein key", {
  skip_if_not(file.exists(dbnsfp_path), "dbNSFP build not on disk yet")
  env <- load_dbnsfp_for_region("3", 52401000L, 52410000L, dbnsfp_path = dbnsfp_path)
  hit <- lookup_dbnsfp_by_aa(env, "3", 280L, "K")
  expect_type(hit, "list")
  expect_named(hit, "dbnsfp")
  # Same fixture values — protein-key path should find the same row
  expect_equal(hit$dbnsfp$revel$score,                   0.19, tolerance = 0.02)
  expect_equal(hit$dbnsfp$alphamissense$am_pathogenicity, 0.12, tolerance = 0.02)
})

test_that("lookup_dbnsfp_by_aa returns NULL for missing aa change", {
  skip_if_not(file.exists(dbnsfp_path), "dbNSFP build not on disk yet")
  env <- load_dbnsfp_for_region("3", 52401000L, 52410000L, dbnsfp_path = dbnsfp_path)
  expect_null(lookup_dbnsfp_by_aa(env, "3", 280L, "Z"))
})
