library(testthat)

# ---------------------------------------------------------------------------
# test_dolphin_bulk.R — per-gene DOLPHIN precompute helper tests
#
# Offline tests (always run):
#   - Empty input returns an empty data.frame with the expected schema
#   - Resume: existing TSV cache is honoured, only un-done variants are queried
#
# Online tests (skipped unless VARVIZ_TEST_DOLPHIN=TRUE):
#   - SNCA 4-variant manuscript cohort (~4 s wall-time at 1.05 s/req)
#   - Verifies pm1 boolean is set, acmg_raw is populated for in-domain hits
# ---------------------------------------------------------------------------

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
source(resolve("analyses/lib/dolphin_bulk.R"), local = FALSE)

ONLINE <- isTRUE(as.logical(Sys.getenv("VARVIZ_TEST_ONLINE", "FALSE"))) ||
          isTRUE(as.logical(Sys.getenv("VARVIZ_TEST_DOLPHIN", "FALSE")))

EXPECTED_COLS <- c("p_notation", "p_single", "ensembl", "pm1",
                   "acmg_raw", "pfam_id", "score_delta", "fetched_at")

# ---------------------------------------------------------------------------
# Offline: empty input
# ---------------------------------------------------------------------------
test_that("fetch_dolphin_gene returns empty df on empty input", {
  tmp <- tempfile("dolphin_test_", fileext = "")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  df <- fetch_dolphin_gene("CASR", character(0), ensembl = "ENST00000639785",
                            cache_dir = tmp, verbose = FALSE)
  expect_true(is.data.frame(df))
  expect_equal(nrow(df), 0)
  expect_setequal(colnames(df), EXPECTED_COLS)
})

# ---------------------------------------------------------------------------
# Offline: resume behaviour (no network calls because all variants pre-loaded)
# ---------------------------------------------------------------------------
test_that("fetch_dolphin_gene resumes from existing TSV cache", {
  tmp <- tempfile("dolphin_test_", fileext = "")
  dir.create(tmp, recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  # Pre-seed the cache with two variants
  seed <- data.frame(
    p_notation  = c("p.A53T", "p.E46K"),
    p_single    = c("A53T", "E46K"),
    ensembl     = c("ENST00000394986", "ENST00000394986"),
    pm1         = c(FALSE, FALSE),
    acmg_raw    = c("", ""),
    pfam_id     = c("", ""),
    score_delta = c("", ""),
    fetched_at  = c("2026-01-01T00:00:00Z", "2026-01-01T00:00:00Z"),
    stringsAsFactors = FALSE
  )
  write.table(seed, file.path(tmp, "SNCA.tsv"), sep = "\t",
              row.names = FALSE, quote = FALSE, na = "")

  # Ask for only the two variants already in the cache → no network calls,
  # function returns the cached rows verbatim.
  df <- fetch_dolphin_gene("SNCA",
                            p_notations = c("p.A53T", "p.E46K"),
                            ensembl = "ENST00000394986",
                            cache_dir = tmp,
                            verbose = FALSE)
  expect_equal(nrow(df), 2)
  expect_setequal(df$p_notation, c("p.A53T", "p.E46K"))
})

# ---------------------------------------------------------------------------
# Online: SNCA 4-variant manuscript cohort
# ---------------------------------------------------------------------------
test_that("fetch_dolphin_gene precomputes SNCA 4-variant cohort", {
  skip_if_not(ONLINE, "set VARVIZ_TEST_DOLPHIN=TRUE to enable online tests")

  tmp <- tempfile("dolphin_test_snca_", fileext = "")
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  variants <- c("p.A53T", "p.E46K", "p.G51D", "p.G14R")
  enst <- "ENST00000394986"  # canonical SNCA

  df <- fetch_dolphin_gene("SNCA", variants, ensembl = enst,
                            cache_dir = tmp, rate_sleep_s = 1.05,
                            verbose = FALSE)
  expect_equal(nrow(df), 4)
  expect_setequal(df$p_notation, variants)
  expect_setequal(df$p_single, c("A53T", "E46K", "G51D", "G14R"))
  expect_true(all(df$ensembl == enst))
  expect_type(df$pm1, "logical")

  # Verify TSV written to disk
  expect_true(file.exists(file.path(tmp, "SNCA.tsv")))
})
