library(testthat)

# ---------------------------------------------------------------------------
# test_dolphin.R — DOLPHIN PM1 pathway helper tests
#
# Offline tests (always run):
#   - 3-letter -> 1-letter amino-acid conversion
#   - dolphin_fires_pm1() parsing on mock responses (results[].acmg schema)
#   - Sentinel "results"-as-string handling (variant outside Pfam domain)
#
# Online tests (skipped unless VARVIZ_TEST_DOLPHIN=TRUE):
#   - Canonical Pfam hotspot positive: FBN1 G560V in EGF_CA domain (PM1+)
#   - Canonical out-of-domain negative: CASR p.Asn1074Asp (C-term disordered)
#   - Gene -> canonical ENST resolution via Ensembl REST
#   - Cache reuse on repeated identical queries
# ---------------------------------------------------------------------------

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src_path <- resolve("analyses/lib/dolphin.R")
source(src_path, local = FALSE)

ONLINE <- isTRUE(as.logical(Sys.getenv("VARVIZ_TEST_ONLINE", "FALSE"))) ||
          isTRUE(as.logical(Sys.getenv("VARVIZ_TEST_DOLPHIN", "FALSE")))

# ---------------------------------------------------------------------------
# Offline: amino-acid notation conversion
# ---------------------------------------------------------------------------
test_that(".dolphin_to_single_letter handles 3-letter form", {
  expect_equal(.dolphin_to_single_letter("p.Ser296Asn"), "S296N")
  expect_equal(.dolphin_to_single_letter("Ser296Asn"),    "S296N")
  expect_equal(.dolphin_to_single_letter("p.Ala150Thr"),  "A150T")
})

test_that(".dolphin_to_single_letter passes through 1-letter form", {
  expect_equal(.dolphin_to_single_letter("p.S296N"), "S296N")
  expect_equal(.dolphin_to_single_letter("S296N"),    "S296N")
  expect_equal(.dolphin_to_single_letter("A150T"),    "A150T")
})

test_that(".dolphin_to_single_letter handles stop codon (Ter / *)", {
  expect_equal(.dolphin_to_single_letter("p.Arg100Ter"), "R100*")
  expect_equal(.dolphin_to_single_letter("R100*"),       "R100*")
})

# ---------------------------------------------------------------------------
# Offline: dolphin_fires_pm1() parsing on mock responses.
#
# Real DOLPHIN response schema:
#   { variation: "G560V", ensembl: "...", results: [ { acmg: "PM1;", ... } ] }
# Or, when out of domain:
#   { variation: "...", results: "This mutation is not within a protein domain ..." }
# ---------------------------------------------------------------------------
test_that("dolphin_fires_pm1 returns FALSE for NULL / empty input", {
  expect_false(dolphin_fires_pm1(NULL))
  expect_false(dolphin_fires_pm1(list()))
  expect_false(dolphin_fires_pm1(list(results = NULL)))
  expect_false(dolphin_fires_pm1(list(results = list())))
})

test_that("dolphin_fires_pm1 returns FALSE on 'not within domain' sentinel", {
  resp <- list(results = "This mutation is not within a protein domain or this mutation doesn't exist in our database")
  expect_false(dolphin_fires_pm1(resp))
})

test_that("dolphin_fires_pm1 detects PM1 in results[].acmg field", {
  expect_true(dolphin_fires_pm1(list(results = list(list(acmg = "PM1;")))))
  expect_true(dolphin_fires_pm1(list(results = list(list(acmg = "PM1; PP3")))))
  expect_true(dolphin_fires_pm1(list(results = list(list(acmg = "PP3; PM1")))))
  expect_true(dolphin_fires_pm1(list(results = list(list(acmg = "PM1")))))
})

test_that("dolphin_fires_pm1 does NOT fire on PM2 (DOLPHIN's freq-based tag)", {
  # We must not consume DOLPHIN's PM2: it would conflict with VarViz's
  # gnomAD-derived PM2. Only the literal PM1 tag is acceptable.
  expect_false(dolphin_fires_pm1(list(results = list(list(acmg = "PM2;")))))
})

test_that("dolphin_fires_pm1 does not match PM10, PM11, BP1, PS1, etc. as PM1", {
  expect_false(dolphin_fires_pm1(list(results = list(list(acmg = "PM10;")))))
  expect_false(dolphin_fires_pm1(list(results = list(list(acmg = "PM11;")))))
  expect_false(dolphin_fires_pm1(list(results = list(list(acmg = "BP1;")))))
  expect_false(dolphin_fires_pm1(list(results = list(list(acmg = "PS1;")))))
})

test_that("dolphin_fires_pm1 walks multi-domain results arrays", {
  multi <- list(results = list(
    list(acmg = "PM2;"),
    list(acmg = "PM1;")  # PM1 found in second domain hit
  ))
  expect_true(dolphin_fires_pm1(multi))

  none <- list(results = list(
    list(acmg = "PM2;"),
    list(acmg = "BP4;")
  ))
  expect_false(dolphin_fires_pm1(none))
})

test_that("dolphin_fires_pm1 handles list-of-strings acmg field", {
  expect_true(dolphin_fires_pm1(list(results = list(list(acmg = list("PM1", "PP3"))))))
  expect_false(dolphin_fires_pm1(list(results = list(list(acmg = list("BP1", "BP4"))))))
})

test_that("dolphin_fires_pm1 ignores results entries without acmg", {
  resp <- list(results = list(list(name_dom = "EGF_CA")))
  expect_false(dolphin_fires_pm1(resp))
})

# ---------------------------------------------------------------------------
# Offline: fetch_dolphin argument validation
# ---------------------------------------------------------------------------
test_that("fetch_dolphin errors when neither gene nor ensembl provided", {
  expect_error(fetch_dolphin(p_notation = "S296N"),
               "provide either `ensembl` or `gene`")
})

test_that("fetch_dolphin returns NULL on empty p_notation", {
  expect_null(fetch_dolphin(ensembl = "ENST00000316623", p_notation = ""))
})

# ---------------------------------------------------------------------------
# Online: real DOLPHIN API calls (skipped unless VARVIZ_TEST_DOLPHIN=TRUE)
# ---------------------------------------------------------------------------
test_that("fetch_dolphin returns parsed JSON for FBN1 G560V (paper example)", {
  skip_if_not(ONLINE, "set VARVIZ_TEST_DOLPHIN=TRUE to enable online tests")
  resp <- fetch_dolphin(ensembl = "ENST00000316623", p_notation = "G560V")
  expect_true(is.list(resp))
  expect_equal(resp$variation, "G560V")
  expect_true(is.list(resp$results))
  expect_true(length(resp$results) >= 1)
})

test_that("dolphin_fires_pm1 fires for FBN1 G560V in EGF_CA domain", {
  skip_if_not(ONLINE, "set VARVIZ_TEST_DOLPHIN=TRUE to enable online tests")
  # FBN1 G560V is the canonical PM1+ example in the DOLPHIN apidoc.
  expect_true(dolphin_pm1_call(ensembl = "ENST00000316623", p_notation = "G560V"))
})

test_that("dolphin_fires_pm1 does NOT fire for CASR p.Asn1074Asp (C-term disordered)", {
  skip_if_not(ONLINE, "set VARVIZ_TEST_DOLPHIN=TRUE to enable online tests")
  # Disordered tail: DOLPHIN returns "results" as a string sentinel, no PM1.
  expect_false(dolphin_pm1_call(ensembl = "ENST00000639785",
                                 p_notation = "p.Asn1074Asp"))
})

test_that("dolphin_canonical_enst resolves CASR symbol", {
  skip_if_not(ONLINE, "set VARVIZ_TEST_DOLPHIN=TRUE to enable online tests")
  enst <- dolphin_canonical_enst("CASR")
  expect_true(is.character(enst))
  expect_true(startsWith(enst, "ENST"))
})

test_that("fetch_dolphin caches repeated calls", {
  skip_if_not(ONLINE, "set VARVIZ_TEST_DOLPHIN=TRUE to enable online tests")
  rm(list = ls(envir = .dolphin_local_cache), envir = .dolphin_local_cache)

  t1 <- system.time(r1 <- fetch_dolphin(ensembl = "ENST00000316623",
                                         p_notation = "G560V"))["elapsed"]
  t2 <- system.time(r2 <- fetch_dolphin(ensembl = "ENST00000316623",
                                         p_notation = "G560V"))["elapsed"]

  expect_identical(r1, r2)
  expect_true(t2 < 0.05 || t2 < t1)
})
