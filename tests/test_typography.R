library(testthat)

# Sources only typography.R (not the heavy server.R) — the typography tier
# system lives in its own small file so unit tests don't pull the full app
# library/data stack. Run from project root: `Rscript -e "testthat::test_file('tests/test_typography.R')"`
typo_path <- if (file.exists("typography.R")) "typography.R" else "../typography.R"
source(typo_path, local = TRUE)

test_that("font tier globals derive from base_font with correct multipliers", {
  expect_true(exists("vv_base_font"))
  expect_true(exists("vv_big"))
  expect_true(exists("vv_medium"))
  expect_true(exists("vv_small"))

  expect_equal(vv_big,    vv_base_font * 1.00)
  expect_equal(vv_medium, vv_base_font * 0.75)
  expect_equal(vv_small,  vv_base_font * 0.55)
})

test_that("default base_font is 12 pt for live app", {
  expect_equal(vv_base_font, 12)
})

