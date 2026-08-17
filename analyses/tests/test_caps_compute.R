library(testthat)

caps_path <- if (file.exists("analyses/lib/caps_compute.R")) {
  "analyses/lib/caps_compute.R"
} else {
  "../analyses/lib/caps_compute.R"
}
source(caps_path, local = TRUE)

test_that("wilson_ci returns p_hat=0.5 with bounds bracketing it", {
  ci <- wilson_ci(50, 100, conf = 0.95)
  expect_equal(ci$p_hat, 0.5)
  expect_lt(ci$lower, 0.5)
  expect_gt(ci$upper, 0.5)
  # approximately symmetric for moderate n
  expect_lt(abs((0.5 - ci$lower) - (ci$upper - 0.5)), 0.01)
})

test_that("wilson_ci handles k=0 and k=n edge cases", {
  ci_zero <- wilson_ci(0, 10, conf = 0.95)
  expect_equal(ci_zero$p_hat, 0)
  expect_equal(ci_zero$lower, 0)
  expect_gt(ci_zero$upper, 0)

  ci_full <- wilson_ci(10, 10, conf = 0.95)
  expect_equal(ci_full$p_hat, 1)
  expect_lt(ci_full$lower, 1)
  expect_equal(ci_full$upper, 1)
})

test_that("wilson_ci returns NA for n=0", {
  ci <- wilson_ci(0, 0, conf = 0.95)
  expect_true(is.na(ci$p_hat))
  expect_true(is.na(ci$lower))
  expect_true(is.na(ci$upper))
})

test_that("wilson_ci bounds within [0,1]", {
  for (n in c(5, 50, 500)) {
    for (k in c(0, ceiling(n/3), floor(2*n/3), n)) {
      ci <- wilson_ci(k, n, conf = 0.95)
      expect_true(ci$lower >= 0 && ci$lower <= 1)
      expect_true(ci$upper >= 0 && ci$upper <= 1)
      expect_true(ci$p_hat >= ci$lower && ci$p_hat <= ci$upper)
    }
  }
})

test_that("caps_score sign reflects observed vs expected", {
  expect_gt(caps_score(observed = 0.6, expected = 0.4), 0)
  expect_lt(caps_score(observed = 0.3, expected = 0.5), 0)
  expect_equal(caps_score(observed = 0.5, expected = 0.5), 0)
})

test_that("caps_score returns NA when expected is 0 or NA", {
  expect_true(is.na(caps_score(observed = 0.5, expected = 0)))
  expect_true(is.na(caps_score(observed = 0.5, expected = NA_real_)))
  expect_true(is.na(caps_score(observed = NA_real_, expected = 0.5)))
})

test_that("compute_caps_table returns one row per bin with required columns", {
  df <- data.frame(
    bin       = c("A", "A", "A", "B", "B", "C"),
    singleton = c(TRUE, FALSE, TRUE, FALSE, FALSE, TRUE)
  )
  expected_per_bin <- c(A = 0.5, B = 0.5, C = 0.5)
  out <- compute_caps_table(df, bin_col = "bin", singleton_col = "singleton",
                            expected_per_bin = expected_per_bin, conf = 0.95)
  expect_equal(nrow(out), 3)
  expect_setequal(
    colnames(out),
    c("bin", "n", "n_singletons", "p_observed", "p_expected",
      "caps", "ci_lower", "ci_upper")
  )
  expect_equal(sort(out$bin), c("A", "B", "C"))
  # n column matches input distribution: A=3, B=2, C=1
  out_sorted <- out[order(out$bin), ]
  expect_equal(out_sorted$n, c(3, 2, 1))
  expect_equal(out_sorted$n_singletons, c(2, 0, 1))
})

test_that("compute_caps_table CI bounds bracket p_observed", {
  df <- data.frame(
    bin       = rep("X", 100),
    singleton = c(rep(TRUE, 30), rep(FALSE, 70))
  )
  out <- compute_caps_table(df, bin_col = "bin", singleton_col = "singleton",
                            expected_per_bin = c(X = 0.2), conf = 0.95)
  expect_equal(out$n, 100)
  expect_equal(out$n_singletons, 30)
  expect_equal(out$p_observed, 0.3)
  expect_true(out$ci_lower < 0.3 && out$ci_upper > 0.3)
  expect_true(out$caps > 0)  # 0.3 > 0.2
})

test_that("compute_caps_table missing expected entry yields NA caps", {
  df <- data.frame(
    bin       = c("X", "X", "Y"),
    singleton = c(TRUE, FALSE, TRUE)
  )
  out <- compute_caps_table(df, bin_col = "bin", singleton_col = "singleton",
                            expected_per_bin = c(X = 0.5), conf = 0.95)
  out_y <- out[out$bin == "Y", ]
  expect_true(is.na(out_y$p_expected))
  expect_true(is.na(out_y$caps))
})
