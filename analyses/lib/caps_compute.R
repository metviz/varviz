# caps_compute.R — CAPS (Consequence-Adjusted Proportion of Singletons)
# computation library, per Gudkov et al. 2025 (Bioinformatics Advances).
#
# Provides:
#   wilson_ci(k, n, conf)               — Wilson score CI for proportion k/n
#   caps_score(observed, expected)      — (observed - expected) / expected
#   compute_caps_table(df, bin_col, singleton_col, expected_per_bin, conf)
#                                        — returns per-bin n, p_observed,
#                                          p_expected, caps, CI bounds
#
# Used by Pass 2 Panel C (analyses/08_panel_c_caps.R) to test whether VarViz's
# pathogenicity bins are consistent with population-genetic singleton-fraction
# expectations on gnomAD-present variants. Panel C provides an orthogonal
# (population-genetic) third truth track alongside Panel A (clinical) and
# Panel B (functional).

suppressMessages({
  library(dplyr)
})

#' Wilson score CI for a binomial proportion k/n.
#' @param k integer count of successes
#' @param n integer total
#' @param conf confidence level (default 0.95)
#' @return list(p_hat, lower, upper); list of NA's if n == 0
wilson_ci <- function(k, n, conf = 0.95) {
  if (n == 0L) return(list(p_hat = NA_real_, lower = NA_real_, upper = NA_real_))
  z      <- qnorm(1 - (1 - conf) / 2)
  p_hat  <- k / n
  denom  <- 1 + z^2 / n
  centre <- (p_hat + z^2 / (2 * n)) / denom
  half   <- (z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))) / denom
  # Clamp into [0, 1] AND ensure lower <= p_hat <= upper. Without the explicit
  # min/max with p_hat, floating-point in centre/half yields lower ~= 3e-17 at
  # k=0 (instead of exactly 0), breaking the bracketing invariant.
  list(
    p_hat = p_hat,
    lower = min(p_hat, max(0, centre - half)),
    upper = max(p_hat, min(1, centre + half))
  )
}

#' CAPS score: (observed - expected) / expected.
#' Positive = enrichment of singletons (suggests purifying selection);
#' negative = depletion (suggests neutrality / age of allele).
#' @return scalar; NA if expected is 0 / NA or observed is NA
caps_score <- function(observed, expected) {
  if (is.na(observed) || is.na(expected) || expected == 0) return(NA_real_)
  (observed - expected) / expected
}

#' Per-bin CAPS table with Wilson CI on observed singleton proportion.
#' @param df data.frame with one row per variant
#' @param bin_col column name (string) containing the bin label per variant
#' @param singleton_col column name (string) containing the per-variant singleton flag (logical)
#' @param expected_per_bin named vector mapping bin label -> expected singleton fraction
#' @param conf confidence level for Wilson CI (default 0.95)
#' @return data.frame: bin, n, n_singletons, p_observed, p_expected, caps, ci_lower, ci_upper
compute_caps_table <- function(df, bin_col, singleton_col, expected_per_bin, conf = 0.95) {
  bins <- as.character(df[[bin_col]])
  sing <- as.logical(df[[singleton_col]])

  rows <- lapply(unique(bins), function(b) {
    mask <- bins == b
    n <- sum(mask)
    k <- sum(sing[mask], na.rm = TRUE)
    ci <- wilson_ci(k, n, conf)
    p_obs <- ci$p_hat
    p_exp <- if (!is.null(expected_per_bin) && b %in% names(expected_per_bin)) {
      expected_per_bin[[b]]
    } else NA_real_
    data.frame(
      bin          = b,
      n            = n,
      n_singletons = k,
      p_observed   = p_obs,
      p_expected   = p_exp,
      caps         = caps_score(p_obs, p_exp),
      ci_lower     = ci$lower,
      ci_upper     = ci$upper,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}
