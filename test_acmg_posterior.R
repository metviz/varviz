# Self-check for acmg_posterior (§4.1). Run: Rscript test_acmg_posterior.R
# ponytail: mirror of the fn in server.R; keep C and prior_p in sync with it.
acmg_posterior <- function(points, prior_p = 0.10, C = 2.0813) {
  if (length(points) != 1 || !is.finite(points)) return(NA_real_)
  prior_odds <- prior_p / (1 - prior_p)
  post_odds  <- prior_odds * C^points
  post_odds / (1 + post_odds)
}

stopifnot(
  abs(acmg_posterior(0)  - 0.10) < 1e-9,   # 0 pts -> posterior == prior
  abs(acmg_posterior(6)  - 0.90) < 5e-3,   # LP cutpoint ~ 0.90 by construction
  abs(acmg_posterior(10) - 0.99) < 1e-2,   # P threshold ~ 0.99
  acmg_posterior(-4) < 0.10,               # benign evidence pulls below prior
  is.na(acmg_posterior(NA)),
  is.na(acmg_posterior(c(1, 2)))           # non-scalar rejected
)
cat("acmg_posterior: all checks pass\n")
