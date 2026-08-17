# Self-check for acmg_posterior (§4.1). Run: Rscript test_acmg_posterior.R
# ponytail: extracts the REAL acmg_posterior out of server.R rather than mirroring
# it, so this fails if C, prior_p, or the formula drifts. server.R is a trusted
# in-repo build artifact (never user input), which makes eval()'ing a fixed range
# safe here — same pattern as test_pp4.R / test_clinvar_star_weighting.R.
src <- readLines("server.R", warn = FALSE)
beg <- grep("^acmg_posterior <- function", src)
end <- beg - 1 + which(src[beg:length(src)] == "}")[1]
eval(parse(text = paste(src[beg:end], collapse = "\n")))

stopifnot(
  abs(acmg_posterior(0)  - 0.10) < 1e-9,   # 0 pts -> posterior == prior
  abs(acmg_posterior(6)  - 0.90) < 5e-3,   # LP cutpoint ~ 0.90 by construction
  abs(acmg_posterior(10) - 0.99) < 1e-2,   # P threshold ~ 0.99
  acmg_posterior(-4) < 0.10,               # benign evidence pulls below prior
  is.na(acmg_posterior(NA)),
  is.na(acmg_posterior(c(1, 2)))           # non-scalar rejected
)
cat("acmg_posterior: all checks pass\n")
