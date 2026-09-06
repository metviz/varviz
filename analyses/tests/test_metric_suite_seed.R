# Self-check for the seeded AUROC bootstrap in analyses/lib/metric_suite.R.
#
# The 95% CI on AUROC comes from a sample.int() bootstrap. Without a seed the
# Pass-Blind CIs quoted in the manuscript moved from run to run. Now
# compute_metric_suite(seed = 1L) is deterministic by default, seed = NULL
# restores the old behaviour, and the caller's RNG state is left untouched
# either way. Plain R, no testthat. Run from project root:
#   Rscript analyses/tests/test_metric_suite_seed.R

source("analyses/lib/metric_suite.R")

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok" else "FAIL", msg))
  if (!isTRUE(ok)) fails <<- fails + 1L
}

set.seed(2024)
n <- 80
truth <- factor(rep(c("Pathogenic", "Benign"), each = n / 2), levels = c("Pathogenic", "Benign"))
score <- c(rnorm(n / 2, 1), rnorm(n / 2, -1))
pred  <- ifelse(score > 1, "Pathogenic", ifelse(score > 0, "VUS-High",
         ifelse(score > -1, "VUS-Low", "Benign")))
auc_ci <- function(...) {
  m <- compute_metric_suite(truth, pred, score, boot_n = 300L, ...)
  unlist(m[m$metric == "AUROC", c("lower", "upper")])
}

# 1. default seed -> identical CI on repeat
a <- auc_ci(); b <- auc_ci()
check(identical(a, b), "default seed: two runs give the same AUROC CI")

# 2. a different seed gives a different CI (bootstrap really re-drew)
c <- auc_ci(seed = 99L)
check(!identical(a, c), "different seed: CI changes")

# 3. seed = NULL is non-deterministic (old behaviour, opt-in)
set.seed(1); d1 <- auc_ci(seed = NULL)
set.seed(2); d2 <- auc_ci(seed = NULL)
check(!identical(d1, d2), "seed = NULL: bootstrap follows the ambient RNG")

# 4. the caller's RNG stream is not disturbed by the internal set.seed
set.seed(7); before <- runif(3)
set.seed(7); invisible(auc_ci()); after <- runif(3)
check(identical(before, after), "caller's RNG state restored after the seeded bootstrap")

# 5. point estimate is unaffected by the seed at all
m1 <- compute_metric_suite(truth, pred, score, boot_n = 300L, seed = 1L)
m2 <- compute_metric_suite(truth, pred, score, boot_n = 300L, seed = 5L)
check(identical(m1$estimate[m1$metric == "AUROC"], m2$estimate[m2$metric == "AUROC"]),
      "AUROC point estimate independent of seed")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
