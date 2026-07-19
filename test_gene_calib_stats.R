# test_gene_calib_stats.R  — run: Rscript test_gene_calib_stats.R
source("gene_calib_stats.R")

# oddspath_tier boundaries (reuse §4.1 C)
stopifnot(
  oddspath_tier(11.83) == "Moderate",   # the Bindra REVEL>=0.7 value
  oddspath_tier(2.5)   == "Supporting",
  oddspath_tier(20)    == "Strong",
  oddspath_tier(400)   == "Very strong",
  oddspath_tier(1.5)   == "None"
)

# perfect separation -> finite LR+ via Haldane, not Inf
r <- gene_lr_plus(pos = c(0.9,0.95,0.99), neg = c(0.1,0.2,0.05), t = 0.5)
stopifnot(is.finite(r$lr), r$lr > 1, r$sens == 1, r$spec == 1)

# hand 2x2: pos>=t: 8/10 (sens .8), neg>=t: 2/10 (spec .8) -> LR+ approx .8/.2 = 4 (Haldane shifts slightly)
pos <- c(rep(1,8), rep(0,2)); neg <- c(rep(1,2), rep(0,8))
r2 <- gene_lr_plus(pos, neg, t = 1)
stopifnot(abs(r2$lr - 4) < 0.61, r2$lo < r2$lr, r2$hi > r2$lr, r2$tp == 8, r2$fp == 2)

# gene_optimal respects the sensitivity floor
pos3 <- c(0.6,0.7,0.8,0.9,0.95); neg3 <- c(0.1,0.2,0.3,0.65,0.72)
opt_hi <- gene_optimal(pos3, neg3, min_sens = 0.90)   # must keep >=90% of pos -> low t
opt_lo <- gene_optimal(pos3, neg3, min_sens = 0.50)   # may pick higher t, higher LR+
stopifnot(!is.null(opt_hi), opt_hi$sens >= 0.90 - 1e-9,
          opt_lo$lr >= opt_hi$lr)

# floor impossible -> NULL
stopifnot(is.null(gene_optimal(pos3, neg3, min_sens = 1.01)))

cat("gene_calib_stats: all checks pass\n")
