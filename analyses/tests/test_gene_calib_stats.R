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

# gene_calibration: assemble per-predictor validation + optimal, with LOO
mk <- function(hg, rv) data.frame(hgvsp = hg, revel = rv,
                                  am = rv, cadd = rv*40, stringsAsFactors = FALSE)
path_arm   <- mk(sprintf("p.A%dV", 1:30), c(rep(0.95, 25), rep(0.60, 5)))
benign_arm <- mk(sprintf("p.G%dS", 1:30), c(rep(0.10, 27), rep(0.80, 3)))

res <- gene_calibration(path_arm, benign_arm, min_sens = 0.80)
stopifnot(
  res$loo_applied == FALSE,
  res$predictors$REVEL$n_pos == 30, res$predictors$REVEL$n_neg == 30,
  res$predictors$REVEL$confident == TRUE,                       # both arms >= 20
  nrow(res$predictors$REVEL$validation) == 3,                   # 3 REVEL cut-points
  !is.null(res$predictors$REVEL$optimal),
  res$predictors$REVEL$optimal$lr > 1
)

# LOO: dropping a P/LP member reduces n_pos by 1
res2 <- gene_calibration(path_arm, benign_arm, query_hgvsp = "p.A1V", min_sens = 0.80)
stopifnot(res2$loo_applied == TRUE, res2$predictors$REVEL$n_pos == 29)

# small benign arm -> not confident, still computes
res3 <- gene_calibration(path_arm, benign_arm[1:8, ], min_sens = 0.80)
stopifnot(res3$predictors$REVEL$confident == FALSE,
          !is.null(res3$predictors$REVEL$validation))

# hgvsp parse: ClinVar name string -> dbNSFP hgvsp (1-letter)
stopifnot(
  parse_clinvar_hgvsp("NM_000314.8(PTEN):c.388C>G (p.Arg130Gly)") == "p.R130G",
  parse_clinvar_hgvsp("NM_000257.4(MYH7):c.1988G>A (p.Arg663His)") == "p.R663H",
  is.na(parse_clinvar_hgvsp("NM_000314.8(PTEN):c.209+1G>A")),        # splice, no p.
  identical(
    parse_clinvar_hgvsp(c("x (p.Gly12Asp)", "no-protein", "y (p.Ter100Cys)")),
    c("p.G12D", NA_character_, NA_character_)),                       # Ter ref not a missense start
  is.na(parse_clinvar_hgvsp("NM_x(GENE):c.394del (p.Lys132del)")),          # in-frame del -> NA, no crash
  is.na(parse_clinvar_hgvsp("NM_x(TP53):c.388del (p.Arg130Glyfs*45)")),     # frameshift -> NA, not p.R130G
  is.na(parse_clinvar_hgvsp("x (p.Xaa130Gly)")),                            # unknown ref token -> NA
  identical(parse_clinvar_hgvsp(c("a (p.Gly12Asp)", "b (p.Lys132del)")),
            c("p.G12D", NA_character_)),                                     # mixed vector: valid survives, no throw
  parse_clinvar_hgvsp("x (p.Arg130Ter)") == "p.R130*"                       # nonsense ALT still maps
)

cat("gene_calib_stats: all checks pass\n")
