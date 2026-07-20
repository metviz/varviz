# gene_calib_stats.R — pure Evidence-Strength / LR+ helpers.
# No network, no shiny. Sourced by server.R (like typography.R).
# Evidence Strength (LR+/OddsPath) is prior-free; prior enters only at acmg_posterior.

# OddsPath (LR+) -> Tavtigian Evidence-Strength tier. Boundaries = prior-free odds,
# C = 2.0813 = 350^(1/8), consistent with acmg_posterior in server.R.
oddspath_tier <- function(lr) {
  if (!is.finite(lr) || lr <= 1) return("None")
  if (lr >= 350)    return("Very strong")
  if (lr >= 18.7)   return("Strong")
  if (lr >= 4.33)   return("Moderate")
  if (lr >= 2.0813) return("Supporting")
  "None"
}

# LR+ at threshold t. pos/neg = numeric score vectors (P/LP arm, LB/B arm),
# higher score = more damaging. Haldane-Anscombe +0.5 on the 2x2 for lr/CI;
# raw (uncorrected) sens/spec returned for the sensitivity floor.
gene_lr_plus <- function(pos, neg, t) {
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  n_pos <- length(pos); n_neg <- length(neg)
  tp <- sum(pos >= t); fn <- n_pos - tp
  fp <- sum(neg >= t); tn <- n_neg - fp
  a <- tp + 0.5; b <- fn + 0.5; c <- fp + 0.5; d <- tn + 0.5   # Haldane cells
  sens_c <- a / (a + b); spec_c <- d / (c + d)
  lr <- sens_c / (1 - spec_c)
  var_ln <- (1 - sens_c)/(sens_c*(a + b)) + spec_c/((1 - spec_c)*(c + d))
  se <- sqrt(var_ln)
  list(lr = lr,
       lo = exp(log(lr) - 1.96*se),
       hi = exp(log(lr) + 1.96*se),
       sens = if (n_pos) tp / n_pos else NA_real_,
       spec = if (n_neg) tn / n_neg else NA_real_,
       tp = tp, fp = fp, n_pos = n_pos, n_neg = n_neg)
}

# gene-optimal threshold: max LR+ subject to raw sens >= min_sens.
# candidate thresholds = distinct pos scores (exact grid). NULL if none clears floor.
gene_optimal <- function(pos, neg, min_sens = 0.90) {
  pos <- pos[!is.na(pos)]; neg <- neg[!is.na(neg)]
  if (length(pos) == 0) return(NULL)
  best <- NULL
  for (t in sort(unique(pos))) {
    r <- gene_lr_plus(pos, neg, t)
    if (is.na(r$sens) || r$sens + 1e-9 < min_sens) next
    if (is.null(best) || r$lr > best$lr) { r$t <- t; best <- r }
  }
  best
}

# Predictor config: column name in the scored arm + Pejaver cut-points (tier -> threshold).
PREDICTORS <- list(
  REVEL = list(col = "revel", cuts = c(Supporting = 0.644, Moderate = 0.773, Strong = 0.932)),
  AM    = list(col = "am",    cuts = c(Supporting = 0.564)),
  CADD  = list(col = "cadd",  cuts = c(Supporting = 28.1, Moderate = 35))
)

# 3-letter -> 1-letter amino acid. Ter included for alt (nonsense), excluded as a
# missense ref (a p.Ter### start is not a missense variant).
.AA3TO1 <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",
             His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",
             Thr="T",Trp="W",Tyr="Y",Val="V",Ter="*")

# ClinVar `name` -> dbNSFP hgvsp like "p.R130G"; NA unless a clean missense p.Xxx###Yyy.
parse_clinvar_hgvsp <- function(name) {
  vapply(name, function(nm) {
    if (is.na(nm)) return(NA_character_)
    g <- regmatches(nm, regexec("p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3})", nm))[[1]]
    if (length(g) != 4) return(NA_character_)
    ref <- .AA3TO1[[g[2]]]; alt <- .AA3TO1[[g[4]]]
    if (is.null(ref) || is.null(alt) || ref == "*") return(NA_character_)
    paste0("p.", ref, g[3], alt)
  }, character(1), USE.NAMES = FALSE)
}
