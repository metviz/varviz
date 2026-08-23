#!/usr/bin/env Rscript
# Recompute every manuscript-reported benchmark statistic from the derived runs.
# Base R only -- no package dependencies, so a reviewer can run this on a bare install.
#
# Input : analyses/derived/varviz_classifications_mds.tsv  (90,701 rows; carries
#         Pass-Full, Pass-Blind, Pass-Blind+DOLPHIN and Pass-Blind+MDS in one file)
# Output: analyses/repro/out/computed.tsv  (id, value)
#
# Run from the repository root:  Rscript analyses/repro/01_recompute.R

args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "analyses/derived/varviz_classifications_mds.tsv"
OUT <- if (length(args) >= 2) args[2] else "analyses/repro/out/computed.tsv"

stopifnot(file.exists(IN))
d <- read.delim(IN, sep = "\t", quote = "", stringsAsFactors = FALSE,
                colClasses = "character", check.names = FALSE)
n <- nrow(d)
cat(sprintf("[repro] loaded %s: %d rows, %d cols\n", IN, n, ncol(d)))

BINS <- c("Pathogenic", "Likely Pathogenic", "VUS-High", "VUS-Mid", "VUS-Low",
          "Likely Benign", "Benign")
VUS  <- c("VUS-High", "VUS-Mid", "VUS-Low")

full  <- d$varviz_classification_full
blind <- d$varviz_classification_blind
has_aug <- all(c("varviz_classification_blind_mds","varviz_classification_blind_dolphin",
                 "mds_pm1","dolphin_pm1") %in% names(d))
bmds  <- if (has_aug) d$varviz_classification_blind_mds else NULL
bdol  <- if (has_aug) d$varviz_classification_blind_dolphin else NULL

res <- list()
put <- function(id, v) res[[id]] <<- v

# ---- distributions over the full universe (Supp S4.2 / Table S3) -------------
pct <- function(x, set) 100 * sum(x %in% set) / n
put("dist_full_plp",   pct(full,  c("Pathogenic", "Likely Pathogenic")))
put("dist_blind_plp",  pct(blind, c("Pathogenic", "Likely Pathogenic")))
put("dist_full_vush",  pct(full,  "VUS-High"))
put("dist_blind_vush", pct(blind, "VUS-High"))
put("dist_full_lbb",   pct(full,  c("Likely Benign", "Benign")))
put("dist_blind_lbb",  pct(blind, c("Likely Benign", "Benign")))
if (has_aug) {
  put("dist_bmds_plp",  pct(bmds, c("Pathogenic", "Likely Pathogenic")))
  put("dist_bmds_vush", pct(bmds, "VUS-High"))
}

# ---- Pass-Full -> Pass-Blind transitions (Supp S4.3 / Table S4) -------------
shifted <- full != blind
put("shift_n",   sum(shifted))
put("shift_pct", 100 * sum(shifted) / n)

tr <- function(a, b) sum(full == a & blind == b)
put("trans_lp_vush_n", tr("Likely Pathogenic", "VUS-High"))
put("trans_p_lp_n",    tr("Pathogenic", "Likely Pathogenic"))
put("trans_lp_vush_pct", 100 * tr("Likely Pathogenic", "VUS-High") /
                            sum(full == "Likely Pathogenic"))
put("trans_p_lp_pct",    100 * tr("Pathogenic", "Likely Pathogenic") /
                            sum(full == "Pathogenic"))
put("shift_lp_vush_share", 100 * tr("Likely Pathogenic", "VUS-High") / sum(shifted))
put("shift_p_lp_share",    100 * tr("Pathogenic", "Likely Pathogenic") / sum(shifted))

if (has_aug) {
# ---- MDS augmentation (MS 3.1 / Supp S3.2) ----------------------------------
# "augmented" = MDS contributed a PM1 tag on the Pass-Blind row.
mds_tag <- trimws(d$mds_pm1)
mds_fires <- nzchar(mds_tag) & mds_tag != "NA" & tolower(mds_tag) != "false"
has_pm1 <- grepl("PM1", d$tags_blind, fixed = TRUE)
put("mds_fires_n",   sum(mds_fires))
put("mds_fires_pct", 100 * sum(mds_fires) / n)
put("mds_aug_n",     sum(mds_fires & !has_pm1))
put("mds_aug_pct",   100 * sum(mds_fires & !has_pm1) / n)
put("mds_scorable_n",   sum(!is.na(suppressWarnings(as.numeric(d$mds)))))
put("mds_scorable_pct", 100 * sum(!is.na(suppressWarnings(as.numeric(d$mds)))) / n)
# "resolved" = was VUS under Pass-Blind, is no longer VUS once MDS is applied.
put("mds_resolved", sum(blind %in% VUS & !(bmds %in% VUS)))
put("mds_changed",  sum(blind != bmds))

# ---- DOLPHIN comparison on the same baseline (Supp S3.2) --------------------
dol_tag <- trimws(d$dolphin_pm1)
dol_fires <- nzchar(dol_tag) & dol_tag != "NA" & tolower(dol_tag) != "false"
put("dol_fires_n", sum(dol_fires))
put("dol_aug_n",   sum(dol_fires & !has_pm1))   # same definition as mds_aug_n
put("dol_cross_n", sum(blind != bdol))
put("dol_vush_lp",   sum(blind == "VUS-High" & bdol == "Likely Pathogenic"))
put("dol_vusm_vush", sum(blind == "VUS-Mid"  & bdol == "VUS-High"))
put("dol_vusl_vusm", sum(blind == "VUS-Low"  & bdol == "VUS-Mid"))
put("dol_resolved",  sum(blind %in% VUS & !(bdol %in% VUS)))
put("mds_vs_dol_ratio", res$mds_resolved / max(res$dol_resolved, 1))

# ---- McNemar: MDS vs DOLPHIN on resolving a Pass-Blind VUS ------------------
vus_rows <- blind %in% VUS
a <- !(bmds[vus_rows] %in% VUS)   # MDS resolved it
b <- !(bdol[vus_rows] %in% VUS)   # DOLPHIN resolved it
b01 <- sum(a & !b); b10 <- sum(!a & b)
put("mcnemar_b01", b01)
put("mcnemar_b10", b10)
put("mcnemar_chi2", if ((b01 + b10) > 0) (abs(b01 - b10))^2 / (b01 + b10) else NA)
# Effect size the reviewers asked for in place of an uninformative p-value.
put("mds_resolve_rate_pct", 100 * sum(a) / sum(vus_rows))
put("dol_resolve_rate_pct", 100 * sum(b) / sum(vus_rows))
put("resolve_rate_diff_pp", 100 * (sum(a) - sum(b)) / sum(vus_rows))

}  # end has_aug

# ---- PM1_strong carried into Pass-Blind (reviewer D2 ceiling) ---------------
put("blind_pm1_strong_n",   sum(grepl("PM1_strong", d$tags_blind, fixed = TRUE)))
put("blind_pm1_strong_pct", 100 * sum(grepl("PM1_strong", d$tags_blind, fixed = TRUE)) / n)

dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
out <- data.frame(id = names(res), value = unlist(res, use.names = FALSE))
write.table(out, OUT, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("[repro] wrote %d computed values -> %s\n", nrow(out), OUT))
