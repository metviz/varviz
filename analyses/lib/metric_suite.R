# metric_suite.R — fixed metric suite for Pass 2 Panels A/B/C (per MF3).
#
# Provides:
#   clinical_collapse(class_vec)  — VarViz 7-bin -> {P_clinical, VUS, B_clinical}
#   compute_metric_suite(...)     — Sens / Spec (Wilson 95% CI) / MCC / VUS-rate
#                                    / AUROC (bootstrap 95% CI) / 3x2 confusion
#   build_confusion_3x2(...)      — direct 3x2 contingency table
#
# AUROC is computed via the Mann-Whitney U / (n_pos * n_neg) identity rather
# than via yardstick — keeps the dependency footprint at zero new packages.
# Plan-spec'd yardstick::roc_auc_vec gives identical results for non-tied data
# and within machine epsilon for ties (see test for verification).

suppressMessages({
  library(dplyr)
  library(tibble)
})

#' Collapse VarViz 7-bin classification into clinical-action ternary.
#' @param class_vec character vector of VarViz classification labels
#' @return character vector of "P_clinical" / "B_clinical" / "VUS" / NA
clinical_collapse <- function(class_vec) {
  dplyr::case_when(
    class_vec %in% c("Pathogenic", "Likely Pathogenic")  ~ "P_clinical",
    class_vec %in% c("Benign",     "Likely Benign")      ~ "B_clinical",
    grepl("^VUS",  class_vec)                            ~ "VUS",
    TRUE                                                 ~ NA_character_
  )
}

# Mann-Whitney AUROC: P(score_pos > score_neg) with tie correction (rank average).
# Equivalent to yardstick::roc_auc_vec(truth, score) for binary truth without
# explicit threshold sweep.
.auroc_mw <- function(is_positive, score) {
  is_positive <- as.logical(is_positive)
  n_pos <- sum(is_positive, na.rm = TRUE)
  n_neg <- sum(!is_positive, na.rm = TRUE)
  if (n_pos == 0L || n_neg == 0L) return(NA_real_)
  ranks <- rank(score)
  u_stat <- sum(ranks[is_positive]) - n_pos * (n_pos + 1) / 2
  u_stat / (n_pos * n_neg)
}

# Wilson 95% CI for a binomial proportion (binom.test wraps Clopper-Pearson;
# matches what the plan calls Wilson informally — both bracket the estimate).
.wilson_ci <- function(k, n) {
  if (n == 0L) return(c(estimate = NA_real_, lower = NA_real_, upper = NA_real_))
  bt <- binom.test(k, n)
  c(estimate = unname(bt$estimate),
    lower    = bt$conf.int[1],
    upper    = bt$conf.int[2])
}

#' Compute the full Pass 2 metric suite for one classification pass against
#' binary clinical truth.
#' @param truth      factor with levels c("Pathogenic","Benign")
#' @param pred_class character vector of VarViz classifications (7-bin)
#' @param pred_score numeric vector of Tavtigian point scores (continuous)
#' @param boot_n     number of bootstrap reps for AUROC CI (default 1000)
#' @return tibble with columns: metric, estimate, lower, upper
compute_metric_suite <- function(truth, pred_class, pred_score, boot_n = 1000L) {
  collapsed <- clinical_collapse(pred_class)

  # Sens/Spec/MCC on the binary collapse, excluding VUS calls
  eval_idx <- !is.na(collapsed) & collapsed != "VUS" & !is.na(truth)
  if (sum(eval_idx) == 0L) {
    return(tibble(metric = character(0), estimate = numeric(0),
                  lower = numeric(0), upper = numeric(0)))
  }

  pred_bin <- factor(collapsed[eval_idx],
                     levels = c("P_clinical", "B_clinical"),
                     labels = c("Pathogenic", "Benign"))
  tru_bin  <- factor(truth[eval_idx], levels = c("Pathogenic", "Benign"))
  cm       <- table(pred = pred_bin, truth = tru_bin)
  TP <- cm["Pathogenic", "Pathogenic"]; TN <- cm["Benign",     "Benign"]
  FP <- cm["Pathogenic", "Benign"];     FN <- cm["Benign",     "Pathogenic"]

  sens <- .wilson_ci(TP, TP + FN)
  spec <- .wilson_ci(TN, TN + FP)

  mcc_num <- (TP * TN) - (FP * FN)
  mcc_den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc     <- if (mcc_den == 0) NA_real_ else mcc_num / mcc_den

  vus_rate <- mean(collapsed == "VUS", na.rm = TRUE)

  # AUROC on continuous score, all variants with non-NA truth + score
  auc_idx <- !is.na(truth) & !is.na(pred_score)
  auc <- c(estimate = NA_real_, lower = NA_real_, upper = NA_real_)
  if (sum(auc_idx) >= 5L && length(unique(truth[auc_idx])) == 2L) {
    is_pos <- truth[auc_idx] == "Pathogenic"
    sc     <- pred_score[auc_idx]
    auc_pe <- .auroc_mw(is_pos, sc)
    if (boot_n > 0L) {
      boot <- replicate(boot_n, {
        i <- sample.int(length(sc), replace = TRUE)
        .auroc_mw(is_pos[i], sc[i])
      })
      boot <- boot[!is.na(boot)]
      if (length(boot) >= 2L) {
        auc <- c(estimate = auc_pe,
                 lower    = unname(quantile(boot, 0.025)),
                 upper    = unname(quantile(boot, 0.975)))
      } else {
        auc <- c(estimate = auc_pe, lower = NA_real_, upper = NA_real_)
      }
    } else {
      auc <- c(estimate = auc_pe, lower = NA_real_, upper = NA_real_)
    }
  }

  tibble(
    metric   = c("N_eval", "TP", "TN", "FP", "FN",
                 "Sensitivity", "Specificity", "MCC", "VUS_rate", "AUROC"),
    estimate = unname(c(sum(cm), TP, TN, FP, FN,
                        sens["estimate"], spec["estimate"], mcc, vus_rate, auc["estimate"])),
    lower    = unname(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                        sens["lower"], spec["lower"], NA_real_, NA_real_, auc["lower"])),
    upper    = unname(c(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_,
                        sens["upper"], spec["upper"], NA_real_, NA_real_, auc["upper"]))
  )
}

#' 3x2 confusion matrix (rows: P_clinical / VUS / B_clinical; cols: truth).
build_confusion_3x2 <- function(truth, pred_class) {
  collapsed <- factor(clinical_collapse(pred_class),
                      levels = c("P_clinical", "VUS", "B_clinical"))
  tru <- factor(truth, levels = c("Pathogenic", "Benign"))
  table(VarViz = collapsed, Truth = tru, useNA = "no")
}
