library(testthat)

src_path <- if (file.exists("analyses/lib/metric_suite.R")) {
  "analyses/lib/metric_suite.R"
} else {
  "../analyses/lib/metric_suite.R"
}
source(src_path, local = TRUE)

# -----------------------------------------------------------------------------
# clinical_collapse
# -----------------------------------------------------------------------------
test_that("clinical_collapse maps {P, LP} -> P_clinical and {B, LB} -> B_clinical", {
  expect_equal(clinical_collapse("Pathogenic"),         "P_clinical")
  expect_equal(clinical_collapse("Likely Pathogenic"),  "P_clinical")
  expect_equal(clinical_collapse("Benign"),             "B_clinical")
  expect_equal(clinical_collapse("Likely Benign"),      "B_clinical")
})

test_that("clinical_collapse maps all VUS-* tiers to VUS", {
  expect_equal(clinical_collapse("VUS-High"), "VUS")
  expect_equal(clinical_collapse("VUS-Mid"),  "VUS")
  expect_equal(clinical_collapse("VUS-Low"),  "VUS")
  expect_equal(clinical_collapse("VUS"),      "VUS")
})

test_that("clinical_collapse returns NA for unknown classifications", {
  expect_true(is.na(clinical_collapse(NA_character_)))
  expect_true(is.na(clinical_collapse("Conflicting")))
})

# -----------------------------------------------------------------------------
# compute_metric_suite
# -----------------------------------------------------------------------------
test_that("perfect classifier yields sensitivity=1, specificity=1, MCC=1, AUROC=1", {
  truth <- factor(c(rep("Pathogenic", 10), rep("Benign", 10)),
                  levels = c("Pathogenic", "Benign"))
  pred  <- c(rep("Pathogenic", 10), rep("Benign", 10))
  score <- c(rep(10, 10), rep(-10, 10))
  res <- compute_metric_suite(truth, pred, score)
  expect_equal(res$estimate[res$metric == "Sensitivity"], 1)
  expect_equal(res$estimate[res$metric == "Specificity"], 1)
  expect_equal(res$estimate[res$metric == "MCC"], 1)
  expect_equal(res$estimate[res$metric == "AUROC"], 1)
  expect_equal(res$estimate[res$metric == "VUS_rate"], 0)
})

test_that("anti-classifier (always wrong) yields MCC=-1, AUROC=0", {
  truth <- factor(c(rep("Pathogenic", 10), rep("Benign", 10)),
                  levels = c("Pathogenic", "Benign"))
  pred  <- c(rep("Benign", 10), rep("Pathogenic", 10))
  score <- c(rep(-10, 10), rep(10, 10))
  res <- compute_metric_suite(truth, pred, score)
  expect_equal(res$estimate[res$metric == "MCC"], -1)
  expect_equal(res$estimate[res$metric == "AUROC"], 0)
})

test_that("VUS_rate counts all VUS-* tiers", {
  truth <- factor(c("Pathogenic","Pathogenic","Benign","Benign"),
                  levels = c("Pathogenic","Benign"))
  pred  <- c("VUS-High","VUS-Mid","VUS-Low","Benign")
  score <- c(2, 1, 0, -5)
  res <- compute_metric_suite(truth, pred, score)
  expect_equal(res$estimate[res$metric == "VUS_rate"], 0.75)
})

test_that("LP and LB get folded into P/B halves of the binary collapse", {
  truth <- factor(c(rep("Pathogenic", 4), rep("Benign", 4)),
                  levels = c("Pathogenic", "Benign"))
  pred  <- c("Pathogenic", "Likely Pathogenic", "Likely Pathogenic", "Pathogenic",
             "Benign", "Likely Benign", "Likely Benign", "Benign")
  score <- c(8, 6, 6, 8, -8, -4, -4, -8)
  res <- compute_metric_suite(truth, pred, score)
  expect_equal(res$estimate[res$metric == "Sensitivity"], 1)
  expect_equal(res$estimate[res$metric == "Specificity"], 1)
})

test_that("Wilson CI: lower<=estimate<=upper, both within [0,1]", {
  truth <- factor(c(rep("Pathogenic", 8), rep("Benign", 8), "Pathogenic", "Benign"),
                  levels = c("Pathogenic","Benign"))
  pred  <- c(rep("Pathogenic", 8), rep("Benign", 8), "Benign", "Pathogenic")
  score <- c(rep(8, 8), rep(-8, 8), -1, 1)
  res <- compute_metric_suite(truth, pred, score)
  for (m in c("Sensitivity", "Specificity")) {
    row <- res[res$metric == m, ]
    expect_true(row$lower    >= 0 && row$lower    <= 1)
    expect_true(row$upper    >= 0 && row$upper    <= 1)
    expect_true(row$estimate >= row$lower && row$estimate <= row$upper)
  }
})

test_that("AUROC bootstrap CI brackets the point estimate", {
  set.seed(42)
  n <- 60
  truth <- factor(rep(c("Pathogenic","Benign"), each = n/2),
                  levels = c("Pathogenic","Benign"))
  score <- c(rnorm(n/2, mean = 1), rnorm(n/2, mean = -1))
  pred  <- ifelse(score > 1, "Pathogenic",
           ifelse(score > 0, "VUS-High",
           ifelse(score > -1, "VUS-Low", "Benign")))
  res <- compute_metric_suite(truth, pred, score)
  auc_row <- res[res$metric == "AUROC", ]
  expect_true(auc_row$estimate > 0.7)
  expect_true(auc_row$lower <= auc_row$estimate && auc_row$estimate <= auc_row$upper)
  expect_true(auc_row$lower >= 0 && auc_row$upper <= 1)
})

# -----------------------------------------------------------------------------
# build_confusion_3x2
# -----------------------------------------------------------------------------
test_that("confusion matrix has 3 rows (P/VUS/B) and 2 cols (truth)", {
  truth <- c("Pathogenic","Pathogenic","Benign","Benign")
  pred  <- c("Pathogenic","VUS-Mid","Benign","VUS-High")
  cm <- build_confusion_3x2(truth, pred)
  expect_equal(nrow(cm), 3)
  expect_equal(ncol(cm), 2)
  expect_equal(rownames(cm), c("P_clinical","VUS","B_clinical"))
  expect_equal(colnames(cm), c("Pathogenic","Benign"))
})

test_that("confusion matrix totals match input", {
  truth <- c(rep("Pathogenic",10), rep("Benign",10))
  pred  <- c(rep("Pathogenic", 7), rep("VUS-Mid", 3),
             rep("Benign", 6), rep("VUS-Low", 4))
  cm <- build_confusion_3x2(truth, pred)
  expect_equal(sum(cm), 20)
  expect_equal(sum(cm[, "Pathogenic"]), 10)
  expect_equal(sum(cm[, "Benign"]), 10)
})
