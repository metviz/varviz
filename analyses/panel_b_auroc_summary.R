#!/usr/bin/env Rscript
# Per-scoreset AUROC for BOTH passes on the MaveDB DMS panel.
#
# panel_b_dms_reval.R computes per-study metrics for Pass-Full only, so the
# Pass-Blind median and range quoted in Supplementary §S2.3 had no generator in
# the repository and could not be re-derived. This supplies both, from the same
# per-study quartile binning the panel uses, so the quoted summary is
# reproducible from a run directory.
#
# Usage:
#   REVAL_DIR=analyses/ps_final_v231/classifications \
#     Rscript analyses/panel_b_auroc_summary.R
suppressMessages({library(dplyr); library(readr); library(tidyr)})
# The panels' own AUROC (Mann-Whitney U, ties at 0.5), not pROC: the repository
# does not depend on pROC and the panel figures are computed this way.
source("analyses/lib/metric_suite.R")

REVAL_DIR <- Sys.getenv("REVAL_DIR", "analyses/ps_final_v231/classifications")
MAVEDB    <- Sys.getenv("MAVEDB", "analyses/derived/mavedb_canonical.tsv")

files <- list.files(REVAL_DIR, pattern = "__dual\\.tsv$", full.names = TRUE)
if (!length(files)) stop("no __dual.tsv files under ", REVAL_DIR)
cls <- bind_rows(lapply(files, read_tsv, show_col_types = FALSE)) |>
  distinct(gene, p_notation, .keep_all = TRUE)
mave <- read_tsv(MAVEDB, show_col_types = FALSE)

score_col <- intersect(c("score_raw", "score"), names(mave))[1]
if (is.na(score_col)) stop("no score column in ", MAVEDB)

joined <- inner_join(mave, cls, by = c("gene", "p_notation"))
joined$score_dms <- joined[[score_col]]

# Per-study quartiles, identical to the panel: bottom 25% LOF, top 25% WT-like.
joined <- joined |>
  group_by(study) |>
  mutate(q25 = quantile(score_dms, .25, na.rm = TRUE),
         q75 = quantile(score_dms, .75, na.rm = TRUE),
         truth = case_when(score_dms <= q25 ~ 1, score_dms >= q75 ~ 0,
                           TRUE ~ NA_real_)) |>
  ungroup() |>
  filter(!is.na(truth))

auroc <- function(truth, score) {
  score <- suppressWarnings(as.numeric(score))
  ok <- !is.na(truth) & !is.na(score)
  truth <- truth[ok]; score <- score[ok]
  if (length(unique(truth)) < 2L) return(NA_real_)
  pos <- score[truth == 1]; neg <- score[truth == 0]
  r <- rank(c(pos, neg))
  (sum(r[seq_along(pos)]) - length(pos) * (length(pos) + 1) / 2) /
    (length(pos) * length(neg))
}

per <- joined |>
  group_by(gene, study) |>
  summarize(n     = n(),
            full  = auroc(truth, as.numeric(varviz_pts_full)),
            blind = auroc(truth, as.numeric(varviz_pts_blind)),
            .groups = "drop")

cat(sprintf("\nscoresets: %d   genes: %d\n", nrow(per), n_distinct(per$gene)))
print(as.data.frame(per), row.names = FALSE, digits = 3)

for (nm in c("full", "blind")) {
  v <- per[[nm]]; v <- v[!is.na(v)]
  cat(sprintf("\nPass-%s: median %.3f  range %.3f-%.3f  (n=%d scoresets)\n",
              if (nm == "full") "Full" else "Blind",
              median(v), min(v), max(v), length(v)))
}

cat("\nPer-gene mean AUROC:\n")
print(as.data.frame(
  per |> group_by(gene) |>
    summarize(n_scoresets = n(),
              mean_full = mean(full, na.rm = TRUE),
              mean_blind = mean(blind, na.rm = TRUE),
              delta = mean(blind, na.rm = TRUE) - mean(full, na.rm = TRUE),
              .groups = "drop") |> arrange(desc(mean_full))),
  row.names = FALSE, digits = 3)
