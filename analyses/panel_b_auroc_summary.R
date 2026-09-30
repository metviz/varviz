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
source("analyses/lib/run_inputs.R")
source("analyses/lib/dms_polarity.R")
DMS_POL <- load_dms_polarity()
DMS_OFF <- load_dms_offsets()

DMS_ELIG <- load_dms_eligibility()
REVAL_DIR <- Sys.getenv("REVAL_DIR", "analyses/ps_final_v231/classifications")
# The UNIVERSE, not mavedb_canonical.tsv. The canonical table is raw MaveDB in each
# scoreset's own numbering; the universe has the reference-numbering correction
# applied, which is the numbering the classifications are keyed on. Reading the raw
# table meant TP53 00000059-a-1 joined on construct positions and matched nothing once
# the universe was corrected. Its sibling scripts 07_panel_b_dms.R and
# panel_b_dms_reval.R already read the universe.
MAVEDB    <- Sys.getenv("MAVEDB", "analyses/derived/variant_universe.tsv")

files <- list.files(REVAL_DIR, pattern = "__dual\\.tsv$", full.names = TRUE)
if (!length(files)) stop("no __dual.tsv files under ", REVAL_DIR)
cls <- bind_rows(lapply(files, read_tsv, show_col_types = FALSE)) |>
  distinct(gene, p_notation, .keep_all = TRUE)
mave <- read_tsv(MAVEDB, show_col_types = FALSE)
if ("source" %in% names(mave)) mave <- mave[mave$source == "MaveDB", , drop = FALSE]

score_col <- intersect(c("score_raw", "score"), names(mave))[1]
if (is.na(score_col)) stop("no score column in ", MAVEDB)

joined <- inner_join(mave, cls, by = c("gene", "p_notation"))
joined$score_dms <- joined[[score_col]]

# Per-study quartiles, identical to the panel, but which tail is loss-of-function
# is per scoreset: nine of the 21 assays score higher for more damage, so a single
# convention reported 1 - the true AUROC for those. See analyses/lib/dms_polarity.R.
# Drop scoresets numbered on a different sequence before anything else: their
# variants were classified at the wrong residues, so neither sign convention makes
# the AUROC mean anything.
# Scoresets whose axis does not measure pathogenicity cannot serve as truth however
# correctly they are read. Recorded per scoreset with its reason in
# analyses/derived/mavedb_scoreset_polarity.tsv; see analyses/lib/dms_polarity.R.
joined$dms_inelig <- vapply(as.character(joined$study),
                      function(s) dms_ineligible_reason(DMS_ELIG, s), character(1))
inelig <- sort(unique(joined$study[nzchar(joined$dms_inelig)]))
if (length(inelig)) {
  # One message per distinct reason: six identical paragraphs is noise, and the
  # scoreset list is the part a reader needs alongside it.
  for (r in unique(vapply(inelig, function(s) dms_ineligible_reason(DMS_ELIG, s),
                          character(1)))) {
    grp <- inelig[vapply(inelig, function(s) identical(dms_ineligible_reason(DMS_ELIG, s), r),
                         logical(1))]
    cat(sprintf("[panel_b] EXCLUDED %d scoreset(s): %s\n", length(grp),
                paste(grp, collapse = ", ")))
    cat(sprintf("[panel_b]   reason: %s\n", r))
  }
  joined <- joined[!nzchar(joined$dms_inelig), , drop = FALSE]
}

joined$dms_off <- dms_pending_offset(joined)
shifted <- sort(unique(joined$study[joined$dms_off != 0L]))
if (length(shifted)) {
  cat(sprintf("\n[panel_b] EXCLUDED, numbered on a different sequence (offset != 0): %s\n",
              paste(shifted, collapse = ", ")))
  cat("[panel_b]   these need the universe rebuilt on reference numbering and re-classified\n")
  joined <- joined[joined$dms_off == 0L, , drop = FALSE]
}

joined$dms_dir <- vapply(as.character(joined$study),
                         function(s) dms_direction_for(DMS_POL, s), character(1))
joined <- joined |>
  group_by(study) |>
  mutate(q25 = quantile(score_dms, .25, na.rm = TRUE),
         q75 = quantile(score_dms, .75, na.rm = TRUE),
         truth = dms_truth(score_dms, q25, q75, dms_dir)) |>
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

# This table is redirected to analyses/reports/panel_b_auroc_*.txt, so the stamp has
# to travel in the stream rather than sit in the terminal that happened to run it.
prov <- run_provenance_lines()
cat("\n", paste0("# ", prov, collapse = "\n"), "\n", sep = "")
