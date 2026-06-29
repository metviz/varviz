#!/usr/bin/env Rscript
# 11_merge_dolphin.R — DOLPHIN-augmented Pass-Blind classifications.
#
# Joins per-gene DOLPHIN PM1 calls (analyses/raw/dolphin/by_gene/{GENE}.tsv,
# produced by analyses/run_dolphin_benchmark.R) into the dual-pass output
# (analyses/derived/varviz_classifications.tsv) and emits a new file
# (analyses/derived/varviz_classifications_dolphin.tsv) with two added
# columns:
#   dolphin_pm1               : TRUE / FALSE / NA (NA = variant not in DOLPHIN cache)
#   varviz_classification_blind_dolphin
#   varviz_pts_blind_dolphin
#
# Rule for adding DOLPHIN PM1 to Pass-Blind:
#   * If `pm1_pathway` is "" or "clinvar_hotspot", the Pass-Blind run currently
#     has NO PM1 (clinvar_hotspot is stripped by `strip_clinvar_tags()`; the
#     empty pathway never fired PM1 at all).
#   * If DOLPHIN reports PM1 for that (gene, p_notation) pair, add +2 Tavtigian
#     points and re-derive the 7-level bin.
#   * If `pm1_pathway` is "ccrs" or "uniprot_domain", Pass-Blind retains its
#     existing PM1 and we leave the row unchanged (no double-counting).
#
# Tavtigian bin thresholds (mirrors server.R):
#   Pathogenic        >= 10
#   Likely Pathogenic  6 .. 9
#   VUS-High          +4 .. +5
#   VUS-Mid           +2 .. +3
#   VUS-Low           -3 .. +1
#   Likely Benign     -4 .. -6
#   Benign            <= -7

suppressPackageStartupMessages({
  library(data.table)
})

UNIVERSE_TSV    <- "analyses/derived/varviz_classifications.tsv"
DOLPHIN_DIR     <- "analyses/raw/dolphin/by_gene"
OUT_TSV         <- "analyses/derived/varviz_classifications_dolphin.tsv"

# ---------------------------------------------------------------------------
# Bin function — single source of truth for re-derivation.
# ---------------------------------------------------------------------------
pts_to_bin <- function(pts) {
  ifelse(pts >= 10,            "Pathogenic",
  ifelse(pts >= 6,             "Likely Pathogenic",
  ifelse(pts >= 4,             "VUS-High",
  ifelse(pts >= 2,             "VUS-Mid",
  ifelse(pts >= -3,            "VUS-Low",
  ifelse(pts >= -6,            "Likely Benign",
                                "Benign"))))))
}

# ---------------------------------------------------------------------------
# Load universe + build DOLPHIN PM1 lookup
# ---------------------------------------------------------------------------
cat("[merge] reading", UNIVERSE_TSV, "...\n")
uni <- fread(UNIVERSE_TSV, sep = "\t", header = TRUE)
cat("[merge] universe:", nrow(uni), "rows,", ncol(uni), "cols\n")

cat("[merge] reading DOLPHIN per-gene TSVs from", DOLPHIN_DIR, "...\n")
dolphin_files <- list.files(DOLPHIN_DIR, pattern = "\\.tsv$", full.names = TRUE)
cat("[merge]   found", length(dolphin_files), "gene caches\n")

dolphin_long <- rbindlist(lapply(dolphin_files, function(f) {
  g <- sub("\\.tsv$", "", basename(f))
  d <- tryCatch(fread(f, sep = "\t", header = TRUE),
                error = function(e) {
                  cat("[merge]   skip", g, "(parse error)\n"); NULL
                })
  if (is.null(d) || nrow(d) == 0) return(NULL)
  data.table(gene = g, p_notation = d$p_notation, dolphin_pm1 = as.logical(d$pm1))
}), fill = TRUE)
cat("[merge] dolphin_long:", nrow(dolphin_long), "rows across",
    length(unique(dolphin_long$gene)), "genes\n")

setkey(dolphin_long, gene, p_notation)
setkey(uni, gene, p_notation)
uni <- dolphin_long[uni]   # left-join (every uni row, dolphin_pm1 may be NA)
setcolorder(uni, c(setdiff(colnames(uni), "dolphin_pm1"), "dolphin_pm1"))

# ---------------------------------------------------------------------------
# Augment Pass-Blind w/ DOLPHIN PM1 where applicable
# ---------------------------------------------------------------------------
# Apply rule: add PM1 (+2 pts) to Pass-Blind iff
#   dolphin_pm1 is TRUE  AND  pm1_pathway is "" or "clinvar_hotspot"
#
# We do NOT touch rows where pm1_pathway is "ccrs" / "uniprot_domain" — those
# already retain PM1 in Pass-Blind (strip_clinvar_tags only strips the
# clinvar_hotspot pathway, never CCRS or uniprot_domain).

augment_mask <- !is.na(uni$dolphin_pm1) & uni$dolphin_pm1 &
                (uni$pm1_pathway == "" | uni$pm1_pathway == "clinvar_hotspot")

uni[, varviz_pts_blind_dolphin       := varviz_pts_blind]
uni[, varviz_classification_blind_dolphin := varviz_classification_blind]
uni[augment_mask, varviz_pts_blind_dolphin       := varviz_pts_blind + 2L]
uni[augment_mask, varviz_classification_blind_dolphin := pts_to_bin(varviz_pts_blind_dolphin)]

# ---------------------------------------------------------------------------
# Diagnostic summary
# ---------------------------------------------------------------------------
cat("\n[merge] dolphin coverage:\n")
print(uni[, .N, by = .(has_dolphin = !is.na(dolphin_pm1), dolphin_pm1 = dolphin_pm1)])

cat("\n[merge] augmentation summary (rows that gained PM1 from DOLPHIN):\n")
cat("  augmented:", sum(augment_mask), "/", nrow(uni),
    sprintf("(%.2f%%)\n", 100 * sum(augment_mask) / nrow(uni)))

cat("\n[merge] Pass-Blind classification distribution (BEFORE -> AFTER DOLPHIN):\n")
before <- table(uni$varviz_classification_blind)
after  <- table(uni$varviz_classification_blind_dolphin)
shifts <- merge(
  as.data.frame(before, stringsAsFactors = FALSE),
  as.data.frame(after,  stringsAsFactors = FALSE),
  by = "Var1", all = TRUE
)
colnames(shifts) <- c("bin", "before", "after")
shifts$delta <- shifts$after - shifts$before
shifts <- shifts[order(-abs(shifts$delta)), ]
print(shifts, row.names = FALSE)

cat("\n[merge] subclass-resolved Pass-Blind transitions (before -> after DOLPHIN):\n")
xt <- table(
  before = uni$varviz_classification_blind,
  after  = uni$varviz_classification_blind_dolphin
)
print(xt)

# ---------------------------------------------------------------------------
# Write augmented TSV
# ---------------------------------------------------------------------------
cat("\n[merge] writing", OUT_TSV, "...\n")
fwrite(uni, OUT_TSV, sep = "\t", quote = FALSE, na = "")
cat("[merge] DONE.\n")
