#!/usr/bin/env Rscript
# ============================================================================
# preprocess_ccrs_to_rdata.R
# 
# Reads CCRStoAAC_output.tsv, selects only VarViz-relevant columns,
# compresses into a compact data.table, and adds to VarViz.RData
#
# Usage:
#   Rscript preprocess_ccrs_to_rdata.R <path_to_CCRStoAAC_output.tsv> <path_to_VarViz.RData>
#
# Example:
#   Rscript preprocess_ccrs_to_rdata.R ./CCRStoAAC_output.tsv ./data/VarViz.RData
# ============================================================================

library(data.table)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript preprocess_ccrs_to_rdata.R <CCRStoAAC_output.tsv> <VarViz.RData>\n")
  cat("\nExample:\n")
  cat("  Rscript preprocess_ccrs_to_rdata.R ./CCRStoAAC_output.tsv ./data/VarViz.RData\n")
  quit(status = 1)
}

ccrs_file  <- args[1]
rdata_file <- args[2]

cat("=== VarViz CCRS Preprocessing ===\n\n")

# ── Step 1: Read the full CCRStoAAC output ──
cat("[1/5] Reading CCRStoAAC_output.tsv ...\n")
cat("      File:", ccrs_file, "\n")

ccrs_full <- fread(ccrs_file, sep = "\t", header = TRUE, na.strings = c("NA", ""))

cat("      Total rows:", format(nrow(ccrs_full), big.mark = ","), "\n")
cat("      Total columns:", ncol(ccrs_full), "\n")
cat("      Unique genes:", length(unique(ccrs_full$ensembl_gene_name)), "\n")
cat("      Unique proteins:", length(unique(ccrs_full$uniprotAcc)), "\n\n")

# ── Step 2: Select only VarViz-relevant columns ──
cat("[2/5] Selecting VarViz-relevant columns ...\n")

keep_cols <- c(
  "ensembl_gene_name",     # gene lookup key
  "uniprotAcc",            # UniProt ID match
  "aac_pos",               # amino acid position
  "aac_weighted_pct",      # CCR percentile (0-100)
  "aac_weighted_pct_cat",  # categorized constraint (0-6, -1/-2/-3)
  "CONSERVATION_score"     # ScoreCons inter-species conservation
)

# Verify all columns exist
missing <- setdiff(keep_cols, names(ccrs_full))
if (length(missing) > 0) {
  stop("Missing columns in input file: ", paste(missing, collapse = ", "))
}

ccrs_slim <- ccrs_full[, ..keep_cols]
cat("      Selected", length(keep_cols), "columns\n\n")

# ── Step 3: Clean and optimize data types ──
cat("[3/5] Cleaning and optimizing data types ...\n")

# Convert position to integer
ccrs_slim[, aac_pos := as.integer(aac_pos)]

# Convert percentile to numeric (already should be, but ensure)
ccrs_slim[, aac_weighted_pct := as.numeric(aac_weighted_pct)]

# Convert category to integer
ccrs_slim[, aac_weighted_pct_cat := as.integer(aac_weighted_pct_cat)]

# Convert conservation score to numeric
ccrs_slim[, CONSERVATION_score := as.numeric(CONSERVATION_score)]

# Remove rows where gene name is NA/empty (shouldn't happen but safety)
ccrs_slim <- ccrs_slim[!is.na(ensembl_gene_name) & ensembl_gene_name != ""]

# Set keys for fast lookup
setkey(ccrs_slim, ensembl_gene_name, uniprotAcc)

cat("      Rows after cleaning:", format(nrow(ccrs_slim), big.mark = ","), "\n")

# Size comparison
full_size <- object.size(ccrs_full)
slim_size <- object.size(ccrs_slim)
cat(sprintf("      Size reduction: %.1f MB → %.1f MB (%.0f%% smaller)\n",
            full_size / 1e6, slim_size / 1e6,
            (1 - as.numeric(slim_size) / as.numeric(full_size)) * 100))
cat("\n")

# ── Step 4: Summary statistics ──
cat("[4/5] Summary statistics ...\n")

# CCR percentile distribution
cat("      CCR percentile distribution:\n")
cat(sprintf("        Unconstrained (pct = 0):    %s positions\n",
            format(sum(ccrs_slim$aac_weighted_pct == 0, na.rm = TRUE), big.mark = ",")))
cat(sprintf("        Constrained (pct > 0):      %s positions\n",
            format(sum(ccrs_slim$aac_weighted_pct > 0, na.rm = TRUE), big.mark = ",")))
cat(sprintf("        Highly constrained (≥95):    %s positions\n",
            format(sum(ccrs_slim$aac_weighted_pct >= 95, na.rm = TRUE), big.mark = ",")))
cat(sprintf("        Top constrained (≥99):       %s positions\n",
            format(sum(ccrs_slim$aac_weighted_pct >= 99, na.rm = TRUE), big.mark = ",")))
cat(sprintf("        Not available (pct < 0):     %s positions\n",
            format(sum(ccrs_slim$aac_weighted_pct < 0, na.rm = TRUE), big.mark = ",")))

# Conservation score coverage
cons_available <- sum(!is.na(ccrs_slim$CONSERVATION_score) & ccrs_slim$CONSERVATION_score >= 0)
cat(sprintf("        Conservation scores available: %s positions\n",
            format(cons_available, big.mark = ",")))
cat("\n")

# ── Step 5: Add to VarViz.RData ──
cat("[5/5] Adding ccrs_data to VarViz.RData ...\n")
cat("      Loading existing:", rdata_file, "\n")

# Load existing RData — this brings in gene_data (and any other objects)
existing_objects <- load(rdata_file)
cat("      Existing objects:", paste(existing_objects, collapse = ", "), "\n")

# Rename for clarity in the RData
ccrs_data <- ccrs_slim

# Save everything back — existing objects + ccrs_data
all_objects <- c(existing_objects, "ccrs_data")
save(list = all_objects, file = rdata_file, compress = "xz")

cat("      Saved objects:", paste(all_objects, collapse = ", "), "\n")

# Final file size
final_size <- file.info(rdata_file)$size
cat(sprintf("      Final RData file size: %.1f MB\n", final_size / 1e6))

cat("\n=== Done! ===\n")
cat("VarViz can now look up CCRS by gene name:\n")
cat('  ccrs_data[ensembl_gene_name == "KCNQ2"]\n')
cat('  ccrs_data[uniprotAcc == "O43526-1"]\n')
cat("\nNo more BED file or bedr dependency needed!\n")
