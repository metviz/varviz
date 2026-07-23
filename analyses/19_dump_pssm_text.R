#!/usr/bin/env Rscript
# 19_dump_pssm_text.R — text companions for pfam_pssm_human.rds.
#
# The rds stores the score matrices as int8 (binary); this writes grep-able text:
#   *_map.tsv       human residue -> (family, column) lookup
#   *_families.tsv  per-family column count + human residues mapped
# Use when the rds already exists (e.g. built by an earlier run of 15p) and you
# just want the text views without rebuilding.
#
# Usage:  Rscript analyses/19_dump_pssm_text.R [path/to/pfam_pssm_human.rds]
RDS <- if (length(commandArgs(trailingOnly = TRUE)))
         commandArgs(trailingOnly = TRUE)[1] else "analyses/derived/pfam_pssm_human.rds"
stopifnot(file.exists(RDS))
out <- readRDS(RDS)
map <- out$map; pssm <- out$pssm

MAP_TSV <- sub("\\.rds$", "_map.tsv", RDS)
FAM_TSV <- sub("\\.rds$", "_families.tsv", RDS)
write.table(map[, c("id", "residue", "family", "column")], MAP_TSV,
            sep = "\t", row.names = FALSE, quote = FALSE)
fam <- data.frame(
  family      = names(pssm),
  columns     = vapply(pssm, nrow, integer(1)),
  human_resid = as.integer(table(factor(map$family, levels = names(pssm)))),
  stringsAsFactors = FALSE)
fam <- fam[order(-fam$human_resid), ]
write.table(fam, FAM_TSV, sep = "\t", row.names = FALSE, quote = FALSE)

cat(sprintf("wrote %s (%s rows)\n      %s (%s families)\n",
            MAP_TSV, format(nrow(map), big.mark = ","),
            FAM_TSV, format(nrow(fam), big.mark = ",")))
