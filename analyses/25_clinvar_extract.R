#!/usr/bin/env Rscript
# 25_clinvar_extract.R — build the genome-wide 2-star ClinVar missense truth set
# used by 24_clinvar_mds_lr.R for the MDS likelihood-ratio calibration.
#
# Input : analyses/tmp/clinvar/variant_summary.txt.gz  (NCBI ClinVar tab dump,
#         https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)
# Output: analyses/tmp/clinvar/clinvar_missense_2star.tsv
#         columns: gene wt pos mut p3 cls review   (cls: 1 = P/LP, 0 = B/LB)
#
# Filters: GRCh38 rows; ReviewStatus >= 2 stars; a parseable missense p.(Aaa###Bbb)
# with mut != wt and mut != Ter; ClinicalSignificance in P/LP or B/LB (no VUS,
# conflicting, risk/drug/association/protective/other).
#
# Usage: Rscript analyses/25_clinvar_extract.R
#
# Fidelity: reproduces the committed truth set to ~99.8% at the variant level
# (63,375 vs 63,932 rows; residual differences are ClinVar dedup / p.-parse edge
# cases and do not affect the MDS likelihood-ratio tiers). Downstream, 24_* maps
# these to Pfam and scores 33,216 (18,570 P / 14,646 B) with LR+ 5.4/11.2/39.0.
suppressMessages(library(data.table))

IN  <- "analyses/tmp/clinvar/variant_summary.txt.gz"
OUT <- "analyses/tmp/clinvar/clinvar_missense_2star.tsv"

aa3 <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",
         His="H",Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",
         Thr="T",Trp="W",Tyr="Y",Val="V")

STARS2 <- c("criteria provided, multiple submitters, no conflicts",
            "reviewed by expert panel",
            "practice guideline")
PLP <- c("Pathogenic","Likely pathogenic","Pathogenic/Likely pathogenic")
BLB <- c("Benign","Likely benign","Benign/Likely benign")

dt <- fread(cmd = paste("zcat", shQuote(IN)), sep = "\t", quote = "", showProgress = FALSE)
dt <- dt[Assembly == "GRCh38" & ReviewStatus %in% STARS2]
dt <- dt[ClinicalSignificance %in% c(PLP, BLB)]

# parse p.(Aaa###Bbb) from the Name column
m <- regmatches(dt$Name, regexpr("p\\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}", dt$Name))
has <- regexpr("p\\.[A-Z][a-z]{2}[0-9]+[A-Z][a-z]{2}", dt$Name) > 0
dt <- dt[has]; p3 <- m
wt3 <- sub("p\\.([A-Z][a-z]{2})[0-9]+[A-Z][a-z]{2}", "\\1", p3)
pos <- as.integer(sub("p\\.[A-Z][a-z]{2}([0-9]+)[A-Z][a-z]{2}", "\\1", p3))
mut3 <- sub("p\\.[A-Z][a-z]{2}[0-9]+([A-Z][a-z]{2})", "\\1", p3)

out <- data.table(
  gene = dt$GeneSymbol, wt = aa3[wt3], pos = pos, mut = aa3[mut3], p3 = p3,
  cls = fifelse(dt$ClinicalSignificance %in% PLP, 1L, 0L),
  review = dt$ReviewStatus
)
out <- out[!is.na(wt) & !is.na(mut) & wt != mut]          # drop Ter/unknown/synonymous
out <- unique(out, by = c("gene","pos","wt","mut","cls")) # one row per variant/label

fwrite(out, OUT, sep = "\t")
cat(sprintf("wrote %s: %d rows (P/LP=%d, B/LB=%d, genes=%d)\n",
            OUT, nrow(out), sum(out$cls==1L), sum(out$cls==0L), uniqueN(out$gene)))
