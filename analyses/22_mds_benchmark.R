#!/usr/bin/env Rscript
# 22_mds_benchmark.R — re-run the dual-pass benchmark's Path-4 augmentation with
# MDS instead of DOLPHIN, and report the clinical metrics side by side.
#
# The domain-alignment pathway was evaluated in two stages in the original
# benchmark: the in-engine classify ran with Path 4 disabled (DOLPHIN's rate
# limit made per-variant calls infeasible), and DOLPHIN PM1 calls were merged in
# afterwards (11_merge_dolphin.R) under the rule "add +2 to Pass-Blind iff the
# pathway score fires AND pm1_pathway is '' or 'clinvar_hotspot'". MDS needs no
# service, so we compute its PM1 call offline for every variant here and apply
# the identical augmentation rule, then recompute Panel A on the MDS-augmented
# Pass-Blind. Reproducing Full / Blind / Blind+DOLPHIN alongside Blind+MDS
# validates the pipeline and yields the true MDS figures for the manuscript.
#
# Usage:  Rscript analyses/22_mds_benchmark.R
suppressMessages({library(data.table); library(readr); library(dplyr)})
source("analyses/lib/pssm_lookup.R")
source("analyses/lib/metric_suite.R")

MDS_THRESHOLD <- -4
IN   <- "analyses/derived/varviz_classifications_dolphin.tsv"
UNIV <- "analyses/derived/variant_universe.tsv"
OUT  <- "analyses/derived/varviz_classifications_mds.tsv"

pts_to_bin <- function(pts) ifelse(pts>=10,"Pathogenic",ifelse(pts>=6,"Likely Pathogenic",
  ifelse(pts>=4,"VUS-High",ifelse(pts>=2,"VUS-Mid",ifelse(pts>=-3,"VUS-Low",
  ifelse(pts>=-6,"Likely Benign","Benign"))))))

cls <- fread(IN, sep = "\t")
tbl <- pssm_table_load("data/pfam_pssm_human.rds")
load("data/VarViz.RData")                                    # gene_data
acc_of   <- setNames(gene_data$uniprot_id, gene_data$gene_name)
ae       <- fread("analyses/tmp/mane/uniprot_acc_entry.tsv") # accession -> entry (from 21)
entry_of <- setNames(ae$entry, ae$accession)

# ── Compute MDS PM1 for every variant (offline keyed join) ─────────────────
aa3 <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",His="H",
         Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",Thr="T",Trp="W",
         Tyr="Y",Val="V")
# p_notation appears in both 3-letter (p.Arg528Trp) and 1-letter (p.N290S) forms
# across the pipeline files; parse either.
p  <- sub("^p\\.", "", cls$p_notation)
m3 <- regmatches(p, regexec("^([A-Za-z]{3})([0-9]+)([A-Za-z]{3})$", p))
m1 <- regmatches(p, regexec("^([A-Z])([0-9]+)([A-Z])$", p))
onelet <- function(three, one, idx3, idx1) {
  if (length(three) == 4) unname(aa3[three[idx3]]) else
  if (length(one) == 4)   one[idx1] else NA_character_
}
cls[, wt  := mapply(function(a,b) onelet(a,b,2L,2L), m3, m1)]
cls[, pos := mapply(function(a,b) if (length(a)==4) as.integer(a[3]) else if (length(b)==4) as.integer(b[3]) else NA_integer_, m3, m1)]
cls[, mut := mapply(function(a,b) onelet(a,b,4L,4L), m3, m1)]
cls[, entry := entry_of[acc_of[gene]]]
cls[, vid := .I]

mp <- as.data.table(tbl$map); setkeyv(mp, c("id","residue"))
Q_MIN <- tbl$quant[["min"]]; Q_MAX <- tbl$quant[["max"]]
deq <- function(q) (as.numeric(q)+127)/254*(Q_MAX-Q_MIN)+Q_MIN
scor <- cls[!is.na(entry) & !is.na(pos) & !is.na(wt) & !is.na(mut)]
hits <- mp[scor[, .(id=entry, residue=pos, vid, wt, mut)],
           on=.(id,residue), nomatch=0L, allow.cartesian=TRUE]
hits <- hits[!is.na(family) & nzchar(family)]
hits[, delta := {
  fam <- as.character(.BY$family)
  M <- if (length(fam) == 1L && !is.na(fam) && fam %in% names(tbl$pssm)) tbl$pssm[[fam]] else NULL
  if (is.null(M)) rep(NA_real_, .N) else {
    mc <- match(mut, colnames(M)); wc <- match(wt, colnames(M))
    ok <- column <= nrow(M) & !is.na(mc) & !is.na(wc); d <- rep(NA_real_, .N)
    d[ok] <- deq(M[cbind(column[ok], mc[ok])]) - deq(M[cbind(column[ok], wc[ok])]); d
  }
}, by = family]
per <- hits[!is.na(delta), .(mds = min(delta)), by = vid]
cls[per, mds := i.mds, on = "vid"]
cls[, mds_pm1 := !is.na(mds) & mds <= MDS_THRESHOLD]

# ── Augment Pass-Blind with MDS (identical rule to DOLPHIN merge) ───────────
aug <- cls$mds_pm1 & (cls$pm1_pathway == "" | cls$pm1_pathway == "clinvar_hotspot")
cls[, varviz_pts_blind_mds := varviz_pts_blind]
cls[, varviz_classification_blind_mds := varviz_classification_blind]
cls[aug, varviz_pts_blind_mds := varviz_pts_blind + 2L]
cls[aug, varviz_classification_blind_mds := pts_to_bin(varviz_pts_blind_mds)]
fwrite(cls, OUT, sep = "\t", quote = FALSE, na = "")

# ── MDS vs DOLPHIN firing agreement (where both are known) ─────────────────
cmp <- cls[!is.na(dolphin_pm1)]
cat(sprintf("\n[bench] MDS vs DOLPHIN PM1 firing on %s variants with a DOLPHIN call:\n",
            format(nrow(cmp), big.mark=",")))
cat(sprintf("  agree %.1f%%  (both fire %d, both no %d, MDS-only %d, DOLPHIN-only %d)\n",
    100*mean(cmp$mds_pm1 == cmp$dolphin_pm1),
    sum(cmp$mds_pm1 & cmp$dolphin_pm1), sum(!cmp$mds_pm1 & !cmp$dolphin_pm1),
    sum(cmp$mds_pm1 & !cmp$dolphin_pm1), sum(!cmp$mds_pm1 & cmp$dolphin_pm1)))
cat(sprintf("  MDS PM1 fires on %s / %s variants overall (%.1f%%); augments %s Pass-Blind rows\n",
    format(sum(cls$mds_pm1), big.mark=","), format(nrow(cls), big.mark=","),
    100*mean(cls$mds_pm1), format(sum(aug), big.mark=",")))

# ── Panel A on VariBench labels: Full / Blind / Blind+DOLPHIN / Blind+MDS ───
univ <- as.data.table(read_tsv(UNIV, show_col_types = FALSE))[source == "VariBench"]
j <- merge(cls, univ[, .(gene, p_notation, label)], by = c("gene","p_notation"))
j[, truth := fifelse(grepl("[Pp]athogenic", label), "Pathogenic",
              fifelse(grepl("[Bb]enign", label), "Benign", NA_character_))]
j <- j[!is.na(truth)]
cat(sprintf("\n[bench] VariBench-labelled with classifications: %d (P=%d, B=%d)\n",
            nrow(j), sum(j$truth=="Pathogenic"), sum(j$truth=="Benign")))

grab <- function(ms, k) ms$estimate[ms$metric == k]
row1 <- function(name, cl, sc) {
  ms <- compute_metric_suite(j$truth, cl, sc)
  cat(sprintf("  %-16s AUROC %.3f  Sens %.3f  Spec %.3f  MCC %.3f  VUS %.3f\n",
      name, grab(ms,"AUROC"), grab(ms,"Sensitivity"), grab(ms,"Specificity"),
      grab(ms,"MCC"), grab(ms,"VUS_rate")))
}
cat("\n===== Panel A (VariBench, dual-pass) =====\n")
row1("Pass-Full",        j$varviz_classification_full,          j$varviz_pts_full)
row1("Pass-Blind",       j$varviz_classification_blind,         j$varviz_pts_blind)
row1("Blind + DOLPHIN",  j$varviz_classification_blind_dolphin, j$varviz_pts_blind_dolphin)
row1("Blind + MDS",      j$varviz_classification_blind_mds,     j$varviz_pts_blind_mds)
cat("\nwritten: ", OUT, "\n")
