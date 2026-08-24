#!/usr/bin/env Rscript
# 23_mds_threshold_heldout.R — is the MDS PM1 threshold held out from Panel A?
#
# analyses/21_assess_mds.R picks the -4 firing threshold on domain-resident
# VariBench (PON-PS_D2) variants and applies no gene filter, so the 744
# variants Panel A scores against — the PON-PS_D2 variants in the 14 benchmark
# genes — sit inside the calibration pool. The MDS increment on Panel A
# (Pass-Blind AUROC 0.933 -> 0.958) is therefore not fully held out.
#
# This quantifies the overlap and re-selects the threshold with the benchmark
# genes removed, reporting the LR+ table three ways:
#   ALL        every scorable PON-PS_D2 variant (what 21_assess_mds.R used)
#   IN-SAMPLE  benchmark-gene variants only (the contaminated portion)
#   HELD OUT   benchmark genes excluded (the honest calibration)
#
# If -4 still delivers Moderate strength (LR+ >= 4.33) on HELD OUT, the
# production threshold needs no change and the overlap is a disclosure item,
# not a correction.
#
# Also runs leave-one-gene-out across the 14 benchmark genes so threshold
# stability can be reported rather than asserted.
#
# Usage:  Rscript analyses/23_mds_threshold_heldout.R

suppressMessages({library(data.table)})
source("analyses/lib/pssm_lookup.R")

BENCH <- c("BAP1","BRCA1","CASR","DDR2","GCK","KCNH2","KCNQ1","KRAS","LDLR",
           "NUDT15","PTEN","SLC13A5","TP53","TSHR")

tbl <- pssm_table_load("data/pfam_pssm_human.rds")
vb  <- fread("analyses/derived/varibench_canonical.tsv")
vb  <- vb[label %in% c("Pathogenic","Benign")]
load("data/VarViz.RData")
acc_of <- setNames(gene_data$uniprot_id, gene_data$gene_name)
ae <- fread("analyses/tmp/mane/uniprot_acc_entry.tsv")
entry_of <- setNames(ae$entry, ae$accession)

aa3 <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",His="H",
         Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",Thr="T",Trp="W",
         Tyr="Y",Val="V")
p <- sub("^p\\.", "", vb$p_notation)
m <- regmatches(p, regexec("^([A-Za-z]{3})([0-9]+)([A-Za-z]{3})$", p))
vb[, wt  := vapply(m, function(x) if (length(x)==4) unname(aa3[x[2]]) else NA_character_, "")]
vb[, pos := vapply(m, function(x) if (length(x)==4) as.integer(x[3]) else NA_integer_, integer(1))]
vb[, mut := vapply(m, function(x) if (length(x)==4) unname(aa3[x[4]]) else NA_character_, "")]
vb[, entry := entry_of[acc_of[gene]]]
vb <- vb[!is.na(entry) & !is.na(pos) & !is.na(wt) & !is.na(mut)]
vb[, vid := .I]

mp <- as.data.table(tbl$map); setkeyv(mp, c("id","residue"))
Q_MIN <- tbl$quant[["min"]]; Q_MAX <- tbl$quant[["max"]]
deq <- function(q) (as.numeric(q) + 127) / 254 * (Q_MAX - Q_MIN) + Q_MIN

hits <- mp[vb[, .(id = entry, residue = pos, vid, wt, mut)],
           on = .(id, residue), nomatch = 0L, allow.cartesian = TRUE]
hits[, delta := {
  # as.character() is load-bearing: pssm_table_load() stores map$family as a
  # factor to save memory, and tbl$pssm[[<factor>]] indexes by integer LEVEL,
  # not by family name -- silently returning the wrong matrix for every group.
  M <- tbl$pssm[[as.character(.BY$family)]]
  if (is.null(M)) rep(NA_real_, .N) else {
    mc <- match(mut, colnames(M)); wc <- match(wt, colnames(M))
    ok <- column <= nrow(M) & !is.na(mc) & !is.na(wc)
    d <- rep(NA_real_, .N)
    d[ok] <- deq(M[cbind(column[ok], mc[ok])]) - deq(M[cbind(column[ok], wc[ok])])
    d
  }
}, by = family]
per <- hits[!is.na(delta), .(mds = min(delta)), by = vid]
vb[per, mds := i.mds, on = "vid"]
scored <- vb[!is.na(mds)]

strength <- function(lr) {
  if (is.infinite(lr) || lr >= 18.7) "STRONG"
  else if (lr >= 4.33) "moderate"
  else if (lr >= 2.08) "supporting"
  else "-"
}

lr_table <- function(d, label) {
  np <- sum(d$label=="Pathogenic"); nb <- sum(d$label=="Benign")
  if (np < 20 || nb < 20) {
    cat(sprintf("\n== %s ==  n=%d (P=%d, B=%d) -- too few to calibrate\n",
                label, nrow(d), np, nb)); return(invisible(NULL))
  }
  s <- -d$mds; y <- as.integer(d$label=="Pathogenic"); r <- rank(s)
  auc <- (sum(r[y==1]) - np*(np+1)/2) / (as.numeric(np)*nb)
  cat(sprintf("\n== %s ==  n=%s (P=%d, B=%d)   AUC = %.3f\n",
              label, format(nrow(d), big.mark=","), np, nb, auc))
  cat(sprintf("  %6s %7s %7s %8s  %s\n","t","sens","spec","LR+","strength"))
  for (t in c(-2,-3,-4,-5,-6,-8,-12)) {
    tp <- sum(d$label=="Pathogenic" & d$mds<=t); fp <- sum(d$label=="Benign" & d$mds<=t)
    sens <- tp/np; fpr <- fp/nb; lr <- if (fpr==0) Inf else sens/fpr
    cat(sprintf("  %6d %6.1f%% %6.1f%% %8s  %s\n", t, 100*sens, 100*(1-fpr),
        if (is.infinite(lr)) "Inf" else sprintf("%.1f", lr), strength(lr)))
  }
  invisible(auc)
}

inb <- scored$gene %in% BENCH
cat("===== MDS threshold: calibration/evaluation overlap =====\n")
cat(sprintf("scorable PON-PS_D2 variants:        %s\n", format(nrow(scored), big.mark=",")))
cat(sprintf("  in the 14 benchmark genes:        %s (%.1f%%)  <- also Panel A's evaluation set\n",
            format(sum(inb), big.mark=","), 100*sum(inb)/nrow(scored)))
cat(sprintf("  outside them:                     %s\n", format(sum(!inb), big.mark=",")))
cat(sprintf("distinct benchmark genes present:   %d of 14\n",
            length(unique(scored$gene[inb]))))

lr_table(scored,        "ALL scorable (what 21_assess_mds.R used)")
lr_table(scored[inb],   "IN-SAMPLE: benchmark genes only")
lr_table(scored[!inb],  "HELD OUT: benchmark genes excluded")

cat("\n===== leave-one-gene-out, threshold -4 =====\n")
cat(sprintf("  %-10s %7s %7s %8s  %s\n","excluded","sens","spec","LR+","strength"))
for (g in c("(none)", BENCH)) {
  d <- if (g == "(none)") scored else scored[gene != g]
  np <- sum(d$label=="Pathogenic"); nb <- sum(d$label=="Benign")
  tp <- sum(d$label=="Pathogenic" & d$mds<=-4); fp <- sum(d$label=="Benign" & d$mds<=-4)
  sens <- tp/np; fpr <- fp/nb; lr <- if (fpr==0) Inf else sens/fpr
  cat(sprintf("  %-10s %6.1f%% %6.1f%% %8s  %s\n", g, 100*sens, 100*(1-fpr),
      if (is.infinite(lr)) "Inf" else sprintf("%.1f", lr), strength(lr)))
}

cat("\nRead the -4 row of HELD OUT: if it still reaches moderate (LR+ >= 4.33),\nthe production threshold is unchanged by the overlap and the finding is a\ndisclosure item rather than a correction.\n")
