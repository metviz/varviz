#!/usr/bin/env Rscript
# Genome-wide ClinVar robustness: MDS LR+ tiers on 2-star missense P/LP vs B/LB.
suppressMessages({library(data.table)})
source("analyses/lib/pssm_lookup.R")
tbl <- pssm_table_load("data/pfam_pssm_human.rds")
deq <- function(q) (as.numeric(q)+127)/254*(attr(tbl,"Q_MAX") %||% 0) # placeholder
`%||%` <- function(a,b) if(is.null(a)) b else a
# dequant constants from 22_mds_benchmark.R
Q_MIN <- -13.2103; Q_MAX <- 3.5686
deq <- function(q) (as.numeric(q)+127)/254*(Q_MAX-Q_MIN)+Q_MIN

ae  <- fread("analyses/tmp/mane/uniprot_acc_entry.tsv")
entry_of <- setNames(ae$entry, ae$accession)
gm  <- fread("analyses/derived/gene_mane.tsv")
acc_of <- setNames(gm$uniprot_id, gm$gene_name)
cv <- fread("analyses/tmp/clinvar/clinvar_missense_2star.tsv")
cv[, entry := entry_of[acc_of[gene]]]
cv <- cv[!is.na(entry) & nzchar(entry)]
cv[, vid := .I]

mp <- as.data.table(tbl$map); setkeyv(mp, c("id","residue"))
hits <- mp[cv, on=.(id=entry, residue=pos), nomatch=0L, allow.cartesian=TRUE]
hits <- hits[!is.na(family) & nzchar(family)]
hits[, delta := {
  fam <- as.character(.BY$family)
  M <- if (fam %in% names(tbl$pssm)) tbl$pssm[[fam]] else NULL
  if (is.null(M)) rep(NA_real_, .N) else {
    mc <- match(mut, colnames(M)); wc <- match(wt, colnames(M))
    ok <- column <= nrow(M) & !is.na(mc) & !is.na(wc); d <- rep(NA_real_, .N)
    d[ok] <- deq(M[cbind(column[ok], mc[ok])]) - deq(M[cbind(column[ok], wc[ok])]); d
  }
}, by=family]
per <- hits[!is.na(delta), .(mds=min(delta)), by=.(vid,cls)]
cat(sprintf("MDS scored: %d (P=%d B=%d)\n", nrow(per), sum(per$cls==1), sum(per$cls==0)))
P <- per[cls==1]$mds; B <- per[cls==0]$mds
wil <- function(k,n){p=k/n;z=1.96;d=1+z^2/n;c=(p+z^2/2/n)/d;h=z*sqrt(p*(1-p)/n+z^2/4/n^2)/d;c(c-h,c+h)}
band <- function(lr) if(!is.finite(lr)) "Strong+" else if(lr<4.3)"Supporting" else if(lr<18.7)"Moderate" else "Strong"
cat(sprintf("%-9s %6s %6s %8s %-13s %s\n","thresh","sens","spec","LR+","LR+ 95%CI","band"))
for(t in c(-4,-6,-8,-10,-12)){
  tp=sum(P<=t);fp=sum(B<=t);sens=tp/length(P);fpr=fp/length(B)
  lr=if(fpr>0)sens/fpr else Inf
  sci=wil(tp,length(P));fci=wil(fp,length(B))
  lo=if(fci[2]>0)sci[1]/fci[2] else NA; hi=if(fci[1]>0)sci[2]/fci[1] else Inf
  cat(sprintf("MDS<=%-4d %6.3f %6.3f %8.1f %5.1f-%-6.1f %s\n",t,sens,1-fpr,lr,lo,hi,band(lr)))
}
saveRDS(per, "analyses/tmp/clinvar/clinvar_mds.rds")
cat("done\n")
