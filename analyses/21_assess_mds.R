#!/usr/bin/env Rscript
# 21_assess_mds.R — does the Missense Disfavour Score separate pathogenic from
# benign, and where should the PM1 threshold sit?
#
# Until now MDS was assessed only two ways: it reproduces DOLPHIN's deltas on 5
# variants, and the -4 PM1 cutoff was eyeballed. This runs it against a labelled
# benchmark (VariBench: 6,057 pathogenic / 10,491 benign missense) restricted to
# the variants MDS can actually score (inside a Pfam domain), and reports:
#   - score distribution by label
#   - discrimination (AUC, Mann-Whitney)
#   - LR+ at each candidate threshold, mapped to the ClinGen/Pejaver evidence
#     strengths (supporting >=2.08, moderate >=4.33, strong >=18.7)
# so the threshold is picked from data, not by eye.
#
# Usage:  Rscript analyses/21_assess_mds.R
suppressMessages({library(httr2); library(data.table)})
source("analyses/lib/pssm_lookup.R")

# ── Load table + labels + the gene->accession->entry bridges ───────────────
tbl <- pssm_table_load("data/pfam_pssm_human.rds")
vb  <- fread("analyses/derived/varibench_canonical.tsv")
vb  <- vb[label %in% c("Pathogenic", "Benign")]
load("data/VarViz.RData")                                   # gene_data
acc_of <- setNames(gene_data$uniprot_id, gene_data$gene_name)

# accession -> UniProt entry name, one cached reviewed-human stream.
ENTRY_CACHE <- "analyses/tmp/mane/uniprot_acc_entry.tsv"
if (file.exists(ENTRY_CACHE)) {
  ae <- fread(ENTRY_CACHE)
} else {
  message("[assess] streaming accession->entry (one-time)")
  rows <- list(); cursor <- NULL
  repeat {
    rq <- request("https://rest.uniprot.org/uniprotkb/search") |>
      req_url_query(query = "reviewed:true AND organism_id:9606",
                    fields = "accession,id", format = "tsv", size = 500)
    if (!is.null(cursor)) rq <- req_url_query(rq, cursor = cursor)
    rp <- req_perform(rq); df <- fread(text = resp_body_string(rp))
    if (!nrow(df)) break; rows[[length(rows)+1]] <- df
    lk <- resp_headers(rp)[["link"]]; if (is.null(lk) || !grepl("next", lk)) break
    cursor <- sub(".*[?&]cursor=([^&>]+).*", "\\1", lk)
  }
  ae <- rbindlist(rows); setnames(ae, c("accession", "entry"))
  fwrite(ae, ENTRY_CACHE, sep = "\t")
}
entry_of <- setNames(ae$entry, ae$accession)

# ── Parse p.Arg528Trp -> pos 528, wt R, mut W ──────────────────────────────
aa3 <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",Gln="Q",Glu="E",Gly="G",His="H",
         Ile="I",Leu="L",Lys="K",Met="M",Phe="F",Pro="P",Ser="S",Thr="T",Trp="W",
         Tyr="Y",Val="V")
p <- sub("^p\\.", "", vb$p_notation)
m <- regmatches(p, regexec("^([A-Za-z]{3})([0-9]+)([A-Za-z]{3})$", p))
vb[, wt  := vapply(m, function(x) if (length(x)==4) unname(aa3[x[2]]) else NA_character_, "")]
vb[, pos := vapply(m, function(x) if (length(x)==4) as.integer(x[3]) else NA_integer_, integer(1))]
vb[, mut := vapply(m, function(x) if (length(x)==4) unname(aa3[x[4]]) else NA_character_, "")]
vb[, entry := entry_of[acc_of[gene]]]

# ── Score every scorable variant (keyed join, not per-variant map scan) ────
vb <- vb[!is.na(entry) & !is.na(pos) & !is.na(wt) & !is.na(mut)]
vb[, vid := .I]
message(sprintf("[assess] scoring %s variants with a resolvable entry", format(nrow(vb), big.mark=",")))

mp <- as.data.table(tbl$map); setkeyv(mp, c("id", "residue"))
Q_MIN <- tbl$quant[["min"]]; Q_MAX <- tbl$quant[["max"]]
deq <- function(q) (as.numeric(q) + 127) / 254 * (Q_MAX - Q_MIN) + Q_MIN

# variant -> every (family, column) its residue maps to
hits <- mp[vb[, .(id = entry, residue = pos, vid, wt, mut)],
           on = .(id, residue), nomatch = 0L, allow.cartesian = TRUE]
# delta per hit, vectorised within each family's matrix
hits[, delta := {
  M <- tbl$pssm[[.BY$family]]
  if (is.null(M)) rep(NA_real_, .N) else {
    mc <- match(mut, colnames(M)); wc <- match(wt, colnames(M))
    ok <- column <= nrow(M) & !is.na(mc) & !is.na(wc)
    d <- rep(NA_real_, .N)
    d[ok] <- deq(M[cbind(column[ok], mc[ok])]) - deq(M[cbind(column[ok], wc[ok])])
    d
  }
}, by = family]
per <- hits[!is.na(delta), .(mds = min(delta)), by = vid]   # most disfavoured
vb[per, mds := i.mds, on = "vid"]
scored <- vb[!is.na(mds)]      # in a Pfam domain
np <- sum(scored$label=="Pathogenic"); nb <- sum(scored$label=="Benign")

cat("\n===== MDS assessment on VariBench =====\n")
cat(sprintf("labelled missense:            %s\n", format(nrow(vb), big.mark=",")))
cat(sprintf("scorable (in a Pfam domain):  %s (%.1f%%)  [%d path / %d benign]\n",
            format(nrow(scored), big.mark=","), 100*nrow(scored)/nrow(vb), np, nb))

cat("\n-- score distribution (more negative = more disfavoured) --\n")
qs <- c(0.1,0.25,0.5,0.75,0.9)
for (lab in c("Pathogenic","Benign"))
  cat(sprintf("  %-11s median %6.2f   IQR [%6.2f, %6.2f]\n", lab,
      median(scored[label==lab]$mds),
      quantile(scored[label==lab]$mds,.25), quantile(scored[label==lab]$mds,.75)))

# AUC via Mann-Whitney (lower MDS = pathogenic, so flip sign for "score")
s <- -scored$mds; y <- as.integer(scored$label=="Pathogenic")
r <- rank(s); auc <- (sum(r[y==1]) - np*(np+1)/2) / (as.numeric(np)*nb)
cat(sprintf("\n-- discrimination --\n  AUC = %.3f  (pathogenic vs benign)\n", auc))

# LR+ at candidate thresholds; predict pathogenic if MDS <= t
cat("\n-- LR+ by threshold (predict pathogenic if MDS <= t) --\n")
cat(sprintf("  %6s  %7s %7s  %8s  %6s  %s\n","t","sens","spec","LR+","%fires","strength"))
strength <- function(lr) if (is.infinite(lr)||lr>=18.7)"STRONG" else if(lr>=4.33)"moderate" else if(lr>=2.08)"supporting" else "-"
for (t in c(-2,-3,-4,-5,-6,-8)) {
  tp <- sum(scored$label=="Pathogenic" & scored$mds<=t); fp <- sum(scored$label=="Benign" & scored$mds<=t)
  sens <- tp/np; fpr <- fp/nb; lr <- if (fpr==0) Inf else sens/fpr
  cat(sprintf("  %6d  %6.1f%% %6.1f%%  %8s  %5.1f%%  %s\n", t, 100*sens, 100*(1-fpr),
      if (is.infinite(lr)) "Inf" else sprintf("%.1f", lr),
      100*(tp+fp)/nrow(scored), strength(lr)))
}
cat("\nThe current server.R threshold is -4. Read the row for it above to see the\nevidence strength it actually delivers on this benchmark.\n")
