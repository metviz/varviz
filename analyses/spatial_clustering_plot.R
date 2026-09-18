# Visual prototype: permutation-significant ClinVar hotspot (gnomAD null) vs CCRS,
# with the saturating +/-15 window and ClinVar P/LP density, along the residue axis.
# ISOLATION: prototype only; sources server.R for data parity, writes PNGs. No app edits.

suppressMessages({ library(dplyr); library(ggplot2) })
source("server.R")

set.seed(42)
W <- 15L; B <- 2000L; ALPHA <- 0.05
# gnomAD-null cutoff = VarViz Max AF (Whiffin/Ware), prevalence-derived. Defaults
# below match the harness (1-in-2000, hetA 0.5, hetG 1.0, pen 1.0, monoallelic).
PREV_1_IN_N <- 2000; HETA <- 0.5; HETG <- 1.0; PEN <- 1.0; INH <- "monoallelic"
AF_RARE <- if (INH == "monoallelic") { 0.5*(1/PREV_1_IN_N)*HETA*HETG*(1/PEN) } else { sqrt(1/PREV_1_IN_N)*HETA*sqrt(HETG)*(1/sqrt(PEN)) }
cat(sprintf("[null] Max AF cutoff = %.3g (prevalence 1 in %d)\n", AF_RARE, PREV_1_IN_N))
GENES <- { a <- commandArgs(trailingOnly = TRUE); if (length(a)) a else c("TP53","KCNQ1") }
OUTDIR <- "analyses/weekly_test/plots"; dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

local_density <- function(pts, L, W) {
  h <- tabulate(pmin(pmax(round(pts),1L),L), nbins=L); cs <- c(0, cumsum(h))
  lo <- pmax(1L,(1:L)-W); hi <- pmin(L,(1:L)+W); cs[hi+1L]-cs[lo]
}

plot_gene <- function(gene) {
  cv <- extract_clinvar(gene, "path"); gd <- extract_gnomad(gene)
  ccr <- tryCatch(extract_ccrs(gene, NULL), error=function(e) NULL)
  grow <- gene_data[gene_data$gene_name==gene,]
  L <- suppressWarnings(as.integer(grow$mane_prot_len[1]))
  if (is.na(L)||L<=0) L <- as.integer(max(c(cv$prot_pos, gd$prot_pos), na.rm=TRUE))
  cv_pos <- cv$prot_pos[!is.na(cv$prot_pos) & cv$prot_pos>=1 & cv$prot_pos<=L]
  bg <- gd$prot_pos[!is.na(gd$prot_pos) & gd$gnomad_allele_freq<AF_RARE]; bg <- bg[bg>=1 & bg<=L]
  n <- length(cv_pos)
  obs <- local_density(cv_pos, L, W)
  ge <- integer(L); for (b in seq_len(B)) ge <- ge + (local_density(sample(bg, n, replace=TRUE), L, W) >= obs)
  pval <- ge / B; sig <- pval < (ALPHA/L)
  win_fire <- obs >= 3L
  ccrs_pos <- if (!is.null(ccr) && nrow(ccr)>0) ccr$prot_pos[ccr$prot_pos>=1 & ccr$prot_pos<=L] else integer(0)

  # everything on one continuous y axis: density scaled into band [2.4, 6], tracks as rows below
  dmax <- max(obs, 1)
  df_den <- data.frame(pos=1:L, y=2.4 + 3.6*obs/dmax, sig=sig)
  yb <- c("+/-15 window (fires>=3)"=1.7, "Permutation peak (gnomAD null)"=1.0, "CCRS >=90th"=0.3)
  cols <- c("+/-15 window (fires>=3)"="#c94b4b",
            "Permutation peak (gnomAD null)"="#2c7fb8","CCRS >=90th"="#41ab5d")
  tracks <- rbind(
    data.frame(pos=which(win_fire), track="+/-15 window (fires>=3)"),
    data.frame(pos=which(sig),      track="Permutation peak (gnomAD null)"),
    if (length(ccrs_pos)) data.frame(pos=ccrs_pos, track="CCRS >=90th") else NULL
  )
  tracks$y <- yb[tracks$track]
  labs_df <- data.frame(x=1, y=unname(yb), lab=names(yb))

  p <- ggplot() +
    geom_ribbon(data=df_den, aes(pos, ymin=2.4, ymax=y), fill="grey82") +
    geom_ribbon(data=subset(df_den, sig), aes(pos, ymin=2.4, ymax=y), fill="#2c7fb8", alpha=.7) +
    geom_point(data=tracks, aes(pos, y, colour=track), shape=15, size=1.4, show.legend=FALSE) +
    geom_hline(yintercept=2.25, colour="grey60", linewidth=.3) +
    geom_text(data=labs_df, aes(x, y, label=lab, colour=lab), hjust=0, vjust=-0.8, size=2.7, show.legend=FALSE) +
    scale_colour_manual(values=cols) +
    labs(title=sprintf("%s (L=%d, ClinVar P/LP n=%d)  |  grey=ClinVar +/-15 density, blue=permutation-significant hotspot (gnomAD null)", gene, L, n),
         x="residue", y=NULL) +
    xlim(1,L) + theme_minimal(base_size=10) +
    theme(plot.title=element_text(size=8.5), axis.text.y=element_blank(),
          panel.grid.major.y=element_blank(), panel.grid.minor.y=element_blank())

  out <- file.path(OUTDIR, sprintf("hotspot_%s.png", gene))
  ggsave(out, p, width=11, height=4, dpi=130)
  cat(sprintf("%s: window fires %d (%.0f%%), permutation peaks %d (%.0f%%), CCRS %d -> %s\n",
              gene, sum(win_fire), 100*mean(win_fire), sum(sig), 100*mean(sig), length(ccrs_pos), out))
}
invisible(lapply(GENES, plot_gene))
