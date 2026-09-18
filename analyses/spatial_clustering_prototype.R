# Prototype: permutation-based ClinVar hotspot PM1 with a gnomAD rare-variant null.
#
# Idea (from SpatialClustering / de Wiel 2017, adapted): instead of firing PM1
# whenever a variant has >=K ClinVar P/LP neighbours within +/-15 residues (a
# fixed density heuristic that SATURATES in densely-curated genes), fire only at
# positions whose local ClinVar density is SIGNIFICANT under a permutation null.
#
# Two nulls compared:
#   uniform  — ClinVar positions permuted uniformly over the protein (naive)
#   gnomad   — ClinVar positions resampled from the observed gnomAD rare-variant
#              positional density (controls for mutability/coverage/mappability;
#              also, since gnomAD rare variants are depleted at functional sites,
#              this sharpens contrast at true hotspots)
#
# Signal set = ClinVar P/LP protein positions (still ClinVar-derived => circular,
# stripped in Pass-Blind, exactly like the current clinvar_hotspot pathway).
#
# ISOLATION: prototype only. Does NOT modify server.R/ui.R. Sources server.R to
# reuse the app's own extract_clinvar / extract_gnomad for data parity.

suppressMessages({ library(dplyr) })
source("server.R")

set.seed(42)
W       <- 15L        # neighbourhood half-window (matches current pathway)
B       <- 2000L      # permutations
AF_RARE <- 1e-4       # gnomAD "rare" cutoff for the null background
ALPHA   <- 0.05
.args   <- commandArgs(trailingOnly = TRUE)
GENES   <- if (length(.args)) .args else c("TP53", "KCNQ1")  # override on CLI

# sliding-window local count of `pts` at every residue 1..L (box filter, O(L))
local_density <- function(pts, L, W) {
  h <- tabulate(pmin(pmax(round(pts), 1L), L), nbins = L)   # per-residue counts
  cs <- c(0, cumsum(h))
  lo <- pmax(1L, (1:L) - W); hi <- pmin(L, (1:L) + W)
  cs[hi + 1L] - cs[lo]
}

analyse_gene <- function(gene) {
  cat(sprintf("\n===== %s =====\n", gene))
  cv <- tryCatch(extract_clinvar(gene, "path"), error = function(e) NULL)
  gd <- tryCatch(extract_gnomad(gene),          error = function(e) NULL)
  grow <- gene_data[gene_data$gene_name == gene, ]
  L <- suppressWarnings(as.integer(grow$mane_prot_len[1]))
  if (is.na(L) || L <= 0) L <- as.integer(max(c(cv$prot_pos, gd$prot_pos), na.rm = TRUE))

  cv_pos <- cv$prot_pos[!is.na(cv$prot_pos)]
  cv_pos <- cv_pos[cv_pos >= 1 & cv_pos <= L]
  null_bg <- gd$prot_pos[!is.na(gd$prot_pos) & gd$gnomad_allele_freq < AF_RARE]
  null_bg <- null_bg[null_bg >= 1 & null_bg <= L]
  n <- length(cv_pos)
  cat(sprintf("L=%d  ClinVar P/LP positions n=%d  gnomAD rare background m=%d\n",
              L, n, length(null_bg)))
  if (n < 3 || length(null_bg) < 10) { cat("  insufficient data — skip\n"); return(invisible()) }

  obs <- local_density(cv_pos, L, W)

  # permutation null density profiles
  perm_ge <- function(sampler) {
    ge <- integer(L)
    for (b in seq_len(B)) {
      d <- local_density(sampler(), L, W)
      ge <- ge + (d >= obs)
    }
    ge / B                      # per-residue empirical p-value P(null>=obs)
  }
  p_unif <- perm_ge(function() sample.int(L, n, replace = TRUE))
  p_gno  <- perm_ge(function() sample(null_bg, n, replace = TRUE))

  bonf <- ALPHA / L
  # current pathway proxy: fires wherever a variant has >=2 ClinVar P/LP in +/-W
  fire_window <- obs >= 2L
  fire_unif   <- p_unif < bonf
  fire_gno    <- p_gno  < bonf

  pct <- function(x) sprintf("%d (%.0f%%)", sum(x), 100 * mean(x))
  cat(sprintf("  PM1 FIRES at residue positions (of %d):\n", L))
  cat(sprintf("    +/-15 window (obs>=2)      : %s\n", pct(fire_window)))
  cat(sprintf("    permutation, uniform null  : %s\n", pct(fire_unif)))
  cat(sprintf("    permutation, gnomAD null   : %s\n", pct(fire_gno)))

  # report contiguous significant peaks under the gnomAD null
  peaks <- which(fire_gno)
  if (length(peaks)) {
    runs <- split(peaks, cumsum(c(1, diff(peaks) != 1)))
    cat("    gnomAD-null significant peak regions:\n")
    for (r in runs) cat(sprintf("      %d-%d (max obs density %d)\n",
                                min(r), max(r), max(obs[r])))
  }
  invisible(data.frame(gene, L, n_clinvar = n, n_gnomad_rare = length(null_bg),
                       fire_window = sum(fire_window), fire_unif = sum(fire_unif),
                       fire_gnomad = sum(fire_gno)))
}

res <- lapply(GENES, analyse_gene)
cat("\n===== summary =====\n")
print(do.call(rbind, res), row.names = FALSE)
