#!/usr/bin/env Rscript
# 13_build_domain_pm1.R — build the offline domain-constraint PM1 table.
#
# Streams the local dbNSFP 4.9a corpus (analyses/tmp/), keeps only residues
# inside an InterPro-annotated domain, and scores each residue by how
# conserved it is relative to its own domain. Output is a compact table
# intended to be bundled alongside ccrs_data in data/VarViz.RData, so the
# criterion needs no network call at runtime.
#
# Usage:
#   Rscript analyses/13_build_domain_pm1.R                 # all chromosomes
#   Rscript analyses/13_build_domain_pm1.R 21 22           # named chromosomes
#   Rscript analyses/13_build_domain_pm1.R --genes CASR,SNCA   # gene subset
#
# dbNSFP columns used (4.9a): 12 aapos, 13 genename, 191 GERP++_RS,
# 195 phyloP100way_vertebrate, 449 Interpro_domain. Verified against the
# header rather than assumed — see check_columns() below.

suppressWarnings(suppressMessages({
  library(data.table)
}))
# Run from the repo root (the dbNSFP paths below are repo-relative anyway).
source("analyses/lib/domain_pm1.R")

TMP     <- "analyses/tmp"
OUT_DIR <- "analyses/derived"
COLS    <- c(aapos = 12L, gene = 13L, gerp = 191L, phylop = 195L, domain = 449L)

# ---- guard: confirm the column indices actually match this dbNSFP build ----
check_columns <- function(f) {
  hdr <- strsplit(readLines(gzfile(f), n = 1L), "\t", fixed = TRUE)[[1]]
  want <- c(aapos = "aapos", gene = "genename", gerp = "GERP++_RS",
            phylop = "phyloP100way_vertebrate", domain = "Interpro_domain")
  bad <- character(0)
  for (k in names(want)) {
    if (!identical(hdr[COLS[[k]]], want[[k]]))
      bad <- c(bad, sprintf("col %d is '%s', expected '%s'",
                            COLS[[k]], hdr[COLS[[k]]], want[[k]]))
  }
  if (length(bad))
    stop("dbNSFP column layout differs from 4.9a:\n  ", paste(bad, collapse = "\n  "))
  invisible(TRUE)
}

# ---- stream one chromosome via awk (orders of magnitude faster than R I/O) --
# Emits gene, aapos, domain, gerp, phylop for domain-resident rows only.
# The domain field is ";"-separated per transcript and frequently ".;.;.;." —
# emitting the first non-"." element is what keeps the row count honest.
extract_chr <- function(f, gene_filter = NULL) {
  awk <- sprintf('
    NR==1 { next }
    {
      d=$%d; t=d; gsub(/[.;[:space:]]/,"",t); if (length(t)==0) next
      split($%d,g,";"); split($%d,a,";")
      gene=g[1]; pos=a[1]
      if (gene=="." || pos==".") next
      %s
      n=split(d,parts,";"); dom=""
      for(i=1;i<=n;i++){ if(parts[i]!="."){ dom=parts[i]; break } }
      if (dom=="") next
      print gene "\t" pos "\t" dom "\t" $%d "\t" $%d
    }',
    COLS[["domain"]], COLS[["gene"]], COLS[["aapos"]],
    if (is.null(gene_filter)) "" else
      sprintf('if (!(gene in KEEP)) next', ""),
    COLS[["gerp"]], COLS[["phylop"]])

  pre <- ""
  if (!is.null(gene_filter)) {
    pre <- sprintf('BEGIN{split("%s",k,","); for(i in k) KEEP[k[i]]=1}',
                   paste(gene_filter, collapse = ","))
  }
  cmd <- sprintf("zcat %s | awk -F'\t' '%s%s'", shQuote(f), pre, awk)
  dt <- tryCatch(
    data.table::fread(cmd = cmd, sep = "\t", header = FALSE, quote = "",
                      col.names = c("gene", "aapos", "domain", "gerp", "phylop"),
                      colClasses = "character", showProgress = FALSE),
    error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0L) return(NULL)
  dt[, aapos := suppressWarnings(as.integer(aapos))]
  dt[, gerp   := suppressWarnings(as.numeric(gerp))]
  dt[, phylop := suppressWarnings(as.numeric(phylop))]
  dt[!is.na(aapos)]
}

# ---- main -----------------------------------------------------------------
args  <- commandArgs(trailingOnly = TRUE)
genes <- NULL
if (length(args) && args[1] == "--genes") {
  genes <- strsplit(args[2], ",", fixed = TRUE)[[1]]
  args  <- args[-(1:2)]
}
chrs <- if (length(args)) args else c(1:22, "X", "Y")

files <- file.path(TMP, sprintf("dbNSFP4.9a_variant.chr%s.gz", chrs))
files <- files[file.exists(files)]
if (!length(files)) stop("no dbNSFP chromosome files found under ", TMP)
check_columns(files[1])

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
message("[domain-pm1] ", length(files), " chromosome file(s)",
        if (!is.null(genes)) paste0(", genes: ", paste(genes, collapse = ",")) else "")

parts <- list()
for (f in files) {
  t0 <- Sys.time()
  dt <- extract_chr(f, genes)
  n  <- if (is.null(dt)) 0L else nrow(dt)
  message(sprintf("  %-38s %9s rows  %5.1fs", basename(f), format(n, big.mark = ","),
                  as.numeric(difftime(Sys.time(), t0, units = "secs"))))
  if (n) parts[[length(parts) + 1L]] <- dt
}
if (!length(parts)) stop("no domain-annotated rows extracted")
raw <- data.table::rbindlist(parts)

# Conservation is a per-position property, so collapse the per-variant rows.
# GERP++ is the primary signal (broader coverage in dbNSFP); phyloP fills gaps.
raw[, cons := fifelse(!is.na(gerp), gerp, phylop)]
pos <- raw[!is.na(cons), .(cons = max(cons, na.rm = TRUE)),
           by = .(gene, domain, aapos)]

scored <- data.table::as.data.table(
  score_domain_table(as.data.frame(pos), min_positions = 20L, threshold = 90))

out <- file.path(OUT_DIR, "domain_pm1_scores.rds")
saveRDS(scored, out, compress = "xz")

message("\n[domain-pm1] summary")
message(sprintf("  positions scored      : %s", format(nrow(scored), big.mark = ",")))
message(sprintf("  distinct genes        : %s", format(uniqueN(scored$gene), big.mark = ",")))
message(sprintf("  distinct domains      : %s", format(uniqueN(scored$domain), big.mark = ",")))
message(sprintf("  scorable (>=20 pos)   : %s (%.1f%%)",
                format(sum(!is.na(scored$pctile)), big.mark = ","),
                100 * mean(!is.na(scored$pctile))))
message(sprintf("  PM1 fires             : %s (%.1f%% of scorable)",
                format(sum(scored$pm1, na.rm = TRUE), big.mark = ","),
                100 * mean(scored$pm1, na.rm = TRUE)))
message(sprintf("  in-RAM                : %s", format(object.size(scored), units = "MB")))
message(sprintf("  on disk (xz)          : %.1f MB", file.size(out) / 1e6))
message("  written to ", out)
