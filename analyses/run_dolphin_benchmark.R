#!/usr/bin/env Rscript
# run_dolphin_benchmark.R — drive DOLPHIN precompute across the 14-gene
# benchmark universe, processing genes smallest-first so the early-fail
# diagnostics (timeouts, rate-limit, schema drift) surface on the cheapest
# work first.
#
# Per-gene output: analyses/raw/dolphin/by_gene/{GENE}.tsv
# Resumable: existing TSVs are read at startup and only un-done variants
# are queried. Safe to kill and restart.
#
# Rate limit: dolphin.mmg-gbit.eu enforces 60 req/min; we use 1.05 s/req
# in the bulk helper to stay under the limit with margin.
#
# Wall-time estimate (rate-limited):
#   ~25 hours total across 90,679 remaining variants (CASR already done).
#
# Usage:
#   Rscript analyses/run_dolphin_benchmark.R [GENE1 GENE2 ...]
#   (no args = process all 13 remaining genes in size-ascending order)

suppressPackageStartupMessages({
  source("analyses/lib/dolphin_bulk.R", local = FALSE)
})

UNIVERSE_TSV <- "analyses/derived/varviz_classifications.tsv"
CACHE_DIR    <- "analyses/raw/dolphin/by_gene"
LOG_PATH     <- "analyses/raw/dolphin/run_log.tsv"

# Hardcoded canonical Ensembl transcript IDs for the 14 benchmark genes.
# Bypasses Ensembl REST (rest.ensembl.org), which returned HTTP 500 during the
# first long-running pass and blocked 6 of 13 genes. These ENSTs were captured
# from successful runs in earlier sessions; they match the canonical / MANE
# Select transcripts that the live VarViz app resolves via Ensembl.
GENE_TO_ENST <- c(
  BAP1     = "ENST00000460680",
  BRCA1    = "ENST00000357654",
  CASR     = "ENST00000639785",
  DDR2     = "ENST00000367921",
  GCK      = "ENST00000403799",
  KCNH2    = "ENST00000262186",
  KCNQ1    = "ENST00000155840",
  KRAS     = "ENST00000311936",
  LDLR     = "ENST00000558518",
  NUDT15   = "ENST00000258662",
  PTEN     = "ENST00000371953",
  SLC13A5  = "ENST00000433363",
  SNCA     = "ENST00000394991",
  TP53     = "ENST00000269305",
  TSHR     = "ENST00000298171"
)

dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

# Read the full universe and aggregate variants per gene.
cat("[driver] reading", UNIVERSE_TSV, "...\n")
uni <- read.table(UNIVERSE_TSV, sep = "\t", header = TRUE,
                  stringsAsFactors = FALSE, quote = "", comment.char = "")
stopifnot(all(c("gene", "p_notation") %in% colnames(uni)))

# Gene size table (variants per gene), ascending.
gene_counts <- as.data.frame(table(uni$gene), stringsAsFactors = FALSE)
colnames(gene_counts) <- c("gene", "n")
gene_counts <- gene_counts[order(gene_counts$n), ]

# Default processing order: every gene except CASR (already done in an
# earlier checkpoint). CLI args override.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 0) {
  process <- args
} else {
  process <- gene_counts$gene[gene_counts$gene != "CASR"]
}
cat("[driver] processing", length(process), "gene(s):\n")
print(gene_counts[gene_counts$gene %in% process, , drop = FALSE], row.names = FALSE)

log_event <- function(gene, status, n_total, n_done, elapsed_s, note = "") {
  row <- data.frame(
    ts          = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"),
    gene        = gene,
    status      = status,
    n_total     = n_total,
    n_done      = n_done,
    elapsed_s   = round(elapsed_s, 2),
    note        = note,
    stringsAsFactors = FALSE
  )
  utils::write.table(row, LOG_PATH, sep = "\t", row.names = FALSE,
                     quote = FALSE, na = "",
                     append = file.exists(LOG_PATH),
                     col.names = !file.exists(LOG_PATH))
}

t_start_total <- Sys.time()
for (gene in process) {
  vars <- unique(uni$p_notation[uni$gene == gene])
  vars <- vars[nchar(vars) > 0]
  n <- length(vars)
  cat(sprintf("\n========== %s (n=%d) ==========\n", gene, n))

  t_gene <- Sys.time()
  ok <- tryCatch({
    enst <- GENE_TO_ENST[gene]
    if (is.na(enst) || nchar(enst) == 0) enst <- NULL  # falls back to Ensembl REST
    df <- fetch_dolphin_gene(gene, vars,
                              ensembl = enst,
                              cache_dir = CACHE_DIR,
                              rate_sleep_s = 1.5,
                              flush_every = 50,
                              verbose = TRUE)
    nrow(df)
  }, error = function(e) {
    msg <- paste0("ERROR: ", conditionMessage(e))
    cat("[driver]", gene, "failed:", msg, "\n")
    log_event(gene, "ERROR", n, NA, as.numeric(difftime(Sys.time(), t_gene, units = "secs")), msg)
    -1L
  })

  if (is.numeric(ok) && ok >= 0) {
    elapsed <- as.numeric(difftime(Sys.time(), t_gene, units = "secs"))
    log_event(gene, "DONE", n, ok, elapsed)
    cat(sprintf("[driver] %s done in %.1f s (%.2f s/variant)\n",
                gene, elapsed, if (ok > 0) elapsed / ok else 0))
  }
}

t_total <- as.numeric(difftime(Sys.time(), t_start_total, units = "secs"))
cat(sprintf("\n[driver] ALL DONE in %.1f s (%.2f h)\n",
            t_total, t_total / 3600))
