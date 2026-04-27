## analyses/01_pull_varibench.R  — Task P2-T1.1
## Download VariBench PON-PS D2 held-out clinical truth VCFs and translate
## ClinVar Variation IDs to HGVS protein notation via NCBI E-utilities esummary.
##
## ID-type decision (tested 2026-04-26):
##   VCF column 3 = ClinVar VARIATION ID (e.g., 225696 → VCV000225696).
##   VCF INFO ALLELEID = ClinVar ALLELE ID (e.g., 227511) — a different concept.
##   NCBI esummary?db=clinvar accepts the VARIATION ID (column 3), NOT ALLELEID.
##   Confirmed: id=225696 returns title="NM_001170535.3(ATAD3A):c.1582C>T (p.Arg528Trp)".
##   ALLELEID is retained in output only for traceability.
##
## Run from project root:
##   Rscript analyses/01_pull_varibench.R
##
## Outputs (all gitignored except manifests):
##   analyses/raw/varibench/phenotype_D2_clinvar_pathogenic.vcf.gz
##   analyses/raw/varibench/phenotype_D2_clinvar_benign.vcf.gz
##   analyses/raw/varibench/clinvar_protein_changes/batch_NNNNN.json  (cached)
##   analyses/derived/varibench_canonical.tsv
##   analyses/manifests/<timestamp>__varibench_*.json

suppressMessages({
  library(httr2)
  library(dplyr)
  library(readr)
  library(purrr)
  library(jsonlite)
  library(stringr)
})

source("analyses/lib/manifest.R")

OUT_RAW      <- "analyses/raw/varibench"
OUT_DERIVED  <- "analyses/derived"
OUT_MANIFEST <- "analyses/manifests"
EUTILS_CACHE <- "analyses/raw/varibench/clinvar_protein_changes"

dir.create(OUT_RAW,      recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DERIVED,  recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MANIFEST, recursive = TRUE, showWarnings = FALSE)
dir.create(EUTILS_CACHE, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# Step A: Download both VCFs (resumable; skip if already cached)
# ---------------------------------------------------------------------------
download_vcf <- function(label, url) {
  out_file <- file.path(OUT_RAW, basename(url))
  if (file.exists(out_file)) {
    cat("CACHED:", basename(url), "\n")
    return(out_file)
  }
  cat("DOWNLOADING:", url, "\n")
  req  <- request(url) |>
    req_user_agent("VarViz-Pass2/1.0 (mailto:agasthyametpally5@gmail.com)") |>
    req_timeout(180)
  resp <- req_perform(req)
  raw  <- resp_body_raw(resp)
  writeBin(raw, out_file)
  write_manifest(OUT_MANIFEST,
                 sprintf("varibench_%s_vcf", label),
                 url, list(), raw, resp_status(resp))
  out_file
}

vcf_paths <- c(
  pathogenic = download_vcf(
    "pathogenic",
    "https://structure.bmc.lu.se/VariBench/data/ponps/phenotype_D2_clinvar_pathogenic.vcf.gz"
  ),
  benign = download_vcf(
    "benign",
    "https://structure.bmc.lu.se/VariBench/data/ponps/phenotype_D2_clinvar_benign.vcf.gz"
  )
)

# ---------------------------------------------------------------------------
# Step B: Parse VCFs
#   col 3 = ClinVar Variation ID  (used for esummary lookup)
#   INFO   = ALLELEID (stored for traceability), GENEINFO, CLNHGVS
# ---------------------------------------------------------------------------
parse_vcf <- function(path, label_collapsed) {
  con  <- gzfile(path, "rt")
  on.exit(close(con))
  rows <- readLines(con)
  rows <- rows[!startsWith(rows, "#")]

  fields <- strsplit(rows, "\t", fixed = TRUE)

  df <- tibble(
    chrom      = vapply(fields, `[`, character(1), 1),
    pos        = vapply(fields, `[`, character(1), 2),
    clinvar_id = vapply(fields, `[`, character(1), 3),   # Variation ID
    ref        = vapply(fields, `[`, character(1), 4),
    alt        = vapply(fields, `[`, character(1), 5),
    info       = vapply(fields, `[`, character(1), 8)
  )

  # Filter to missense only
  df <- df |> filter(grepl("missense_variant", info, fixed = TRUE))

  # Extract INFO fields
  df$gene             <- str_match(df$info, "GENEINFO=([^:|;]+)")[, 2]
  df$hgvs_g           <- str_match(df$info, "CLNHGVS=([^;]+)")[, 2]
  df$clinvar_alleleid <- str_match(df$info, "ALLELEID=([0-9]+)")[, 2]
  df$label            <- label_collapsed
  df$source           <- "VariBench"
  df$set              <- "PON-PS_D2"

  df |> select(-info)
}

vb <- bind_rows(
  parse_vcf(vcf_paths["pathogenic"], "Pathogenic"),
  parse_vcf(vcf_paths["benign"],     "Benign")
)

cat("Parsed:", nrow(vb), "missense rows\n")
cat("Unique Variation IDs:", length(unique(vb$clinvar_id)), "\n")

# ---------------------------------------------------------------------------
# Step C: Translate ClinVar Variation IDs → HGVS protein via NCBI esummary
#   Endpoint: esummary.fcgi?db=clinvar&id=<csv_ids>&retmode=json
#   Batch size: 200 IDs (safe under NCBI's ~500 limit)
#   Rate limit: Sys.sleep(0.34) ≈ 2.9 req/sec (no API key → max 3 req/sec)
#   Per-batch JSON cached to EUTILS_CACHE/batch_NNNNN.json for resumability
# ---------------------------------------------------------------------------

ESUM_URL <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

`%||%` <- function(a, b) if (is.null(a)) b else a

fetch_batch <- function(ids, batch_idx) {
  cache_file <- file.path(EUTILS_CACHE, sprintf("batch_%05d.json", batch_idx))
  if (file.exists(cache_file)) {
    return(jsonlite::fromJSON(cache_file, simplifyVector = FALSE))
  }
  Sys.sleep(0.34)  # NCBI rate limit: 3 req/sec without API key
  req <- request(ESUM_URL) |>
    req_url_query(db = "clinvar", id = paste(ids, collapse = ","), retmode = "json") |>
    req_user_agent("VarViz-Pass2/1.0 (mailto:agasthyametpally5@gmail.com)") |>
    req_timeout(120) |>
    req_retry(max_tries = 3, backoff = function(i) 2^i)
  resp <- req_perform(req)
  raw  <- resp_body_raw(resp)
  writeBin(raw, cache_file)
  write_manifest(
    OUT_MANIFEST,
    sprintf("clinvar_esummary_batch%05d", batch_idx),
    ESUM_URL,
    list(db = "clinvar", n_ids = length(ids), retmode = "json"),
    raw,
    resp_status(resp)
  )
  jsonlite::fromJSON(cache_file, simplifyVector = FALSE)
}

extract_p_notation <- function(title) {
  # Title example: "NM_000059.4(BRCA2):c.5946delT (p.Ser1982Argfs*22)"
  # Extract the p.Xxx part from within parentheses
  m <- str_match(title, "\\(p\\.([A-Za-z0-9*?=]+)\\)")
  if (is.na(m[1, 1])) NA_character_ else paste0("p.", m[, 2])
}

unique_ids <- unique(na.omit(vb$clinvar_id))
cat("Translating", length(unique_ids), "ClinVar Variation IDs in batches of 200\n")
batches <- split(unique_ids, ceiling(seq_along(unique_ids) / 200L))

p_lookup <- list()
for (i in seq_along(batches)) {
  if (i %% 10 == 0 || i == 1 || i == length(batches)) {
    cat(sprintf("[%3d/%d] batch with %d ids\n", i, length(batches), length(batches[[i]])))
  }
  resp_data <- tryCatch(
    fetch_batch(batches[[i]], i),
    error = function(e) {
      cat("  ERR batch", i, ":", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(resp_data)) next

  uids <- resp_data$result$uids
  if (is.null(uids)) next

  for (uid in uids) {
    rec   <- resp_data$result[[uid]]
    if (is.null(rec)) next
    title <- rec$title %||%
             tryCatch(rec$variation_set[[1]]$variation_name, error = function(e) "") %||%
             ""
    p_lookup[[as.character(uid)]] <- extract_p_notation(title)
  }
}

cat(sprintf("p_lookup populated for %d unique IDs\n", length(p_lookup)))

# Map back to the full data frame (use clinvar_id — Variation ID)
vb$p_notation <- unlist(p_lookup[vb$clinvar_id], use.names = FALSE)

# Coerce to character; NULL entries become "NULL" string after unlist
vb$p_notation <- as.character(vb$p_notation)
vb$p_notation[vb$p_notation == "NULL"] <- NA_character_

unmatched <- sum(is.na(vb$p_notation))
cat(sprintf(
  "Translated: %d / %d (%.1f%% unmatched)\n",
  sum(!is.na(vb$p_notation)), nrow(vb), 100 * unmatched / nrow(vb)
))

# ---------------------------------------------------------------------------
# Step D: Write canonical TSV
# ---------------------------------------------------------------------------
canonical <- vb |>
  filter(!is.na(p_notation), !is.na(gene)) |>
  select(gene, p_notation, hgvs_g, clinvar_alleleid, label, source, set)

write_tsv(canonical, file.path(OUT_DERIVED, "varibench_canonical.tsv"))
cat("Wrote", nrow(canonical), "canonical rows.\n")
cat("Label breakdown:\n")
print(table(canonical$label))
cat("Top genes:\n")
print(head(sort(table(canonical$gene), decreasing = TRUE), 10))
