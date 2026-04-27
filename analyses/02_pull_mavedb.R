# Pulls MaveDB DMS score sets for BENCHMARK_GENES into analyses/raw/mavedb/.
# Multi-gene ingest per MF2 (15-gene benchmark scope). Writes manifests and a
# canonical TSV of all variant-level scores normalized to p_notation.
#
# Discovery basis (cached at /tmp/mavedb_genes.json from earlier query of
# https://api.mavedb.org/api/v1/score-sets/mapped-genes): only 9 of the 15
# BENCHMARK_GENES have any mapped MaveDB scoresets. URNs hard-pinned below to
# avoid relying on text search.
#
# Run from project root: `Rscript analyses/02_pull_mavedb.R`
#
# Outputs:
#   analyses/raw/mavedb/<gene>__<urn-suffix>.csv  (raw download — gitignored)
#   analyses/derived/mavedb_canonical.tsv         (normalized union — gitignored)
#   analyses/manifests/<timestamp>__mavedb_<urn>.json  (committed)

suppressMessages({
  library(httr2)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(stringr)
})
source("analyses/lib/manifest.R")

OUT_RAW       <- "analyses/raw/mavedb"
OUT_DERIVED   <- "analyses/derived"
OUT_MANIFEST  <- "analyses/manifests"
META_CACHE    <- "analyses/raw/mavedb/_metadata"

dir.create(OUT_RAW,    recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DERIVED, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MANIFEST, recursive = TRUE, showWarnings = FALSE)
dir.create(META_CACHE,  recursive = TRUE, showWarnings = FALSE)

API_BASE <- "https://api.mavedb.org/api/v1"

# ---------------------------------------------------------------------------
# URN-per-gene pinning (from mapped-genes ∩ BENCHMARK_GENES, 2026-04-26).
# For genes with >3 scoresets (BRCA1=30, KRAS=22) we restrict to the principal
# study family (Findlay 2018 BRCA1 HDR; Bandaru 2024-style KRAS drug screens).
# ---------------------------------------------------------------------------
GENE_TO_URNS <- list(
  TP53   = c("urn:mavedb:00000059-a-1"),                      # Kotler 2018 (mapped)
  BRCA1  = c("urn:mavedb:00000097-a-1", "urn:mavedb:00000097-b-1",
             "urn:mavedb:00000097-c-1", "urn:mavedb:00000097-d-1",
             "urn:mavedb:00000097-e-1", "urn:mavedb:00000097-f-1",
             "urn:mavedb:00000097-g-1", "urn:mavedb:00000097-h-1",
             "urn:mavedb:00000097-i-1", "urn:mavedb:00000097-j-1",
             "urn:mavedb:00000097-k-1", "urn:mavedb:00000097-l-1",
             "urn:mavedb:00000097-m-1", "urn:mavedb:00000097-n-1",
             "urn:mavedb:00000097-o-1", "urn:mavedb:00000097-p-1",
             "urn:mavedb:00000097-q-1", "urn:mavedb:00000097-r-1",
             "urn:mavedb:00000097-s-1", "urn:mavedb:00000097-t-1",
             "urn:mavedb:00000097-u-1", "urn:mavedb:00000097-v-1",
             "urn:mavedb:00000097-w-1", "urn:mavedb:00000097-x-1",
             "urn:mavedb:00000097-y-1", "urn:mavedb:00000097-z-1",
             "urn:mavedb:00000081-a-1", "urn:mavedb:00000081-a-2",
             "urn:mavedb:00000093-a-1", "urn:mavedb:00000003-a-2"),
  LDLR   = c("urn:mavedb:00000023-a-2"),                      # Matreyek 2021
  PTEN   = c("urn:mavedb:00000101-a-1", "urn:mavedb:00000102-0-1"),  # Mighell + Matreyek
  KRAS   = c("urn:mavedb:00000115-a-1", "urn:mavedb:00000115-a-2",
             "urn:mavedb:00000115-a-3", "urn:mavedb:00000115-a-4",
             "urn:mavedb:00000115-a-5", "urn:mavedb:00000115-a-6"),  # First 6 of 22 (drug-response screens; rest similar)
  KCNQ1  = c("urn:mavedb:00000094-a-14", "urn:mavedb:00000094-a-15"),
  KCNH2  = c("urn:mavedb:00001216-a-1"),
  GCK    = c("urn:mavedb:00000096-a-1", "urn:mavedb:00000096-b-1"),
  NUDT15 = c("urn:mavedb:00000055-a-1", "urn:mavedb:00000055-b-1",
             "urn:mavedb:00000055-0-1")
)

cat(sprintf("[mavedb] Will fetch %d scoresets across %d genes\n",
            sum(lengths(GENE_TO_URNS)), length(GENE_TO_URNS)))

# ---------------------------------------------------------------------------
# fetch_metadata(urn) -- cache + return list with title, pmid_list, n_variants
# ---------------------------------------------------------------------------
fetch_metadata <- function(urn) {
  meta_file <- file.path(META_CACHE, paste0(URLencode(urn, reserved = TRUE), ".json"))
  if (file.exists(meta_file) && file.size(meta_file) > 0) {
    return(fromJSON(meta_file, simplifyVector = FALSE))
  }
  url <- sprintf("%s/score-sets/%s", API_BASE, URLencode(urn, reserved = TRUE))
  resp <- request(url) |>
    req_user_agent("VarViz-Pass2/1.0 (research)") |>
    req_timeout(60) |> req_retry(max_tries = 3) |> req_perform()
  if (resp_status(resp) != 200) {
    warning(sprintf("[mavedb] metadata HTTP %d for %s", resp_status(resp), urn))
    return(NULL)
  }
  raw <- resp_body_raw(resp)
  writeBin(raw, meta_file)
  fromJSON(meta_file, simplifyVector = FALSE)
}

# ---------------------------------------------------------------------------
# fetch_scores(urn) -- download scores CSV; cache + manifest
# ---------------------------------------------------------------------------
fetch_scores <- function(urn, gene) {
  urn_suffix <- gsub("^urn:mavedb:", "", urn)
  out_file   <- file.path(OUT_RAW, sprintf("%s__%s.csv", gene, urn_suffix))
  if (file.exists(out_file) && file.size(out_file) > 0) {
    cat(sprintf("[mavedb] CACHED %s\n", basename(out_file)))
    return(out_file)
  }
  url <- sprintf("%s/score-sets/%s/scores", API_BASE, URLencode(urn, reserved = TRUE))
  cat(sprintf("[mavedb] GET %s\n", url))
  resp <- request(url) |>
    req_user_agent("VarViz-Pass2/1.0 (research)") |>
    req_timeout(120) |> req_retry(max_tries = 3) |> req_perform()
  if (resp_status(resp) != 200) {
    warning(sprintf("[mavedb] scores HTTP %d for %s", resp_status(resp), urn))
    return(NULL)
  }
  raw <- resp_body_raw(resp)
  writeBin(raw, out_file)
  write_manifest(
    out_dir       = OUT_MANIFEST,
    label         = sprintf("mavedb_%s", urn_suffix),
    url           = url,
    params        = list(urn = urn, gene = gene),
    response_raw  = raw,
    http_status   = resp_status(resp)
  )
  out_file
}

# ---------------------------------------------------------------------------
# AA three-letter -> one-letter helpers (shared with VariBench)
# ---------------------------------------------------------------------------
AA_3to1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G",
  His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S",
  Thr="T", Trp="W", Tyr="Y", Val="V", Ter="*", Sec="U", Pyl="O"
)
normalize_p_notation <- function(hgvs_p) {
  if (is.na(hgvs_p) || !nzchar(hgvs_p)) return(NA_character_)
  s <- sub("^[^:]*:", "", hgvs_p)
  s <- gsub("[()]", "", s)
  # Already one-letter form (p.R175H)?
  m1 <- regmatches(s, regexec("^p\\.([A-Z*])([0-9]+)([A-Z*])$", s))[[1]]
  if (length(m1) == 4) return(paste0("p.", m1[2], m1[3], m1[4]))
  # Three-letter form (p.Arg175His)
  m3 <- regmatches(s, regexec("^p\\.([A-Z][a-z]{2}|\\*)([0-9]+)([A-Z][a-z]{2}|\\*)$", s))[[1]]
  if (length(m3) != 4) return(NA_character_)
  ref1 <- if (m3[2] == "*") "*" else AA_3to1[m3[2]]
  alt1 <- if (m3[4] == "*") "*" else AA_3to1[m3[4]]
  if (is.na(ref1) || is.na(alt1)) return(NA_character_)
  paste0("p.", ref1, m3[3], alt1)
}
normalize_p_notation_v <- Vectorize(normalize_p_notation, USE.NAMES = FALSE)

# ---------------------------------------------------------------------------
# normalize_scoreset(csv_path, urn, gene) -> canonical tibble
# ---------------------------------------------------------------------------
normalize_scoreset <- function(csv_path, urn, gene) {
  # Note: do NOT pass comment="#" — MaveDB accessions contain '#' as a separator
  # (e.g., "urn:mavedb:00000059-a-1#1") and read_csv would truncate them.
  df <- tryCatch(read_csv(csv_path, show_col_types = FALSE),
                 error = function(e) NULL)
  if (is.null(df) || nrow(df) == 0) return(NULL)
  if (!"hgvs_pro" %in% colnames(df)) {
    warning(sprintf("[mavedb] no hgvs_pro in %s (cols: %s)",
                    basename(csv_path), paste(colnames(df), collapse=",")))
    return(NULL)
  }
  if (!"score" %in% colnames(df)) {
    warning(sprintf("[mavedb] no score col in %s", basename(csv_path)))
    return(NULL)
  }
  df |>
    filter(grepl("^p\\.", hgvs_pro), !is.na(score)) |>
    transmute(
      gene       = gene,
      p_notation = normalize_p_notation_v(hgvs_pro),
      hgvs_pro_raw = hgvs_pro,
      score_raw  = as.numeric(score),
      study      = gsub("^urn:mavedb:", "", urn),
      source     = "MaveDB"
    ) |>
    filter(!is.na(p_notation))
}

# ---------------------------------------------------------------------------
# Main loop
# ---------------------------------------------------------------------------
all_canonical <- list()
total_scoresets <- 0
total_failed    <- 0

for (gene in names(GENE_TO_URNS)) {
  for (urn in GENE_TO_URNS[[gene]]) {
    total_scoresets <- total_scoresets + 1
    csv_path <- tryCatch(fetch_scores(urn, gene), error = function(e) {
      message(sprintf("[mavedb] FETCH ERROR %s: %s", urn, e$message))
      NULL
    })
    if (is.null(csv_path) || !file.exists(csv_path)) {
      total_failed <- total_failed + 1
      next
    }
    canon <- normalize_scoreset(csv_path, urn, gene)
    if (is.null(canon) || nrow(canon) == 0) {
      total_failed <- total_failed + 1
      next
    }
    all_canonical[[length(all_canonical) + 1L]] <- canon
    cat(sprintf("[mavedb] %s %s -> %d missense rows\n", gene, urn, nrow(canon)))
    Sys.sleep(0.2)  # polite to MaveDB
  }
}

if (length(all_canonical) == 0L) {
  stop("[mavedb] no scoresets normalized successfully")
}

canonical <- bind_rows(all_canonical) |>
  # Deduplicate within (gene, p_notation, study) -- one row per scoreset
  distinct(gene, p_notation, study, .keep_all = TRUE) |>
  arrange(gene, study, p_notation)

write_tsv(canonical, file.path(OUT_DERIVED, "mavedb_canonical.tsv"))

cat("\n[mavedb] Summary\n")
cat(sprintf("  Scoresets attempted:   %d\n", total_scoresets))
cat(sprintf("  Scoresets failed:      %d\n", total_failed))
cat(sprintf("  Total missense rows:   %d\n", nrow(canonical)))
cat(sprintf("  Unique (gene,p_notation): %d\n",
            n_distinct(canonical$gene, canonical$p_notation)))
cat("\nPer-gene breakdown:\n")
print(canonical |> count(gene, study) |> arrange(gene, desc(n)), n = Inf)
