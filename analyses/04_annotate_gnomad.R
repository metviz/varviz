# Annotates the variant universe with gnomAD v4.1 presence + allele count.
# For each unique gene in the universe, queries gnomAD GraphQL once and caches
# the raw JSON response under analyses/raw/gnomad/<gene>.json. Resumable: a
# re-run skips cached genes and avoids re-hitting the API.
#
# Output schema (universe + 5 new columns):
#   ..., gnomad_ac, gnomad_an, gnomad_af, gnomad_singleton, gnomad_present
#
# The CAPS analysis (Task 3.3, Gudkov 2025) restricts to gnomad_present == TRUE.
#
# Run from project root: `Rscript analyses/04_annotate_gnomad.R`

suppressMessages({
  library(httr2)
  library(dplyr)
  library(readr)
  library(jsonlite)
  library(purrr)
  library(tibble)
})
source("analyses/lib/manifest.R")

# ---------------------------------------------------------------------------
# 3-letter → 1-letter HGVS-p normalizer. gnomAD's hgvsp is typically
# `p.Arg175His` (three-letter); our universe p_notation is one-letter
# (`p.R175H`) per analyses/03_build_universe.R. We must convert before joining.
# Mirrors the helper in 03_build_universe.R; intentionally inlined to keep
# Phase 1.4 self-contained (no cross-script lib refactor in scope).
# ---------------------------------------------------------------------------
AA_3to1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G",
  His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S",
  Thr="T", Trp="W", Tyr="Y", Val="V", Ter="*", Sec="U", Pyl="O"
)
to_one_letter <- function(p) {
  if (is.null(p) || is.na(p) || !nzchar(p)) return(NA_character_)
  s <- sub("^[^:]*:", "", p)
  s <- gsub("[()]", "", s)
  if (grepl("^p\\.[A-Z*][0-9]+[A-Z*]$", s)) return(s)
  m <- regmatches(s, regexec("^p\\.([A-Z][a-z]{2}|\\*)([0-9]+)([A-Z][a-z]{2}|\\*)$", s))[[1]]
  if (length(m) != 4) return(NA_character_)
  ref1 <- if (m[2] == "*") "*" else AA_3to1[m[2]]
  alt1 <- if (m[4] == "*") "*" else AA_3to1[m[4]]
  if (is.na(ref1) || is.na(alt1)) return(NA_character_)
  paste0("p.", ref1, m[3], alt1)
}
to_one_letter_v <- Vectorize(to_one_letter, USE.NAMES = FALSE)

UNIVERSE_IN  <- "analyses/derived/variant_universe.tsv"
OUT_FILE     <- "analyses/derived/variant_universe_gnomad.tsv"
CACHE_DIR    <- "analyses/raw/gnomad"
MANIFEST_DIR <- "analyses/manifests"
dir.create(CACHE_DIR, recursive = TRUE, showWarnings = FALSE)

GQL_ENDPOINT <- "https://gnomad.broadinstitute.org/api"

# gnomAD v4.1, GRCh38, MANE Select transcripts. exome+genome AC/AN summed.
GQL_QUERY <- '
query GeneVariants($symbol: String!) {
  gene(gene_symbol: $symbol, reference_genome: GRCh38) {
    variants(dataset: gnomad_r4) {
      variant_id
      hgvsp
      consequence
      transcript_id
      exome { ac an }
      genome { ac an }
    }
  }
}
'

`%||%` <- function(a, b) if (is.null(a)) b else a

fetch_gnomad_gene <- function(gene) {
  cache_file <- file.path(CACHE_DIR, paste0(gene, ".json"))
  if (file.exists(cache_file)) {
    return(fromJSON(cache_file, simplifyVector = FALSE))
  }

  req <- request(GQL_ENDPOINT) |>
    req_method("POST") |>
    req_user_agent("VarViz-Pass2/1.0 (mailto:agasthyametpally5@gmail.com)") |>
    req_timeout(180) |>
    req_retry(max_tries = 3, backoff = function(i) 2^i) |>
    req_body_json(list(query = GQL_QUERY, variables = list(symbol = gene)))

  resp <- tryCatch(req_perform(req), error = function(e) {
    cat(sprintf("  [ERROR] gnomAD fetch failed for %s: %s\n", gene, conditionMessage(e)))
    NULL
  })
  if (is.null(resp)) return(NULL)

  raw <- resp_body_raw(resp)
  writeBin(raw, cache_file)

  write_manifest(
    out_dir       = MANIFEST_DIR,
    label         = sprintf("gnomad_gene_%s", gene),
    url           = GQL_ENDPOINT,
    params        = list(symbol = gene, dataset = "gnomad_r4", reference_genome = "GRCh38"),
    response_raw  = raw,
    http_status   = resp_status(resp)
  )

  fromJSON(rawToChar(raw), simplifyVector = FALSE)
}

# Extract per-variant AC/AN/AF rows from a gnomAD GraphQL response.
extract_variant_rows <- function(gene, gnomad_resp) {
  if (is.null(gnomad_resp$data$gene)) return(tibble())
  vars <- gnomad_resp$data$gene$variants
  if (is.null(vars) || length(vars) == 0) return(tibble())

  map_dfr(vars, function(v) {
    if (is.null(v$hgvsp)) return(NULL)
    ex_ac <- v$exome$ac  %||% 0L
    ex_an <- v$exome$an  %||% 0L
    gn_ac <- v$genome$ac %||% 0L
    gn_an <- v$genome$an %||% 0L
    total_ac <- ex_ac + gn_ac
    total_an <- ex_an + gn_an
    tibble(
      gene             = gene,
      hgvsp            = v$hgvsp,
      consequence      = v$consequence %||% NA_character_,
      gnomad_ac        = as.integer(total_ac),
      gnomad_an        = as.integer(total_an),
      gnomad_af        = if (total_an > 0) total_ac / total_an else NA_real_,
      gnomad_singleton = (total_ac == 1L)
    )
  })
}

# --- Main loop ---
universe <- read_tsv(UNIVERSE_IN, show_col_types = FALSE)
genes <- sort(unique(universe$gene))
cat(sprintf("[gnomad] Annotating %d genes\n", length(genes)))

annotations <- map_dfr(seq_along(genes), function(i) {
  gene <- genes[i]
  cached <- file.exists(file.path(CACHE_DIR, paste0(gene, ".json")))
  cat(sprintf("[gnomad] [%2d/%d] %s%s\n",
              i, length(genes), gene, if (cached) " (cached)" else ""))
  resp <- fetch_gnomad_gene(gene)
  if (is.null(resp)) return(tibble())
  extract_variant_rows(gene, resp)
})

# Normalize gnomAD hgvsp to one-letter to match universe$p_notation
annotations <- annotations |>
  mutate(p_notation = to_one_letter_v(hgvsp)) |>
  filter(!is.na(p_notation)) |>
  # Same (gene, p_notation) can repeat across nucleotide variants; pick the one
  # with the highest total AC so we annotate the variant that's actually present.
  group_by(gene, p_notation) |>
  slice_max(order_by = gnomad_ac, n = 1, with_ties = FALSE) |>
  ungroup()

joined <- universe |>
  left_join(
    annotations |>
      select(gene, p_notation, gnomad_ac, gnomad_an, gnomad_af, gnomad_singleton),
    by = c("gene", "p_notation")
  ) |>
  mutate(
    gnomad_present   = !is.na(gnomad_ac),
    gnomad_singleton = ifelse(is.na(gnomad_singleton), FALSE, gnomad_singleton)
  )

write_tsv(joined, OUT_FILE)

cat("\n[gnomad] Universe with gnomAD annotation:\n")
cat(sprintf("  Total rows:       %d\n", nrow(joined)))
cat(sprintf("  gnomAD-present:   %d  (%.1f%%)\n",
            sum(joined$gnomad_present),
            100 * mean(joined$gnomad_present)))
cat(sprintf("  Singletons (AC=1):%d\n", sum(joined$gnomad_singleton, na.rm = TRUE)))

cat("\n[gnomad] gnomAD-present by gene:\n")
print(joined |> group_by(gene) |>
        summarize(n_total = n(),
                  n_present = sum(gnomad_present),
                  n_singleton = sum(gnomad_singleton, na.rm = TRUE),
                  .groups = "drop") |>
        arrange(desc(n_present)),
      n = Inf)
