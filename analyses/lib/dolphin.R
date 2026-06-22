# dolphin.R — DOLPHIN REST API client for the Pfam-based PM1 pathway
#
# Implements the 4th PM1 pathway in the VarViz ACMG/AMP engine, alongside the
# existing CCRS-percentile, UniProt-annotated-domain, and ClinVar-15-residue-
# hotspot pathways. DOLPHIN derives PM1 from eukaryotic Pfam position-specific
# scoring matrices (Corcuff et al., 2023, Front. Bioinform. 3:1127341),
# strictly orthogonal to gnomAD-derived population evidence and to ClinVar-
# derived clinical evidence.
#
# Endpoint:    GET https://dolphin.mmg-gbit.eu/api/VariantAnnotation
# Required:    variation={A123V} (single-letter), ensembl={ENST...},
#              frequencies=1  (the apidoc claims optional, but the server
#                              returns HTTP 500 without it — verified 2026-06-20)
# Response:    JSON with `results` array; each element has its own `acmg`
#              field (e.g. "PM1;" or "PM2;").
#              When the variant is outside any annotated Pfam domain,
#              `results` is returned as a plain string instead of an array.
# Auth:        none (public academic API)
# Rate limit:  60 requests / minute (X-RateLimit-Limit header).
#
# Caveats:
#   - Pass `gene` symbol → resolve to canonical ENST via Ensembl REST first.
#   - DOLPHIN also exposes `Dolphin_AF` (extrapolated cross-species frequency)
#     in the response; we ignore it. VarViz consumes only the ACMG PM1 tag,
#     which is derived from the Pfam alignment scoring matrix and is strictly
#     orthogonal to gnomAD allele frequency (used for VarViz PM2/BA1/BS1/BS2).
#   - DOLPHIN may also emit "PM2;" tags based on its own frequency analysis;
#     these MUST NOT be consumed by VarViz (would conflict with gnomAD-based
#     PM2). Our parser only matches the literal "PM1" tag.
#
# Dependencies (loaded by sourcing this file before use):
#   httr, jsonlite
#
# Cache: in-session via cache_get/cache_set("dolphin", ...) when sourced into
#        server.R. For standalone (harness) use, falls back to a local
#        environment-scoped cache.

# ---------------------------------------------------------------------------
# Standalone cache fallback (only used if cache_get/cache_set are not defined,
# e.g. when called from the dual-pass classification harness).
# ---------------------------------------------------------------------------
.dolphin_local_cache <- new.env(parent = emptyenv())

.dolphin_cache_get <- function(key) {
  if (exists("cache_get", mode = "function")) {
    return(cache_get("dolphin", key))
  }
  if (exists(key, envir = .dolphin_local_cache, inherits = FALSE)) {
    return(get(key, envir = .dolphin_local_cache))
  }
  NULL
}

.dolphin_cache_set <- function(key, value) {
  if (exists("cache_set", mode = "function")) {
    cache_set("dolphin", key, value)
    return(invisible(NULL))
  }
  assign(key, value, envir = .dolphin_local_cache)
  invisible(NULL)
}

# ---------------------------------------------------------------------------
# Amino-acid 3-letter -> 1-letter conversion.
# ---------------------------------------------------------------------------
.dolphin_aa3to1_map <- c(
  Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
  Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
  Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
  Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
  Ter = "*", Sec = "U", Pyl = "O"
)

.dolphin_to_single_letter <- function(p_notation) {
  s <- trimws(as.character(p_notation))
  s <- sub("^p\\.", "", s)
  if (grepl("^[A-Z*][0-9]+[A-Z*]$", s)) return(s)
  m <- regmatches(s, regexec("^([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|\\*)$", s))[[1]]
  if (length(m) == 4) {
    ref3 <- .dolphin_aa3to1_map[tools::toTitleCase(tolower(m[2]))]
    alt3 <- .dolphin_aa3to1_map[tools::toTitleCase(tolower(m[4]))]
    if (!is.na(ref3) && !is.na(alt3)) {
      return(paste0(ref3, m[3], alt3))
    }
  }
  s
}

# ---------------------------------------------------------------------------
# Resolve gene symbol -> canonical Ensembl transcript ID via Ensembl REST.
# Cached per gene symbol. Returns ENST string on success, NA_character_ on
# failure (network error, unknown gene, no canonical/longest transcript).
# ---------------------------------------------------------------------------
dolphin_canonical_enst <- function(gene_symbol, timeout_s = 20) {
  gene_symbol <- trimws(as.character(gene_symbol))
  if (nchar(gene_symbol) == 0) return(NA_character_)
  cache_key <- paste0("ENST:", gene_symbol)

  cached <- .dolphin_cache_get(cache_key)
  if (!is.null(cached)) {
    if (identical(cached, NA)) return(NA_character_)
    return(cached)
  }

  result <- tryCatch({
    url <- paste0(
      "https://rest.ensembl.org/lookup/symbol/homo_sapiens/",
      utils::URLencode(gene_symbol, reserved = TRUE),
      "?content-type=application/json&expand=1"
    )
    resp <- httr::GET(url, httr::timeout(timeout_s),
                      httr::add_headers(Accept = "application/json"))
    if (httr::status_code(resp) != 200) {
      message("[DOLPHIN] Ensembl lookup HTTP ", httr::status_code(resp),
              " for ", gene_symbol)
      return(NA)
    }
    gene_json <- jsonlite::fromJSON(
      httr::content(resp, "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )
    transcripts <- gene_json$Transcript
    if (is.null(transcripts) || length(transcripts) == 0) return(NA)

    canon_tx <- NULL
    for (tx in transcripts) {
      if (!is.null(tx$is_canonical) && tx$is_canonical == 1) {
        canon_tx <- tx
        break
      }
    }
    if (is.null(canon_tx)) {
      cds_lengths <- vapply(transcripts, function(tx) {
        if (!is.null(tx$Translation$length)) as.integer(tx$Translation$length) else 0L
      }, integer(1))
      canon_tx <- transcripts[[which.max(cds_lengths)]]
    }
    if (is.null(canon_tx) || is.null(canon_tx$id)) return(NA)
    canon_tx$id
  }, error = function(e) {
    message("[DOLPHIN] Ensembl lookup error for ", gene_symbol, ": ", e$message)
    NA
  })

  .dolphin_cache_set(cache_key, result)
  if (identical(result, NA)) return(NA_character_)
  result
}

# ---------------------------------------------------------------------------
# Main fetch function. Returns parsed JSON list on success, NULL on miss/error.
# Caller may pass either gene symbol (resolved internally to ENST) or an
# explicit ENST via the `ensembl` argument.
# ---------------------------------------------------------------------------
fetch_dolphin <- function(gene = NULL, p_notation, ensembl = NULL,
                          timeout_s = 15) {
  p_notation <- as.character(p_notation)
  if (nchar(p_notation) == 0) return(NULL)

  variation_1 <- .dolphin_to_single_letter(p_notation)

  if (is.null(ensembl) || nchar(as.character(ensembl)) == 0) {
    if (is.null(gene) || nchar(as.character(gene)) == 0) {
      stop("fetch_dolphin: provide either `ensembl` or `gene`.")
    }
    ensembl <- dolphin_canonical_enst(gene, timeout_s = timeout_s)
    if (is.na(ensembl)) return(NULL)
  }
  ensembl <- trimws(as.character(ensembl))

  cache_key <- paste0(ensembl, ":", variation_1)
  cached <- .dolphin_cache_get(cache_key)
  if (!is.null(cached)) {
    if (identical(cached, NA)) return(NULL)
    return(cached)
  }

  url <- paste0(
    "https://dolphin.mmg-gbit.eu/api/VariantAnnotation",
    "?variation=", utils::URLencode(variation_1, reserved = TRUE),
    "&ensembl=", utils::URLencode(ensembl, reserved = TRUE),
    "&frequencies=1"
  )

  result <- tryCatch({
    resp <- httr::GET(url, httr::timeout(timeout_s),
                      httr::add_headers(Accept = "application/json"))
    sc <- httr::status_code(resp)
    if (sc != 200) {
      message("[DOLPHIN] HTTP ", sc, " for ", ensembl, " ", variation_1)
      return(NA)
    }
    txt <- httr::content(resp, "text", encoding = "UTF-8")
    if (nchar(txt) == 0) return(NA)
    jsonlite::fromJSON(txt, simplifyVector = FALSE)
  }, error = function(e) {
    message("[DOLPHIN] Error for ", ensembl, " ", variation_1, ": ", e$message)
    NA
  })

  .dolphin_cache_set(cache_key, result)
  if (identical(result, NA)) return(NULL)
  result
}

# ---------------------------------------------------------------------------
# PM1 predicate. Walks the `results` array of a DOLPHIN response and returns
# TRUE iff any per-domain record contains a PM1 tag in its `acmg` field.
# Conservative: any parsing failure, absent annotation, or `results` returned
# as the string sentinel ("This mutation is not within a protein domain ...")
# returns FALSE.
#
# Notes:
#   - Only the literal "PM1" tag is matched. DOLPHIN may also emit "PM2;"
#     (from its own frequency analysis), which we intentionally ignore since
#     VarViz already derives PM2 from gnomAD.
#   - The regex avoids spuriously matching "PM10", "PM11", etc.
# ---------------------------------------------------------------------------
dolphin_fires_pm1 <- function(dolphin_resp) {
  if (is.null(dolphin_resp)) return(FALSE)
  results <- dolphin_resp$results
  if (is.null(results)) return(FALSE)
  if (is.character(results)) return(FALSE)  # "mutation not within domain ..."
  if (!is.list(results) || length(results) == 0) return(FALSE)

  pm1_rx <- "(^|[;[:space:]])PM1([[:space:]]|;|$)"
  for (r in results) {
    if (!is.list(r)) next
    acmg_field <- r$acmg
    if (is.null(acmg_field) || length(acmg_field) == 0) next
    acmg_str <- paste(unlist(acmg_field), collapse = ";")
    if (nchar(acmg_str) == 0) next
    if (grepl(pm1_rx, acmg_str)) return(TRUE)
  }
  FALSE
}

# ---------------------------------------------------------------------------
# Convenience: one-call wrapper returning the boolean PM1 decision.
# Either `ensembl` or `gene` must be provided.
# ---------------------------------------------------------------------------
dolphin_pm1_call <- function(gene = NULL, p_notation, ensembl = NULL,
                             timeout_s = 15) {
  resp <- fetch_dolphin(gene = gene, p_notation = p_notation,
                        ensembl = ensembl, timeout_s = timeout_s)
  dolphin_fires_pm1(resp)
}
