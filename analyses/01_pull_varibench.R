# Pulls VariBench missense pathogenicity test sets into analyses/raw/varibench/
# and writes canonical TSVs to analyses/derived/.
#
# VariBench host: https://structure.bmc.lu.se/VariBench/ (HTTPS, after redirect from HTTP).
# Two test sets are downloaded:
#   * F1  — Dataset 1 F1: ClinVar one-star variants (14,819 rows; 7,346 benign + 7,473 pathogenic)
#   * F11 — Dataset 1 F11: Exclude LP/LB, ClinVar Sept 2016 (7,766 rows; only "pathogenic"/"benign"
#           assertions; higher stringency)
# Files are tab-separated despite the .csv extension. Columns differ between F1 and F11; both are
# normalized to the canonical schema (gene, p_notation, hgvs_p_raw, label, review_status, source, set).
#
# The canonical TSVs are filtered to the seven VarViz-supported genes
# (TP53, PPP2R5C, BAP1, DDR2, TSHR, SLC13A5, CASR). Per-gene/label counts in F1
# (probed 2026-04-26): TP53 62P+6B, CASR 18P+6B, BAP1 0P+2B, TSHR 1P+0B; others 0.
# Filtering happens here so downstream `03_build_universe.R` works on a small, on-panel set.
#
# Run from project root:
#   Rscript analyses/01_pull_varibench.R
#
# Outputs:
#   analyses/raw/varibench/D1_F1_clinvar_one_star.tsv      (raw download — gitignored)
#   analyses/raw/varibench/D1_F11_clinvar_exclude_LP_LB.tsv
#   analyses/derived/varibench_canonical.tsv               (panel-filtered union — gitignored)
#   analyses/manifests/<timestamp>__varibench_*.json       (one per fetch — committed)

suppressMessages({
  library(httr2)
  library(dplyr)
  library(readr)
  library(stringr)
})
source("analyses/lib/manifest.R")

OUT_RAW       <- "analyses/raw/varibench"
OUT_DERIVED   <- "analyses/derived"
OUT_MANIFEST  <- "analyses/manifests"
dir.create(OUT_RAW,     recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_DERIVED, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_MANIFEST, recursive = TRUE, showWarnings = FALSE)

VARVIZ_GENES <- c("TP53", "PPP2R5C", "BAP1", "DDR2", "TSHR", "SLC13A5", "CASR")

# ---------------------------------------------------------------------------
# fetch_varibench_file(): download one VariBench file, write raw + manifest.
# Returns the local file path.
# ---------------------------------------------------------------------------
fetch_varibench_file <- function(url, out_name, label) {
  out_file <- file.path(OUT_RAW, out_name)
  message("[varibench] GET ", url)
  resp <- request(url) |>
    req_timeout(60) |>
    req_user_agent("VarViz-Pass2/1.0 (research)") |>
    req_perform()
  if (resp_status(resp) != 200) stop("[varibench] HTTP ", resp_status(resp), " for ", url)
  raw_bytes <- resp_body_raw(resp)
  writeBin(raw_bytes, out_file)

  write_manifest(
    out_dir       = OUT_MANIFEST,
    label         = label,
    url           = url,
    params        = list(),
    response_raw  = raw_bytes,
    http_status   = resp_status(resp)
  )
  out_file
}

# ---------------------------------------------------------------------------
# Three-letter ↔ one-letter helpers — convert NP_xxxxx.x:p.Arg387Gln to p.R387Q
# ---------------------------------------------------------------------------
AA_3to1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G",
  His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S",
  Thr="T", Trp="W", Tyr="Y", Val="V", Ter="*", Stop="*", Sec="U", Pyl="O"
)

normalize_p_notation <- function(hgvs_p) {
  # Examples to handle:
  #   "NP_000633.2:p.Arg387Gln"          -> "p.R387Q"
  #   "p.Arg387Gln"                       -> "p.R387Q"
  #   "p.R387Q"                           -> "p.R387Q"
  #   "NP_000633.2:p.(Arg387Gln)"         -> "p.R387Q"
  #   "NP_000633.2:p.Arg387Ter"           -> "p.R387*"
  s <- sub("^[^:]*:", "", hgvs_p)             # strip protein-id prefix
  s <- gsub("[()]", "", s)                    # strip parentheses
  m <- regmatches(s, regexec("^p\\.([A-Za-z]{3}|\\*)([0-9]+)([A-Za-z]{3}|\\*)$", s))[[1]]
  if (length(m) != 4) return(NA_character_)   # not a clean missense substitution
  ref3 <- m[2]; pos <- m[3]; alt3 <- m[4]
  ref1 <- if (ref3 == "*") "*" else AA_3to1[ref3]
  alt1 <- if (alt3 == "*") "*" else AA_3to1[alt3]
  if (is.na(ref1) || is.na(alt1)) return(NA_character_)
  paste0("p.", ref1, pos, alt1)
}
normalize_p_notation_v <- Vectorize(normalize_p_notation, USE.NAMES = FALSE)

# ---------------------------------------------------------------------------
# canonicalize(): map a VariBench frame to the shared schema.
# ---------------------------------------------------------------------------
canonicalize <- function(df, source_label) {
  df |>
    transmute(
      gene          = symbol,
      hgvs_p_raw    = hgvs_p,
      p_notation    = normalize_p_notation_v(hgvs_p),
      label         = clinical_significance,
      review_status = review_status,
      source        = source_label,
      set           = "test"
    ) |>
    filter(!is.na(p_notation), gene %in% VARVIZ_GENES) |>
    distinct(gene, p_notation, .keep_all = TRUE)
}

# ---------------------------------------------------------------------------
# Pull F1 (Dataset 1, one-star ClinVar)
# ---------------------------------------------------------------------------
F1_URL  <- "https://structure.bmc.lu.se/VariBench/data/variationtype/substitutions/test/Dataset1/VTD_SCR_TESTD_D1_F1_ClinVar_one_Star.csv"
f1_path <- fetch_varibench_file(F1_URL, "D1_F1_clinvar_one_star.tsv", "varibench_D1_F1_clinvar_one_star")
f1_df   <- read_tsv(f1_path, show_col_types = FALSE)
cat("[varibench] F1 rows: ", nrow(f1_df),
    "; cols: ", paste(colnames(f1_df), collapse = ","), "\n", sep = "")
f1_canonical <- canonicalize(f1_df, "VariBench-D1-F1")
cat("[varibench] F1 panel-filtered: ", nrow(f1_canonical), "\n", sep = "")

# ---------------------------------------------------------------------------
# Pull F11 (Dataset 1, Exclude LP/LB; higher stringency)
# ---------------------------------------------------------------------------
F11_URL  <- "https://structure.bmc.lu.se/VariBench/data/variationtype/substitutions/test/Dataset1/VTD_SCR_TESTD_D1_F11_Exclude_Lp_LB_ClinVar_Sep_2016.csv"
f11_path <- fetch_varibench_file(F11_URL, "D1_F11_clinvar_exclude_LP_LB.tsv", "varibench_D1_F11_exclude_LP_LB")
f11_df   <- read_tsv(f11_path, show_col_types = FALSE)
cat("[varibench] F11 rows: ", nrow(f11_df),
    "; cols: ", paste(colnames(f11_df), collapse = ","), "\n", sep = "")
f11_canonical <- canonicalize(f11_df, "VariBench-D1-F11")
cat("[varibench] F11 panel-filtered: ", nrow(f11_canonical), "\n", sep = "")

# ---------------------------------------------------------------------------
# Union, deduplicate (prefer higher-stringency F11 row when same variant in both)
# ---------------------------------------------------------------------------
canonical <- bind_rows(f11_canonical, f1_canonical) |>
  distinct(gene, p_notation, .keep_all = TRUE) |>
  arrange(gene, p_notation)

write_tsv(canonical, file.path(OUT_DERIVED, "varibench_canonical.tsv"))

cat("\n[varibench] canonical written to ", file.path(OUT_DERIVED, "varibench_canonical.tsv"), "\n", sep = "")
cat("[varibench] union rows: ", nrow(canonical), "\n", sep = "")
cat("[varibench] per-gene/label breakdown:\n")
print(canonical |> count(gene, label) |> arrange(gene, label), n = Inf)
