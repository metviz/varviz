# Builds the unified Pass 2 variant universe by unioning VariBench (clinical
# truth) and MaveDB (functional truth) canonical TSVs, then defensively
# filtering to BENCHMARK_GENES per MF2.
#
# Output schema (one row per (gene, p_notation, source)):
#   gene, p_notation, source, label, score_raw, study, set, hgvs_g,
#   clinvar_alleleid, hgvs_pro_raw
#
# Two source groups:
#   - VariBench rows have label (Pathogenic / Benign), set ("PON-PS_D2"),
#     hgvs_g, clinvar_alleleid; score_raw and study are NA.
#   - MaveDB rows have score_raw and study (URN suffix); label, set,
#     hgvs_g, clinvar_alleleid are NA.
#
# Run from project root: `Rscript analyses/03_build_universe.R`
#
# Outputs:
#   analyses/derived/variant_universe.tsv            (gitignored)

suppressMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
})
source("analyses/lib/benchmark_genes.R")

# ---------------------------------------------------------------------------
# p-notation normalizer: convert three-letter (p.Leu344Pro) → one-letter
# (p.L344P). VariBench uses three-letter, MaveDB uses one-letter; this brings
# them into a common form so (gene, p_notation) joins work across sources.
# ---------------------------------------------------------------------------
AA_3to1 <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G",
  His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S",
  Thr="T", Trp="W", Tyr="Y", Val="V", Ter="*", Sec="U", Pyl="O"
)
to_one_letter <- function(p) {
  if (is.na(p) || !nzchar(p)) return(NA_character_)
  s <- sub("^[^:]*:", "", p)
  s <- gsub("[()]", "", s)
  # Already one-letter form (p.R175H or p.R175*)?
  if (grepl("^p\\.[A-Z*][0-9]+[A-Z*]$", s)) return(s)
  # Three-letter form (p.Arg175His)
  m <- regmatches(s, regexec("^p\\.([A-Z][a-z]{2}|\\*)([0-9]+)([A-Z][a-z]{2}|\\*)$", s))[[1]]
  if (length(m) != 4) return(NA_character_)
  ref1 <- if (m[2] == "*") "*" else AA_3to1[m[2]]
  alt1 <- if (m[4] == "*") "*" else AA_3to1[m[4]]
  if (is.na(ref1) || is.na(alt1)) return(NA_character_)
  paste0("p.", ref1, m[3], alt1)
}
to_one_letter_v <- Vectorize(to_one_letter, USE.NAMES = FALSE)

OUT       <- "analyses/derived/variant_universe.tsv"
VB_PATH   <- "analyses/derived/varibench_canonical.tsv"
MV_PATH   <- "analyses/derived/mavedb_canonical.tsv"

stopifnot(file.exists(VB_PATH), file.exists(MV_PATH))

vb <- read_tsv(VB_PATH, show_col_types = FALSE)
mv <- read_tsv(MV_PATH, show_col_types = FALSE)

cat(sprintf("[universe] VariBench in:   %d rows\n", nrow(vb)))
cat(sprintf("[universe] MaveDB in:      %d rows\n", nrow(mv)))

# --- Normalize each source to a common wide schema ---
vb_canonical <- vb |>
  transmute(
    gene             = as.character(gene),
    p_notation       = to_one_letter_v(as.character(p_notation)),  # 3-letter → 1-letter
    p_notation_raw   = as.character(p_notation),                   # keep three-letter for traceability
    source           = as.character(source),
    label            = as.character(label),
    score_raw        = NA_real_,
    study            = NA_character_,
    set              = as.character(set),
    hgvs_g           = as.character(hgvs_g),
    clinvar_alleleid = as.character(clinvar_alleleid),
    hgvs_pro_raw     = NA_character_
  ) |>
  filter(!is.na(p_notation))  # drop frameshifts, deletions, synonymous (=) etc.

mv_canonical <- mv |>
  transmute(
    gene             = gene,
    p_notation       = to_one_letter_v(p_notation),  # idempotent for already-one-letter form
    p_notation_raw   = p_notation,
    source           = source,
    label            = NA_character_,
    score_raw        = score_raw,
    study            = study,
    set              = NA_character_,
    hgvs_g           = NA_character_,
    clinvar_alleleid = NA_character_,
    hgvs_pro_raw     = hgvs_pro_raw
  ) |>
  filter(!is.na(p_notation))

# --- Union, defensive filter to BENCHMARK_GENES, deduplicate within source ---
universe <- bind_rows(vb_canonical, mv_canonical) |>
  filter(gene %in% BENCHMARK_GENES) |>
  distinct(gene, p_notation, source, study, .keep_all = TRUE) |>
  arrange(gene, p_notation, source)

# Hard sanity check
stopifnot(all(universe$gene %in% BENCHMARK_GENES))

write_tsv(universe, OUT)

cat(sprintf("\n[universe] Output:         %d rows  (%d unique gene+p_notation pairs)\n",
            nrow(universe), n_distinct(universe$gene, universe$p_notation)))
cat(sprintf("[universe] Unique genes:   %d / %d BENCHMARK_GENES\n",
            n_distinct(universe$gene), length(BENCHMARK_GENES)))

absent <- setdiff(BENCHMARK_GENES, universe$gene)
if (length(absent) > 0) {
  cat(sprintf("[universe] Genes with no rows: %s\n", paste(absent, collapse = ", ")))
}

cat("\n[universe] By source:\n")
print(table(universe$source))

cat("\n[universe] By gene × source:\n")
gene_x_source <- universe |> count(gene, source) |>
  tidyr::pivot_wider(names_from = source, values_from = n, values_fill = 0L)
# Compute totals safely (some columns may not exist if all rows are one source)
gene_x_source$VariBench <- if ("VariBench" %in% colnames(gene_x_source)) gene_x_source$VariBench else 0L
gene_x_source$MaveDB    <- if ("MaveDB"    %in% colnames(gene_x_source)) gene_x_source$MaveDB    else 0L
gene_x_source$total     <- gene_x_source$VariBench + gene_x_source$MaveDB
print(gene_x_source |> arrange(desc(total)), n = Inf)

# Variants with BOTH clinical and functional evidence (intersection diagnostic)
both_idx <- universe |>
  group_by(gene, p_notation) |>
  summarize(has_vb = any(source == "VariBench"),
            has_mv = any(source == "MaveDB"),
            .groups = "drop") |>
  filter(has_vb & has_mv)

cat(sprintf("\n[universe] Variants with BOTH VariBench label AND MaveDB score: %d\n",
            nrow(both_idx)))
cat("           (these support both Panel A clinical concordance and Panel B DMS concordance)\n")
cat("           Per-gene intersection breakdown:\n")
print(both_idx |> count(gene, sort = TRUE), n = Inf)
