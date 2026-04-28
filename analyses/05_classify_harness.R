# Phase 2.1b — Dual-pass classification harness (MF1)
#
# For each variant in the universe, drive server.R's existing fetch + classify
# functions to produce TWO classifications:
#   Pass-Full  — uses all ACMG criteria including ClinVar-derived ones
#                (PS1*, PM5, PM1-hotspot, PP5*, BP6*)
#   Pass-Blind — same engine, but with ClinVar-derived tags withheld via
#                analyses/lib/clinvar_blind.R::strip_clinvar_tags()
#
# Single source of truth: this script CALLS build_variant_table() — it does
# not re-implement tag derivation. Pass-Full classification is bit-identical
# to what the live VarViz app produces. Per server.R commit 0146843,
# build_variant_table() exposes ACMG_PM1_Pathway as a column so the helper
# can selectively strip PM1 only when its firing pathway was the ClinVar
# 15aa hotspot.
#
# Resumability: per-gene checkpoint files at analyses/classifications/<gene>__dual.tsv.
# Re-running skips genes whose checkpoint already exists.
#
# Run from project root: `Rscript analyses/05_classify_harness.R`
#
# Output (gitignored per Pass 2 policy):
#   analyses/classifications/<gene>__dual.tsv  (per-gene checkpoint)
#   analyses/derived/varviz_classifications.tsv (concatenated summary)

suppressMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

cat("[harness] Sourcing server.R...\n")
t0 <- Sys.time()
suppressMessages(source("server.R"))
cat(sprintf("[harness] Sourced server.R in %.1f sec\n",
            as.numeric(Sys.time() - t0, units = "secs")))

source("analyses/lib/clinvar_blind.R")

UNIVERSE_IN    <- "analyses/derived/variant_universe_gnomad.tsv"
CHECKPOINT_DIR <- "analyses/classifications"
SUMMARY_OUT    <- "analyses/derived/varviz_classifications.tsv"
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

universe <- read_tsv(UNIVERSE_IN, show_col_types = FALSE)
genes <- sort(unique(universe$gene))
cat(sprintf("[harness] Universe: %d rows / %d genes\n", nrow(universe), length(genes)))

`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# Extract integer position from one-letter HGVS-p (p.X175Y -> 175).
extract_pos <- function(p) suppressWarnings(as.integer(sub("^p\\.[A-Z*]([0-9]+)[A-Z*]$", "\\1", p)))

# -----------------------------------------------------------------------------
# Per-gene driver: fetches gene-level data once, then calls build_variant_table
# on the gene's variants in a single batch. Returns a tibble with one wide row
# per variant: tags_full + tags_blind + Pass-Full + Pass-Blind classifications.
# -----------------------------------------------------------------------------
classify_gene <- function(gene_name) {
  gene_attrib_row <- gene_data[gene_data$gene_name == gene_name, ]
  if (nrow(gene_attrib_row) == 0) {
    cat(sprintf("  [%s] NOT in gene_data — skipping\n", gene_name))
    return(NULL)
  }
  uid <- as.character(gene_attrib_row$uniprot_id[1])

  # Gene-level fetches (one call each)
  pfam_d    <- tryCatch(extract_pfam(uid),                                 error = function(e) NULL)
  uniprot_d <- tryCatch(extract_uniprot_feature_data(uid),                 error = function(e) NULL)
  gnomad_d  <- tryCatch(extract_gnomad(gene_name),                         error = function(e) NULL)
  clinvar_d <- tryCatch(extract_clinvar(gene_name),                        error = function(e) NULL)
  ccrs_d    <- tryCatch(extract_ccrs(gene_name, pfam_d$primaryAccession),  error = function(e) NULL)
  af_d      <- tryCatch(extract_alphafold_plddt(uid),                      error = function(e) NULL)
  mean_d    <- tryCatch(get_mean_pathogenicity(uid),                       error = function(e) NULL)
  gi_d      <- tryCatch(extract_gene_info_uniprot(uid, gene_name),         error = function(e) NULL)

  hgnc_for_clingen <- if (!is.null(gi_d) && !is.null(gi_d$hgnc_id) && nchar(gi_d$hgnc_id) > 0) gi_d$hgnc_id else NULL
  clingen_d <- tryCatch(
    fetch_clingen_validity(gene_name, hgnc_id = hgnc_for_clingen),
    error = function(e) list(classification = NA_character_, moi = "", disease = "", source = "ClinGen LDH")
  )

  afs_dest <- file.path(getwd(), paste0(uid, "-F1-aa-substitutions.csv"))
  afs_d <- if (file.exists(afs_dest)) {
    tryCatch(data.table::fread(afs_dest), error = function(e) NULL)
  } else NULL

  # Build highlight_df for this gene's universe variants
  rows <- universe[universe$gene == gene_name, ]
  hl <- data.frame(
    Mutation = rows$p_notation,
    prot_pos = extract_pos(rows$p_notation),
    gene     = rows$gene,
    stringsAsFactors = FALSE
  )
  hl <- hl[!is.na(hl$prot_pos), , drop = FALSE]
  if (nrow(hl) == 0) return(NULL)

  cat(sprintf("  [%s] uid=%s, %d variants, calling build_variant_table()...\n",
              gene_name, uid, nrow(hl)))
  t1 <- Sys.time()
  vtbl <- tryCatch(
    build_variant_table(
      hl, af_d, mean_d, afs_d, gnomad_d, clinvar_d,
      pfam_d, uniprot_d, ccrs_d,
      af_cutoff = 0.0001, ac_cutoff = 13,
      clinvar_missense = NULL, consurf_data = NULL,
      denovo_status     = "not_denovo",
      inh_param         = "monoallelic",
      cutoff_method     = "calc_af",
      prevalence_1_in_n = 2000,
      allelic_het = 0.5, genetic_het = 1.0, penetrance = 1.0,
      pop_size = 125748, conf_interval = 0.95,
      clingen_disease_param = clingen_d$disease %||% "",
      clingen_moi_param     = clingen_d$moi     %||% "",
      consurf_file_name     = ""
    ),
    error = function(e) { cat(sprintf("    ERROR: %s\n", conditionMessage(e))); NULL }
  )
  if (is.null(vtbl) || nrow(vtbl) == 0) return(NULL)
  cat(sprintf("    build_variant_table %.1fs (%d rows)\n",
              as.numeric(Sys.time() - t1, units = "secs"), nrow(vtbl)))

  # Vectorized dual-pass over the gene's rows
  tags_full_str_vec <- as.character(vtbl$ACMG_Tags)
  pm1_path_vec      <- as.character(vtbl$ACMG_PM1_Pathway)

  parse_tags <- function(s) if (nchar(s) > 0) trimws(strsplit(s, ",")[[1]]) else character(0)

  out_rows <- map_dfr(seq_len(nrow(vtbl)), function(i) {
    tags_full  <- parse_tags(tags_full_str_vec[i])
    tags_blind <- strip_clinvar_tags(tags_full, pm1_pathway = pm1_path_vec[i])
    res_full   <- tryCatch(classify_acmg(tags_full),  error = function(e) list(classification = NA_character_, pts = NA_real_, rule = NA_character_))
    res_blind  <- tryCatch(classify_acmg(tags_blind), error = function(e) list(classification = NA_character_, pts = NA_real_, rule = NA_character_))
    tibble(
      gene                          = gene_name,
      p_notation                    = as.character(vtbl$Variant[i]),
      pm1_pathway                   = pm1_path_vec[i],
      tags_full                     = tags_full_str_vec[i],
      tags_blind                    = paste(tags_blind, collapse = ","),
      varviz_classification_full    = res_full$classification,
      varviz_pts_full               = res_full$pts,
      varviz_rule_full              = res_full$rule,
      varviz_classification_blind   = res_blind$classification,
      varviz_pts_blind              = res_blind$pts,
      varviz_rule_blind             = res_blind$rule
    )
  })
  out_rows
}

# -----------------------------------------------------------------------------
# Main loop with per-gene checkpointing
# -----------------------------------------------------------------------------
run_one <- function(gene_name) {
  ckpt <- file.path(CHECKPOINT_DIR, paste0(gene_name, "__dual.tsv"))
  if (file.exists(ckpt)) {
    cat(sprintf("  [%s] cached -> %s\n", gene_name, ckpt))
    return(read_tsv(ckpt, show_col_types = FALSE))
  }
  result <- classify_gene(gene_name)
  if (!is.null(result) && nrow(result) > 0) write_tsv(result, ckpt)
  result
}

cat(sprintf("\n[harness] Running dual-pass on %d genes\n", length(genes)))
all_results <- list()
for (i in seq_along(genes)) {
  g <- genes[i]
  cat(sprintf("[%2d/%d] %s\n", i, length(genes), g))
  r <- run_one(g)
  if (!is.null(r)) all_results[[g]] <- r
}

if (length(all_results) == 0) {
  cat("[harness] FAIL: no gene yielded results\n"); quit(status = 1)
}

summary_df <- bind_rows(all_results)
write_tsv(summary_df, SUMMARY_OUT)

# -----------------------------------------------------------------------------
# Summary diagnostics
# -----------------------------------------------------------------------------
cat(sprintf("\n[harness] Dual-pass classification complete: %d variants\n", nrow(summary_df)))
cat("\n[harness] Pass-Full classification distribution:\n")
print(table(summary_df$varviz_classification_full, useNA = "ifany"))
cat("\n[harness] Pass-Blind classification distribution:\n")
print(table(summary_df$varviz_classification_blind, useNA = "ifany"))

n_disagree <- sum(summary_df$varviz_classification_full != summary_df$varviz_classification_blind, na.rm = TRUE)
cat(sprintf("\n[harness] Variants where Full vs Blind classifications disagree: %d / %d (%.1f%%)\n",
            n_disagree, nrow(summary_df), 100 * n_disagree / nrow(summary_df)))

cat("\n[harness] PM1 pathway distribution (Pass-Full):\n")
print(table(summary_df$pm1_pathway[nchar(summary_df$pm1_pathway) > 0], useNA = "ifany"))

cat("\n[harness] Variants with clinvar_hotspot PM1 (where Pass-Blind strips PM1):\n")
hot <- summary_df[summary_df$pm1_pathway == "clinvar_hotspot", ]
cat(sprintf("  count: %d\n", nrow(hot)))
if (nrow(hot) > 0) {
  shifts <- table(hot$varviz_classification_full, hot$varviz_classification_blind)
  cat("  Full -> Blind shift among clinvar_hotspot PM1 variants:\n")
  print(shifts)
}
