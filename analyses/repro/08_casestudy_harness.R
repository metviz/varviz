# Dual-pass classification harness — MDS corroborate build.
#
# Differs from analyses/05_classify_harness.R only in configuration: server.R
# Path 4 now scores MDS for every variant, originating PM1 where no other path
# fired and adding its LR tier to the base strength (capped at Strong) where
# one did. The ClinVar-hotspot upgrade records its pre-upgrade tag so
# strip_clinvar_tags() can demote it under Pass-Blind.
#
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
# Resumability: per-gene checkpoint files at analyses/repro/casestudy/classifications/<gene>__dual.tsv.
# Re-running skips genes whose checkpoint already exists.
#
# Run from project root: `Rscript analyses/05_classify_harness.R`
#
# Output (gitignored per Pass 2 policy):
#   analyses/repro/casestudy/classifications/<gene>__dual.tsv  (per-gene checkpoint)
#   analyses/repro/casestudy/summary.tsv (concatenated summary)

suppressMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

cat("[harness] Sourcing server.R...\n")
t0 <- Sys.time()
suppressMessages(source("server.R"))
options(varviz.mds_tiered = TRUE)   # LR tiers: +2 / +3 / +4 at MDS <= -4 / -8 / -12
options(varviz.mds_pm1    = TRUE)   # MDS evaluated for every variant (originate + corroborate)
cat("[harness] varviz.mds_tiered = TRUE ; varviz.mds_pm1 = TRUE (ungated Path 4)\n")
cat(sprintf("[harness] Sourced server.R in %.1f sec\n",
            as.numeric(Sys.time() - t0, units = "secs")))

source("analyses/lib/clinvar_blind.R")
source("analyses/lib/local_predictors.R")
source("analyses/lib/harness_guard.R")

UNIVERSE_IN    <- "analyses/repro/casestudy_variants.tsv"
CHECKPOINT_DIR <- "analyses/repro/casestudy/classifications"
SUMMARY_OUT    <- "analyses/repro/casestudy/summary.tsv"
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
  # A failed fetch is recorded in options(varviz.harness_fetch_failed) so
  # run_one() refuses to checkpoint the gene (see analyses/lib/harness_guard.R).
  pfam_d    <- harness_fetch("pfam",    gene_name, extract_pfam(uid))
  uniprot_d <- harness_fetch("uniprot", gene_name, extract_uniprot_feature_data(uid))
  gnomad_d  <- harness_fetch("gnomad",  gene_name, extract_gnomad(gene_name))
  clinvar_d <- harness_fetch("clinvar", gene_name, extract_clinvar(gene_name))
  ccrs_d    <- harness_fetch("ccrs",    gene_name, extract_ccrs(gene_name, pfam_d$primaryAccession))
  af_d      <- harness_fetch("alphafold", gene_name, extract_alphafold_plddt(uid))
  mean_d    <- harness_fetch("mean_path", gene_name, get_mean_pathogenicity(uid))
  gi_d      <- harness_fetch("gene_info", gene_name, extract_gene_info_uniprot(uid, gene_name))

  hgnc_for_clingen <- if (!is.null(gi_d) && !is.null(gi_d$hgnc_id) && nchar(gi_d$hgnc_id) > 0) gi_d$hgnc_id else NULL
  clingen_d <- tryCatch(
    fetch_clingen_validity(gene_name, hgnc_id = hgnc_for_clingen),
    error = function(e) list(classification = NA_character_, moi = "", disease = "", source = "ClinGen LDH")
  )

  afs_dest <- file.path(getwd(), paste0(uid, "-F1-aa-substitutions.csv"))
  afs_d <- if (file.exists(afs_dest)) {
    tryCatch(data.table::fread(afs_dest), error = function(e) NULL)
  } else NULL

  # ---------------------------------------------------------------------------
  # Local-prediction pre-loads + fetch_dbnsfp() monkey-patch
  #
  # PRIMARY: single-source local dbNSFP 4.9a (analyses/tmp/...custombuild.bgz,
  # tabix-indexed). One per-gene tabix call; protein-key index resolves hgvsp
  # directly. Full predictor parity with the live app's MyVariant->dbNSFP path
  # (SIFT, PP2 HDIV/HVAR, LRT, MutationTaster, FATHMM, PROVEAN, MetaSVM/LR/RNN,
  # CADD, DANN, GERP_RS/NR, PhyloP/PhastCons multi-way + REVEL + AlphaMissense).
  #
  # FALLBACK: legacy REVEL+AM+UCSC synthesis preserved for graceful degradation
  # when dbNSFP misses a variant. UCSC fetch is the slow link (per-gene API
  # call); kept warm so the fallback works without cold-cache hits.
  #
  # The patched fetch_dbnsfp returns the dbNSFP hit when found; otherwise
  # synthesizes from REVEL+AM+UCSC; otherwise defers to the remote MyVariant
  # call (saved as fetch_dbnsfp_remote). on.exit restores the original binding
  # so the patch is scoped to this classify_gene() call.
  # ---------------------------------------------------------------------------
  exon_df_local <- harness_fetch("ensembl_exons", gene_name, fetch_ensembl_exons(gene_name))
  chrom_for_gene <- if (!is.null(exon_df_local) && nrow(exon_df_local) > 0) {
    as.character(exon_df_local$chr[1])
  } else NA_character_

  revel_dt_local <- if (!is.na(chrom_for_gene)) {
    tryCatch(load_revel_for_chrom(chrom_for_gene), error = function(e) NULL)
  } else NULL

  am_dt_local <- tryCatch(load_alphamissense_csv(uid), error = function(e) NULL)

  # dbNSFP 4.9a single-source pre-load (preferred path; covers SIFT, PP2,
  # MetaSVM/LR/RNN, CADD, DANN, GERP_RS, plus REVEL/AM/PhyloP/PhastCons —
  # full predictor parity with the live VarViz app via MyVariant->dbNSFP).
  # Per-gene region pulled by tabix once. The protein-key index built inside
  # load_dbnsfp_for_region lets the patched fetch_dbnsfp resolve hgvsp
  # directly without aa->genomic translation.
  dbnsfp_env_local <- if (!is.na(chrom_for_gene) && !is.null(exon_df_local) &&
                          nrow(exon_df_local) > 0L) {
    g_start <- min(exon_df_local$genomic_start, na.rm = TRUE)
    g_end   <- max(exon_df_local$genomic_end,   na.rm = TRUE)
    harness_fetch("dbnsfp_region", gene_name,
                  load_dbnsfp_for_region(chrom_for_gene, g_start, g_end))
  } else NULL
  # chrom_for_gene = NA (Ensembl down) silently disables the dbNSFP and REVEL
  # preloads: conservation reads GERP=NA / PhyloP~0, cons_strong flips, PM1/PP3
  # drop, row count stays correct. Record it so run_one() refuses the checkpoint.
  if (file.exists(DBNSFP_DEFAULT_PATH) && is.null(dbnsfp_env_local)) {
    options(varviz.localpred_failed =
              union(getOption("varviz.localpred_failed", character(0)), gene_name))
    cat(sprintf("    [%s] LOCAL PREDICTORS FAILED TO LOAD (chrom=%s) -- gene will not be checkpointed\n",
                gene_name, chrom_for_gene))
  }

  # Protein length: prefer AM CSV max position, fall back to universe max.
  prot_length_for_gene <- if (!is.null(am_dt_local) && nrow(am_dt_local) > 0) {
    suppressWarnings(max(as.integer(sub("^[A-Z*]([0-9]+).*$", "\\1",
                                        am_dt_local$protein_variant)),
                         na.rm = TRUE))
  } else {
    universe_pos <- extract_pos(universe$p_notation[universe$gene == gene_name])
    if (length(universe_pos) > 0L) max(universe_pos, na.rm = TRUE) else NA_integer_
  }

  ucsc_cons_df_local <- if (is.finite(prot_length_for_gene) && prot_length_for_gene > 0) {
    tryCatch(fetch_conservation_scores(gene_name, prot_length_for_gene),
             error = function(e) NULL)
  } else NULL

  # Hash UCSC conservation by aa_pos for O(1) per-variant lookup.
  ucsc_cons_env <- new.env(hash = TRUE, parent = emptyenv())
  if (!is.null(ucsc_cons_df_local) && nrow(ucsc_cons_df_local) > 0 &&
      "phylop100_raw" %in% colnames(ucsc_cons_df_local)) {
    for (i in seq_len(nrow(ucsc_cons_df_local))) {
      assign(as.character(ucsc_cons_df_local$aa_pos[i]),
             list(PhyloP_100V    = ucsc_cons_df_local$phylop100_raw[i],
                  PhyloP_470M    = ucsc_cons_df_local$phylop470_raw[i],
                  PhastCons_100V = ucsc_cons_df_local$phastcons_raw[i]),
             envir = ucsc_cons_env)
    }
  }

  cat(sprintf("    [%s] local caches: dbNSFP=%s, REVEL=%s, AM=%s, UCSC=%d aa\n",
              gene_name,
              if (!is.null(dbnsfp_env_local))
                sprintf("%d rows", length(ls(dbnsfp_env_local)))
              else "NA",
              if (!is.null(revel_dt_local))
                sprintf("chr%s/%.1fM", chrom_for_gene, nrow(revel_dt_local) / 1e6)
              else "NA",
              if (!is.null(am_dt_local)) sprintf("%d", nrow(am_dt_local)) else "NA",
              length(ls(ucsc_cons_env))))

  fetch_dbnsfp_remote <- get("fetch_dbnsfp", envir = globalenv())
  patched_fetch_dbnsfp <- function(gn, hgvsp) {
    aa_pos_int <- suppressWarnings(as.integer(sub("^p\\.[A-Z*]([0-9]+).*$", "\\1", hgvsp)))
    alt_aa_chr <- sub("^p\\.[A-Z*][0-9]+([A-Z*])$", "\\1", hgvsp)
    if (!nzchar(alt_aa_chr) || identical(alt_aa_chr, hgvsp)) alt_aa_chr <- NA_character_
    # Reference residue disambiguates rows that share (aapos, aaalt) across
    # transcripts -- without it a G51D query can be served an N51D row.
    ref_aa_chr <- sub("^p\\.([A-Z*])[0-9]+[A-Z*]$", "\\1", hgvsp)
    if (!nzchar(ref_aa_chr) || identical(ref_aa_chr, hgvsp)) ref_aa_chr <- NA_character_

    # PRIMARY PATH: single-source dbNSFP 4.9a lookup. Returns the full raw
    # MyVariant-shaped hit (all 458 dbNSFP cols → SIFT, PP2 HDIV/HVAR, LRT,
    # MutationTaster, FATHMM, PROVEAN, MetaSVM/LR/RNN, REVEL, AlphaMissense,
    # CADD, DANN, GERP_RS/NR, PhyloP/PhastCons multi-way). Apples-to-apples
    # with the live app's MyVariant→dbNSFP path.
    if (!is.null(dbnsfp_env_local) && !is.na(aa_pos_int) && !is.na(alt_aa_chr)) {
      hit <- tryCatch(
        lookup_dbnsfp_by_aa(dbnsfp_env_local, chrom_for_gene, aa_pos_int, alt_aa_chr,
                            aa_ref = ref_aa_chr),
        error = function(e) NULL
      )
      if (!is.null(hit)) return(hit)
    }

    # FALLBACK PATH: legacy REVEL+AM+UCSC synthesis (kept for graceful
    # degradation when dbNSFP missed a variant — e.g. region-edge cases or
    # variants outside the pre-loaded window).
    revel_score <- NA_real_
    if (!is.null(revel_dt_local) && !is.null(exon_df_local) &&
        !is.na(aa_pos_int) && !is.na(alt_aa_chr)) {
      gpos <- tryCatch(aa_to_genomic(aa_pos_int, exon_df_local),
                       error = function(e) NA_integer_)
      if (!is.na(gpos)) {
        strand <- exon_df_local$strand[1]
        cps <- if (isTRUE(strand == 1)) c(gpos, gpos + 1L, gpos + 2L)
                                         else c(gpos, gpos - 1L, gpos - 2L)
        revel_score <- tryCatch(
          lookup_revel_by_aa(revel_dt_local, chrom_for_gene, cps, alt_aa_chr),
          error = function(e) NA_real_
        )
      }
    }

    am_result <- if (!is.null(am_dt_local)) lookup_alphamissense(am_dt_local, hgvsp)
                 else list(score = NA_real_, class = NA_character_)

    cons <- list(PhyloP_100V = NA_real_, PhastCons_100V = NA_real_,
                 GERP_RS = NA_real_, PhyloP_470M = NA_real_)
    if (!is.na(aa_pos_int)) {
      k <- as.character(aa_pos_int)
      if (exists(k, envir = ucsc_cons_env, inherits = FALSE)) {
        hit <- get(k, envir = ucsc_cons_env, inherits = FALSE)
        cons$PhyloP_100V    <- hit$PhyloP_100V
        cons$PhastCons_100V <- hit$PhastCons_100V
        cons$PhyloP_470M    <- hit$PhyloP_470M
        # GERP_RS stays NA — UCSC tracks don't expose GERP. PM1_strong logic
        # has 4 other conservation gates so this is a tolerable miss.
      }
    }

    any_local <- !is.na(revel_score) || !is.na(am_result$score) ||
                 !is.na(cons$PhyloP_100V) || !is.na(cons$PhastCons_100V) ||
                 !is.na(cons$PhyloP_470M)
    if (!any_local) return(fetch_dbnsfp_remote(gn, hgvsp))

    # Return RAW MyVariant.info shape — parse_dbnsfp_scores() (server.R:2507)
    # consumes hit$dbnsfp via safe_extract_num(d, "<key>", "<sub>", "score") and
    # converts to the structured shape that build_variant_table() reads. Synthesizing
    # the parsed shape directly would bypass parse_dbnsfp_scores and break the rest
    # of the pipeline (verdict fields, classify_pred wrappers, etc.).
    list(
      dbnsfp = list(
        revel         = list(score = revel_score),
        alphamissense = list(am_pathogenicity = am_result$score,
                             am_class         = am_result$class),
        phylop = list(
          `100way_vertebrate`  = list(score = cons$PhyloP_100V),
          `470way_mammalian`   = list(score = cons$PhyloP_470M)
        ),
        phastcons = list(
          `100way_vertebrate` = list(score = cons$PhastCons_100V)
        )
        # gerp++ omitted — UCSC tracks don't expose it; safe_extract returns NA.
      )
    )
  }

  assign("fetch_dbnsfp", patched_fetch_dbnsfp, envir = globalenv())
  on.exit(assign("fetch_dbnsfp", fetch_dbnsfp_remote, envir = globalenv()), add = TRUE)

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

  # Per-gene analysis parameters. Defaults are the app defaults the manuscript
  # reports (Supplementary S1): monoallelic, prevalence 1 in 2,000, allelic
  # heterogeneity 0.2, genetic heterogeneity 1.0, penetrance 0.5 -> maximum
  # credible allele frequency 1.0e-4 (Whiffin 2017) and, over 251,496 alleles
  # at the 95% Poisson bound, maximum credible allele count 34. pop_size is an
  # ALLELE number (2N), not a count of individuals. A universe may override any
  # of these per gene with a column of the same name; SNCA, for example, is run
  # at prevalence 1 in 10,000 (max AF 2.0e-5, max AC 9).
  .col <- function(name, default) {
    if (!name %in% names(universe)) return(default)
    v <- universe[[name]][universe$gene == gene_name]
    v <- v[!is.na(v) & nzchar(as.character(v))]
    if (length(v) == 0) default else v[1]
  }
  .p <- list(
    inh_param         = as.character(.col("inh_param", "monoallelic")),
    af_cutoff         = as.numeric(.col("af_cutoff", 0.0001)),
    ac_cutoff         = as.numeric(.col("ac_cutoff", 34)),
    prevalence_1_in_n = as.numeric(.col("prevalence_1_in_n", 2000)),
    allelic_het       = as.numeric(.col("allelic_het", 0.2)),
    genetic_het       = as.numeric(.col("genetic_het", 1.0)),
    penetrance        = as.numeric(.col("penetrance", 0.5))
  )
  cat(sprintf("  [%s] params: inh=%s af_cutoff=%.3g ac_cutoff=%g prev=1/%s hetA=%s hetG=%s pen=%s\n",
              gene_name, .p$inh_param, .p$af_cutoff, .p$ac_cutoff,
              .p$prevalence_1_in_n, .p$allelic_het, .p$genetic_het, .p$penetrance))

  cat(sprintf("  [%s] uid=%s, %d variants, calling build_variant_table()...\n",
              gene_name, uid, nrow(hl)))
  t1 <- Sys.time()
  vtbl <- tryCatch(
    build_variant_table(
      hl, af_d, mean_d, afs_d, gnomad_d, clinvar_d,
      pfam_d, uniprot_d, ccrs_d,
      af_cutoff = .p$af_cutoff, ac_cutoff = .p$ac_cutoff,
      clinvar_missense = NULL, consurf_data = NULL,
      denovo_status     = "not_denovo",
      inh_param         = .p$inh_param,
      cutoff_method     = "calc_af",
      prevalence_1_in_n = .p$prevalence_1_in_n,
      allelic_het = .p$allelic_het, genetic_het = .p$genetic_het,
      penetrance = .p$penetrance,
      pop_size = 251496, conf_interval = 0.95,
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
  before <- sentinel_snapshot()
  result <- classify_gene(gene_name)
  failed <- new_sentinels(before, sentinel_snapshot())
  # A source outage during this gene drops PM1 pathways or predictors while
  # leaving the row count intact. Refuse the checkpoint so a re-run redoes the
  # gene instead of caching bad numbers.
  if (length(failed)) {
    cat(sprintf("  [%s] ABORT: required data source failed (%s); checkpoint NOT written\n",
                gene_name, paste(failed, collapse = ", ")))
    return(NULL)
  }
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

if (length(all_results) < length(genes)) {
  cat(sprintf("[harness] FAIL: %d of %d genes aborted or skipped; summary NOT written (re-run to retry them)\n",
              length(genes) - length(all_results), length(genes)))
  quit(status = 1)
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
