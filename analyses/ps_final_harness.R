# Dual-pass classification harness — FINAL shipping build.
#
# Run on server.R patched 2026-08-21 12:27 (req_timeout(30) + req_retry(3) on
# the load-bearing UniProt fetch). Earlier runs -- ps_mds_corroborate and
# ps_mds_consfree -- are contaminated by silent failures of that call and must
# not be used; see analyses/repro/README.md "OPEN ISSUE".
#
# Config (submission, 2026-09-06): MDS Moderate-only, scored for every variant
# (originate + corroborate). Corroboration takes the STRONGER of pathway PM1
# and MDS, never the sum; the LR tiers (+3/+4) stay off unless
# VARVIZ_MDS_TIERED=TRUE is set for the exploratory calibration run.
#
# Differs from analyses/05_classify_harness.R only in configuration and in the
# sharding / VARVIZ_* env interface. The ClinVar-hotspot upgrade records its
# pre-upgrade tag so strip_clinvar_tags() can demote it under Pass-Blind.
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
# Resumability: per-gene checkpoint files at analyses/ps_final/classifications/<gene>__dual.tsv.
# Re-running skips genes whose checkpoint already exists.
#
# Run from project root: `Rscript analyses/05_classify_harness.R`
#
# Output (gitignored per Pass 2 policy):
#   analyses/ps_final/classifications/<gene>__dual.tsv  (per-gene checkpoint)
#   analyses/ps_final/summary.tsv (concatenated summary)

suppressMessages({
  library(dplyr)
  library(readr)
  library(purrr)
})

cat("[harness] Sourcing server.R...\n")
t0 <- Sys.time()
suppressMessages(source("server.R"))
# VARVIZ_MDS_PM1=FALSE runs the MDS-off ablation through THIS harness rather
# than ps_nomds_harness.R, so the ablation gets the same four data-source
# sentinels. A comparator scored without them is not comparable: the stale
# analyses/ps_nomds run (Jul 30) has no run.log and cannot be audited at all.
.env_flag <- function(name, default) {
  v <- toupper(Sys.getenv(name, default))
  if (!v %in% c("TRUE", "FALSE")) stop(name, " must be TRUE or FALSE, got: ", Sys.getenv(name))
  v == "TRUE"
}
# PM2 weight sensitivity: 1 = ClinGen SVI Supporting (default), 2 = Richards
# 2015 Moderate. Lets the SVI change be measured with everything else fixed.
.pm2_pts <- suppressWarnings(as.integer(Sys.getenv("VARVIZ_PM2_POINTS", "1")))
if (is.na(.pm2_pts) || !.pm2_pts %in% c(1L, 2L))
  stop("VARVIZ_PM2_POINTS must be 1 or 2, got: ", Sys.getenv("VARVIZ_PM2_POINTS"))
options(varviz.pm2_points = .pm2_pts)

# Restores the retired AlphaMissense+REVEL PP3->Strong upgrade for sensitivity
# analysis. FALSE in the shipped configuration.
options(varviz.pp3_ps3_proxy = .env_flag("VARVIZ_PP3_PROXY", "FALSE"))
# VARVIZ_META_CONSENSUS=FALSE disables the 2-of-3 MetaSVM/MetaLR/MetaRNN branch
# that raises PP3 to moderate, so its contribution can be measured against a run
# that is otherwise identical.
options(varviz.pp3_meta_consensus = .env_flag("VARVIZ_META_CONSENSUS", "TRUE"))
# VARVIZ_AM_CALIBRATED=TRUE replaces AlphaMissense's developer threshold (0.564
# supporting, 0.34 benign) with the Bergquist 2025 ClinGen SVI calibration
# (supporting 0.792, moderate 0.906, strong 0.990; BP4 supporting 0.169).
options(varviz.am_calibrated = .env_flag("VARVIZ_AM_CALIBRATED", "FALSE"))
# VARVIZ_AM_NO_OVERRIDE=TRUE keeps the calibrated thresholds but restores the
# old precedence: AlphaMissense speaks only when no other tool has given PP3.
# The difference against the plain calibrated run is the override's share.
options(varviz.am_no_override = .env_flag("VARVIZ_AM_NO_OVERRIDE", "FALSE"))
# VARVIZ_AM_HYBRID=TRUE keeps the developer threshold at the supporting rung and
# lets only the calibrated moderate-and-above rungs raise a level another tool
# already set - the half of the calibration that gains true positives.
options(varviz.am_hybrid = .env_flag("VARVIZ_AM_HYBRID", "TRUE"))
# VARVIZ_PP3_3PT=TRUE enables the 3-point PP3 rung between Moderate and Strong,
# reported as a calibrated interval by Pejaver 2022 and Bergquist 2025 and
# expected to enter a future edition of the guidelines.
options(varviz.pp3_3pt = .env_flag("VARVIZ_PP3_3PT", "FALSE"))

.mds_on     <- .env_flag("VARVIZ_MDS_PM1",    "TRUE")   # MDS scored for every variant
.mds_tiered <- .env_flag("VARVIZ_MDS_TIERED", "FALSE")  # exploratory +3/+4 tiers; off for submission
options(varviz.mds_pm1    = .mds_on)
options(varviz.mds_tiered = .mds_on && .mds_tiered)
cat(sprintf("[harness] varviz.mds_pm1 = %s ; varviz.mds_tiered = %s ; PM2 = %d pt\n",
            .mds_on, .mds_on && .mds_tiered, .pm2_pts))
cat(sprintf("[harness] Sourced server.R in %.1f sec\n",
            as.numeric(Sys.time() - t0, units = "secs")))

source("analyses/lib/clinvar_blind.R")
source("analyses/lib/local_predictors.R")
source("analyses/lib/harness_guard.R")

UNIVERSE_IN    <- Sys.getenv("VARVIZ_UNIVERSE", "analyses/derived/variant_universe_gnomad.tsv")
.out_dir <- Sys.getenv("VARVIZ_OUT_DIR", "analyses/ps_final")
CHECKPOINT_DIR <- file.path(.out_dir, "classifications")
SUMMARY_OUT    <- file.path(.out_dir, "summary.tsv")
# summary.tsv is the FINAL shipping build the manuscript numbers were read
# from. Overwriting it requires the explicit opt-in VARVIZ_FORCE=1.
.force <- identical(Sys.getenv("VARVIZ_FORCE", ""), "1")
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

universe <- read_tsv(UNIVERSE_IN, show_col_types = FALSE)
genes <- sort(unique(universe$gene))
# Shard support: VARVIZ_GENES="A,B" restricts this process to those genes.
# Per-gene checkpoints are disjoint files, so shards never collide.
.shard <- Sys.getenv("VARVIZ_GENES", "")
if (nzchar(.shard)) {
  want <- trimws(strsplit(.shard, ",")[[1]])
  miss <- setdiff(want, genes)
  if (length(miss)) stop("VARVIZ_GENES names genes absent from universe: ", paste(miss, collapse=", "))
  genes <- want
  cat(sprintf("[harness] SHARD: %s\n", paste(genes, collapse=", ")))
}
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
  # Every failure is recorded in options(varviz.harness_fetch_failed) so
  # run_one() refuses the checkpoint (analyses/lib/harness_guard.R). The old
  # sentinel list covered only uniprot / localpred / clinvar_batch.
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

  afs_dest <- am_csv_path(uid)  # server.R helper: app dir if present, else scratch
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
    tryCatch(load_dbnsfp_for_region(chrom_for_gene, g_start, g_end),
             error = function(e) NULL)
  } else NULL

  # Protein length. UniProt's sequence length is authoritative and pfam_d is
  # already fetched above; the AlphaMissense CSV agrees with it because it holds
  # every substitution of every residue. The old order tried AM first and then
  # fell back to the highest variant position in the universe, which equals the
  # protein length only by accident. On the 226-variant RASopathy cohort, where
  # 15 of 16 genes had no AlphaMissense file, that truncated every one of them
  # (SOS1 to 1237 aa instead of 1333, MRAS to 71 instead of 208), and
  # fetch_conservation_scores() then built its conservation array to the wrong
  # length. Fall back to the universe maximum only when both real sources are
  # gone, and say so loudly, because everything downstream is position-indexed.
  prot_length_for_gene <- {
    len_uniprot <- suppressWarnings(as.integer(
      if (is.null(pfam_d) || is.null(pfam_d$length)) NA else pfam_d$length))[1]
    len_am <- if (!is.null(am_dt_local) && nrow(am_dt_local) > 0) {
      suppressWarnings(max(as.integer(sub("^[A-Z*]([0-9]+).*$", "\\1",
                                          am_dt_local$protein_variant)),
                           na.rm = TRUE))
    } else NA_integer_
    if (isTRUE(!is.na(len_uniprot) && len_uniprot > 0)) {
      len_uniprot
    } else if (isTRUE(is.finite(len_am) && len_am > 0)) {
      len_am
    } else {
      universe_pos <- extract_pos(universe$p_notation[universe$gene == gene_name])
      len_universe <- if (length(universe_pos) > 0L) max(universe_pos, na.rm = TRUE) else NA_integer_
      cat(sprintf("    [%s] WARNING: no UniProt or AlphaMissense protein length; using the universe maximum (%s aa), which truncates every position-indexed annotation past it\n",
                  gene_name, as.character(len_universe)))
      len_universe
    }
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
  record_cache_degradation(gene_name, dbnsfp_env_local, revel_dt_local, am_dt_local,
                           length(ls(ucsc_cons_env)), prot_length_for_gene)

  # Silent-degradation guard, third vector (after the two UniProt fetches).
  # fetch_ensembl_exons() failing -- rest.ensembl.org times out under parallel
  # load -- leaves chrom_for_gene = NA, which silently disables BOTH the dbNSFP
  # and REVEL preloads. Conservation then reads as GERP=NA / PhyloP~0, flipping
  # cons_strong and dropping PM1/PP3, while the row count stays correct. This is
  # what corrupted BRCA1 in the first sharded ps_final run (1,326 classification
  # changes, no failure line anywhere). Record it so run_one() refuses the
  # checkpoint rather than caching a gene scored without its predictors.
  if (file.exists(DBNSFP_DEFAULT_PATH) && is.null(dbnsfp_env_local)) {
    options(varviz.localpred_failed =
              union(getOption("varviz.localpred_failed", character(0)), gene_name))
    cat(sprintf("    [%s] LOCAL PREDICTORS FAILED TO LOAD (chrom=%s) -- gene will not be checkpointed\n",
                gene_name, chrom_for_gene))
  }

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

  # Disease-model parameters. A universe may carry per-gene values in optional
  # columns (inh_param, af_cutoff, prevalence_1_in_n, allelic_het, genetic_het,
  # penetrance); where a column is absent the historical constant is used, so
  # universes without them classify exactly as before. These feed the maximum
  # credible allele frequency and therefore PM2, and inh_param additionally
  # switches the BS1/PM2 thresholds (server.R:5118-5144), so a recessive gene
  # run as monoallelic over-fires PM2.
  .col <- function(name, default) {
    if (!name %in% names(universe)) return(default)
    v <- universe[[name]][universe$gene == gene_name]
    v <- v[!is.na(v) & nzchar(as.character(v))]
    if (length(v) == 0) default else v[1]
  }
  .p <- list(
    inh_param         = as.character(.col("inh_param", "monoallelic")),
    af_cutoff         = as.numeric(.col("af_cutoff", 0.0001)),
    prevalence_1_in_n = as.numeric(.col("prevalence_1_in_n", 2000)),
    allelic_het       = as.numeric(.col("allelic_het", 0.2)),
    genetic_het       = as.numeric(.col("genetic_het", 1.0)),
    penetrance        = as.numeric(.col("penetrance", 0.5))
  )
  cat(sprintf("  [%s] params: inh=%s af_cutoff=%.3g prev=1/%s hetA=%s hetG=%s pen=%s\n",
              gene_name, .p$inh_param, .p$af_cutoff, .p$prevalence_1_in_n,
              .p$allelic_het, .p$genetic_het, .p$penetrance))

  cat(sprintf("  [%s] uid=%s, %d variants, calling build_variant_table()...\n",
              gene_name, uid, nrow(hl)))
  t1 <- Sys.time()
  vtbl <- tryCatch(
    build_variant_table(
      hl, af_d, mean_d, afs_d, gnomad_d, clinvar_d,
      pfam_d, uniprot_d, ccrs_d,
      af_cutoff = .p$af_cutoff, ac_cutoff = 34,
      clinvar_missense = NULL, consurf_data = NULL,
      denovo_status     = "not_denovo",
      inh_param         = .p$inh_param,
      cutoff_method     = "calc_af",
      prevalence_1_in_n = .p$prevalence_1_in_n,
      allelic_het = .p$allelic_het, genetic_het = .p$genetic_het,
      penetrance = .p$penetrance,
      # pop_size is an ALLELE number (2N) over 125,748 gnomAD exomes, not a
      # count of individuals; with af_cutoff 1.0e-4 it gives max credible AC 34.
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
    # na = character(): a cached gene must round-trip unchanged. read_tsv's
    # default turns the empty pm1_pathway of a no-PM1 variant into NA, and
    # write_tsv then serialises it as the literal "NA" -- so a concatenation
    # pass silently rewrote 5,714 empty cells and broke the nchar() > 0 filter
    # every PM1-pathway tally uses. Keeping "" as "" makes the pass idempotent.
    df <- read_tsv(ckpt, show_col_types = FALSE, na = character())
    # A column that is blank for every row of this gene (pm1_pathway when the
    # gene has no PM1 evidence anywhere) is type-guessed as logical, and
    # write_tsv then serialises it as the literal "NA". Every column written
    # here is character or integer and integers are never blank, so any logical
    # column is that artefact: restore the empty strings that are on disk.
    lg <- vapply(df, is.logical, logical(1))
    if (any(lg)) df[lg] <- lapply(df[lg], function(x) rep("", length(x)))
    return(df)
  }
  before <- sentinel_snapshot()          # all four sentinel families
  result <- classify_gene(gene_name)
  after  <- sentinel_snapshot()
  # A UniProt outage during this gene drops the domain/site PM1 pathways while
  # leaving the row count intact -- the exact silent corruption that spoiled
  # ps_mds_corroborate and ps_mds_consfree. Refuse the checkpoint so a re-run
  # redoes the gene instead of caching bad numbers.
  if (length(setdiff(after, before))) {
    cat(sprintf("  [%s] ABORT: required data source failed (%s); checkpoint NOT written\n",
                gene_name, paste(setdiff(after, before), collapse = ", ")))
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
  cat(sprintf("[harness] FAIL: %d of %d genes aborted or skipped; nothing written (re-run to retry them)\n",
              length(genes) - length(all_results), length(genes)))
  quit(status = 1)
}

summary_df <- bind_rows(all_results)
# A shard holds only its own genes, so writing SUMMARY_OUT here would clobber
# the file with a partial run. Only the final unsharded pass -- which finds
# every gene cached and concatenates them -- is allowed to write it.
if (nzchar(.shard)) {
  cat(sprintf("[harness] SHARD done (%s); summary.tsv left for the final pass\n", .shard))
  quit(status = 0)
}
# A degraded predictor cache lowers PP3 strength without changing the row count,
# so it must never pass unnoticed. VARVIZ_ALLOW_DEGRADED=1 is the explicit opt-in.
.degraded <- cache_degradations()
if (length(.degraded)) {
  cat("[harness] DEGRADED PREDICTOR CACHES:\n")
  for (d in .degraded) cat("  - ", d, "\n", sep = "")
  if (!identical(Sys.getenv("VARVIZ_ALLOW_DEGRADED", ""), "1")) {
    cat("[harness] REFUSE: predictor caches were incomplete for the genes above;\n",
        "          set VARVIZ_ALLOW_DEGRADED=1 to accept and record this.\n", sep = "")
    quit(status = 1)
  }
  cat("[harness] VARVIZ_ALLOW_DEGRADED=1: proceeding with the gaps above recorded.\n")
}

if (file.exists(SUMMARY_OUT) && !.force) {
  cat(sprintf("[harness] REFUSE: %s already exists; set VARVIZ_FORCE=1 to overwrite\n", SUMMARY_OUT))
  quit(status = 1)
}
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
