#!/usr/bin/env Rscript
# run_single_variant.R — classify variants headlessly, exactly as the app does.
#
# Drives build_variant_table(), the same function server.R calls, so the output
# is what the UI would render. Useful for checking a fix without clicking
# through the browser, and for diffing a variant across releases.
#
# Usage:
#   Rscript analyses/run_single_variant.R SLC6A1 p.Ala305Thr,p.Ala305Val
#   Rscript analyses/run_single_variant.R CASR p.Arg185Gln --inh=biallelic --prev=100000
#   Rscript analyses/run_single_variant.R POLR3A p.Gly854Arg --out=/tmp/polr3a.tsv
#   Rscript analyses/run_single_variant.R SLC6A1 p.Ala305Thr --plot=/tmp/slc6a1.png
#   Rscript analyses/run_single_variant.R --self-test      # arg parsing only, no network
#
# Options
#   --inh=       monoallelic | biallelic          (default monoallelic)
#   --prev=      prevalence, 1 in N               (2000)
#   --hetA=      allelic heterogeneity            (0.2)
#   --hetG=      genetic heterogeneity            (1)
#   --pen=       penetrance                       (0.5)
#   --pop=       population size                  (251496)
#   --ci=        confidence interval              (0.95)
#   --af_cutoff= max credible allele frequency    (1e-04)
#   --ac_cutoff= max credible allele count        (34)
#   --out=       write the variant table as TSV, with the same commented
#                run-parameter header the app export writes
#   --plot=      write the track stack as .png or .pdf (pLDDT, AF pathogenicity,
#                gnomAD frequency, mutation density, ClinVar/PTM/CCRS/Hot,
#                domain bar) -- the same panels and builders the app uses
#   --width=     plot width in inches             (11)
#   --height=    plot height in inches            (14)
#   --dpi=       raster resolution                (150)
#   --label=     variant label style passed to pfamplot()
#
# Defaults match the app's own: monoallelic, prevalence 1 in 2000, hetA 0.2,
# hetG 1, penetrance 0.5, pop 251496, CI 0.95, af_cutoff 1e-4, ac_cutoff 34.
# Run from the project root — server.R sources by relative path.

DEFAULTS <- list(inh = "monoallelic", prev = 2000, hetA = 0.2, hetG = 1,
                 pen = 0.5, pop = 251496, ci = 0.95,
                 af_cutoff = 1e-04, ac_cutoff = 34,
                 out = "", plot = "", width = 11, height = 14, dpi = 150,
                 label = "")

# Flags are --key=value; the two positionals are gene and a comma-separated
# variant list. Numeric defaults decide which values get coerced, so a typo in a
# flag name is an error rather than a silently ignored argument.
parse_args <- function(argv) {
  p <- DEFAULTS
  flags <- grepl("^--", argv)
  for (a in argv[flags]) {
    kv <- sub("^--", "", a)
    if (!grepl("=", kv)) next
    k <- sub("=.*$", "", kv); v <- sub("^[^=]*=", "", kv)
    if (!k %in% names(p)) stop("unknown flag --", k, " (known: ",
                               paste(names(p), collapse = ", "), ")")
    p[[k]] <- if (is.numeric(DEFAULTS[[k]])) as.numeric(v) else v
  }
  pos <- argv[!flags]
  if (length(pos) < 2) stop("need a gene and at least one variant")
  p$gene <- pos[1]
  p$variants <- trimws(strsplit(paste(pos[-1], collapse = ","), ",")[[1]])
  p$variants <- p$variants[nzchar(p$variants)]
  p
}

argv <- commandArgs(trailingOnly = TRUE)

if ("--self-test" %in% argv) {
  a <- parse_args(c("SLC6A1", "p.Ala305Thr,p.Ala305Val"))
  stopifnot(a$gene == "SLC6A1", length(a$variants) == 2, a$inh == "monoallelic", a$prev == 2000)
  b <- parse_args(c("CASR", "p.Arg185Gln", "--inh=biallelic", "--prev=1e5"))
  stopifnot(b$inh == "biallelic", b$prev == 1e5, is.numeric(b$prev))
  c3 <- parse_args(c("X", "p.A1B", "p.C2D"))          # variants as separate words
  stopifnot(length(c3$variants) == 2)
  stopifnot(inherits(try(parse_args(c("X", "p.A1B", "--nope=1")), silent = TRUE), "try-error"))
  d <- parse_args(c("X", "p.A1B", "--plot=/tmp/x.png", "--dpi=300"))
  stopifnot(d$plot == "/tmp/x.png", d$dpi == 300, is.numeric(d$dpi))
  cat("self-test OK\n"); quit(status = 0)
}

p <- parse_args(argv)
suppressMessages(source("server.R"))
`%||%` <- function(a, b) if (is.null(a)) b else a

uid <- as.character(gene_data$uniprot_id[gene_data$gene_name == p$gene][1])
if (is.na(uid) || !nzchar(uid)) stop("no UniProt accession for gene ", p$gene)
cat(sprintf("VarViz %s | %s (%s) | %d variant(s) | %s, 1 in %s, penetrance %s\n\n",
            VARVIZ_VERSION, p$gene, uid, length(p$variants), p$inh,
            format(p$prev, big.mark = ","), p$pen))

pfam_d    <- extract_pfam(uid)
if (is.null(pfam_d)) stop("UniProt fetch failed for ", uid, " — cannot score without prot_len")

# Same gate the app applies: a variant whose stated reference residue is not
# what the canonical sequence carries is scored at the wrong position.
mm <- check_reference_residues(p$variants, pfam_d$sequence$value)
if (!is.null(mm)) {
  for (i in seq_len(nrow(mm)))
    cat(sprintf("SKIPPED %s — position %d carries %s, not %s\n",
                mm$variant[i], mm$pos[i], mm$found[i], mm$stated[i]))
  p$variants <- setdiff(p$variants, mm$variant)
  if (!length(p$variants)) stop("no variants left after the reference check")
  cat("\n")
}

hl <- data.frame(Mutation = p$variants,
                 prot_pos = sapply(p$variants, extract_protein_position),
                 gene = p$gene, stringsAsFactors = FALSE)
afs_p <- am_csv_path(uid)
cg    <- fetch_clingen_validity(p$gene)

af_d      <- extract_alphafold_plddt(uid)
mean_d    <- get_mean_pathogenicity(uid)
afs_d     <- if (file.exists(afs_p)) data.table::fread(afs_p) else NULL
gnomad_d  <- extract_gnomad(p$gene)
clinvar_d <- extract_clinvar(p$gene)
uniprot_d <- extract_uniprot_feature_data(uid)
ccrs_d    <- extract_ccrs(p$gene, pfam_d$primaryAccession)

vtbl <- build_variant_table(
  hl, af_d, mean_d, afs_d, gnomad_d, clinvar_d,
  pfam_d, uniprot_d, ccrs_d,
  af_cutoff = p$af_cutoff, ac_cutoff = p$ac_cutoff,
  clinvar_missense = NULL, consurf_data = NULL, denovo_status = "not_denovo",
  inh_param = p$inh, cutoff_method = "calc_af",
  prevalence_1_in_n = p$prev, allelic_het = p$hetA, genetic_het = p$hetG,
  penetrance = p$pen, pop_size = p$pop, conf_interval = p$ci,
  clingen_disease_param = cg$disease %||% "", clingen_moi_param = cg$moi %||% "",
  consurf_file_name = "")

SHOW <- c("ClinVar", "ClinVar_Stars", "ClinVar_Name", "ClinVar_Trait", "ClinVar_Match",
          "gnomAD_AF", "gnomAD_AC", "CCRS", "Domain", "MDS_Score",
          "ACMG_Tags", "ACMG_PM1_Pathway", "PM1_Derivation",
          "ClinGen_Disease", "ClinGen_MOI", "ClinGen_Class",
          "Final_ACMG_Tags", "Final_Classification", "Final_Points")
cat("\n")
for (i in seq_len(nrow(vtbl))) {
  cat("=== ", as.character(vtbl$Variant[i]), " ===\n", sep = "")
  for (f in intersect(SHOW, colnames(vtbl))) {
    v <- as.character(vtbl[[f]][i])
    if (length(v) && !is.na(v) && nzchar(v)) cat(sprintf("  %-20s %s\n", f, v))
  }
}

# ── TSV, with the app's own run-parameter header ────────────────────────────
# Same shape the download button writes, so a headless run and an app export can
# be diffed directly. Empty parameters are dropped, as they are in the app.
if (nzchar(p$out)) {
  hdr_val <- function(col) {
    if (!col %in% colnames(vtbl) || nrow(vtbl) == 0) return("")
    v <- vtbl[[col]][1]
    if (is.null(v) || is.na(v)) "" else as.character(v)
  }
  meta <- c(
    "varviz_version"        = VARVIZ_VERSION,
    "run"                   = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"),
    "source"                = "analyses/run_single_variant.R",
    "gene"                  = p$gene,
    "variants_requested"    = paste(p$variants, collapse = ", "),
    "variants"              = as.character(nrow(vtbl)),
    "inheritance"           = hdr_val("Analysis_Inheritance"),
    "prevalence"            = hdr_val("Analysis_Prevalence"),
    "allelic_heterogeneity" = hdr_val("Analysis_Allelic_Het"),
    "genetic_heterogeneity" = hdr_val("Analysis_Genetic_Het"),
    "penetrance"            = hdr_val("Analysis_Penetrance"),
    "af_cutoff"             = hdr_val("Analysis_AF_Cutoff"),
    "ac_cutoff"             = hdr_val("Analysis_AC_Cutoff"),
    "clingen_disease"       = hdr_val("ClinGen_Disease"),
    "clingen_moi"           = hdr_val("ClinGen_MOI"),
    "clingen_classification"= hdr_val("ClinGen_Class"))
  meta <- meta[nzchar(meta)]
  con <- file(p$out, open = "wt"); on.exit(close(con), add = TRUE)
  writeLines(c("# VarViz variant summary",
               if (length(meta)) paste0("# ", names(meta), ": ", unname(meta)), "#"), con)
  write.table(vtbl, con, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat("\nwrote ", p$out, " (", nrow(vtbl), " rows x ", ncol(vtbl), " cols)\n", sep = "")
}

# ── Plots, using the app's own builders ─────────────────────────────────────
if (nzchar(p$plot)) {
  L <- pfam_data_len <- pfam_d$sequence$length
  ud <- uniprot_d
  ptm_types <- c("mod_res", "lipid", "carbohyd", "crosslnk", "disulfid")
  gene_ptm_data <- if (!is.null(ud) && nrow(ud) > 0) {
    pr <- ud[ud$type %in% ptm_types, , drop = FALSE]
    if (nrow(pr) > 0)
      data.frame(uniprot_id = pr$uniprot_id, location = pr$start,
                 final_ptm_group = ifelse(!is.na(pr$mod_res_group), pr$mod_res_group, pr$type),
                 stringsAsFactors = FALSE)
    else data.frame(uniprot_id = character(), location = numeric(),
                    final_ptm_group = character(), stringsAsFactors = FALSE)
  } else data.frame(uniprot_id = character(), location = numeric(),
                    final_ptm_group = character(), stringsAsFactors = FALSE)

  # Same "Hot" row the app draws: permutation-significant ClinVar hotspots.
  ps_positions <- tryCatch({
    L_ps <- suppressWarnings(as.integer(gene_data$mane_prot_len[gene_data$gene_name == p$gene][1]))
    if (is.na(L_ps) || L_ps < 1) L_ps <- L
    prof <- clinvar_hotspot_profile(p$gene, clinvar_d, gnomad_d, L_ps, p$af_cutoff)
    if (!is.null(prof)) which(prof) else integer(0)
  }, error = function(e) integer(0))

  panels <- list(
    tryCatch(plot_pLDDT(af_d, hl, prot_length = L), error = function(e) NULL),
    tryCatch(plot_afmps(mean_d, hl, prot_length = L), error = function(e) NULL),
    if (!is.null(gnomad_d) && nrow(gnomad_d) > 0)
      tryCatch(gnomad_freqplot(p$gene, p$af_cutoff, hl, prot_length = L), error = function(e) NULL),
    tryCatch(densityplot(gnomad_d, clinvar_d, L, p$ac_cutoff, hl, gene_ptm_data, ccrs_d, NULL),
             error = function(e) NULL),
    tryCatch(clinvar_ccrsplot(pfam_d, uniprot_d, clinvar_d, gene_ptm_data, ccrs_d,
                              highlight = hl, ps_positions = ps_positions),
             error = function(e) NULL),
    tryCatch(pfamplot(pfam_d, uniprot_d, clinvar_d, hl, p$label, for_plotly = FALSE),
             error = function(e) NULL))
  panels <- panels[!sapply(panels, is.null)]
  if (!length(panels)) stop("every panel failed to build")

  aligned <- cowplot::align_plots(plotlist = panels, align = "v", axis = "lr")
  gg <- cowplot::plot_grid(plotlist = aligned, ncol = 1)
  ggplot2::ggsave(p$plot, gg, width = p$width, height = p$height, dpi = p$dpi,
                  limitsize = FALSE)
  cat("wrote ", p$plot, " (", length(panels), " panels, ", p$width, "x", p$height,
      " in @ ", p$dpi, " dpi)\n", sep = "")
}
