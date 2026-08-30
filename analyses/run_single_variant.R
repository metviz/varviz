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
#   Rscript analyses/run_single_variant.R --self-test      # arg parsing only, no network
#
# Defaults match the app's own: monoallelic, prevalence 1 in 2000, hetA 0.2,
# hetG 1, penetrance 0.5, pop 251496, CI 0.95, af_cutoff 1e-4, ac_cutoff 34.
# Run from the project root — server.R sources by relative path.

DEFAULTS <- list(inh = "monoallelic", prev = 2000, hetA = 0.2, hetG = 1,
                 pen = 0.5, pop = 251496, ci = 0.95,
                 af_cutoff = 1e-04, ac_cutoff = 34, out = "")

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

vtbl <- build_variant_table(
  hl,
  extract_alphafold_plddt(uid), get_mean_pathogenicity(uid),
  if (file.exists(afs_p)) data.table::fread(afs_p) else NULL,
  extract_gnomad(p$gene), extract_clinvar(p$gene),
  pfam_d, extract_uniprot_feature_data(uid), extract_ccrs(p$gene, pfam_d$primaryAccession),
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

if (nzchar(p$out)) {
  write.table(vtbl, p$out, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
  cat("\nwrote ", p$out, " (", nrow(vtbl), " rows x ", ncol(vtbl), " cols)\n", sep = "")
}
