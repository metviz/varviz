# Validate the permutation-significant clinvar_hotspot PM1 change on TP53.
# Sources server.R, fetches TP53 gene-level data, classifies DBD-hotspot vs
# outside-DBD positions, prints PM1 tag + pathway. Isolation: read-only test.
suppressMessages(source("server.R"))
source("analyses/lib/clinvar_blind.R")

gene <- "TP53"
uid  <- gene_data$uniprot_id[gene_data$gene_name == gene][1]
L    <- gene_data$mane_prot_len[gene_data$gene_name == gene][1]

pfam_d    <- tryCatch(extract_pfam(uid), error=function(e) NULL)
uniprot_d <- tryCatch(extract_uniprot_feature_data(uid), error=function(e) NULL)
gnomad_d  <- tryCatch(extract_gnomad(gene), error=function(e) NULL)
clinvar_d <- tryCatch(extract_clinvar(gene), error=function(e) NULL)
af_d      <- tryCatch(extract_alphafold_plddt(uid), error=function(e) NULL)
mean_d    <- tryCatch(get_mean_pathogenicity(uid), error=function(e) NULL)

# quick look at the profile itself
prof <- clinvar_hotspot_profile(gene, clinvar_d, gnomad_d, L, 1e-4)
cat(sprintf("profile: %s, significant residues = %d / %d\n",
            if (is.null(prof)) "NULL(fallback)" else "computed",
            if (is.null(prof)) NA else sum(prof), L))
for (p in c(175, 248, 36, 334))
  cat(sprintf("  pos %d significant? %s\n", p, if (is.null(prof)) "n/a" else prof[p]))

variants <- c("p.R175H","p.R248Q","p.P36S","p.G334V")
hl <- data.frame(Mutation=variants,
                 prot_pos=as.integer(sub("^p\\.[A-Z]([0-9]+)[A-Z]$","\\1",variants)),
                 gene=gene, stringsAsFactors=FALSE)
vt <- build_variant_table(hl, af_d, mean_d, NULL, gnomad_d, clinvar_d,
        pfam_d, uniprot_d, NULL, af_cutoff=1e-4, ac_cutoff=13,
        clinvar_missense=NULL, consurf_data=NULL, denovo_status="not_denovo",
        inh_param="monoallelic", cutoff_method="calc_af", prevalence_1_in_n=2000,
        allelic_het=0.5, genetic_het=1.0, penetrance=1.0, pop_size=125748,
        conf_interval=0.95, clingen_disease_param="", clingen_moi_param="",
        consurf_file_name="")
cat("\n=== classifications ===\n")
for (i in seq_len(nrow(vt)))
  cat(sprintf("  %-9s  tags=[%s]  pm1_pathway=%s\n",
      vt$Variant[i], vt$ACMG_Tags[i], vt$ACMG_PM1_Pathway[i]))
