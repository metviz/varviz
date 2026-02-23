suppressMessages(library(dplyr))
suppressMessages(library(shiny))
suppressMessages(library(Cairo)) 
suppressMessages(library(shinyjs))
suppressMessages(library(optparse))
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))
suppressMessages(library(RCurl))
suppressMessages(library(XML))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(zoo))
suppressMessages(library(ggrepel))
suppressMessages(library(data.table))
suppressMessages(library(stringi))
suppressMessages(library(reshape2))
suppressMessages(library(bedr))
suppressMessages(library(plotly))
suppressMessages(library(httr))
suppressMessages(library(httr2))
suppressMessages(library(readr))
suppressMessages(library(webshot))
suppressMessages(library(htmlwidgets))
theme_set(theme_cowplot(font_size=12)) 

#rm(list=ls())

file_path="/data"

load(paste(file_path,"MutEViz.RData",sep="/"))

gene_list = sort(gene_data$gene_name)

plotsetA <- c("clinvar","density","freq","depth")
plotsetB <- c("clinvar","freq","depth")
plotsetC <- c("clinvar","density","depth") 
plotsetD <- c("clinvar","density","freq") 
plotsetE <- c("density","freq","depth") 
plotsetF <- c("clinvar","density")
plotsetG <- c("clinvar","freq") 
plotsetH <- c("clinvar","depth") 
plotsetI <- c("density","freq")
plotsetJ <- c("density","depth") 
plotsetK <- c("freq","depth")

curl.cainfo = "/data/cacert.pem"

#function to extract the protein position
extract_protein_position <- function(x){ as.numeric(lapply(regmatches(x , gregexpr("[[:digit:]]+", x)), `[`, 1)) }

print(x)

# Function to fetch UniProt features for a given UniProt ID
extract_pfam <- function(x) { 
  tryCatch({
    # Build the request
    resp <- request("https://rest.uniprot.org/uniprotkb/") %>%
      req_url_path_append(paste0(x, ".json")) %>%
      req_url_query(fields = "accession,id,protein_name,gene_primary,ft_repeat,ft_region,ft_domain,ft_compbias,ft_coiled,ft_motif,ft_zn_fing,length") %>%
      req_error(body = function(resp) {
        paste("Error:", resp_status(resp), resp_status_desc(resp))
      }) %>%
      req_perform()
    
    # Parse Uniprot JSON response
    pfam_data <- resp %>%
      resp_body_string() %>%
      fromJSON()
    # Define pfam domain color palette
    color_palette <- c("#2dcf00", "#ff5353", "#5b5bff", "#ebd61d", "#ba21e0",
                       "#ff9c42", "#ff7dff", "#b9264f", "#baba21", "#c48484",
                       "#1f88a7", "#cafeb8", "#4a9586", "#ceb86c", "#0e180f")
    # Create a list to store used colors for each feature type
    used_colors <- list()
    
    # Function to process features of a each type 
    process_feature_type <- function(features, feature_type) {
      type_features <- features[features$type == feature_type, ]
      if (nrow(type_features) == 0) return(NULL)
      
      # Get available colors (excluding those used by other types)
      available_colors <- setdiff(color_palette, unlist(used_colors))
      
      # If we've used all colors, reset the available colors
      if (length(available_colors) == 0) {
        available_colors <- color_palette
      }
      
      # Assign colors, cycling through the available colors if needed
      colors <- rep_len(available_colors, nrow(type_features))
      
      # Update the used colors for this feature type
      used_colors[[feature_type]] <<- unique(colors)
      
      data.frame(
        type = type_features$type,
        description = type_features$description,
        start = type_features$location$start$value,
        end = type_features$location$end$value,
        color = colors,
        stringsAsFactors = FALSE
      )
    }
    
    # Get unique feature types
    feature_types <- unique(pfam_data$features$type)
    
    # Process each feature type and add to pfam_data
    for (feature_type in feature_types) {
      processed_df <- process_feature_type(pfam_data$features, feature_type)
      if (!is.null(processed_df)) {
        pfam_data[[feature_type]] <- processed_df
      }
    }
    #protein Length
    pfam_data$length <- pfam_data$sequence$length
    return(pfam_data)
  }, error = function(e) {
    stop(paste("Request failed:", e$message))
  })
}


# --- Alphafold ----#

extract_alphafold_plddt_stream <- function(x) {
  url <- sprintf("https://alphafold.ebi.ac.uk/files/AF-%s-F1-confidence_v4.json", x)
  message("Streaming ", url)
  print(url)
  con <- curl(url, handle = new_handle(timeout = 120))
  on.exit(close(con), add = TRUE)
  
  # stream_in will read & parse incrementally
  json_data <- stream_in(con, verbose = FALSE)
  
  if (nrow(json_data) == 0) {
    warning("No data in stream for ", url)
    return(NULL)
  }
  
  # same downstream processing as above…
  af <- data.frame(
    residueNumber     = json_data$residueNumber,
    confidenceScore   = json_data$confidenceScore,
    confidenceCategory= json_data$confidenceCategory,
    stringsAsFactors  = FALSE
  )
  af$confidenceScore[is.na(af$confidenceScore)] <- 0
  af$Confidence_Level <- factor(
    classify_confidence(af$confidenceScore),
    levels = names(confidence_colors)
  )
  af
}

# Example usage with Uniprot ID
#uniprot_id <- "Q9NZC2"
af <- extract_alphafold_plddt(x)

if (is.null(af) || any(is.na(af$confidenceScore))) {
  stop("Invalid or missing values in AlphaFold data.")
}

# Function to create the plot
plot_pLDDT <- function(af) {
  # Define confidence colors
  confidence_colors <- c(
    "D" = "orange",
    "L" = "yellow",
    "M" = "lightblue",
    "H" = "blue"
  )

  # Ensure confidenceCategory is a factor with predefined levels
  af$confidenceCategory <- factor(af$confidenceCategory, levels = names(confidence_colors))
  
  # Create the plot
  p <- ggplot(af, aes(
    x = residueNumber, y = confidenceScore,
    fill = confidenceCategory,  # Use confidenceCategory for coloring
    text = paste(
      "Residue:", residueNumber,
      "<br>Confidence Score:", confidenceScore,
      "<br>Confidence:", confidenceCategory
    )
  )) +
    geom_hline(yintercept = 10, colour = "#BB0000", linetype = "dashed") +
    geom_bar(stat = "identity", width = 1) +
    scale_fill_manual(values = confidence_colors) +
    labs(
      x = "",
      y = "AlphaFold Confidence (pLDDT)",  element_text(size = rel(0.45)) 
    ) +
    theme_minimal() +
    theme(legend.position = "none")
  p <- p + ggplot2::theme(
    axis.title.y = element_text(size = rel(0.5)),
    legend.position="none"
  )
  return(p)
}

# --- Alphafold ----#

# --- Alphafold Mean Pathogenicity Scores ----#

# Function to get mean pathogenicity data
get_mean_pathogenicity <- function(x) {
  url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", x, "-F1-aa-substitutions.csv")
  afs_data <- read_csv(url)
  
  # Compute mean pathogenicity per position
  mean_data <- afs_data %>%
    mutate(Position = as.numeric(substr(protein_variant, 2, nchar(protein_variant)-1))) %>%
    group_by(Position) %>%
    summarise(mean_pathogenicity = mean(am_pathogenicity, na.rm = TRUE),
              variants = paste(protein_variant[am_pathogenicity >= 0.564], collapse = ", ")) %>%
    mutate(
      category = case_when(
        mean_pathogenicity < 0.34 ~ "Benign",
        mean_pathogenicity >= 0.34 & mean_pathogenicity < 0.564 ~ "Uncertain",
        mean_pathogenicity >= 0.564 ~ "Pathogenic"
      )
    )
  
  return(mean_data)
}

# Load ClinVar variant data
inputclinvar_data <- fread("/home/raghu/Downloads/Simple-ClinVar/Path_clinvar.feb19.goldstars.pf", sep = "\t", header=TRUE, stringsAsFactors=FALSE)

# gene_data_pathogenic <- function(gene_symbol) {
#   gene_data <- inputclinvar_data %>%
#     filter(GeneSymbol == gene_symbol, ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic")) %>%
#     mutate(
#       ClinicalSignificance = factor(ClinicalSignificance, levels = c("Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic")),
#       y_position = case_when(
#         ClinicalSignificance == "Likely pathogenic" ~ 0.5,
#         ClinicalSignificance == "Pathogenic/Likely pathogenic" ~ 1,
#         ClinicalSignificance == "Pathogenic" ~ 1.5
#       )
#     ) %>%
#     filter(!is.na(pos_aa)) %>%
#     mutate(pos_aa = as.numeric(pos_aa))
#   
#   return(gene_data)
# }

# Modified function integrated into server.R using existing pfam_data
create_gene_plot_pathogenic <- function(gene_symbol, inputclinvar_data, pfam_data) {
  gene_data <- inputclinvar_data %>%
    filter(GeneSymbol == gene_symbol, ClinicalSignificance %in% c("Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic")) %>%
    mutate(
      ClinicalSignificance = factor(ClinicalSignificance, levels = c("Likely pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic")),
      y_position = case_when(
        ClinicalSignificance == "Likely pathogenic" ~ 0.5,
        ClinicalSignificance == "Pathogenic/Likely pathogenic" ~ 1,
        ClinicalSignificance == "Pathogenic" ~ 1.5
      )
    ) %>%
    filter(!is.na(pos_aa)) %>%
    mutate(pos_aa = as.numeric(pos_aa))
  
  if (nrow(gene_data) == 0) {
    stop(paste("No pathogenic data found for gene:", gene_symbol))
  }
  
  max_pos <- max(c(gene_data$pos_aa, pfam_data$end), na.rm = TRUE)
  
  p <- ggplot() +
    geom_rect(aes(xmin = 1, xmax = max_pos, ymin = -0.1, ymax = 0), 
              fill = "grey80", alpha = 0.8, color = "black", size = 0.3) +
    
    geom_rect(data = pfam_data,
              aes(xmin = start, xmax = end, ymin = -0.1, ymax = 0.1, fill = description),
              alpha = 1, color = "black", size = 0.3) +
    geom_text(data = pfam_data,
              aes(x = (start + end)/2, y = -0.05, label = description),
              vjust = 0.5, size = 3, color = "black") +
    
    geom_segment(data = gene_data,
                 aes(x = pos_aa, xend = pos_aa, y = 0, yend = y_position),
                 color = "gray70") +
    geom_point(data = gene_data,
               aes(x = pos_aa, y = y_position, color = ClinicalSignificance),
               size = 3) +
    
    scale_color_manual(values = c(
      "Likely pathogenic" = "red",
      "Pathogenic/Likely pathogenic" = "#CD5B45",
      "Pathogenic" = "darkred"
    )) +
    
    scale_fill_manual(values = setNames(pfam_data$color, pfam_data$description)) +
    
    scale_x_continuous(breaks = seq(0, max_pos, by = 50), minor_breaks = NULL) +
    scale_y_continuous(breaks = c(0, 0.5, 1, 1.5), 
                       labels = c("Domain/Region", "(Likely pathogenic) 0.5",
                                  "(Pathogenic/Likely pathogenic) 1", "(Pathogenic) 1.5"),
                       limits = c(-0.5, 2)) +
    labs(title = paste("Pathogenicity Along", gene_symbol, "Gene"),
         x = "Amino Acid Position", y = "",
         color = "Clinical Significance", fill = "Domain/Region") +
    theme_minimal() +
    theme(legend.position = "right",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.y = element_text(size = 10, angle = 0, hjust = 1))
  
  p = ggplotly(p, tooltip = c("pos_aa", "ClinicalSignificance"))
  return(p)
}



# Example UniProt ID
#uniprot_id <- "Q8NBP7"
mean_data <- get_mean_pathogenicity(x)

# Function to create the mean Path scores plot
plot_afmps <- function(mean_data) {
# Define color scale matching the image
color_scale <- c(
  "Benign" = "#1E90FF",      # Dodger blue
  "Uncertain" = "#FFD700",   # Gold
  "Pathogenic" = "#FF4500"   # OrangeRed
)

# Create bar chart with gradient coloring
p <- ggplot(mean_data, aes(x = Position, y = mean_pathogenicity, fill = category,
                                 text = paste("Position:", Position,
                                              "<br>AF Mean Pathogenicity:", round(mean_pathogenicity, 3),
                                              "<br>Category:", category,
                                              "<br>Pathogenic Variants:", variants))) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = color_scale) +
  geom_hline(yintercept = 0.34, linetype = "dashed", color = "gray50") +
  geom_hline(yintercept = 0.564, linetype = "dashed", color = "gray50") +
  annotate("text", x = max(mean_data$Position)*0.05, y = 0.17, label = "Benign", color = "gray30") +
  annotate("text", x = max(mean_data$Position)*0.05, y = 0.45, label = "Uncertain", color = "gray30") +
  annotate("text", x = max(mean_data$Position)*0.05, y = 0.78, label = "Pathogenic", color = "gray30") +
  labs(x = "", y = "AF Mean Pathogenicity Scores", element_text(size = rel(0.5)) )+
  theme_minimal() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.position = "none")

  p <- p + ggplot2::theme(
  axis.title.y = element_text(size = rel(0.5)),
  legend.position="none"
)

  return(p)
}
# --- Alphafold Mean Pathogenicity Scores ----#

#https://www.ebi.ac.uk/interpro/protein/UniProt/Q8NBP7/
#https://www.ebi.ac.uk/proteins/api/proteins?accession=Q8NBP7
#Extract Uniprot info from EBI data
extract_ebi <- function(x){
  Ebi_url <- paste("https://www.ebi.ac.uk/proteins/api/features?offset=0&size=100&accession=",x,sep = "")
  Ebi_data <- jsonlite::fromJSON(Ebi_url)
  EBi_features=as.data.frame(Ebi_data$features)
  Ebi_variants=EBi_features[EBi_features$type == "VARIANT",]
  Ebi_variants$begin=as.numeric(Ebi_variants$begin)
  Ebi_variants=Ebi_variants[order(as.numeric(Ebi_variants$begin)),]
  
  Ebi_domains=EBi_features[EBi_features$type == "DOMAIN",]
  Ebi_domains$begin=as.numeric(Ebi_domains$begin)
  Ebi_domains$end=as.numeric(Ebi_domains$end)
  Ebi_domains=Ebi_domains[order(as.numeric(Ebi_domains$begin)),]
  
  return(Ebi_variants,EBi_features, Ebi_domains)
}


#Extract Clinvar data
extract_clinvar <- function(gene_name){
  gene_clinvar_data=clinvar_data[clinvar_data$genename==gene_name,]
  gene_clinvar_data$clinvar_goldstar = as.character(gene_clinvar_data$clinvar_goldstar) 
  return(gene_clinvar_data)
}


#Extract Uniprot feature data
extract_uniprot_feature_data <- function(uniprotID){
  gene_uniprot_data= uniprot_feature_data[uniprot_feature_data$uniprot_id==uniprotID, ]
  gene_uniprot_data$start = as.numeric(gene_uniprot_data$start)
  gene_uniprot_data$end = as.numeric(gene_uniprot_data$end)
  return(gene_uniprot_data)
}


#Extract gnomad data
extract_gnomad <- function(chr,start,stop){
 query.regions <- paste(chr,":",start,"-",stop ,sep="")
 gnomad_vcf <- paste(file_path,"gnomad.exomes.ucon.r2.0.1.sites.GRCh38.noVEP.ANN.processed.xls.gz",sep = "/")
gene_gnomad_data <- tabix(query.regions, gnomad_vcf, check.chr = FALSE)
 if(is.null(gene_gnomad_data)) {
   gene_gnomad_data <- data.frame(matrix(ncol = 20, nrow = 0))
 }
 colnames(gene_gnomad_data) <- c("chrom","pos","ref","alt","hom","gene","transcript","cMut","pMutv","prot_pos","gnomad_allele1_count","gnomad_allele2_count","gnomad_allele3_count","gnomad_AN","gnomad_allele1_freq","gnomad_allele2_freq","gnomad_allele3_freq","gnomad_allele_count","gnomad_allele_freq","dominant")
 gene_gnomad_data$prot_pos <- as.numeric(gene_gnomad_data$prot_pos)
 gene_gnomad_data$gnomad_allele_count <- as.numeric(gene_gnomad_data$gnomad_allele_count)
 gene_gnomad_data$gnomad_allele_freq <- as.numeric(gene_gnomad_data$gnomad_allele_freq)
 return(gene_gnomad_data)
}


# Extract gnomAD data using embedded GraphQL query
# Function to query gnomAD API with gene symbol
query_gnomad_api <- function(gene_name) {
  gene_symbol=gene_name
  query <- "
  query GeneVariants($geneId: String!, $referenceGenome: ReferenceGenomeId!, $datasetId: DatasetId!) {
    gene(gene_symbol: $geneId, reference_genome: $referenceGenome) {
      variants(dataset: $datasetId) {
        variant_id
        hgvsp
        exome {
          af
          ac_hom
        }
      }
    }
  }"
  
  url <- "https://gnomad.broadinstitute.org/api"
  
  response <- POST(
    url = url,
    body = list(
      query = query,
      variables = list(
        geneId = gene_name,
        referenceGenome = "GRCh38",
        datasetId = "gnomad_r4"
      )
    ),
    encode = "json",
    content_type_json()
  )
  
  if (status_code(response) != 200) {
    stop("API request failed: ", content(response, "text"))
  }
  
  fromJSON(content(response, "text"), simplifyDataFrame = FALSE)
}

# Function to parse variant JSON data
parse_gnomad_variants <- function(json_parsed) {
  variants <- json_parsed$data$gene$variants
  
  three_to_one <- c(
    Ala = "A", Arg = "R", Asn = "N", Asp = "D", Cys = "C",
    Gln = "Q", Glu = "E", Gly = "G", His = "H", Ile = "I",
    Leu = "L", Lys = "K", Met = "M", Phe = "F", Pro = "P",
    Ser = "S", Thr = "T", Trp = "W", Tyr = "Y", Val = "V",
    Ter = "*", Sec = "U"
  )
  
  parse_hgvsp <- function(hgvsp_string) {
    if (is.null(hgvsp_string) || length(hgvsp_string) == 0 || is.na(hgvsp_string)) {
      return(c(NA, NA, NA))
    }
    
    hgvsp_string <- str_extract(hgvsp_string, "p\\.[^ ]+")
    if (is.na(hgvsp_string) || length(hgvsp_string) == 0) {
      return(c(NA, NA, NA))
    }
    
    if (str_detect(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|X|Ter).*")) {
      parts <- str_match(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|X|Ter)")
      return(c(
        three_to_one[parts[2]],
        parts[3],
        ifelse(parts[4] %in% c("X", "Ter"), "*", three_to_one[parts[4]])
      ))
    }
    
    if (str_detect(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)del")) {
      parts <- str_match(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)del")
      return(c(three_to_one[parts[2]], parts[3], "del"))
    }
    
    if (str_detect(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)[A-Za-z]{3}?fs")) {
      parts <- str_match(hgvsp_string, "p\\.([A-Za-z]{3})([0-9]+)[A-Za-z]{3}?fs")
      return(c(three_to_one[parts[2]], parts[3], "fs"))
    }
    
    if (str_detect(hgvsp_string, "p\\.Met1\\?")) return(c("M", 1, "?"))
    
    return(c(NA, NA, NA))
  }
  
  
  map_df(variants, function(v) {
    hgvsp_fields <- tryCatch({
      parsed <- parse_hgvsp(v$hgvsp)
      if (length(parsed) != 3) c(NA, NA, NA) else parsed
    }, error = function(e) c(NA, NA, NA))
    
    data.frame(
      aa_pos = as.integer(hgvsp_fields[2]),
      aa_ref = hgvsp_fields[1],
      aa_alt = hgvsp_fields[3],
      hgvsp = ifelse(is.null(v$hgvsp), NA, v$hgvsp),
      exome_af = ifelse(is.null(v$exome$af), 0, v$exome$af),
      exome_ac_hom = ifelse(is.null(v$exome$ac_hom), 0, v$exome$ac_hom),
      stringsAsFactors = FALSE
    )
  })
}  


extract_gnomad_coverage <- function(chr,start,stop){
  query.regions <- paste(chr,":",start,"-",stop ,sep="")
  gnomad_coverage_vcf <- paste(file_path,"gnomad.exomes.coverage.median_protcoord_ucon.tsv.gz",sep = "/")
  gene_gnomad_coverage_data <- tabix(query.regions, gnomad_coverage_vcf, check.chr = FALSE)
  if(is.null(gene_gnomad_coverage_data)) {
    gene_gnomad_coverage_data <- data.frame(matrix(ncol = 13, nrow = 0))
  }
  colnames(gene_gnomad_coverage_data) <- c("chr","pos","prot_pos","ensemble_transcript","median")
  gene_gnomad_coverage_data$prot_pos <- as.numeric(gene_gnomad_coverage_data$prot_pos)
  gene_gnomad_coverage_data$median <- as.numeric(gene_gnomad_coverage_data$median)
  return(gene_gnomad_coverage_data)
}

#Extract ccrs data
extract_ccrs <- function(chr,start,stop){
  query.regions <- paste(chr,":",start,"-",stop ,sep="")
  ccrs_vcf <- paste(file_path,"ccrs_b38_ucon.bed.gz",sep = "/")
  gene_ccrs_data <- tabix(query.regions, ccrs_vcf, check.chr = FALSE)
  if(is.null(gene_ccrs_data)) {
    gene_ccrs_data <- data.frame(matrix(ncol = 4, nrow = 0)) 
  }
  colnames(gene_ccrs_data) <- c("chr","pos","prot_pos","ensemble_transcript")
  gene_ccrs_data$prot_pos <- as.numeric(gene_ccrs_data$prot_pos)
  gene_ccrs_data = unique(gene_ccrs_data[,c("prot_pos","ensemble_transcript")])
  return(gene_ccrs_data)
}


#Extract DBNSFP data
extract_dbnsfp <- function(chr,start,stop){
  query.regions <- paste(chr,":",start,"-",stop ,sep="")
  dbnsfp_vcf <- paste(file_path,"dbNSFP4.0b1c_ucon.gnomad.exomes.coverage.median.xls.gz",sep = "/")
  gene_dbnsfp_data <- tabix(query.regions, dbnsfp_vcf, check.chr = FALSE)
  if(is.null(gene_dbnsfp_data)) {
    gene_dbnsfp_data <- data.frame(matrix(ncol = 69, nrow = 0)) 
  }
  colnames(gene_dbnsfp_data) <- c("hg19_chr","hg19_pos","hg38_chr","hg38_pos","ref","alt","aaref","aaalt","rs_dbSNP151","prot_pos","genename","Ensembl_geneid","Ensembl_transcriptid","Ensembl_proteinid","Uniprot_acc","APPRIS","GENCODE_basic","TSL","VEP_canonical","refcodon","codonpos","GERP.._RS","X1000Gp3_AC","X1000Gp3_AF","TWINSUK_AC","TWINSUK_AF","ALSPAC_AC","ALSPAC_AF","UK10K_AC","UK10K_AF","ESP6500_AA_AC","ESP6500_AA_AF","ESP6500_EA_AC","ESP6500_EA_AF","ExAC_AC","ExAC_AF","gnomAD_exomes_AC","gnomAD_exomes_AN","gnomAD_exomes_AF","gnomAD_exomes_nhomalt","gnomAD_exomes_controls_AC","gnomAD_exomes_controls_AN","gnomAD_exomes_controls_AF","gnomAD_exomes_controls_nhomalt","gnomAD_exomes_controls_SAS_nhomalt","gnomAD_exomes_controls_POPMAX_AC","gnomAD_exomes_controls_POPMAX_AN","gnomAD_exomes_controls_POPMAX_AF","gnomAD_exomes_controls_POPMAX_nhomalt","gnomAD_genomes_AC","gnomAD_genomes_AN","gnomAD_genomes_AF","gnomAD_genomes_nhomalt","gnomAD_genomes_POPMAX_AC","gnomAD_genomes_POPMAX_AN","gnomAD_genomes_POPMAX_AF","gnomAD_genomes_POPMAX_nhomalt","gnomAD_genomes_controls_AC","gnomAD_genomes_controls_AN","gnomAD_genomes_controls_AF","gnomAD_genomes_controls_nhomalt","gnomAD_genomes_controls_POPMAX_AC","gnomAD_genomes_controls_POPMAX_AN","gnomAD_genomes_controls_POPMAX_AF","gnomAD_genomes_controls_POPMAX_nhomalt","GTEx_V7_gene","GTEx_V7_tissue","Geuvadis_eQTL_target_gene","median")
  gene_dbnsfp_data$prot_pos <- as.numeric(gene_dbnsfp_data$prot_pos)
  return(gene_dbnsfp_data)
}

theme_Publication <- function(base_size=14, base_family="Times") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.position="none",
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

scale_fill_Publication <- function(...){
  library(scales)
  discrete_scale("fill","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}

scale_colour_Publication <- function(...){
  library(scales)
  discrete_scale("colour","Publication",manual_pal(values = c("#386cb0","#fdb462","#7fc97f","#ef3b2c","#662506","#a6cee3","#fb9a99","#984ea3","#ffff33")), ...)
  
}



# Updated function to use extract_gnomad_combined()
gnomad_depthplot <- function(gene_name, af_cutoff = 0.01, highlight_positions = NULL) {
  gnomad_data <- extract_gnomad_combined(gene_name, af_cutoff)

  ymaxlim <- af_cutoff + af_cutoff * 0.1

  gnomad_data$freq_status <- ifelse(gnomad_data$gnomad_allele_freq > af_cutoff, "AboveCutoff", "BelowCutoff")
  gnomad_data$plot_freq <- ifelse(gnomad_data$gnomad_allele_freq > af_cutoff, af_cutoff, gnomad_data$gnomad_allele_freq)
  gnomad_data$color <- ifelse(gnomad_data$dominant == "dominant", "orange", "#2CA25F")
  gnomad_data$shape <- ifelse(gnomad_data$freq_status == "AboveCutoff", 17, 16)

  prot_length <- max(gnomad_data$prot_pos, na.rm = TRUE)

  p <- ggplot() +
    geom_point(data = gnomad_data, aes(x = prot_pos, y = plot_freq, color = color, shape = factor(shape)), size = 2) +
    scale_shape_manual(values = c(`16` = 16, `17` = 17)) +
    scale_color_identity() +
    scale_y_reverse(limits = c(ymaxlim, -0.000005), expand = c(0, 0)) +
    xlim(-prot_length * 0.01, prot_length + prot_length * 0.01) +
    labs(y = "gnomAD freq", x = "Protein Position", title = paste("gnomAD Rainfall Plot:", gene_symbol)) +
    theme_bw(base_size = 10) +
    theme(
      axis.title.x = element_text(size = 10),
      axis.title.y = element_text(size = 10),
      legend.position = "none"
    )

  if (!is.null(highlight_positions) && length(highlight_positions) > 0) {
    highlight <- data.frame(prot_pos = highlight_positions)
    p <- p + geom_vline(data = highlight, aes(xintercept = prot_pos), linetype = "dotted", color = "blue")
  }

  return(p)
}



gnomad_freqplot <- function(gene_name, af_cutoff) {
  json_parsed <- query_gnomad_api(gene_name)
  gnomad_data <- parse_gnomad_variants(json_parsed)
  prot_length <- max(gnomad_data$aa_pos, na.rm = TRUE)
  ymaxlim <- af_cutoff + af_cutoff * 0.1
  
  gnomad_data <- gnomad_data %>%
    mutate(
      freq_status = ifelse(exome_af > af_cutoff, "AboveCutoff", "BelowCutoff"),
      plot_freq = ifelse(exome_af > af_cutoff, af_cutoff, exome_af),
      shape = ifelse(freq_status == "AboveCutoff", 17, 16),
      color = case_when(
        freq_status == "AboveCutoff" & exome_ac_hom > 0 ~ "orange",
        freq_status == "AboveCutoff" ~ "#2CA25F",
        exome_ac_hom > 0 ~ "orange",
        TRUE ~ "#66C2A4"
      )
    )
  
  p <- ggplot(gnomad_data, aes(x = aa_pos, y = plot_freq)) +
    geom_point(aes(shape = factor(shape), color = color), size = 2) +
    scale_shape_manual(values = c(`16` = 16, `17` = 17)) +
    scale_color_identity() +
    scale_y_reverse(limits = c(ymaxlim, -0.0000001), expand = c(0, 0)) +
    xlim(-prot_length * 0.01, prot_length + prot_length * 0.01) +
    labs(y = "gnomAD freq") +
    theme_bw(base_size = 10) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = rel(0.5)),
      plot.margin = unit(c(0, 0, 0, 0), "cm"),
      legend.position = "none"
    )
  
  return(p)
}


densityplot <- function(gnomad_data,clinvar_data,prot_length,allele_count,highlight,gene_ptm_data,gene_ccrs_data,user_path_variants) 
{
  kde.adjust<-0.05
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_density(data=gnomad_data[gnomad_data$gnomad_allele_count>as.numeric(allele_count),], aes(x=prot_pos),adjust=kde.adjust,fill='blue',alpha = 0.5, color="blue") 
  p <- p + ggplot2::xlim(-prot_length * 0.01, prot_length + prot_length * 0.01) 
  p <- p + ggplot2::labs(y = "Mutation Density") 
  p <- p + ggplot2::labs(x = "") 
  p <- p + ggplot2::theme(
    axis.title.y = element_text(size = rel(0.5)),
    legend.position="none"
  )
  if(nrow(clinvar_data) > 2 ){
    p = p + ggplot2::geom_density(aes(x=unique(clinvar_data$prot_pos)),adjust=kde.adjust, fill='red', alpha = 0.5,color="red") 
  }
  if(nrow(user_path_variants) > 2){
    p = p + ggplot2::geom_density(aes(x=unique(user_path_variants$V1)),adjust=kde.adjust, fill='green', alpha = 0.5,color="green") 
  }
  
  if(nrow(highlight) > 0 ){
    p = p + geom_vline(data=highlight,aes(xintercept=prot_pos), linetype="dotted")
  }
  return(p)
}

clinvar_ccrsplot <- function(pfam_data,uniprot_data,gene_clinvar_data,gene_ptm_data,gene_ccrs_data) 
{
  begin = end = NULL
  p <- ggplot2::ggplot()
  p <- p + ggplot2::ylim(0, 2.1) 
  p <- p + xlim(-pfam_data$sequence$length * 0.01, pfam_data$sequence$length + pfam_data$sequence$length * 0.01)
  p <- p + ggplot2::labs(y = "Pathogenic report") 
  p <- p + ggplot2::theme(
    axis.title.y = element_text(size = rel(0.5)),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position="none")
  p <- p + ggplot2::annotate("rect", xmin=0, xmax=pfam_data$sequence$length, ymin=1.7, ymax=1.85, fill="#C0C0C0", alpha=0.2)
  p <- p + annotate("text", label = "Clinvar Mis", x = 0, y = 2, size = 3, colour = "red")
  p <- p + ggplot2::annotate("rect", xmin=0, xmax=pfam_data$sequence$length, ymin=1.2, ymax=1.35, fill="#C0C0C0", alpha=0.2)
  p <- p + annotate("text", label = "Clinvar LOF", x = 0, y = 1.5, size = 3, colour = "red")
  p <- p + ggplot2::annotate("rect", xmin=0, xmax=pfam_data$sequence$length, ymin=0.7, ymax=0.85, fill="#C0C0C0", alpha=0.2)
  p <- p + annotate("text", label = "PTM", x = 0, y = 1, size = 3, colour = "red")
  p <- p + ggplot2::annotate("rect", xmin=0, xmax=pfam_data$sequence$length, ymin=0.1, ymax=0.25, fill="#C0C0C0", alpha=0.2)
  p <- p + annotate("text", label = "CCRS", x = 0, y = 0.5, size = 3, colour = "red")
  p <- p + ggplot2::scale_color_manual(values=c('0'='#808080','1'='#FF8C00','2'='#FF00FF','distal'='#999999','site'='#E69F00','proximal'='#56B4E9'))
  
  if( nrow(gene_clinvar_data[gene_clinvar_data$type=="missense_variant",]) > 0){
    p <- p + ggplot2::geom_segment(data = gene_clinvar_data[gene_clinvar_data$type=="missense_variant",],aes(x=prot_pos,xend=prot_pos,y=1.7,yend=1.85,color=clinvar_goldstar),size=0.5) 
  }
  if( nrow(gene_clinvar_data[gene_clinvar_data$type!="missense_variant",]) > 0){
    p <- p + ggplot2::geom_segment(data = gene_clinvar_data[gene_clinvar_data$type!="missense_variant",],aes(x=prot_pos,xend=prot_pos,y=1.2,yend=1.35,color=clinvar_goldstar),size=0.5) 
  }
  if(nrow(gene_ptm_data) > 0 ){
    p <- p + ggplot2::geom_segment(data=gene_ptm_data,aes(x=location,xend=location,y=0.7,yend=0.85,color=final_ptm_group),size=0.5) 
  }
  if(nrow(gene_ccrs_data) > 0 ){
    p <- p + ggplot2::geom_segment(data=gene_ccrs_data,aes(x=prot_pos,xend=prot_pos,y=0.1,yend=0.25),colour="#ADFF2F",size=0.5) 
  }
  return(p)
}


conservplot <- function(consurf_score,prot_length) 
{
  consurf_colors <- as.factor(c( '#10c8d1',  '#8cffff',  '#d7ffff',  '#eaffff',  '#ffffff',  '#fcedf4',  '#fac9de',  '#f07dab',  '#a02560'))
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_line(data=consurf_score,aes(x=POS,y=COLOR, colour=COLOR))
  p <- p + ggplot2::scale_colour_gradientn(na.value = '#ffff96', colours = consurf_colors ) 
  p <- p + ggplot2::ylim(0, 9)
  p <- p + ggplot2::xlim(-prot_length * 0.01, prot_length + prot_length * 0.01) 
  p <- p + ggplot2::labs(y = "Conservation Scores") 
  p <- p + ggplot2::theme(
    #axis.ticks.y = element_blank(), 
    #axis.text.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = rel(0.5)),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position="none")
  return(p)
}

pfamplot <- function(pfam_data,uniprot_data,gene_clinvar_data,highlight,label) 
{
  begin = end = NULL
  p <- ggplot2::ggplot()
  p <- p + ggplot2::ylim(0, 3) 
  p <- p + xlim(-pfam_data$sequence$length * 0.01, pfam_data$sequence$length + pfam_data$sequence$length * 0.01)
  p <- p + ggplot2::labs(y = "Uniprot Domains") 
  p <- p + ggplot2::labs(x = "Amino acid number") 
  p <- p + ggplot2::geom_segment(aes(x=0,xend=pfam_data$sequence$length,y=1.2,yend=1.2),size=2,color="grey") 
  if( length(pfam_data$Domain) > 0){
    p <- p + ggplot2::geom_rect(data = pfam_data$Domain, aes(xmin = start, xmax = end, ymin = 0.9, ymax = 1.5), fill = pfam_data$Domain$color, size = 2) 
    p <- p + ggplot2::geom_text(data=pfam_data$Domain , aes(x=end -((end - start)/2) , y = 1.2, label = description), size=2 ) 
    #    if(nrow(pfam_data$motifs[pfam_data$motifs$type=="low_complexity",])>0){
    #      p <- p + ggplot2::geom_rect(data = pfam_data$motifs[pfam_data$motifs$type=="low_complexity",], aes(xmin = metadata$start, xmax = metadata$end, ymin = 1.5, ymax = 1.6,alpha=0.2), color="black") 
    #      p <- p + ggplot2::geom_rect(data = pfam_data$motifs[pfam_data$motifs$type=="low_complexity",], aes(xmin = metadata$start, xmax = metadata$end, ymin = 0.8, ymax = 0.9,alpha=0.2), color="black") 
    #    }
  } 
  if( length(pfam_data$Region) > 0){
    p <- p + ggplot2::geom_rect(data = pfam_data$Region, aes(xmin = start, xmax = end, ymin = 0.9, ymax = 1.5,fill=description), fill = pfam_data$Region$color, color="black") 
    p <- p + ggplot2::geom_text(data=pfam_data$Region , aes(x=end -((end - start)/2) , y = 1.2, label = description), size=2 ) 
  }
  
  p <- p + ggplot2::annotate(geom="text", x=pfam_data$sequence$length/2, y=0.2, label=paste(pfam_data$uniProtkbId,"-",pfam_data$primaryAccession,"-",pfam_data$proteinDescription$recommendedName$fullName$value,"(", pfam_data$sequence$length,"aa)",sep=" "), color="black",size=3) 
  
  p <- p + ggplot2::theme( 
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.y = element_text(size = rel(0.5)), 
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    legend.position="none")
  
  if(nrow(highlight) > 0 ){
    p <- p + geom_segment(data=highlight,aes(x=prot_pos,y=1.5,xend=prot_pos,yend=1.9), linetype = "solid", color="black")
    p <- p + geom_point(data=highlight,aes(x=prot_pos,y=1.9), color="red",size=2)
    if (label == "Yes"){
      p <- p + geom_text(data=highlight, aes(x =prot_pos, y = 2.1, label = Mutation), size = 2, hjust = 0, nudge_y = 0.05, angle=45, family = 'Times',fontface = 'bold') 
    }
  }
  if(nrow(uniprot_data[uniprot_data$type=="signal",])>0 ){
    p <- p + ggplot2::geom_segment(data=uniprot_data[uniprot_data$type=="signal",],aes(x=start,xend=end,y=1.2,yend=1.2),size=2,color="red") 
  }
  if(nrow(uniprot_data[uniprot_data$type=="disulfid",])>0 ){
    values <- uniprot_data[uniprot_data$type=="disulfid" & !is.na(uniprot_data$end),]
    values <- values[order(values$start),]
    values$y_end <- as.numeric("2")
    k = 0
    for (i in 1:nrow(values)) {
      if( i %% 4 == 0 ) {
        k=0
      }
      values[i,]$y_end = as.numeric(2+k)
      k = k + 0.1
    }
    p <- p + geom_segment(data=values,aes(x =start, y = 1.5, xend = start, yend = y_end),color="blue") 
    p <- p + geom_segment(data=values,aes(x = end, y = 1.5, xend = end, yend = y_end),color="blue") 
    p <- p + geom_segment(data=values,aes(x = start, y = y_end, xend = end, yend = y_end),color="blue")
    
  }
  if(nrow(uniprot_data[uniprot_data$type=="mod_res",])>0 ){
    p <- p + geom_segment(data=uniprot_data[uniprot_data$type=="mod_res",],aes(x = start, y = 0.7, xend = start, yend =0.9),color="black") 
    p <- p + geom_point(data=uniprot_data[uniprot_data$type=="mod_res",],aes(x=start,y=0.7,fill=mod_res_group),shape=23,size =2)
  }
  if(nrow(uniprot_data[uniprot_data$type=="act_site",])>0 ){
    p <- p + geom_segment(data=uniprot_data[uniprot_data$type=="act_site",],aes(x = start, y = 0.7, xend = start, yend =0.9),color="black") 
    p <- p + geom_point(data=uniprot_data[uniprot_data$type=="act_site",],aes(x=start,y=0.7,fill="pink"),shape=22,size =2)
  }
  if(nrow(uniprot_data[uniprot_data$type=="mutagen",])>0 ){
    p <- p + geom_point(data=uniprot_data[uniprot_data$type=="mutagen",],aes(x=start,y=0.7),shape=10,size =2)
  }
  
  return(p)
}




# Define server logic for slider examples
shinyServer(function(input, output) {
  
  # first select input
  output$selects1 <- renderUI({
    selectInput("gene_name", "Choose a gene", gene_list, selected = "PCSK9")
    
  })
  
  shinyjs::onclick("toggleAdvanced", shinyjs::toggle(id = "advanced", anim = TRUE))    
  
  consurf_score <- eventReactive(input$gene_name,{ 
    if( !is.null(input$consurf_file) ) {
      consurf_score=fread(input$consurf_file$datapath, sep= "\t", header=T,quote="")
    }else{
      consurf_score=data.frame() 
    }
  })
  
  user_path_variants <-eventReactive(input$gene_name,{ 
    if( !is.null(input$user_file) ) {
      user_path_variants=fread(input$user_file$datapath, sep= "\t", header=F,quote="")
    }else{
      user_path_variants=data.frame() 
    }
  })
  
  # reset is a shinyjs function
  observeEvent(input$gene_name, { shinyjs::reset("consurf_file")  })
  observeEvent(input$gene_name, { shinyjs::reset("user_file")   })
  
  variants <- eventReactive(input$goButton,{ 
    variants=as.data.frame(strsplit(input$variants,"[,]" )) 
    colnames(variants) <- c("x")
    variants
  })
  highlight = eventReactive(input$goButton,{ 
    highlight=data.frame(variants(), sapply(variants(), extract_protein_position) ) 
    colnames(highlight) <- c("Mutation","prot_pos")
    highlight
  })
  
  data <- reactive({
    myPrev = 1/input$prev
    if(input$inh=="monoallelic"){
      myMaxAF = (1/2) * myPrev * input$hetA * input$hetG * (1/input$pen)
    }
    if(input$inh=="biallelic"){
      myMaxAF = sqrt(myPrev) * input$hetA * sqrt(input$hetG) * (1/sqrt(input$pen))
    }
    myMaxAC = qpois(p=as.numeric(input$CI),
                    lambda=(input$popSize)*(myMaxAF))
    return(list(myMaxAF,myMaxAC))
  }) 
  
  output$maxAC_F <- renderText({  paste("<b>Maximum credible population AF:<font color=\"#FF0000\"><b>", signif(data()[[1]],3), "</b></font><br>", "Maximum tolerated reference AC: <font color=\"#FF0000\"><b>", data()[[2]], "</b></font></b>") })
  
  gene_attrib <- eventReactive(input$goButton,{ gene_data[gene_data$gene_name==input$gene_name, ] })
  uniprotID <- eventReactive(input$goButton,{ as.character(gene_attrib()$uniprot_id)  })
  uniprot_data <- eventReactive(input$goButton,{ extract_uniprot_feature_data(uniprotID())  })
  
  pfam_data <- eventReactive(input$goButton,{ extract_pfam(uniprotID()) })
  ##af<- eventReactive(input$goButton,{ extract_alphafold_plddt(uniprotID()) })
  
  gene_clinvar_data <- eventReactive(input$goButton,{ extract_clinvar(input$gene_name) })
  gene_dbnsfp_data <- eventReactive(input$goButton,{ extract_dbnsfp(gene_attrib()$chr,gene_attrib()$hg19_start,gene_attrib()$hg19_stop) })
  gene_gnomad_data <- eventReactive(input$goButton,{ extract_gnomad(gene_attrib()$chr,gene_attrib()$hg38_start,gene_attrib()$hg38_stop) })
  gene_gnomad_coverage_data <- eventReactive(input$goButton,{ extract_gnomad_coverage(gene_attrib()$chr,gene_attrib()$hg38_start,gene_attrib()$hg38_stop) })
  gene_ccrs_data <- eventReactive(input$goButton,{ extract_ccrs(gene_attrib()$chr,gene_attrib()$hg38_start,gene_attrib()$hg38_stop) })

  output$geneinfo <- DT::renderDataTable( gene_attrib(), options = list(scrollX = "70%"))
  
  final_plot <- reactive({
    
    if (input$goButton == 0)
      return()
    
    gene_clinvar_data = gene_clinvar_data()
    if (input$clinvar_filter == "1 Star or more"){
      gene_clinvar_data_filtered=subset(gene_clinvar_data,clinvar_goldstar  > 0)
    }else if ( input$clinvar_filter == "2 Stars or more") {
      gene_clinvar_data_filtered=subset(gene_clinvar_data, clinvar_goldstar > 1)
    }else {
      gene_clinvar_data_filtered=gene_clinvar_data
    }
    
    gene_gnomad_data=gene_gnomad_data()
    gene_dbnsfp_data=gene_dbnsfp_data()
    pfam_data=pfam_data()
    gene_ptm_data=ptm_data[ptm_data$uniprot_id==gene_attrib()$uniprot_id,]
    
    
    p1 <- pfamplot(pfam_data(),uniprot_data(),gene_clinvar_data,highlight(),input$label)
    #p1 <- p1 + scale_colour_Publication()+ theme_Publication()
    p2 <-plot_afmps(mean_data) 
    #p2 <- clinvar_ccrsplot(pfam_data(),uniprot_data(),gene_clinvar_data,gene_ptm_data,gene_ccrs_data())
    #p2 <- p2 + scale_colour_Publication()+ theme_Publication()
    if( (nrow(gene_gnomad_data)>0) ) {
      p6 <- plot_pLDDT(af)
      #p6 <- gnomad_depthplot(gene_gnomad_coverage_data(),pfam_data$sequence$length)
      #p6 <- p6 + scale_colour_Publication()+ theme_Publication()
      #p5 <- gnomad_freqplot(gene_gnomad_data ,pfam_data$sequence$length,signif(data()[[1]],3))
      #p5 <-gnomad_freqplot (gene_name, max_af)
      p5 <- gnomad_freqplot(input$gene_name, data()[[1]])
      print(input$gene_name)
      #p5 <- p5 + scale_colour_Publication()+ theme_Publication()
      p4 <- densityplot(gene_gnomad_data,gene_clinvar_data_filtered[gene_clinvar_data_filtered$type=="missense_variant",],pfam_data$sequence$length,data()[[2]],highlight(),gene_ptm_data,gene_ccrs_data(),user_path_variants())
      #p4 <- p4 + scale_colour_Publication()+ theme_Publication()
    }
    if (!is.null(input$consurf_file) ) {
      p3 <- conservplot(consurf_score(),pfam_data$sequence$length)
      #p3 <- p3 + scale_colour_Publication()+ theme_Publication()
    }

    #if all boxes are checked
    if(all(plotsetA %in% input$plotselection)){
      #plotsetA <- c("clinvar","density","freq","depth")
            if (!is.null(input$consurf_file) ) {
              if( (nrow(gene_gnomad_data)>0) ) {
                 plotlist <- list(p6,p5,p4,p3,p2,p1)
              }else{
                 plotlist <- list(p4,p3,p2,p1)
              }
           }else{
              if( (nrow(gene_gnomad_data)>0) ) {
                 plotlist <- list(p6,p5,p4,p2,p1)
              }else{
                 plotlist <- list(p4,p2,p1)
              }
           }
     }else if (all(plotsetB %in% input$plotselection)){
            #plotsetB <- c("clinvar","freq","depth") (p2,p5,p6)
            if (!is.null(input$consurf_file) ) {
              if( (nrow(gene_gnomad_data)>0) ) {
                plotlist <- list(p6,p5,p3,p2,p1)
              }else{
                plotlist <- list(p3,p2,p1)
              }
            }else{
              if( (nrow(gene_gnomad_data)>0) ) {
                plotlist <- list(p6,p5,p2,p1)
              }else{
                plotlist <- list(p2,p1)
              }
            }
     }else if (all(plotsetC %in% input$plotselection)){
       #plotsetC <- c("Clinvar","density","depth") (p2,p4,p6)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p4,p3,p2,p1)
         }else{
           plotlist <- list(p4,p3,p2,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p4,p2,p1)
         }else{
           plotlist <- list(p4,p2,p1)
         }
       }
     }else if (all(plotsetD %in% input$plotselection)){
       #plotsetD <- c("clinvar","density","freq") (p2,p4,p5)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p5,p4,p3,p2,p1)
         }else{
           plotlist <- list(p4,p3,p2,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p5,p4,p2,p1)
         }else{
           plotlist <- list(p4,p2,p1)
         }
       }
     }else if (all(plotsetE %in% input$plotselection)){
       #plotsetE <- c("density","freq","depth") (p4,p5,p6)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p5,p4,p3,p1)
         }else{
           plotlist <- list(p4,p3,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p5,p4,p1)
         }else{
           plotlist <- list(p4,p1)
         }
       }
     }else if (all(plotsetF %in% input$plotselection)){
       #plotsetF <- c("clinvar","density") (p2,p4)
       if (!is.null(input$consurf_file) ) {
           plotlist <- list(p4,p3,p2,p1)
       }else{
           plotlist <- list(p4,p2,p1)
       }
     }else if (all(plotsetG %in% input$plotselection)){
       #plotsetG <- c("clinvar","freq") (p2,p5)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p5,p3,p2,p1)
         }else{
           plotlist <- list(p3,p2,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p5,p2,p1)
         }else{
           plotlist <- list(p2,p1)
         }
       }
     }else if (all(plotsetH %in% input$plotselection)){
       #plotsetH <- c("clinvar","depth") (p2,p6)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p3,p2,p1)
         }else{
           plotlist <- list(p3,p2,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p2,p1)
         }else{
           plotlist <- list(p2,p1)
         }
       }
     }else if (all(plotsetI %in% input$plotselection)){
       #plotsetI <- c("density","freq") (p4,p5)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p4,p3,p2,p1)
         }else{
           plotlist <- list(p4,p3,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p5,p4,p1)
         }else{
           plotlist <- list(p4,p1)
         }
       }
     }else if (all(plotsetJ %in% input$plotselection)){
       #plotsetJ <- c("density","depth") (p4,p6)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p4,p3,p1)
         }else{
           plotlist <- list(p4,p3,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p4,p1)
         }else{
           plotlist <- list(p4,p1)
         }
       }
     }else if (all(plotsetJ %in% input$plotselection)){
       #plotsetK <- c("freq","depth")(p5,p6)
       if (!is.null(input$consurf_file) ) {
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p5,p3,p1)
         }else{
           plotlist <- list(p3,p1)
         }
       }else{
         if( (nrow(gene_gnomad_data)>0) ) {
           plotlist <- list(p6,p5,p1)
         }else{
           plotlist <- list(p1)
         }
       }
     }else{
       #Final else part - Uniprot only   
       if (!is.null(input$consurf_file) ) {
         plotlist <- list(p3,p1)
       }else{
         plotlist <- list(p1)
       }
     }
    
    HPlot = subplot(plotlist, nrows = length(plotlist), shareX = TRUE,margin = c(-0.1,0,0,0))
    SPlot = plot_grid(plotlist = plotlist, ncol = 1, align = 'v')
    list(HTMLversion=HPlot,GGversion=SPlot)
    
  })
  
  output$mplot  <- renderPlotly({
    final_plot()$HTMLversion
  })
  
  highlight_data <- eventReactive(input$goButton,{ 
    if( (!is.null(gene_dbnsfp_data())) ) {
      df=merge(highlight(),gene_dbnsfp_data(),by="prot_pos") 
    }else{
      df=data.frame()
    }
  })
  
  output$highlight <- DT::renderDataTable( highlight_data(), options = list(scrollX = "70%"))
  output$clinvar <- DT::renderDataTable( gene_clinvar_data(), options = list(scrollX = "70%"))
  output$dbnsfp <- DT::renderDataTable( gene_dbnsfp_data(), options = list(scrollX = "70%"))
  output$gnomad <- DT::renderDataTable( gene_gnomad_data(), options = list(scrollX = "70%"))
  output$uniprot <- DT::renderDataTable( uniprot_data(), options = list(scrollX = "70%"))
  output$ccrs <- DT::renderDataTable( gene_ccrs_data(), options = list(scrollX = "70%"))
  
  # downloadHandler 
  
  output$down <- downloadHandler(
    filename =  function() {
      paste("MutEviz",input$gene_name, input$format, sep=".")
    },
    content = function(file) {
      if(input$format == "png")
        png(file) # open the png device
      else
        pdf(file) # open the pdf device
      plot(final_plot()$GGversion)
      dev.off()  # turn the device off
    } 
  )
  
  
})
