suppressMessages(library(dplyr))
suppressMessages(library(stringr))
suppressMessages(library(purrr))
suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinycssloaders))
suppressMessages(library(data.table))
suppressMessages(library(jsonlite))
suppressMessages(library(curl))
suppressMessages(library(grid))
suppressMessages(library(cowplot))
suppressMessages(library(ggplot2))
suppressMessages(library(stringi))
suppressMessages(library(plotly))
suppressMessages(library(httr))
suppressMessages(library(httr2))
suppressMessages(library(readr))
suppressMessages(library(wordcloud))
suppressMessages(library(tm))
suppressMessages(tryCatch(library(ggrepel), error=function(e) NULL))
# trackViewer/GenomicRanges/IRanges no longer needed — lollipops drawn natively in ggplot
theme_set(theme_cowplot(font_size=12))

# --- Unified VarViz plot theme for consistent axis fonts ---
# Typography tiers (vv_base_font / vv_big / vv_medium / vv_small) live in
# typography.R so they can be unit-tested without sourcing the full server.
source("typography.R", local = FALSE)

theme_vv <- function() {
  theme(
    text         = element_text(family = "Helvetica"),
    axis.title   = element_text(size = vv_medium, face = "bold", family = "Helvetica"),
    axis.text    = element_text(size = vv_small, family = "Helvetica"),
    plot.margin  = unit(c(0.15, 0.3, 0.15, 0.3), "cm"),
    legend.position = "none"
  )
}

# Helper: add variant dotted lines to a plot (filtered to x-range)
add_variant_vlines <- function(p, highlight, xmax = Inf, xmin = 0, ymax = NULL) {
  if (is.null(highlight) || nrow(highlight) == 0 || !"prot_pos" %in% colnames(highlight)) return(p)
  h <- highlight[highlight$prot_pos >= xmin & highlight$prot_pos <= xmax, , drop = FALSE]
  if (nrow(h) == 0) return(p)
  p <- p + geom_vline(data = h, aes(xintercept = prot_pos),
                       linetype = "dashed", color = "black", linewidth = 0.45, alpha = 0.7)
  p
}

current_directory <- getwd()

load("data/VarViz.RData")
gene_list <- sort(gene_data$gene_name)  # safety sort; run presort_gene_data.R once to pre-sort in RData

# GeVIR percentiles — merged into gene_data by add_gevir_to_rdata.R
# GeVIR_pct: low = missense intolerant (supports PP2), high = missense tolerant (supports BP1)
# Falls back to NA gracefully if column absent (old RData without GeVIR merge)
if (!"GeVIR_pct" %in% colnames(gene_data)) {
  gene_data$GeVIR_pct <- NA_real_
  message("[GeVIR] Column not found in gene_data — run add_gevir_to_rdata.R to enable BP1/PP2 GeVIR support")
}

# ---- API Response Cache ----
# Caches API responses in memory so re-querying the same gene is instant.
# Each cache is a named list keyed by gene_name or uniprotID.
# Cache persists for the lifetime of the R session (cleared on app restart).
api_cache <- new.env(hash = TRUE, parent = emptyenv())
api_cache$pfam        <- list()   # uniprotID -> pfam data
api_cache$uniprot_ft  <- list()   # uniprotID -> uniprot feature data
api_cache$gnomad_api  <- list()   # gene_name -> gnomAD API JSON
api_cache$gnomad_data <- list()   # gene_name -> parsed gnomAD data (for density)
api_cache$clinvar     <- list()   # gene_name -> clinvar data
api_cache$alphafold   <- list()   # uniprotID -> alphafold CSV data
api_cache$af_plddt    <- list()   # uniprotID -> alphafold pLDDT JSON
api_cache$gene_info   <- list()   # gene_name -> UniProt gene info
api_cache$ncbi_summary <- list()  # gene_name -> NCBI gene summary text
api_cache$ccrs         <- list()  # gene_name -> CCRS protein-level data
api_cache$clingen      <- list()  # gene_name -> ClinGen gene validity (classification, MOI, disease)
api_cache$dbnsfp       <- list()  # gene:variant -> dbNSFP MyVariant.info data
api_cache$consurf      <- list()  # uniprotID -> ConSurf-DB parsed data.frame
api_cache$ucsc_cons    <- list()  # gene_name -> data.frame of per-AA conservation scores (PhyloP/PhastCons)

cache_get <- function(cache_name, key) {
  if (key %in% names(api_cache[[cache_name]])) {
    message("[Cache HIT] ", cache_name, " : ", key)
    return(api_cache[[cache_name]][[key]])
  }
  return(NULL)
}

cache_set <- function(cache_name, key, value) {
  api_cache[[cache_name]][[key]] <- value
  invisible(value)
}
#curl.cainfo = "/data/cacert.pem"

#function to extract the protein position
extract_protein_position <- function(uniprotID){ as.numeric(lapply(regmatches(uniprotID , gregexpr("[[:digit:]]+", uniprotID)), `[`, 1)) }

# Function to fetch UniProt features for a given UniProt ID
extract_pfam <- function(uniprotID) {
  cached <- cache_get("pfam", uniprotID)
  if (!is.null(cached)) {
    message("[Pfam] Cache hit for: ", uniprotID)
    return(cached)
  }
  tryCatch({
    # Build the request
    resp <- request("https://rest.uniprot.org/uniprotkb/") %>%
      req_url_path_append(paste0(uniprotID, ".json")) %>%
      req_url_query(fields = "accession,id,protein_name,gene_primary,ft_repeat,ft_region,ft_domain,ft_compbias,ft_coiled,ft_motif,ft_zn_fing,ft_transmem,ft_topo_dom,ft_intramem,length,sequence") %>%
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
    cache_set("pfam", uniprotID, pfam_data)
    return(pfam_data)
  }, error = function(e) {
    message("UniProt API request failed: ", e$message)
    return(NULL) # Return NULL on error
  })
}


# --- Gene Info from UniProt (for GeneInfo tab) ----#
extract_gene_info_uniprot <- function(uniprotID, gene_name) {
  cached <- cache_get("gene_info", uniprotID)
  if (!is.null(cached)) {
    message("[GeneInfo] Cache hit for: ", uniprotID)
    return(cached)
  }
  tryCatch({
    message("[GeneInfo] Fetching UniProt data for: ", uniprotID)
    
    resp <- request("https://rest.uniprot.org/uniprotkb/") %>%
      req_url_path_append(paste0(uniprotID, ".json")) %>%
      req_timeout(30) %>%
      req_error(body = function(resp) {
        paste("Error:", resp_status(resp), resp_status_desc(resp))
      }) %>%
      req_perform()
    
    raw_json <- resp %>% resp_body_string()
    jdata <- fromJSON(raw_json, simplifyVector = FALSE)
    
    message("[GeneInfo] JSON parsed successfully. Top-level keys: ", 
            paste(names(jdata), collapse = ", "))
    
    # --- Gene Name ---
    primary_gene <- tryCatch({
      g <- jdata$genes[[1]]$geneName$value
      if (is.null(g)) gene_name else g
    }, error = function(e) gene_name)
    
    message("[GeneInfo] Gene: ", primary_gene)
    
    # --- Synonyms ---
    synonyms <- tryCatch({
      syns <- jdata$genes[[1]]$synonyms
      if (!is.null(syns) && length(syns) > 0) {
        paste(sapply(syns, function(s) s$value), collapse = ", ")
      } else {
        # Also check alternative names in protein description
        alt <- jdata$proteinDescription$alternativeNames
        if (!is.null(alt) && length(alt) > 0) {
          paste(sapply(alt, function(a) {
            if (!is.null(a$fullName$value)) a$fullName$value else ""
          }), collapse = ", ")
        } else "—"
      }
    }, error = function(e) "—")
    
    # --- Protein Name ---
    protein_name <- tryCatch({
      rec <- jdata$proteinDescription$recommendedName$fullName$value
      if (is.null(rec)) {
        # Fallback to submitted name
        sub <- jdata$proteinDescription$submissionNames[[1]]$fullName$value
        if (is.null(sub)) "Unknown protein" else sub
      } else rec
    }, error = function(e) "Unknown protein")
    
    # --- Protein Length ---
    prot_length <- tryCatch({
      l <- jdata$sequence$length
      if (is.null(l)) "—" else as.character(l)
    }, error = function(e) "—")
    
    # --- Function ---
    func_text <- tryCatch({
      comments <- jdata$comments
      if (is.null(comments) || length(comments) == 0) return("No function annotation available.")
      
      func_comments <- Filter(function(c) {
        !is.null(c$commentType) && c$commentType == "FUNCTION"
      }, comments)
      
      if (length(func_comments) > 0) {
        all_texts <- lapply(func_comments, function(fc) {
          if (!is.null(fc$texts) && length(fc$texts) > 0) {
            paste(sapply(fc$texts, function(t) {
              if (!is.null(t$value)) t$value else ""
            }), collapse = " ")
          } else ""
        })
        result <- paste(unlist(all_texts), collapse = " ")
        if (nchar(trimws(result)) == 0) "No function annotation available." else trimws(result)
      } else "No function annotation available."
    }, error = function(e) {
      message("[GeneInfo] Function extraction error: ", e$message)
      "No function annotation available."
    })
    
    message("[GeneInfo] Function text length: ", nchar(func_text))
    
    # --- Disease Involvement ---
    disease_info <- tryCatch({
      comments <- jdata$comments
      if (is.null(comments) || length(comments) == 0) {
        return(data.frame(Disease = character(), Description = character(),
                          MIM = character(), Notes = character(), stringsAsFactors = FALSE))
      }
      
      disease_comments <- Filter(function(c) {
        !is.null(c$commentType) && c$commentType == "DISEASE"
      }, comments)
      
      message("[GeneInfo] Found ", length(disease_comments), " disease comments")
      
      if (length(disease_comments) > 0) {
        disease_rows <- lapply(disease_comments, function(dc) {
          d_name <- tryCatch({
            v <- dc$disease$diseaseId
            if (is.null(v)) "—" else v
          }, error = function(e) "—")
          
          d_desc <- tryCatch({
            v <- dc$disease$description
            if (is.null(v)) "—" else v
          }, error = function(e) "—")
          
          d_mim <- tryCatch({
            ref <- dc$disease$diseaseCrossReference
            if (!is.null(ref) && !is.null(ref$database) && ref$database == "MIM") {
              ref$id
            } else "—"
          }, error = function(e) "—")
          
          d_note <- tryCatch({
            if (!is.null(dc$texts) && length(dc$texts) > 0) {
              paste(sapply(dc$texts, function(t) {
                if (!is.null(t$value)) t$value else ""
              }), collapse = " ")
            } else if (!is.null(dc$note) && !is.null(dc$note$texts) && length(dc$note$texts) > 0) {
              paste(sapply(dc$note$texts, function(t) {
                if (!is.null(t$value)) t$value else ""
              }), collapse = " ")
            } else "—"
          }, error = function(e) "—")
          
          data.frame(Disease = d_name, Description = d_desc,
                     MIM = d_mim, Notes = d_note, stringsAsFactors = FALSE)
        })
        do.call(rbind, disease_rows)
      } else {
        data.frame(Disease = character(), Description = character(),
                   MIM = character(), Notes = character(), stringsAsFactors = FALSE)
      }
    }, error = function(e) {
      message("[GeneInfo] Disease extraction error: ", e$message)
      data.frame(Disease = character(), Description = character(),
                 MIM = character(), Notes = character(), stringsAsFactors = FALSE)
    })
    
    # --- Cross-references (OMIM, Ensembl, HGNC) ---
    xrefs <- tryCatch({
      x <- jdata$uniProtKBCrossReferences
      if (is.null(x)) list() else x
    }, error = function(e) list())
    
    extract_xref_ids <- function(xrefs, db_name) {
      hits <- Filter(function(x) !is.null(x$database) && x$database == db_name, xrefs)
      if (length(hits) > 0) {
        paste(sapply(hits, function(h) h$id), collapse = ", ")
      } else NULL
    }
    
    omim_ids    <- extract_xref_ids(xrefs, "MIM")
    ensembl_ids <- extract_xref_ids(xrefs, "Ensembl")
    hgnc_ids    <- extract_xref_ids(xrefs, "HGNC")
    
    message("[GeneInfo] OMIM: ", omim_ids, " | Ensembl: ", 
            if(is.null(ensembl_ids)) "NULL" else substr(ensembl_ids, 1, 30),
            " | HGNC: ", hgnc_ids)
    
    # --- Build external links as list of list(name, url, label) ---
    links <- list()
    
    # DECIPHER
    links[[length(links)+1]] <- list(
      name = "DECIPHER",
      url  = paste0("https://www.deciphergenomics.org/gene/", primary_gene, "/overview"),
      label = primary_gene
    )
    
    # AlphaFold
    links[[length(links)+1]] <- list(
      name = "AlphaFold",
      url  = paste0("https://alphafold.ebi.ac.uk/search/text/", uniprotID),
      label = uniprotID
    )
    
    # OMIM
    if (!is.null(omim_ids)) {
      omim_first <- trimws(strsplit(omim_ids, ",")[[1]][1])
      links[[length(links)+1]] <- list(
        name = "OMIM",
        url  = paste0("https://www.omim.org/entry/", omim_first),
        label = omim_first
      )
    }
    
    # Ensembl
    if (!is.null(ensembl_ids)) {
      ens_first <- trimws(strsplit(ensembl_ids, ",")[[1]][1])
      links[[length(links)+1]] <- list(
        name = "Ensembl",
        url  = paste0("https://www.ensembl.org/Homo_sapiens/Transcript/Summary?t=", ens_first),
        label = ens_first
      )
    }
    
    # HGNC
    if (!is.null(hgnc_ids)) {
      hgnc_first <- trimws(strsplit(hgnc_ids, ",")[[1]][1])
      hgnc_num <- gsub("HGNC:", "", hgnc_first)
      links[[length(links)+1]] <- list(
        name = "HGNC",
        url  = paste0("https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/", hgnc_first),
        label = hgnc_first
      )
    }
    
    # UniProt
    links[[length(links)+1]] <- list(
      name = "UniProt",
      url  = paste0("https://www.uniprot.org/uniprotkb/", uniprotID, "/entry"),
      label = uniprotID
    )
    
    # PanelApp
    links[[length(links)+1]] <- list(
      name = "PanelApp",
      url  = paste0("https://panelapp.genomicsengland.co.uk/panels/entities/", primary_gene),
      label = primary_gene
    )
    
    # GenCC
    if (!is.null(hgnc_ids)) {
      hgnc_first <- trimws(strsplit(hgnc_ids, ",")[[1]][1])
      links[[length(links)+1]] <- list(
        name = "GenCC",
        url  = paste0("https://search.thegencc.org/genes/", hgnc_first),
        label = primary_gene
      )
    }
    
    # ClinGen
    if (!is.null(hgnc_ids)) {
      hgnc_first <- trimws(strsplit(hgnc_ids, ",")[[1]][1])
      links[[length(links)+1]] <- list(
        name = "ClinGen",
        url  = paste0("https://search.clinicalgenome.org/kb/genes/", hgnc_first),
        label = primary_gene
      )
    }
    
    # Medline
    links[[length(links)+1]] <- list(
      name = "Medline",
      url  = paste0("https://medlineplus.gov/genetics/gene/", tolower(primary_gene), "/"),
      label = primary_gene
    )
    
    # TRAP-Score
    links[[length(links)+1]] <- list(
      name = "TRAP-Score",
      url  = paste0("https://trap-score.org/Search?query=", primary_gene),
      label = primary_gene
    )
    
    # Gene2Phenotype
    links[[length(links)+1]] <- list(
      name = "Gene2Phenotype",
      url  = paste0("https://www.ebi.ac.uk/gene2phenotype/gene/", primary_gene),
      label = primary_gene
    )
    
    # GTEx
    links[[length(links)+1]] <- list(
      name = "GTEx",
      url  = paste0("https://gtexportal.org/home/gene/", primary_gene),
      label = primary_gene
    )
    
    # ConSurf-DB
    links[[length(links)+1]] <- list(
      name = "ConSurf-DB",
      url  = paste0("https://consurfdb.tau.ac.il/scripts/display_db.cgi?pdb_ID=", uniprotID, "&view_chain=A"),
      label = uniprotID
    )
    
    # Allele Frequency App (Ware et al.)
    links[[length(links)+1]] <- list(
      name = "AF Calculator",
      url  = "https://cardiodb.org/allelefrequencyapp/",
      label = "CardioDB"
    )
    
    result <- list(
      gene_name    = primary_gene,
      uniprot_id   = uniprotID,
      synonyms     = synonyms,
      protein_name = protein_name,
      length       = prot_length,
      function_text = func_text,
      diseases     = disease_info,
      links        = links,
      hgnc_id      = if (!is.null(hgnc_ids)) trimws(strsplit(hgnc_ids, ",")[[1]][1]) else NULL
    )
    
    message("[GeneInfo] Successfully extracted gene info for ", primary_gene)
    cache_set("gene_info", uniprotID, result)
    return(result)
    
  }, error = function(e) {
    message("[GeneInfo] FAILED for ", uniprotID, ": ", e$message)
    return(NULL)
  })
}

# --- Alphafold ----#
# --- Confidence Assignment ---
confidence_colors <- c(
  "Very low (pLDDT 0-50)" = "orange",
  "Low (pLDDT 50-70)" = "yellow",
  "Confident (pLDDT 70-90)" = "lightblue",
  "Very high (pLDDT 90-100)" = "blue"
)

# Function to classify confidence
classify_confidence <- function(pLDDT) {
  sapply(pLDDT, function(x) {
    if (x >= 90) {
      return("Very high (pLDDT 90-100)")
    } else if (x >= 70) {
      return("Confident (pLDDT 70-90)")
    } else if (x >= 50) {
      return("Low (pLDDT 50-70)")
    } else {
      return("Very low (pLDDT 0-50)")
    }
  })
}

extract_alphafold_plddt <- function(uniprotID) {
  cached <- cache_get("af_plddt", uniprotID)
  if (!is.null(cached)) { message("[AlphaFold pLDDT] Cache hit for: ", uniprotID); return(cached) }
  
  # Try multiple AlphaFold versions (v6 newest → v2 oldest)
  # Also try fragment suffixes (-F1, -F2) for multi-fragment predictions
  versions <- c("v6", "v4", "v3", "v2")
  json_raw <- NULL
  
  for (ver in versions) {
    url <- paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprotID, "-F1-confidence_", ver, ".json")
    message("[AlphaFold pLDDT] Trying: ", url)
    
    json_raw <- tryCatch({
      lines <- readLines(url, warn = FALSE)
      if (length(lines) > 0 && nchar(paste(lines, collapse = "")) > 10) lines else NULL
    }, error = function(e) NULL)
    
    if (!is.null(json_raw)) {
      message("[AlphaFold pLDDT] Found data with version ", ver)
      break
    }
  }
  
  # If all versions failed, return NULL
  if (is.null(json_raw) || length(json_raw) == 0) {
    message("[AlphaFold pLDDT] No data for ", uniprotID, " across versions ", paste(versions, collapse=", "))
    return(NULL)
  }
  
  # Try to parse JSON safely
  json_data <- tryCatch({
    fromJSON(paste(json_raw, collapse = ""))
  }, error = function(e) {
    message("[AlphaFold pLDDT] JSON parsing failed. Skipping.")
    return(NULL)
  })
  
  # Validate required fields
  if (is.null(json_data$residueNumber) ||
      is.null(json_data$confidenceScore) ||
      length(json_data$residueNumber) == 0) {
    
    message("[AlphaFold pLDDT] Missing or empty fields. Skipping.")
    return(NULL)
  }
  
  # Build dataframe
  af <- data.frame(
    residueNumber = json_data$residueNumber,
    confidenceScore = ifelse(is.na(json_data$confidenceScore), 0, json_data$confidenceScore)
  )
  
  # Confidence category
  af$confidenceCategory <- classify_confidence(af$confidenceScore)
  
  # Factor matching your colors
  af$Confidence_Level <- factor(
    af$confidenceCategory,
    levels = names(confidence_colors)
  )
  
  cache_set("af_plddt", uniprotID, af)
  return(af)
}


# Function to create the AlphaFold plot
plot_pLDDT <- function(af, highlight = data.frame(), prot_length = NULL) {
  if (is.null(af) || nrow(af) == 0) {
    p <- ggplot() + 
      ggplot2::annotate(geom = "text", x = 0.5, y = 0.5, 
                        label = "AlphaFold pLDDT data not available for this protein", 
                        size = 4.5, color = "#94a3b8", fontface = "italic") +
      labs(y = "AlphaFold pLDDT") +
      theme_bw(base_size = vv_medium) + theme_vv() +
      theme(axis.title.x = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
    return(p)
  }
  xmax <- if (!is.null(prot_length) && prot_length > 0) prot_length else max(af$residueNumber, na.rm = TRUE)
  
  # Define confidence colors
  confidence_colors <- c(
    "Very low (pLDDT 0-50)" = "orange",
    "Low (pLDDT 50-70)" = "yellow", 
    "Confident (pLDDT 70-90)" = "lightblue",
    "Very high (pLDDT 90-100)" = "blue"
  )
  
  # Create the plot
  p <- ggplot(af, aes(x = residueNumber, y = confidenceScore)) +
    geom_bar(aes(fill = Confidence_Level), stat = "identity", width = 1) +
    scale_fill_manual(values = confidence_colors) +
    geom_hline(yintercept = 50, colour = "#BB0000", linetype = "dashed", alpha = 0.7) +
    labs(x = "", y = "AlphaFold pLDDT") +
    scale_x_continuous(limits = c(-xmax * 0.01, xmax + xmax * 0.01), expand = c(0,0)) +
    theme_minimal(base_size = vv_medium) +
    theme_vv()
  
  p <- add_variant_vlines(p, highlight, xmax = xmax)
  
  return(p)
}

# Function to get mean pathogenicity data
get_mean_pathogenicity <- function(uniprotID) {
  cached <- cache_get("alphafold", uniprotID)
  if (!is.null(cached)) { message("[AlphaFold Path] Cache hit for: ", uniprotID); return(cached) }
  
  destination_path <- file.path(current_directory, paste0(uniprotID, "-F1-aa-substitutions.csv"))
  
  # Try multiple filename patterns (AlphaMissense data)
  urls <- c(
    paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprotID, "-F1-aa-substitutions.csv"),
    paste0("https://alphafold.ebi.ac.uk/files/AF-", uniprotID, "-F1-aa-substitutions-v1.csv")
  )
  
  ok <- FALSE
  for (url in urls) {
    message("[AlphaFold] Trying: ", url)
    ok <- tryCatch({
      download.file(url, destination_path, mode = "wb", quiet = TRUE)
      file.exists(destination_path) && file.size(destination_path) > 0
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok) {
      message("[AlphaFold] Found substitution data at: ", url)
      break
    }
  }
  
  # If all downloads failed or file is empty, skip
  if (!ok) {
    message("[AlphaFold] No substitution data found for ", uniprotID, ". Skipping.")
    return(NULL)
  }
  
  # Read CSV safely
  afs_data <- tryCatch({
    fread(destination_path)
  }, error = function(e) {
    message("[AlphaFold] Could not read CSV. Skipping.")
    return(NULL)
  })
  
  if (is.null(afs_data) || nrow(afs_data) == 0) {
    message("[AlphaFold] File contained no rows. Skipping.")
    return(NULL)
  }
  
  # Compute mean pathogenicity
  mean_data <- afs_data %>%
    mutate(Position = as.numeric(substr(protein_variant, 2, nchar(protein_variant) - 1))) %>%
    group_by(Position) %>%
    summarise(
      mean_pathogenicity = mean(am_pathogenicity, na.rm = TRUE),
      variants = paste(protein_variant[am_pathogenicity >= 0.564], collapse = ", ")
    ) %>%
    mutate(
      category = case_when(
        mean_pathogenicity < 0.34 ~ "Benign",
        mean_pathogenicity < 0.564 ~ "Uncertain",
        TRUE ~ "Pathogenic"
      )
    )
  
  cache_set("alphafold", uniprotID, mean_data)
  return(mean_data)
}


# Function to create the mean Path scores plot
plot_afmps <- function(mean_data, highlight = data.frame(), prot_length = NULL) {
  if (is.null(mean_data) || nrow(mean_data) == 0) {
    p <- ggplot() + 
      ggplot2::annotate(geom = "text", x = 0.5, y = 0.5, 
                        label = "AlphaFold Mean Pathogenicity data not available for this protein", 
                        size = 4.5, color = "#94a3b8", fontface = "italic") +
      labs(y = "AF Mean Pathogenicity") +
      theme_bw(base_size = vv_medium) + theme_vv() +
      theme(axis.title.x = element_blank(),
            axis.text = element_blank(), axis.ticks = element_blank(),
            panel.grid = element_blank())
    return(p)
  }
  xmax <- if (!is.null(prot_length) && prot_length > 0) prot_length else max(mean_data$Position, na.rm = TRUE)
  
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
    ggplot2::annotate(geom = "text", x = xmax * 0.05, y = 0.17, label = "Benign", color = "gray30", size = 3) +
    ggplot2::annotate(geom = "text", x = xmax * 0.05, y = 0.45, label = "Uncertain", color = "gray30", size = 3) +
    ggplot2::annotate(geom = "text", x = xmax * 0.05, y = 0.78, label = "Pathogenic", color = "gray30", size = 3) +
    labs(x = "", y = "AF Mean Pathogenicity") +
    scale_x_continuous(limits = c(-xmax * 0.01, xmax + xmax * 0.01), expand = c(0,0)) +
    theme_minimal(base_size = vv_medium) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank()) +
    theme_vv()

  p <- add_variant_vlines(p, highlight, xmax = xmax)

  return(p)
}

 #Extract Clinvar data
 # --- ClinVar data via NCBI E-utilities API ----#
 extract_clinvar <- function(gene_name){
   cached <- cache_get("clinvar", gene_name)
   if (!is.null(cached)) { message("[ClinVar] Cache hit for: ", gene_name); return(cached) }
   tryCatch({
     message("[ClinVar] Fetching variants for gene: ", gene_name)
     
     # Step 1: Search ClinVar for pathogenic/likely pathogenic variants in this gene
     # Try multiple search strategies since ClinVar search syntax varies
     ids <- NULL
     
     # Strategy 1: Standard clinical significance filter
     search_queries <- c(
       paste0(gene_name, '[gene] AND ("pathogenic"[clinsig] OR "likely pathogenic"[clinsig]) AND "homo sapiens"[orgn]'),
       paste0(gene_name, '[gene] AND (pathogenic[clinsig]) AND human[orgn]'),
       paste0(gene_name, '[gene] AND "clinsig pathogenic"[Properties]'),
       paste0(gene_name, '[gene]')
     )
     
     for (qi in seq_along(search_queries)) {
       search_term <- URLencode(search_queries[qi])
       search_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=clinvar&term=",
                            search_term, "&retmax=5000&retmode=json")
       
       search_resp <- httr::GET(search_url)
       search_json <- jsonlite::fromJSON(httr::content(search_resp, "text", encoding = "UTF-8"))
       
       ids <- search_json$esearchresult$idlist
       if (!is.null(ids) && length(ids) > 0) {
         message("[ClinVar] Query ", qi, " returned ", length(ids), " IDs: ", substr(search_queries[qi], 1, 80))
         break
       }
       message("[ClinVar] Query ", qi, " returned 0 results, trying next...")
       Sys.sleep(0.35)
     }
     
     if (is.null(ids) || length(ids) == 0) {
       message("[ClinVar] No variants found for ", gene_name)
       return(data.frame(chr = character(), pos = numeric(), prot_pos = numeric(),
                         name = character(), genename = character(), type = character(),
                         ClinicalSignificance = character(), clinvar_goldstar = character(),
                         stringsAsFactors = FALSE))
     }
     
     message("[ClinVar] Found ", length(ids), " variant IDs, fetching summaries...")
     
     # Step 2: Fetch summaries in batches of 200
     all_rows <- list()
     batch_size <- 200
     
     for (start in seq(1, length(ids), by = batch_size)) {
       batch_ids <- ids[start:min(start + batch_size - 1, length(ids))]
       id_str <- paste(batch_ids, collapse = ",")
       
       summary_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=clinvar&id=",
                             id_str, "&retmode=json")
       summary_resp <- httr::GET(summary_url)
       summary_json <- jsonlite::fromJSON(httr::content(summary_resp, "text", encoding = "UTF-8"),
                                          simplifyVector = FALSE)
       
       results <- summary_json$result
       if (is.null(results)) next
       
       uid_list <- results$uids
       if (is.null(uid_list)) next
       
       # Debug: save raw API response for first batch
       if (start == 1) {
         tryCatch({
           debug_file <- file.path(getwd(), "clinvar_debug_raw.json")
           writeLines(httr::content(summary_resp, "text", encoding = "UTF-8"), debug_file)
           message("[ClinVar DEBUG] Raw API response saved to: ", debug_file)
         }, error = function(e) message("[ClinVar DEBUG] Could not save debug file: ", e$message))
       }
       
       for (uid in uid_list) {
         rec <- results[[uid]]
         if (is.null(rec)) next
         
         # Debug: dump first record structure
         if (length(all_rows) == 0 && uid == uid_list[1]) {
           message("[ClinVar DEBUG] First record keys: ", paste(names(rec), collapse = ", "))
           message("[ClinVar DEBUG] title: ", rec$title)
           message("[ClinVar DEBUG] obj_type: ", rec$obj_type)
           message("[ClinVar DEBUG] protein_change class: ", class(rec$protein_change), 
                   " length: ", length(rec$protein_change),
                   " value: ", paste(unlist(rec$protein_change), collapse="|"))
           # Dump clinical_significance structure
           cs <- rec$clinical_significance
           if (!is.null(cs)) {
             if (is.list(cs)) {
               message("[ClinVar DEBUG] clinical_significance keys: ", paste(names(cs), collapse = ", "))
               message("[ClinVar DEBUG] clinical_significance$description: ", cs$description)
             } else {
               message("[ClinVar DEBUG] clinical_significance (non-list): ", cs)
             }
           } else {
             message("[ClinVar DEBUG] clinical_significance is NULL")
             # Try alternative field names
             gc_desc <- tryCatch(rec$germline_classification$description, error = function(e) "ERROR")
             message("[ClinVar DEBUG] germline_classification$description: ", gc_desc)
           }
         }
         
         # Extract clinical significance — handle multiple possible structures
         clin_sig <- tryCatch({
           # Try standard path first
           desc <- rec$clinical_significance$description
           if (!is.null(desc) && desc != "") desc else {
             # Try germline_classification (newer ClinVar API format)
             gc <- rec$germline_classification$description
             if (!is.null(gc) && gc != "") gc else {
               # Try top-level description
               d2 <- rec$description
               if (!is.null(d2) && d2 != "") d2 else "NA"
             }
           }
         }, error = function(e) "NA")
         
         # Track filtering stats
         if (!exists("debug_counts", inherits = FALSE)) {
           debug_counts <- list(total = 0, skip_not_pathogenic = 0, skip_conflicting = 0, 
                                skip_no_prot_pos = 0, added = 0)
         }
         debug_counts$total <- debug_counts$total + 1
         
         # Skip if not pathogenic/likely pathogenic
         if (!grepl("pathogenic", clin_sig, ignore.case = TRUE)) {
           debug_counts$skip_not_pathogenic <- debug_counts$skip_not_pathogenic + 1
           if (debug_counts$skip_not_pathogenic <= 3) {
             message("[ClinVar SKIP] uid=", uid, " not pathogenic: '", clin_sig, "'")
           }
           next
         }
         # Skip conflicting
         if (grepl("conflicting", clin_sig, ignore.case = TRUE)) {
           debug_counts$skip_conflicting <- debug_counts$skip_conflicting + 1
           next
         }
         
         # Extract review status — handle multiple possible structures
         review_status <- tryCatch({
           rs <- rec$clinical_significance$review_status
           if (!is.null(rs) && rs != "") rs else {
             rs2 <- rec$germline_classification$review_status
             if (!is.null(rs2) && rs2 != "") rs2 else "no assertion criteria provided"
           }
         }, error = function(e) "no assertion criteria provided")
         
         gold_star <- if (grepl("practice guideline", review_status, ignore.case = TRUE)) "4" else if (grepl("expert panel", review_status, ignore.case = TRUE)) "3" else if (grepl("multiple submitters", review_status, ignore.case = TRUE)) "2" else if (grepl("single submitter|criteria provided", review_status, ignore.case = TRUE)) "1" else "0"
         
         # Extract variant title/name
         var_title <- tryCatch({
           t <- rec$title
           if (is.null(t)) "NA" else t
         }, error = function(e) "NA")
         
         # Use protein_change field first (most reliable), fallback to title
         prot_change <- tryCatch({
           pc <- rec$protein_change
           if (is.null(pc) || (is.character(pc) && all(nchar(pc) == 0))) {
             # Try extracting from genes list
             genes_list <- rec$genes
             if (!is.null(genes_list) && length(genes_list) > 0) {
               for (g in genes_list) {
                 gpc <- g$protein_change
                 if (!is.null(gpc) && is.character(gpc) && any(nchar(gpc) > 0)) {
                   pc <- gpc
                   break
                 }
               }
             }
           }
           if (is.null(pc)) { NULL } else {
             # protein_change can be a list of strings (one per transcript)
             if (is.list(pc)) pc <- unlist(pc)
             if (is.character(pc)) pc <- pc[!is.na(pc) & nchar(pc) > 0]
             if (length(pc) == 0) NULL else pc[1]
           }
         }, error = function(e) NULL)
         
         prot_pos <- NA
         var_type <- "other"
         
         # Source for amino acid change: protein_change field, title, or variation_set_name
         aa_source <- NULL
         if (!is.null(prot_change) && nchar(prot_change) > 0) {
           aa_source <- prot_change
         } else if (grepl("\\(p\\.", var_title)) {
           aa_source <- sub(".*\\(p\\.([^)]+)\\).*", "\\1", var_title)
         } else if (grepl("p\\.", var_title)) {
           aa_source <- sub(".*p\\.([^ ,;]+).*", "\\1", var_title)
         } else {
           # Try variation_set_name
           vsn <- tryCatch(rec$variation_set_name, error = function(e) NULL)
           if (!is.null(vsn) && grepl("p\\.", vsn)) {
             aa_source <- sub(".*p\\.([^ ,;)]+).*", "\\1", vsn)
           }
         }
         
         if (!is.null(aa_source) && nchar(aa_source) > 0) {
           # Extract numeric position from the amino acid change
           pos_match <- regmatches(aa_source, regexpr("\\d+", aa_source))
           if (length(pos_match) > 0) {
             prot_pos <- as.numeric(pos_match)
           }
           
           # Classify variant type
           if (grepl("(del|ins|dup|fs|Ter|\\*|=)", aa_source, ignore.case = TRUE)) {
             var_type <- "lof_variant"
           } else {
             var_type <- "missense_variant"
           }
         }
         
         # Also check molecular_consequence_list for variant type
         mc <- tryCatch(rec$molecular_consequence_list, error = function(e) NULL)
         if (!is.null(mc) && length(mc) > 0) {
           mc_str <- paste(unlist(mc), collapse = " ")
           if (grepl("missense", mc_str, ignore.case = TRUE)) var_type <- "missense_variant" else if (grepl("frameshift|nonsense|splice|stop_gained", mc_str, ignore.case = TRUE)) var_type <- "lof_variant"
         }
         
         # Debug first few records
         if (length(all_rows) < 5) {
           message("[ClinVar DEBUG] uid=", uid, " title=", substr(var_title, 1, 60),
                   " clin_sig=", clin_sig, " prot_pos=", prot_pos, " type=", var_type,
                   " protein_change=", if(!is.null(prot_change)) prot_change else "NULL",
                   " aa_source=", if(!is.null(aa_source)) aa_source else "NULL")
         }
         
         # Extract chromosomal location
         chr_val <- "NA"
         pos_val <- NA
         
         variation_set <- tryCatch(rec$variation_set, error = function(e) NULL)
         if (!is.null(variation_set) && length(variation_set) > 0) {
           vs <- variation_set[[1]]
           if (!is.null(vs$variation_loc)) {
             locs <- vs$variation_loc
             for (loc in locs) {
               if (!is.null(loc$assembly_name) && grepl("GRCh37", loc$assembly_name)) {
                 chr_val <- if (!is.null(loc$chr)) loc$chr else "NA"
                 pos_val <- if (!is.null(loc$start)) as.numeric(loc$start) else NA
                 break
               }
             }
           }
         }
         
         # Extract VCV accession (e.g. VCV000017865)
         vcv_id <- tryCatch({
           v <- rec$accession
           if (!is.null(v) && nchar(v) > 0) v else {
             v2 <- rec$variation_set[[1]]$variation_name
             if (!is.null(v2)) paste0("uid:", uid) else paste0("uid:", uid)
           }
         }, error = function(e) paste0("uid:", uid))

         # Extract trait / disease name
         trait_name <- tryCatch({
           # Try trait_set → trait → name
           ts <- rec$trait_set
           if (!is.null(ts) && length(ts) > 0) {
             names_vec <- unlist(lapply(ts, function(t) {
               lapply(t$trait, function(tr) {
                 nl <- tr$name
                 if (!is.null(nl)) {
                   for (n in nl) {
                     if (!is.null(n$`element-value`) && nchar(n$`element-value`) > 0)
                       return(n$`element-value`)
                   }
                 }
                 NULL
               })
             }))
             names_vec <- names_vec[!sapply(names_vec, is.null)]
             if (length(names_vec) > 0) names_vec[[1]] else ""
           } else {
             # Fallback: use title before the gene symbol
             t2 <- sub("^(.*?)\\s*[A-Z0-9]+\\(.*$", "\\1", var_title)
             if (nchar(t2) > 0 && t2 != var_title) trimws(t2) else ""
           }
         }, error = function(e) "")

         if (!is.na(prot_pos)) {
           all_rows[[length(all_rows) + 1]] <- data.frame(
             chr = chr_val,
             pos = pos_val,
             prot_pos = prot_pos,
             name = var_title,
             genename = gene_name,
             type = var_type,
             ClinicalSignificance = clin_sig,
             clinvar_goldstar = gold_star,
             vcv_id = vcv_id,
             trait_name = trait_name,
             stringsAsFactors = FALSE
           )
           if (exists("debug_counts", inherits = FALSE)) debug_counts$added <- debug_counts$added + 1
         } else {
           if (exists("debug_counts", inherits = FALSE)) {
             debug_counts$skip_no_prot_pos <- debug_counts$skip_no_prot_pos + 1
             if (debug_counts$skip_no_prot_pos <= 5) {
               message("[ClinVar SKIP] uid=", uid, " no prot_pos. title=", substr(var_title, 1, 60),
                       " protein_change=", if(!is.null(prot_change)) prot_change else "NULL",
                       " aa_source=", if(!is.null(aa_source)) aa_source else "NULL")
             }
           }
         }
       }
       
       # Rate limit: NCBI allows 3 requests/sec without API key
       Sys.sleep(0.35)
     }
     
     # Debug: print filtering stats
     if (exists("debug_counts", inherits = FALSE)) {
       message("[ClinVar STATS] Total records: ", debug_counts$total,
               " | Skipped (not pathogenic): ", debug_counts$skip_not_pathogenic,
               " | Skipped (conflicting): ", debug_counts$skip_conflicting,
               " | Skipped (no prot_pos): ", debug_counts$skip_no_prot_pos,
               " | Added: ", debug_counts$added)
     }
     
     if (length(all_rows) == 0) {
       message("[ClinVar] No parsed variants for ", gene_name)
       return(data.frame(chr = character(), pos = numeric(), prot_pos = numeric(),
                         name = character(), genename = character(), type = character(),
                         ClinicalSignificance = character(), clinvar_goldstar = character(),
                         vcv_id = character(), trait_name = character(),
                         stringsAsFactors = FALSE))
     }
     
     result <- do.call(rbind, all_rows)
     result$prot_pos <- as.numeric(result$prot_pos)
     result$clinvar_goldstar <- as.character(result$clinvar_goldstar)
     
     message("[ClinVar] Extracted ", nrow(result), " variants (",
             sum(result$type == "missense_variant"), " missense, ",
             sum(result$type == "lof_variant"), " LoF)")
     cache_set("clinvar", gene_name, result)
     return(result)
     
   }, error = function(e) {
     message("[ClinVar] Error: ", e$message)
     data.frame(chr = character(), pos = numeric(), prot_pos = numeric(),
                name = character(), genename = character(), type = character(),
                ClinicalSignificance = character(), clinvar_goldstar = character(),
                vcv_id = character(), trait_name = character(),
                stringsAsFactors = FALSE)
   })
 }

 #Extract Uniprot feature data
 # Extract UniProt feature data via API (replaces local uniprot_feature_data)
 extract_uniprot_feature_data <- function(uniprotID) {
   cached <- cache_get("uniprot_ft", uniprotID)
   if (!is.null(cached)) { message("[UniProt Features] Cache hit for: ", uniprotID); return(cached) }
   tryCatch({
     message("[UniProt Features] Fetching for: ", uniprotID)
     resp <- request("https://rest.uniprot.org/uniprotkb/") %>%
       req_url_path_append(paste0(uniprotID, ".json")) %>%
       req_url_query(fields = "ft_signal,ft_disulfid,ft_mod_res,ft_act_site,ft_mutagen,ft_lipid,ft_carbohyd,ft_crosslnk,ft_transmem,ft_topo_dom,ft_intramem") %>%
       req_timeout(20) %>%
       req_error(body = function(resp) paste("Error:", resp_status(resp))) %>%
       req_perform()
     
     jdata <- fromJSON(resp %>% resp_body_string(), simplifyVector = FALSE)
     features <- jdata$features
     
     if (is.null(features) || length(features) == 0) {
       message("[UniProt Features] No features found")
       return(data.frame(uniprot_id = character(), type = character(),
                         start = numeric(), end = numeric(), description = character(), mod_res_group = character(),
                         stringsAsFactors = FALSE))
     }
     
     # Map UniProt feature type names to short codes used by pfamplot
     type_map <- c(
       "Signal peptide"    = "signal",
       "Disulfide bond"    = "disulfid",
       "Modified residue"  = "mod_res",
       "Active site"       = "act_site",
       "Mutagenesis"       = "mutagen",
       "Lipidation"        = "lipid",
       "Glycosylation"     = "carbohyd",
       "Cross-link"        = "crosslnk",
       "Transmembrane"     = "transmem",
       "Topological domain" = "topo_dom",
       "Intramembrane"     = "intramem"
     )
     
     rows <- lapply(features, function(f) {
       ftype <- type_map[f$type]
       if (is.na(ftype)) return(NULL)
       
       loc_start <- f$location$start$value
       loc_end   <- f$location$end$value
       if (is.null(loc_start)) return(NULL)
       
       # Classify mod_res_group from description
       desc <- if (!is.null(f$description)) f$description else ""
       mod_group <- if (ftype == "mod_res") {
         if (grepl("Phospho", desc, ignore.case = TRUE)) "Phosphorylation" else if (grepl("Acetyl", desc, ignore.case = TRUE)) "Acetylation" else if (grepl("Methyl", desc, ignore.case = TRUE)) "Methylation" else if (grepl("Ubiquitin", desc, ignore.case = TRUE)) "Ubiquitination" else if (grepl("Hydroxy", desc, ignore.case = TRUE)) "Hydroxylation" else if (grepl("Glycyl", desc, ignore.case = TRUE)) "Glycylation" else "Other modification"
       } else NA_character_
       
       data.frame(
         uniprot_id = uniprotID,
         type = ftype,
         start = as.numeric(loc_start),
         end = as.numeric(ifelse(is.null(loc_end), loc_start, loc_end)),
         description = if (!is.null(desc)) desc else "",
         mod_res_group = mod_group,
         stringsAsFactors = FALSE
       )
     })
     
     result <- do.call(rbind, Filter(Negate(is.null), rows))
     if (is.null(result) || nrow(result) == 0) {
       return(data.frame(uniprot_id = character(), type = character(),
                         start = numeric(), end = numeric(), description = character(), mod_res_group = character(),
                         stringsAsFactors = FALSE))
     }
     
     message("[UniProt Features] Found: ",
             sum(result$type == "signal"), " signal, ",
             sum(result$type == "disulfid"), " disulfid, ",
             sum(result$type == "mod_res"), " mod_res, ",
             sum(result$type == "act_site"), " act_site, ",
             sum(result$type == "mutagen"), " mutagen, ",
             sum(result$type == "transmem"), " transmem, ",
             sum(result$type == "topo_dom"), " topo_dom")
     cache_set("uniprot_ft", uniprotID, result)
     return(result)
     
   }, error = function(e) {
     message("[UniProt Features] Error: ", e$message)
     data.frame(uniprot_id = character(), type = character(),
                start = numeric(), end = numeric(), description = character(), mod_res_group = character(),
                stringsAsFactors = FALSE)
   })
 }


#Extract gnomad data
# --- gnomAD data via API (replaces local tabix file) ----#
extract_gnomad <- function(gene_name) {
  cached <- cache_get("gnomad_data", gene_name)
  if (!is.null(cached)) { message("[gnomAD Data] Cache hit for: ", gene_name); return(cached) }
  tryCatch({
    message("[gnomAD] Fetching variant data via API for: ", gene_name)
    json_parsed <- query_gnomad_api(gene_name)
    gnomad_data <- parse_gnomad_variants(json_parsed)
    
    if (is.null(gnomad_data) || nrow(gnomad_data) == 0) {
      message("[gnomAD] No variants returned from API")
      return(data.frame(prot_pos = numeric(), gnomad_allele_count = numeric(),
                        gnomad_allele_freq = numeric(),
                        af_afr = numeric(), af_nfe = numeric(),
                        af_eas = numeric(), af_sas = numeric(),
                        af_fin = numeric(), stringsAsFactors = FALSE))
    }
    
    # Map API columns to legacy column names used by densityplot
    # Also carry through ancestry AFs for the variant table population columns
    anc_cols <- intersect(c("af_afr","af_nfe","af_eas","af_sas","af_amr","af_fin"),
                          colnames(gnomad_data))
    result <- data.frame(
      prot_pos           = gnomad_data$aa_pos,
      aa_ref             = gnomad_data$aa_ref,
      aa_alt             = gnomad_data$aa_alt,
      gnomad_allele_freq = gnomad_data$exome_af,
      gnomad_allele_count = gnomad_data$exome_ac,
      gnomad_data[, anc_cols, drop = FALSE],
      stringsAsFactors = FALSE
    )
    result <- result[!is.na(result$prot_pos), ]
    
    message("[gnomAD] API returned ", nrow(result), " variants with protein positions")
    cache_set("gnomad_data", gene_name, result)
    return(result)
    
  }, error = function(e) {
    message("[gnomAD] API error: ", e$message)
    data.frame(prot_pos = numeric(), gnomad_allele_count = numeric(),
               gnomad_allele_freq = numeric(),
               af_afr = numeric(), af_nfe = numeric(),
               af_eas = numeric(), af_sas = numeric(),
               af_fin = numeric(), stringsAsFactors = FALSE)
  })
}


# Extract gnomAD data using embedded GraphQL query
# Function to query gnomAD API with gene symbol
query_gnomad_api <- function(gene_name) {
  cached <- cache_get("gnomad_api", gene_name)
  if (!is.null(cached)) { message("[gnomAD API] Cache hit for: ", gene_name); return(cached) }
  gene_symbol=gene_name
  # Query both exome and genome data
  query <- "
  query GeneVariants($geneId: String!, $referenceGenome: ReferenceGenomeId!, $datasetId: DatasetId!) {
    gene(gene_symbol: $geneId, reference_genome: $referenceGenome) {
      variants(dataset: $datasetId) {
        variant_id
        hgvsp
        exome {
          ac
          an
          af
          ac_hom
          populations { id ac an ac_hom }
        }
        genome {
          ac
          an
          af
          ac_hom
          populations { id ac an ac_hom }
        }
      }
    }
  }"
  
  url <- "https://gnomad.broadinstitute.org/api"
  
  # gnomAD v4 only (GRCh38, MANE Select transcript)
  # No v2.1.1 fallback — avoids cross-version transcript numbering discrepancies
  # where the same genomic variant has different residue numbers in old vs new transcripts
  for (attempt in 1:3) {
    message("[gnomAD API] Querying gnomad_r4 for ", gene_name,
            if (attempt > 1) paste0(" (retry ", attempt, ")") else "")
    
    response <- tryCatch({
      POST(
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
        content_type_json(),
        timeout(30)
      )
    }, error = function(e) {
      message("[gnomAD API] Request error: ", e$message)
      return(NULL)
    })
    
    if (is.null(response)) { Sys.sleep(2); next }
    
    resp_text <- httr::content(response, "text", encoding = "UTF-8")
    
    if (status_code(response) != 200) {
      message("[gnomAD API] HTTP ", status_code(response))
      Sys.sleep(2)
      next
    }
    
    json_parsed <- fromJSON(resp_text, simplifyDataFrame = FALSE)
    
    # Check for GraphQL errors — retry if transient
    if (!is.null(json_parsed$errors)) {
      err_msg <- substr(toJSON(json_parsed$errors, auto_unbox = TRUE), 1, 300)
      message("[gnomAD API] GraphQL errors: ", err_msg)
      if (grepl("overloaded|timeout|rate", err_msg, ignore.case = TRUE) && attempt < 3) {
        message("[gnomAD API] Transient error, retrying in 3s...")
        Sys.sleep(3)
        next
      }
      break
    }
    
    if (is.null(json_parsed$data$gene)) {
      message("[gnomAD API] No gene data in gnomad_r4 for ", gene_name)
      break
    }
    
    n_variants <- length(json_parsed$data$gene$variants)
    message("[gnomAD API] Got ", n_variants, " variants from gnomad_r4")
    
    if (n_variants > 0) {
      cache_set("gnomad_api", gene_name, json_parsed)
      return(json_parsed)
    }
    break
  }  # end retry loop
  
  message("[gnomAD API] No gnomAD v4 data for ", gene_name)
  return(list(data = list(gene = list(variants = list()))))
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
    
    # Compute per-population AF = ac/an (gnomAD v4 populations have ac+an but no af)
    pop_af <- function(pops, pop_id) {
      if (is.null(pops) || length(pops) == 0) return(NA_real_)
      for (p in pops) {
        if (identical(p$id, pop_id)) {
          ac <- if (!is.null(p$ac)) as.numeric(p$ac) else 0
          an <- if (!is.null(p$an)) as.numeric(p$an) else 0
          return(if (an > 0) ac / an else NA_real_)
        }
      }
      NA_real_
    }
    pops <- if (!is.null(v$exome$populations)) v$exome$populations else if (!is.null(v$genome$populations)) v$genome$populations else NULL

    data.frame(
      aa_pos = as.integer(hgvsp_fields[2]),
      aa_ref = hgvsp_fields[1],
      aa_alt = hgvsp_fields[3],
      hgvsp = ifelse(is.null(v$hgvsp), NA, v$hgvsp),
      exome_ac = if (!is.null(v$exome$ac)) v$exome$ac else if (!is.null(v$genome$ac)) v$genome$ac else 0,
      exome_an = if (!is.null(v$exome$an)) v$exome$an else if (!is.null(v$genome$an)) v$genome$an else 0,
      exome_af = if (!is.null(v$exome$af)) v$exome$af else if (!is.null(v$genome$af)) v$genome$af else 0,
      exome_ac_hom = if (!is.null(v$exome$ac_hom)) v$exome$ac_hom else if (!is.null(v$genome$ac_hom)) v$genome$ac_hom else 0,
      af_afr = pop_af(pops, "afr"),
      af_nfe = pop_af(pops, "nfe"),
      af_eas = pop_af(pops, "eas"),
      af_sas = pop_af(pops, "sas"),
      af_amr = pop_af(pops, "amr"),
      af_fin = pop_af(pops, "fin"),
      stringsAsFactors = FALSE
    )
  })
}  


# ── CCRS: Protein-level constraint from CCRStoAAC (Hasenahuer et al., 2022) ──
# Data source: marciaah/CCRStoAAC-output — CCR percentiles pre-mapped to UniProt
# amino acid positions using gnomAD 3.0 (GRCh38). Loaded as ccrs_data in VarViz.RData.
# No BED file, no bedr/tabix, no genomic-to-protein conversion needed.

extract_ccrs <- function(gene_symbol, uniprot_id = NULL) {
  cache_key <- gene_symbol
  cached <- cache_get("ccrs", cache_key)
  if (!is.null(cached)) {
    message("[CCRS] Cache hit for ", cache_key)
    return(cached)
  }
  
  empty_ccrs <- data.frame(prot_pos = integer(0), ccr_percentile = numeric(0),
                           stringsAsFactors = FALSE)
  
  # Check if ccrs_data is loaded from VarViz.RData
  if (!exists("ccrs_data", envir = .GlobalEnv) && !exists("ccrs_data")) {
    message("[CCRS] ccrs_data not found — was VarViz.RData built with preprocess_ccrs_to_rdata.R?")
    return(empty_ccrs)
  }
  
  tryCatch({
    # Lookup by gene name (primary)
    result <- ccrs_data[ccrs_data$ensembl_gene_name == gene_symbol, ]
    
    # Fallback: try UniProt ID if gene name returns nothing
    if (nrow(result) == 0 && !is.null(uniprot_id)) {
      result <- ccrs_data[ccrs_data$uniprotAcc == uniprot_id | 
                          ccrs_data$uniprotAcc == paste0(uniprot_id, "-1"), ]
    }
    
    if (nrow(result) == 0) {
      message("[CCRS] No CCRStoAAC data for gene: ", gene_symbol)
      return(empty_ccrs)
    }
    
    # Filter to constrained positions (percentile > 0, i.e. not unconstrained or unavailable)
    # Keep all positions ≥ 90th percentile for plotting (matches old threshold)
    constrained <- result[!is.na(result$aac_weighted_pct) & result$aac_weighted_pct >= 90, ]
    
    if (nrow(constrained) == 0) {
      message("[CCRS] Gene ", gene_symbol, " has ", nrow(result), 
              " mapped positions but none ≥ 90th percentile")
      # Return empty but cache the full result for variant table lookups
      cache_set("ccrs", paste0(cache_key, "_full"), result)
      return(empty_ccrs)
    }
    
    # Format for compatibility with existing plotting code (prot_pos column)
    gene_ccrs_data <- data.frame(
      prot_pos = constrained$aac_pos,
      ccr_percentile = constrained$aac_weighted_pct,
      stringsAsFactors = FALSE
    )
    gene_ccrs_data <- unique(gene_ccrs_data)
    gene_ccrs_data <- gene_ccrs_data[order(gene_ccrs_data$prot_pos), ]
    
    message("[CCRS] Found ", nrow(gene_ccrs_data), " constrained positions (≥90th pct) for ", gene_symbol)
    
    # Cache both the plot-ready data AND the full data for variant table
    cache_set("ccrs", cache_key, gene_ccrs_data)
    cache_set("ccrs", paste0(cache_key, "_full"), result)
    
    return(gene_ccrs_data)
  }, error = function(e) {
    message("[CCRS] Error: ", e$message)
    return(empty_ccrs)
  })
}

# Helper: Look up CCRS percentile at a specific amino acid position
# Used by build_variant_table for the CCRS column
lookup_ccrs_at_position <- function(gene_symbol, position) {
  # Try to get full (all positions) data from cache
  full_key <- paste0(gene_symbol, "_full")
  full_data <- cache_get("ccrs", full_key)
  
  if (is.null(full_data)) {
    # If full data not cached, try loading it
    if (exists("ccrs_data")) {
      full_data <- ccrs_data[ccrs_data$ensembl_gene_name == gene_symbol, ]
      if (nrow(full_data) > 0) {
        cache_set("ccrs", full_key, full_data)
      }
    }
  }
  
  if (is.null(full_data) || nrow(full_data) == 0) {
    return(list(in_ccrs = FALSE, percentile = NA, conservation = NA))
  }
  
  row <- full_data[full_data$aac_pos == position, ]
  if (nrow(row) == 0) {
    return(list(in_ccrs = FALSE, percentile = NA, conservation = NA))
  }
  
  row <- row[1, ]  # take first if duplicates
  pct <- row$aac_weighted_pct
  cons <- if ("CONSERVATION_score" %in% names(row)) row$CONSERVATION_score else NA
  
  # Negative percentiles mean "not available" (-1, -2, -3)
  if (is.na(pct) || pct < 0) {
    return(list(in_ccrs = FALSE, percentile = NA, conservation = cons))
  }
  
  return(list(
    in_ccrs = (pct >= 90),
    percentile = round(pct, 1),
    conservation = cons
  ))
}


# ══════════════════════════════════════════════════════════════════════════════
# dbNSFP via MyVariant.info API — pathogenicity scores, ensemble predictors,
# conservation, and population frequencies for each user-input variant
# ══════════════════════════════════════════════════════════════════════════════

# --- Amino acid 3-letter ↔ 1-letter conversion ---
aa3to1_map <- c(
  Ala="A", Arg="R", Asn="N", Asp="D", Cys="C", Gln="Q", Glu="E", Gly="G",
  His="H", Ile="I", Leu="L", Lys="K", Met="M", Phe="F", Pro="P", Ser="S",
  Thr="T", Trp="W", Tyr="Y", Val="V", Ter="*", Xaa="X"
)
aa1to3_map <- setNames(names(aa3to1_map), aa3to1_map)

# Convert p.Leu57Pro → p.L57P (or vice versa)
aa3to1 <- function(hgvsp) {
  # Already 1-letter? (e.g., p.L57P — letters around digits are single uppercase)
  if (grepl("^p\\.[A-Z]\\d+[A-Z*]$", hgvsp)) return(hgvsp)
  s <- hgvsp
  for (aa3 in names(aa3to1_map)) {
    s <- gsub(aa3, aa3to1_map[[aa3]], s)
  }
  s
}

aa1to3 <- function(hgvsp) {
  # Already 3-letter? (e.g., p.Leu57Pro)
  if (grepl("^p\\.[A-Z][a-z]{2}\\d+[A-Z][a-z]{2}$", hgvsp)) return(hgvsp)
  # Extract parts: p.X123Y
  m <- regmatches(hgvsp, regexec("^p\\.([A-Z*])(\\d+)([A-Z*])$", hgvsp))[[1]]
  if (length(m) < 4) return(hgvsp)
  ref_1 <- m[2]; pos <- m[3]; alt_1 <- m[4]
  ref_3 <- if (ref_1 %in% names(aa1to3_map)) aa1to3_map[[ref_1]] else ref_1
  alt_3 <- if (alt_1 %in% names(aa1to3_map)) aa1to3_map[[alt_1]] else alt_1
  paste0("p.", ref_3, pos, alt_3)
}

# Add cache slot for dbNSFP
api_cache$dbnsfp <- list()

# =============================================================================
# ClinGen Gene-Disease Validity
# =============================================================================
# Queries the ClinGen Linked Data Hub (LDH) for gene validity assertions.
# Returns a list with:
#   classification : "Definitive" | "Strong" | "Moderate" | "Limited" |
#                    "Disputed" | "Refuted" | "No Known Disease Relationship"
#   moi            : mode of inheritance string (or "")
#   disease        : curated disease name (or "")
#   source         : "ClinGen LDH"
#
# API endpoint: https://ldh.genome.network/ldh/Gene/id/{HGNC_symbol}/ld
# Returns JSON with ldFor/ld arrays; we extract GeneValidityAssertion entities.
# Requires no authentication. Rate limit is lenient for interactive use.
# =============================================================================

fetch_clingen_validity <- function(gene_symbol, hgnc_id = NULL) {
  # Use HGNC ID as cache key if available (more stable than symbol)
  cache_key <- if (!is.null(hgnc_id) && nchar(hgnc_id) > 0) hgnc_id else gene_symbol
  cached <- cache_get("clingen", cache_key)
  if (!is.null(cached)) return(cached)

  tier_order <- c("Definitive" = 6, "Strong" = 5, "Moderate" = 4,
                  "Limited" = 3, "Disputed" = 2, "Refuted" = 1,
                  "No Known Disease Relationship" = 0)

  result <- list(
    classification = NA_character_,  # best single tier
    moi            = "",
    disease        = "",             # best single disease name
    source         = "ClinGen LDH",
    all_assertions = NULL            # data.frame: disease, classification, moi, submitter
  )

  # ── Helper: parse tier string ─────────────────────────────────────────────
  match_tier <- function(cls_str) {
    cls_str <- as.character(cls_str)
    matched <- names(tier_order)[sapply(names(tier_order),
                 function(t) grepl(t, cls_str, ignore.case = TRUE))]
    if (length(matched) > 0) matched[1] else NA_character_
  }

  # ── Helper: try LDH URL and return parsed assertions or NULL ──────────────
  try_ldh_url <- function(url_suffix) {
    tryCatch({
      url  <- paste0("https://ldh.genome.network/ldh/Gene/id/",
                     URLencode(url_suffix, reserved = TRUE), "/ld")
      message("[ClinGen LDH] Trying: ", url)
      resp <- httr::GET(url, httr::timeout(12),
                        httr::add_headers(Accept = "application/json"))
      if (httr::status_code(resp) != 200) return(NULL)
      json    <- jsonlite::fromJSON(httr::content(resp, "text", encoding = "UTF-8"),
                                   simplifyVector = FALSE)
      ld_list <- json$ld
      if (is.null(ld_list) || length(ld_list) == 0) return(NULL)
      gva_items <- list()
      for (entry in ld_list) {
        gva <- entry$GeneValidityAssertion
        if (!is.null(gva) && length(gva) > 0) gva_items <- c(gva_items, gva)
      }
      if (length(gva_items) == 0) return(NULL)
      gva_items
    }, error = function(e) { message("[ClinGen LDH] Error: ", e$message); NULL })
  }

  # ── Source 1: ClinGen LDH ─────────────────────────────────────────────────
  # Try HGNC ID first (e.g. "HGNC:950"), then fall back to gene symbol
  # Some genes are indexed by HGNC ID in LDH, not by symbol
  ldh_ok <- tryCatch({
    gva_items <- if (!is.null(hgnc_id) && nchar(hgnc_id) > 0)
                   try_ldh_url(hgnc_id) else NULL
    if (is.null(gva_items))
      gva_items <- try_ldh_url(gene_symbol)
    if (is.null(gva_items)) { message("[ClinGen LDH] No data for: ", gene_symbol); FALSE
    } else {
      rows <- lapply(gva_items, function(item) {
            ec  <- item$entContent
            if (is.null(ec)) return(NULL)
            cls_raw <- if (!is.null(ec$classification)) as.character(ec$classification) else
                       if (!is.null(ec$classificationDescription)) as.character(ec$classificationDescription) else ""
            tier <- match_tier(cls_raw)
            if (is.na(tier)) return(NULL)
            dis <- if (!is.null(ec$disease))     as.character(ec$disease) else
                   if (!is.null(ec$diseaseName)) as.character(ec$diseaseName) else ""
            moi <- if (!is.null(ec$modeOfInheritance)) as.character(ec$modeOfInheritance) else ""
            aff <- if (!is.null(ec$affiliation))  as.character(ec$affiliation)  else
                   if (!is.null(ec$submitter))    as.character(ec$submitter)    else "ClinGen"
            data.frame(disease = dis, classification = tier, moi = moi,
                       submitter = aff, stringsAsFactors = FALSE)
          })
          rows <- Filter(Negate(is.null), rows)

          if (length(rows) > 0) {
            df <- do.call(rbind, rows)
            result$all_assertions <<- df
            df$tier_num <- tier_order[df$classification]
            best <- df[which.max(df$tier_num), ]
            result$classification <<- best$classification
            result$disease        <<- best$disease
            result$moi            <<- best$moi
            message("[ClinGen LDH] ", gene_symbol, " -> ", best$classification,
                    " (", nrow(df), " assertions)")
            TRUE
          } else FALSE
    }
  }, error = function(e) {
    message("[ClinGen LDH] Error for ", gene_symbol, ": ", e$message); FALSE
  })

  # ── Source 2: ClinGen Evidence Repo API (VCEP / affiliation-level fallback) ─
  # Hits the erepo REST API which backs clinicalgenome.org/affiliation/50013/
  # Returns per-disease assertions from all VCEPs for this gene
  if (!ldh_ok || is.null(result$all_assertions)) {
    tryCatch({
      # Use passed-in hgnc_id if available; otherwise look up from genenames.org
      erepo_hgnc <- if (!is.null(hgnc_id) && nchar(hgnc_id) > 0) {
        hgnc_id
      } else {
        hgnc_resp <- httr::GET(
          paste0("https://rest.genenames.org/fetch/symbol/", URLencode(gene_symbol, reserved = TRUE)),
          httr::add_headers(Accept = "application/json"), httr::timeout(8)
        )
        tryCatch({
          hj <- jsonlite::fromJSON(httr::content(hgnc_resp, "text", encoding = "UTF-8"),
                                   simplifyVector = FALSE)
          docs <- hj$response$docs
          if (!is.null(docs) && length(docs) > 0) docs[[1]]$hgnc_id else NULL
        }, error = function(e) NULL)
      }

      if (!is.null(erepo_hgnc)) {
        erepo_resp <- httr::GET(
          paste0("https://erepo.clinicalgenome.org/evrepo/api/classifications?hgnc_id=",
                 URLencode(erepo_hgnc, reserved = TRUE), "&limit=50"),
          httr::add_headers(Accept = "application/json"), httr::timeout(12)
        )
        if (httr::status_code(erepo_resp) == 200) {
          ej <- jsonlite::fromJSON(httr::content(erepo_resp, "text", encoding = "UTF-8"),
                                   simplifyVector = FALSE)
          items <- if (!is.null(ej$results)) ej$results else if (is.list(ej)) ej else list()
          rows <- lapply(items, function(it) {
            tier <- match_tier(if (!is.null(it$classification)) it$classification else "")
            if (is.na(tier)) return(NULL)
            dis <- if (!is.null(it$disease)) {
              if (is.list(it$disease) && !is.null(it$disease$label)) it$disease$label
              else as.character(it$disease)
            } else ""
            moi <- if (!is.null(it$modeOfInheritance)) {
              if (is.list(it$modeOfInheritance) && !is.null(it$modeOfInheritance$label))
                it$modeOfInheritance$label
              else as.character(it$modeOfInheritance)
            } else ""
            aff <- if (!is.null(it$affiliation)) {
              if (is.list(it$affiliation) && !is.null(it$affiliation$affiliation_fullname))
                it$affiliation$affiliation_fullname
              else as.character(it$affiliation)
            } else "ClinGen VCEP"
            data.frame(disease = dis, classification = tier, moi = moi,
                       submitter = aff, stringsAsFactors = FALSE)
          })
          rows <- Filter(Negate(is.null), rows)
          if (length(rows) > 0) {
            df <- do.call(rbind, rows)
            result$all_assertions <<- df
            result$source <<- "ClinGen eRepo"
            df$tier_num <- tier_order[df$classification]
            best <- df[which.max(df$tier_num), ]
            result$classification <<- best$classification
            result$disease        <<- best$disease
            result$moi            <<- best$moi
            message("[ClinGen eRepo] ", gene_symbol, " -> ", best$classification,
                    " (", nrow(df), " assertions)")
          }
        }
      }
    }, error = function(e) {
      message("[ClinGen eRepo] Error for ", gene_symbol, ": ", e$message)
    })
  }

  # ── Source 3: GenCC API ───────────────────────────────────────────────────
  # https://thegencc.org/api/v1/validity-prop?hgnc_id=HGNC:XXXXX
  # Supplements or replaces ClinGen when no assertions found above
  gencc_hgnc <- if (!is.null(hgnc_id) && nchar(hgnc_id) > 0) hgnc_id else NULL
  if (is.null(result$all_assertions) || is.na(result$classification)) {
    tryCatch({
      # Look up HGNC ID if not provided
      if (is.null(gencc_hgnc)) {
        hj <- tryCatch({
          hr <- httr::GET(
            paste0("https://rest.genenames.org/fetch/symbol/",
                   URLencode(gene_symbol, reserved = TRUE)),
            httr::add_headers(Accept = "application/json"), httr::timeout(8)
          )
          jsonlite::fromJSON(httr::content(hr, "text", encoding = "UTF-8"),
                             simplifyVector = FALSE)
        }, error = function(e) NULL)
        docs <- if (!is.null(hj)) hj$response$docs else NULL
        gencc_hgnc <- if (!is.null(docs) && length(docs) > 0) docs[[1]]$hgnc_id else NULL
      }

      if (!is.null(gencc_hgnc)) {
        hgnc_num <- sub("HGNC:", "", gencc_hgnc)
        gc_resp <- httr::GET(
          paste0("https://thegencc.org/api/v1/validity-prop?hgnc_id=HGNC:", hgnc_num),
          httr::add_headers(Accept = "application/json"), httr::timeout(12)
        )
        if (httr::status_code(gc_resp) == 200) {
          gj <- jsonlite::fromJSON(httr::content(gc_resp, "text", encoding = "UTF-8"),
                                   simplifyVector = FALSE)
          items <- if (!is.null(gj$data)) gj$data else list()
          rows <- lapply(items, function(it) {
            cls_raw <- if (!is.null(it$classification_title)) it$classification_title else
                       if (!is.null(it$classification))       as.character(it$classification) else ""
            # Map GenCC classifications to ClinGen tier names
            tier <- if (grepl("definitive",         cls_raw, ignore.case=TRUE)) "Definitive"   else
                    if (grepl("strong",              cls_raw, ignore.case=TRUE)) "Strong"       else
                    if (grepl("moderate",            cls_raw, ignore.case=TRUE)) "Moderate"     else
                    if (grepl("limited",             cls_raw, ignore.case=TRUE)) "Limited"      else
                    if (grepl("disputed",            cls_raw, ignore.case=TRUE)) "Disputed"     else
                    if (grepl("refuted",             cls_raw, ignore.case=TRUE)) "Refuted"      else
                    if (grepl("no known|animal|cel", cls_raw, ignore.case=TRUE)) "No Known Disease Relationship" else
                    NA_character_
            if (is.na(tier)) return(NULL)
            dis <- if (!is.null(it$disease_title)) it$disease_title else
                   if (!is.null(it$disease))       as.character(it$disease) else ""
            moi <- if (!is.null(it$moi_title)) it$moi_title else ""
            sub <- if (!is.null(it$submitter_title)) it$submitter_title else "GenCC"
            data.frame(disease=dis, classification=tier, moi=moi,
                       submitter=sub, stringsAsFactors=FALSE)
          })
          rows <- Filter(Negate(is.null), rows)
          if (length(rows) > 0) {
            df <- do.call(rbind, rows)
            result$all_assertions <<- df
            result$source <<- "GenCC"
            df$tier_num <- tier_order[df$classification]
            best <- df[which.max(df$tier_num), ]
            result$classification <<- best$classification
            result$disease        <<- best$disease
            result$moi            <<- best$moi
            message("[GenCC] ", gene_symbol, " -> ", best$classification,
                    " (", nrow(df), " submitters)")
          }
        }
      }
    }, error = function(e) {
      message("[GenCC] Error for ", gene_symbol, ": ", e$message)
    })
  }

  # Only cache if we got real data — don't cache failures/unknowns
  if (!is.na(result$classification)) {
    cache_set("clingen", cache_key, result)
    # Also cache under gene symbol as secondary key so symbol-based lookups hit cache too
    if (cache_key != gene_symbol) cache_set("clingen", gene_symbol, result)
  }
  result
}

# --- Fetch dbNSFP data from MyVariant.info ---
fetch_dbnsfp <- function(gene_name, hgvsp) {
  # Clean inputs: trim whitespace, remove leading "p." duplication
  gene_name <- trimws(gene_name)
  hgvsp <- trimws(hgvsp)
  if (nchar(gene_name) == 0 || nchar(hgvsp) == 0) return(NULL)
  
  cache_key <- paste0(gene_name, ":", hgvsp)
  cached <- cache_get("dbnsfp", cache_key)
  if (!is.null(cached)) {
    if (identical(cached, NA)) return(NULL)  # cached miss
    return(cached)
  }
  
  # Try both 3-letter and 1-letter forms
  hgvsp_1 <- aa3to1(hgvsp)
  hgvsp_3 <- aa1to3(hgvsp)
  variants_to_try <- unique(c(hgvsp, hgvsp_3, hgvsp_1))
  
  for (v in variants_to_try) {
    tryCatch({
      # Build the q= parameter as a single string, then encode the whole thing
      q_param <- paste0("dbnsfp.genename:", gene_name, " AND dbnsfp.hgvsp:", trimws(v))
      url <- paste0("https://myvariant.info/v1/query?",
                     "q=", URLencode(q_param, reserved = FALSE),
                     "&fields=dbnsfp,clinvar,dbsnp,cadd&size=5")
      message("[dbNSFP] Trying: ", gene_name, " ", v, " -> ", url)
      resp <- httr::GET(url, httr::timeout(15))
      if (httr::status_code(resp) != 200) {
        message("[dbNSFP] HTTP ", httr::status_code(resp), " for ", v)
        next
      }
      
      json_text <- httr::content(resp, "text", encoding = "UTF-8")
      parsed <- jsonlite::fromJSON(json_text, simplifyVector = FALSE)
      
      hits <- parsed$hits
      if (!is.null(hits) && length(hits) > 0) {
        result <- hits[[1]]
        message("[dbNSFP] HIT for ", gene_name, " ", v, " — keys: ", 
                paste(names(result$dbnsfp)[1:min(8, length(names(result$dbnsfp)))], collapse = ", "))
        cache_set("dbnsfp", cache_key, result)
        return(result)
      }
    }, error = function(e) {
      message("[dbNSFP] Error querying ", v, ": ", e$message)
    })
    Sys.sleep(0.2)  # rate limit courtesy
  }
  
  message("[dbNSFP] No hits for ", gene_name, " ", hgvsp)
  cache_set("dbnsfp", cache_key, NA)  # cache miss to avoid re-querying
  return(NULL)
}

# =============================================================================
# MULTI-CONSERVATION TRACK — UCSC REST API
# =============================================================================
# Fetches per-base PhyloP 100V, PhyloP 470M, PhastCons 100V
# from the UCSC Genome Browser REST API for the full CDS of a gene.
# Genomic coords are obtained from Ensembl REST API (/lookup/symbol).
# Scores are then mapped back to amino acid positions (codon midpoint).
#
# Track names (hg38):
#   PhyloP 100-way vertebrate : phyloP100wayAll
#   PhyloP 470-way mammalian  : phyloP470way
#   PhastCons 100-way         : phastCons100way
# =============================================================================

#' Get canonical transcript CDS exon structure from Ensembl REST API
#' Returns a data.frame: exon_rank, chr, genomic_start, genomic_end, strand
#' Coordinates are 1-based, inclusive (converted from Ensembl 0-based half-open)
fetch_ensembl_exons <- function(gene_symbol) {
  cache_key <- paste0("exons_", gene_symbol)
  cached <- cache_get("ucsc_cons", cache_key)
  if (!is.null(cached)) return(cached)

  tryCatch({
    # Step 1: Get canonical transcript ID
    url_gene <- paste0(
      "https://rest.ensembl.org/lookup/symbol/homo_sapiens/", gene_symbol,
      "?content-type=application/json&expand=1"
    )
    resp <- httr::GET(url_gene, httr::timeout(20),
                      httr::add_headers(Accept = "application/json"))
    if (httr::status_code(resp) != 200) {
      message("[UCSC Cons] Ensembl lookup failed for ", gene_symbol,
              " (HTTP ", httr::status_code(resp), ")")
      return(NULL)
    }
    gene_json <- jsonlite::fromJSON(
      httr::content(resp, "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )

    # Pick canonical transcript (is_canonical == 1) or longest CDS
    transcripts <- gene_json$Transcript
    if (is.null(transcripts) || length(transcripts) == 0) {
      message("[UCSC Cons] No transcripts found for ", gene_symbol)
      return(NULL)
    }

    canon_tx <- NULL
    for (tx in transcripts) {
      if (!is.null(tx$is_canonical) && tx$is_canonical == 1) {
        canon_tx <- tx
        break
      }
    }
    # Fall back to longest CDS if no canonical flag
    if (is.null(canon_tx)) {
      cds_lengths <- sapply(transcripts, function(tx) {
        if (!is.null(tx$Translation$length)) tx$Translation$length else 0L
      })
      canon_tx <- transcripts[[which.max(cds_lengths)]]
    }

    tx_id  <- canon_tx$id
    strand <- if (!is.null(canon_tx$strand)) canon_tx$strand else 1
    chr    <- if (!is.null(canon_tx$seq_region_name)) canon_tx$seq_region_name else gene_json$seq_region_name

    message("[UCSC Cons] Canonical transcript: ", tx_id, " | chr", chr,
            " | strand: ", strand)

    # Step 2: Get exon structure with CDS phase info
    url_tx <- paste0(
      "https://rest.ensembl.org/lookup/id/", tx_id,
      "?content-type=application/json&expand=1"
    )
    resp2 <- httr::GET(url_tx, httr::timeout(20),
                       httr::add_headers(Accept = "application/json"))
    if (httr::status_code(resp2) != 200) {
      message("[UCSC Cons] Ensembl transcript lookup failed for ", tx_id)
      return(NULL)
    }
    tx_json <- jsonlite::fromJSON(
      httr::content(resp2, "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )

    exons <- tx_json$Exon
    if (is.null(exons) || length(exons) == 0) {
      message("[UCSC Cons] No exons found for ", tx_id)
      return(NULL)
    }

    # Build exon data.frame — convert Ensembl 0-based half-open to 1-based inclusive
    exon_df <- do.call(rbind, lapply(seq_along(exons), function(i) {
      ex <- exons[[i]]
      data.frame(
        exon_rank    = i,
        chr          = as.character(chr),
        genomic_start = as.integer(ex$start),   # Ensembl already 1-based
        genomic_end   = as.integer(ex$end),
        strand        = as.integer(strand),
        stringsAsFactors = FALSE
      )
    }))

    # Sort exons by genomic position on correct strand
    exon_df <- exon_df[order(exon_df$genomic_start), , drop = FALSE]

    # Get CDS start/end from Translation object to trim UTRs
    if (!is.null(tx_json$Translation)) {
      cds_start_genomic <- as.integer(tx_json$Translation$start)
      cds_end_genomic   <- as.integer(tx_json$Translation$end)
      exon_df <- exon_df[
        exon_df$genomic_end   >= cds_start_genomic &
        exon_df$genomic_start <= cds_end_genomic, , drop = FALSE
      ]
      # Trim partial first/last exon to CDS boundary
      exon_df$genomic_start[1] <- max(exon_df$genomic_start[1], cds_start_genomic)
      exon_df$genomic_end[nrow(exon_df)] <- min(
        exon_df$genomic_end[nrow(exon_df)], cds_end_genomic
      )
    }

    # Reverse for minus-strand genes so exons are CDS-order (AA position 1 first)
    if (strand == -1) exon_df <- exon_df[nrow(exon_df):1, , drop = FALSE]
    exon_df$exon_rank <- seq_len(nrow(exon_df))

    message("[UCSC Cons] Exons (CDS trimmed): ", nrow(exon_df),
            " | Total CDS nt: ",
            sum(exon_df$genomic_end - exon_df$genomic_start + 1))

    cache_set("ucsc_cons", cache_key, exon_df)
    return(exon_df)

  }, error = function(e) {
    message("[UCSC Cons] Ensembl exon fetch error: ", e$message)
    return(NULL)
  })
}


#' Convert amino acid position to genomic coordinate (hg38, 1-based)
#' Returns the genomic position of the FIRST base of the codon for aa_pos
#' @param aa_pos  Integer amino acid position (1-based)
#' @param exon_df data.frame from fetch_ensembl_exons()
#' @return integer genomic position, or NA
aa_to_genomic <- function(aa_pos, exon_df) {
  if (is.null(exon_df) || nrow(exon_df) == 0) return(NA_integer_)
  nt_pos <- (as.integer(aa_pos) - 1L) * 3L + 1L  # 1st nt of codon

  strand    <- exon_df$strand[1]
  cumulative <- 0L

  if (strand == 1) {
    for (i in seq_len(nrow(exon_df))) {
      exon_len <- exon_df$genomic_end[i] - exon_df$genomic_start[i] + 1L
      if (cumulative + exon_len >= nt_pos) {
        offset <- nt_pos - cumulative - 1L
        return(as.integer(exon_df$genomic_start[i] + offset))
      }
      cumulative <- cumulative + exon_len
    }
  } else {
    # Minus strand: exon_df already reversed so exon 1 is highest genomic coords
    for (i in seq_len(nrow(exon_df))) {
      exon_len <- exon_df$genomic_end[i] - exon_df$genomic_start[i] + 1L
      if (cumulative + exon_len >= nt_pos) {
        offset <- nt_pos - cumulative - 1L
        return(as.integer(exon_df$genomic_end[i] - offset))
      }
      cumulative <- cumulative + exon_len
    }
  }
  return(NA_integer_)
}


#' Fetch conservation scores from UCSC REST API for a genomic window
#' Returns a list: phylop100, phylop470, phastcons
#' Each element is a named numeric vector: names = genomic positions (as character)

fetch_ucsc_window <- function(chr, start, end) {
  # Ensure chr has "chr" prefix
  chr_ucsc <- if (grepl("^chr", chr)) chr else paste0("chr", chr)

  base_url <- "https://api.genome.ucsc.edu/getData/track"

  # Helper: try a track name, return named numeric vector pos->value or NULL
  fetch_track <- function(track_name) {
    url <- paste0(base_url,
                  "?genome=hg38;track=", track_name,
                  ";chrom=", chr_ucsc,
                  ";start=", start - 1L,   # UCSC REST uses 0-based half-open
                  ";end=",   end)
    tryCatch({
      resp <- httr::GET(url, httr::timeout(30))
      if (httr::status_code(resp) != 200) {
        message("[UCSC Cons] HTTP ", httr::status_code(resp), " for track='", track_name, "'")
        return(NULL)
      }
      raw_text <- httr::content(resp, "text", encoding = "UTF-8")
      parsed <- jsonlite::fromJSON(raw_text, simplifyVector = TRUE)
      top_keys <- names(parsed)

      # UCSC REST API wiggle/bigWig response structure (confirmed from phastCons log):
      # parsed$<chromName>  e.g. parsed$chr14  — this IS the data container
      # It can be:
      #   A) data.frame with cols: start, end, value  (fixedStep / variableStep wig)
      #   B) list with $start (int vec) and $value (num vec)
      #   C) data.frame with chromStart, dataValue (BED-like)
      # Other top-level keys: downloadTime, genome, track, start, end, chrom, itemsReturned

      extract_pv <- function(obj) {
        if (is.null(obj)) return(NULL)
        # data.frame with start + value
        if (is.data.frame(obj)) {
          if ("start" %in% names(obj) && "value" %in% names(obj)) {
            pos <- as.integer(obj$start) + 1L  # 0-based -> 1-based
            val <- as.numeric(obj$value)
            return(list(pos = pos, val = val))
          }
          if ("chromStart" %in% names(obj)) {
            vc <- intersect(c("value","score","dataValue"), names(obj))[1]
            if (!is.na(vc)) {
              pos <- as.integer(obj$chromStart) + 1L
              val <- as.numeric(obj[[vc]])
              return(list(pos = pos, val = val))
            }
          }
        }
        # list with $start and $value vectors
        if (is.list(obj) && !is.data.frame(obj) &&
            !is.null(obj$start) && !is.null(obj$value)) {
          pos <- as.integer(unlist(obj$start)) + 1L
          val <- as.numeric(unlist(obj$value))
          return(list(pos = pos, val = val))
        }
        return(NULL)
      }

      # Try: chrom key (e.g. "chr14"), then track_name key, then "data", then top-level
      candidate_keys <- c(chr_ucsc, track_name, "data")
      for (key in candidate_keys) {
        obj <- parsed[[key]]
        pv  <- extract_pv(obj)
        if (!is.null(pv) && length(pv$pos) > 0) {
          message("[UCSC Cons] ", track_name, ": ", length(pv$pos), " positions via $",
                  key, " | range ", min(pv$pos), "-", max(pv$pos))
          return(setNames(pv$val, as.character(pv$pos)))
        }
      }
      # top-level start/value fallback
      pv <- extract_pv(parsed)
      if (!is.null(pv) && length(pv$pos) > 0) {
        message("[UCSC Cons] ", track_name, ": ", length(pv$pos), " positions via top-level")
        return(setNames(pv$val, as.character(pv$pos)))
      }

      # Dump first 500 chars of response for debugging unknown track structures
      message("[UCSC Cons] ", track_name, ": no data extracted. Top keys: ",
              paste(top_keys, collapse=", "),
              " | Response preview: ",
              substr(gsub("\\s+", " ", raw_text), 1, 300))
      return(NULL)
    }, error = function(e) {
      message("[UCSC Cons] Error fetching ", track_name, ": ", e$message)
      return(NULL)
    })
  }

  # Helper: try primary track name, then fallback names
  fetch_track_with_fallbacks <- function(...) {
    names_to_try <- c(...)
    for (nm in names_to_try) {
      result <- fetch_track(nm)
      if (!is.null(result) && length(result) > 0) return(result)
    }
    return(NULL)
  }

  list(
    phylop100 = fetch_track_with_fallbacks("phyloP100way", "phyloP100wayAll"),
    phylop470 = fetch_track_with_fallbacks("phyloP30way",  "phyloP30wayAll", "phyloP17wayAll"),
    phastcons = fetch_track_with_fallbacks("phastCons100way", "phastCons100wayAll")
  )
}


#' Main function: fetch all 4 conservation scores for a full gene
#' Returns data.frame: aa_pos, phylop100, phylop470, phastcons
#' All scores normalized to [0,1]; raw scores kept in *_raw columns
#' @param gene_symbol  Gene name, e.g. "PPP2R5C"
#' @param prot_length  Integer protein length in amino acids
fetch_conservation_scores <- function(gene_symbol, prot_length) {
  cache_key <- gene_symbol
  cached <- cache_get("ucsc_cons", cache_key)
  if (!is.null(cached)) {
    message("[UCSC Cons] Cache hit for ", gene_symbol)
    return(cached)
  }

  message("[UCSC Cons] Starting conservation fetch for ", gene_symbol,
          " (", prot_length, " aa)")

  # Step 1: Get exon coords
  exon_df <- fetch_ensembl_exons(gene_symbol)
  if (is.null(exon_df)) {
    message("[UCSC Cons] Cannot proceed — no exon data")
    return(NULL)
  }

  chr    <- exon_df$chr[1]
  strand <- exon_df$strand[1]

  gene_start <- min(exon_df$genomic_start)
  gene_end   <- max(exon_df$genomic_end)

  # Step 2: Fetch UCSC data for the full genomic span of CDS
  # Split into windows ≤ 50,000 bp to respect UCSC API limits
  win_size   <- 50000L
  windows    <- seq(gene_start, gene_end, by = win_size)

  # Accumulate windows as lists first, then combine once — avoids O(n^2) c() growth
  raw_windows <- list(phylop100 = list(), phylop470 = list(), phastcons = list())

  for (win_start in windows) {
    win_end <- min(win_start + win_size - 1L, gene_end)
    message("[UCSC Cons] Fetching window chr", chr, ":", win_start, "-", win_end)
    win_data <- fetch_ucsc_window(chr, win_start, win_end)
    Sys.sleep(0.3)  # rate limit courtesy
    for (nm in names(raw_windows)) {
      if (!is.null(win_data[[nm]]) && length(win_data[[nm]]) > 0)
        raw_windows[[nm]] <- c(raw_windows[[nm]], list(win_data[[nm]]))
    }
  }
  # Combine all windows into single named vectors in one pass
  all_tracks <- lapply(raw_windows, function(wlist) {
    if (length(wlist) == 0) return(NULL)
    do.call(c, wlist)
  })

  # Step 3: Map genomic positions -> amino acid positions
  # Convert named vectors to environments for O(1) lookup (avoids %in% names() per AA)
  make_env <- function(named_vec) {
    if (is.null(named_vec) || length(named_vec) == 0) return(NULL)
    e <- new.env(hash = TRUE, parent = emptyenv(), size = length(named_vec))
    for (i in seq_along(named_vec)) assign(names(named_vec)[i], named_vec[[i]], envir = e)
    e
  }
  env100 <- make_env(all_tracks$phylop100)
  env470 <- make_env(all_tracks$phylop470)
  envcst <- make_env(all_tracks$phastcons)

  message("[UCSC Cons] Track sizes: phylop100=", length(all_tracks$phylop100),
          " phylop470=", length(all_tracks$phylop470),
          " phastcons=", length(all_tracks$phastcons))

  aa_positions  <- seq_len(prot_length)
  p100 <- rep(NA_real_, prot_length)
  p470 <- rep(NA_real_, prot_length)
  pcst <- rep(NA_real_, prot_length)

  for (aa in aa_positions) {
    gpos <- aa_to_genomic(aa, exon_df)
    if (is.na(gpos)) next
    gk <- as.character(gpos)
    if (!is.null(env100) && exists(gk, envir = env100, inherits = FALSE))
      p100[aa] <- get(gk, envir = env100, inherits = FALSE)
    if (!is.null(env470) && exists(gk, envir = env470, inherits = FALSE))
      p470[aa] <- get(gk, envir = env470, inherits = FALSE)
    if (!is.null(envcst) && exists(gk, envir = envcst, inherits = FALSE))
      pcst[aa] <- get(gk, envir = envcst, inherits = FALSE)
  }

  result <- data.frame(
    aa_pos        = aa_positions,
    phylop100_raw = p100,
    phylop470_raw = p470,
    phastcons_raw = pcst,
    stringsAsFactors = FALSE
  )

  # Step 4: Normalize each score to [0, 1] (min-max, NA-aware)
  # Raw values preserve the real biological scale in tooltips
  normalize_01 <- function(x) {
    rng <- range(x, na.rm = TRUE)
    if (diff(rng) == 0) return(ifelse(is.na(x), NA_real_, 0.5))
    (x - rng[1]) / diff(rng)
  }

  result$phylop100  <- normalize_01(result$phylop100_raw)
  result$phylop470  <- normalize_01(result$phylop470_raw)
  result$phastcons  <- normalize_01(result$phastcons_raw)


  n_covered <- sum(!is.na(result$phylop100))
  message("[UCSC Cons] Done. ", n_covered, "/", prot_length,
          " AA positions covered by PhyloP 100V")

  cache_set("ucsc_cons", cache_key, result)
  return(result)
}


#' Plot the multi-conservation track
#' 3 overlapping line plots: PhyloP 100V, PhyloP 470M, PhastCons
#' Y-axis: normalized 0–1; raw score in hover tooltip
#' @param cons_data   data.frame from fetch_conservation_scores()
#' @param prot_length Integer protein length
#' @param highlight   data.frame with prot_pos column (user variants for vlines)
plot_multiconservation <- function(cons_data, prot_length, highlight = data.frame()) {

  # Colour palette — distinct, accessible, matches VarViz style
  cons_colors <- c(
    "PhyloP 100V" = "#2196F3",   # blue
    "PhyloP 470M" = "#9C27B0",   # purple
    "PhastCons"   = "#4CAF50"    # green
  )

  # Empty / loading state — return a native plotly (not ggplotly) so xaxis range
  # is set correctly and subplot() can align it with the other tracks
  if (is.null(cons_data) || nrow(cons_data) == 0) {
    msg <- if (is.null(cons_data)) "Fetching conservation scores from UCSC…" else
                                   "No conservation data returned"
    return(
      plotly::plot_ly() %>%
        plotly::add_annotations(
          x = prot_length / 2, y = 0.5, text = msg,
          showarrow = FALSE,
          font = list(size = 11, color = "#94a3b8")
        ) %>%
        plotly::layout(
          yaxis = list(
            title = list(text = "<b>Multi-<br>Conservation</b>",
                         font = list(size = vv_medium, family = "Helvetica")),
            range = c(0, 1), showticklabels = FALSE, showgrid = FALSE
          ),
          xaxis = list(range = c(0, prot_length), showgrid = FALSE),
          margin = list(l = 55, r = 10, t = 5, b = 5)
        )
    )
  }

  # Build long-format for tooltip construction
  # We'll add traces one-by-one in plotly directly for full tooltip control
  p <- plotly::plot_ly() %>%
    plotly::layout(
      yaxis = list(
        title = list(text = "<b>Multi-\nConservation</b>",
                     font = list(size = vv_medium, family = "Helvetica")),
        range      = c(-0.05, 1.1),
        tickvals   = c(0, 0.5, 1),
        ticktext   = c("0", "0.5", "1"),
        tickfont   = list(size = vv_small, family = "Helvetica"),
        showgrid   = TRUE,
        gridcolor  = "#f0f0f0",
        zeroline   = TRUE,
        zerolinecolor = "#dddddd"
      ),
      xaxis = list(
        range    = c(0, prot_length),
        showgrid = FALSE
      ),
      hovermode = "x unified",
      margin    = list(l = 55, r = 10, t = 5, b = 5),
      legend    = list(
        orientation = "h",
        x = 0, y = 1.15,
        font = list(size = 8),
        bgcolor = "rgba(255,255,255,0.7)"
      ),
      showlegend = TRUE
    )

  # Helper: build hover text with raw score
  hover_fmt <- function(label, norm_col, raw_col, unit = "") {
    paste0(label, ": %{customdata:.3f}", unit,
           " (norm: %{y:.3f})<extra></extra>")
  }

  add_cons_trace <- function(plt, col_norm, col_raw, name, color) {
    valid <- !is.na(cons_data[[col_norm]])
    if (sum(valid) == 0) return(plt)
    plt %>% plotly::add_trace(
      x          = cons_data$aa_pos[valid],
      y          = cons_data[[col_norm]][valid],
      customdata = cons_data[[col_raw]][valid],
      type       = "scatter",
      mode       = "lines",
      name       = name,
      line       = list(color = color, width = 1.6, shape = "spline", smoothing = 0.5),
      hovertemplate = paste0(
        "<b>", name, "</b><br>",
        "Pos: %{x}<br>",
        "Raw: %{customdata:.3f}<br>",
        "Norm: %{y:.3f}<extra></extra>"
      ),
      connectgaps = TRUE
    )
  }

  p <- add_cons_trace(p, "phylop100", "phylop100_raw", "PhyloP 100V", cons_colors["PhyloP 100V"])
  p <- add_cons_trace(p, "phylop470", "phylop470_raw", "PhyloP 470M", cons_colors["PhyloP 470M"])
  p <- add_cons_trace(p, "phastcons", "phastcons_raw", "PhastCons",   cons_colors["PhastCons"])

  # Variant vlines
  if (!is.null(highlight) && nrow(highlight) > 0 && "prot_pos" %in% colnames(highlight)) {
    for (vpos in highlight$prot_pos) {
      p <- p %>% plotly::add_segments(
        x = vpos, xend = vpos, y = 0, yend = 1.05,
        line = list(color = "rgba(0,0,0,0.35)", width = 0.8, dash = "dot"),
        showlegend = FALSE, hoverinfo = "skip"
      )
    }
  }

  p
}

#' ggplot static version of multi-conservation track (for PDF/PNG download via cowplot)
#' Mirrors plot_multiconservation() visually but returns a ggplot object, not a plotly widget.
plot_multiconservation_static <- function(cons_data, prot_length, highlight = data.frame()) {

  cons_colors <- c(
    "PhyloP 100V" = "#2196F3",
    "PhyloP 470M" = "#9C27B0",
    "PhastCons"   = "#4CAF50"
  )

  if (is.null(cons_data) || nrow(cons_data) == 0) {
    p <- ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = "Conservation scores unavailable",
                        size = 4, color = "#94a3b8", fontface = "italic") +
      labs(y = "Multi-\nConservation") +
      theme_bw(base_size = vv_medium) + theme_vv() +
      theme(axis.title.x = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), panel.grid = element_blank())
    return(p)
  }

  # Reshape to long format for ggplot
  score_cols   <- c("phylop100", "phylop470", "phastcons")
  score_labels <- c("PhyloP 100V", "PhyloP 470M", "PhastCons")

  long_df <- do.call(rbind, lapply(seq_along(score_cols), function(i) {
    col <- score_cols[i]
    if (!col %in% colnames(cons_data)) return(NULL)
    rows <- !is.na(cons_data[[col]])
    if (sum(rows) == 0) return(NULL)
    data.frame(
      aa_pos = cons_data$aa_pos[rows],
      score  = cons_data[[col]][rows],
      metric = score_labels[i],
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(long_df) || nrow(long_df) == 0) {
    p <- ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = "No conservation data available",
                        size = 4, color = "#94a3b8", fontface = "italic") +
      labs(y = "Multi-\nConservation") +
      theme_bw(base_size = vv_medium) + theme_vv() +
      theme(axis.title.x = element_blank(), axis.text = element_blank(),
            axis.ticks = element_blank(), panel.grid = element_blank())
    return(p)
  }

  long_df$metric <- factor(long_df$metric, levels = score_labels)

  p <- ggplot(long_df, aes(x = aa_pos, y = score, color = metric)) +
    geom_line(linewidth = 0.8, na.rm = TRUE) +
    scale_color_manual(values = cons_colors, name = NULL) +
    scale_x_continuous(limits = c(0, prot_length), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-0.05, 1.1),
                       breaks = c(0, 0.5, 1), labels = c("0", "0.5", "1")) +
    labs(y = "Multi-\nConservation") +
    theme_bw(base_size = vv_medium) + theme_vv() +
    theme(
      axis.title.x = element_blank(),
      legend.position = "top",
      legend.key.size = unit(0.4, "cm"),
      legend.text = element_text(size = 7)
    )

  p <- add_variant_vlines(p, highlight, xmax = prot_length)
  p
}


# --- Safe extraction helpers for deeply nested dbNSFP JSON ---
# MyVariant.info returns dbNSFP fields as lists or vectors; pick first non-NA value
safe_extract <- function(obj, ...) {
  val <- tryCatch({
    keys <- list(...)
    result <- obj
    for (k in keys) {
      if (is.null(result)) return(NA)
      result <- result[[k]]
    }
    if (is.null(result)) return(NA)
    if (is.list(result)) result <- unlist(result)
    result <- result[!is.na(result) & result != "." & result != ""]
    if (length(result) == 0) return(NA)
    result[1]
  }, error = function(e) NA)
  val
}

safe_extract_num <- function(obj, ...) {
  val <- safe_extract(obj, ...)
  suppressWarnings(as.numeric(val))
}

# --- Parse dbNSFP hit into structured list ---
parse_dbnsfp_scores <- function(hit) {
  if (is.null(hit) || identical(hit, NA)) {
    return(list(
      has_data = FALSE,
      pathogenicity = list(), ensemble = list(), conservation = list(), population = list(),
      acmg_tags = character(0)
    ))
  }
  
  d <- hit$dbnsfp
  cv <- hit$clinvar
  
  # ── Pathogenicity Predictors ──
  sift_score <- safe_extract_num(d, "sift", "score")
  sift_pred  <- safe_extract(d, "sift", "pred")
  
  pp2_hdiv_score <- safe_extract_num(d, "polyphen2", "hdiv", "score")
  pp2_hdiv_pred  <- safe_extract(d, "polyphen2", "hdiv", "pred")
  pp2_hvar_score <- safe_extract_num(d, "polyphen2", "hvar", "score")
  pp2_hvar_pred  <- safe_extract(d, "polyphen2", "hvar", "pred")
  
  lrt_score <- safe_extract_num(d, "lrt", "score")
  lrt_pred  <- safe_extract(d, "lrt", "pred")
  
  fathmm_score <- safe_extract_num(d, "fathmm", "score")
  fathmm_pred  <- safe_extract(d, "fathmm", "pred")
  
  provean_score <- safe_extract_num(d, "provean", "score")
  provean_pred  <- safe_extract(d, "provean", "pred")
  
  mt_score <- safe_extract_num(d, "mutationtaster", "score")
  mt_pred  <- safe_extract(d, "mutationtaster", "pred")
  
  # ── Ensemble / Meta Scores ──
  revel     <- safe_extract_num(d, "revel", "score")
  metasvm   <- safe_extract_num(d, "metasvm", "score")
  metasvm_p <- safe_extract(d, "metasvm", "pred")
  metalr    <- safe_extract_num(d, "metalr", "score")
  metalr_p  <- safe_extract(d, "metalr", "pred")
  metarnn   <- safe_extract_num(d, "metarnn", "score")
  metarnn_p <- safe_extract(d, "metarnn", "pred")
  cadd_raw  <- safe_extract_num(d, "cadd", "raw_rankscore")
  cadd_phred <- safe_extract_num(d, "cadd", "phred")
  # CADD is often a separate top-level field in MyVariant.info, not under dbnsfp
  if (is.na(cadd_phred) && !is.null(hit$cadd)) {
    cadd_phred <- safe_extract_num(hit$cadd, "phred")
    cadd_raw   <- safe_extract_num(hit$cadd, "rawscore")
  }
  dann      <- safe_extract_num(d, "dann", "score")
  am_score  <- safe_extract_num(d, "alphamissense", "am_pathogenicity")
  am_class  <- safe_extract(d, "alphamissense", "am_class")

  # ── dbscSNV splice site scores (Jian et al. 2014) ──
  # ada_score: AdaBoost classifier; rf_score: Random Forest classifier
  # Threshold > 0.6 = predicted to affect splicing
  # Note: these are meaningful for any variant near a splice site;
  # for deep-exon missense they will typically be low (~0), which is informative.
  dbscsnv_ada <- safe_extract_num(d, "dbscsnv", "ada_score")
  dbscsnv_rf  <- safe_extract_num(d, "dbscsnv", "rf_score")
  
  # ── Conservation ──
  # Actual structure confirmed from debug:
  #   d$phylop$100way_vertebrate$score
  #   d$phastcons$100way_vertebrate$score
  #   d$gerp++$rs (already works)
  all_keys <- names(d)
  
  gerp_rs    <- safe_extract_num(d, "gerp++", "rs")
  gerp_nr    <- safe_extract_num(d, "gerp++", "nr")
  phylop_100 <- safe_extract_num(d, "phylop", "100way_vertebrate", "score")
  phylop_30  <- safe_extract_num(d, "phylop", "470way_mammalian", "score")  # dbNSFP v4 uses 470-way
  if (is.na(phylop_30)) phylop_30 <- safe_extract_num(d, "phylop", "30way_mammalian", "score")
  phastcons  <- safe_extract_num(d, "phastcons", "100way_vertebrate", "score")
  siphy      <- safe_extract_num(d, "siphy_29way", "pi", "a")
  
  # Also extract primate conservation (bonus)
  phylop_pri <- safe_extract_num(d, "phylop", "17way_primate", "score")
  phastcons_pri <- safe_extract_num(d, "phastcons", "17way_primate", "score")
  
  message("[dbNSFP] Conservation: GERP=", gerp_rs, " PhyloP100=", phylop_100, 
          " PhyloP470M=", phylop_30, " PhastCons=", phastcons)
  
  # ── Population Frequencies ──
  # Population frequencies confirmed absent from dbNSFP object (all keys scanned).
  # Ancestry AFs fetched separately via fetch_pop_by_rsid() using the rsid field.
  # Initialise to NA; overwritten after parse_dbnsfp_scores() returns.
  gnomad_af  <- NA_real_
  gnomad_afr <- NA_real_
  gnomad_eur <- NA_real_
  gnomad_eas <- NA_real_
  gnomad_sas <- NA_real_
  onekg_af   <- NA_real_
  exac_af    <- NA_real_
  message("[dbNSFP] Population fields absent from dbnsfp — will fetch via rsID")
  
  # ── ClinVar from MyVariant ──
  cv_sig <- safe_extract(cv, "rcv", "clinical_significance")
  if (is.na(cv_sig)) cv_sig <- safe_extract(cv, "variant", "clinical_significance")
  cv_review <- safe_extract(cv, "rcv", "review_status")
  
  # Build structured output
  pathogenicity <- list(
    SIFT = list(score = sift_score, pred = sift_pred, 
                verdict = classify_pred("sift", sift_score, sift_pred)),
    PolyPhen2_HDIV = list(score = pp2_hdiv_score, pred = pp2_hdiv_pred,
                          verdict = classify_pred("polyphen2", pp2_hdiv_score, pp2_hdiv_pred)),
    PolyPhen2_HVAR = list(score = pp2_hvar_score, pred = pp2_hvar_pred,
                          verdict = classify_pred("polyphen2", pp2_hvar_score, pp2_hvar_pred)),
    LRT = list(score = lrt_score, pred = lrt_pred,
               verdict = classify_pred("lrt", lrt_score, lrt_pred)),
    FATHMM = list(score = fathmm_score, pred = fathmm_pred,
                  verdict = classify_pred("fathmm", fathmm_score, fathmm_pred)),
    PROVEAN = list(score = provean_score, pred = provean_pred,
                   verdict = classify_pred("provean", provean_score, provean_pred)),
    MutationTaster = list(score = mt_score, pred = mt_pred,
                          verdict = classify_pred("mutationtaster", mt_score, mt_pred))
  )
  
  ensemble <- list(
    REVEL = list(score = revel, verdict = classify_pred("revel", revel, NA)),
    MetaSVM = list(score = metasvm, pred = metasvm_p,
                   verdict = classify_pred("metasvm", metasvm, metasvm_p)),
    MetaLR = list(score = metalr, pred = metalr_p,
                  verdict = classify_pred("metalr", metalr, metalr_p)),
    MetaRNN = list(score = metarnn, pred = metarnn_p,
                   verdict = classify_pred("metarnn", metarnn, metarnn_p)),
    CADD = list(score = cadd_phred, raw = cadd_raw,
                verdict = classify_pred("cadd", cadd_phred, NA)),
    DANN = list(score = dann, verdict = classify_pred("dann", dann, NA)),
    AlphaMissense = list(score = am_score, pred = am_class,
                         verdict = classify_pred("alphamissense", am_score, am_class)),
    dbscSNV_ADA = list(score = dbscsnv_ada,
                       verdict = if (!is.na(dbscsnv_ada)) ifelse(dbscsnv_ada > 0.6, "damaging", "benign") else "unknown"),
    dbscSNV_RF  = list(score = dbscsnv_rf,
                       verdict = if (!is.na(dbscsnv_rf))  ifelse(dbscsnv_rf  > 0.6, "damaging", "benign") else "unknown")
  )
  
  conservation <- list(
    GERP_RS = gerp_rs,
    GERP_NR = gerp_nr,
    PhyloP_100V = phylop_100,
    PhyloP_470M = phylop_30,  # dbNSFP v4 uses 470-way mammalian
    PhastCons_100V = phastcons
  )
  
  population <- list(
    gnomAD_AF = gnomad_af,
    gnomAD_AFR = gnomad_afr,
    gnomAD_NFE = gnomad_eur,
    gnomAD_EAS = gnomad_eas,
    gnomAD_SAS = gnomad_sas,
    OneKG_AF = onekg_af,
    ExAC_AF = exac_af
  )
  
  list(
    has_data = TRUE,
    pathogenicity = pathogenicity,
    ensemble = ensemble,
    conservation = conservation,
    population = population,
    clinvar_sig = if (!is.na(cv_sig)) cv_sig else "",
    clinvar_review = if (!is.na(cv_review)) cv_review else ""
  )
}

# --- Classify prediction as damaging/benign/ambiguous ---
classify_pred <- function(tool, score, pred) {
  # Returns: "damaging", "benign", or "ambiguous"
  pred_lower <- tolower(as.character(pred))
  
  result <- tryCatch({
    switch(tool,
      "sift" = {
        if (!is.na(score) && score < 0.05) "damaging" else if (!is.na(score) && score >= 0.05) "benign" else if (grepl("^d", pred_lower)) "damaging" else if (grepl("^t", pred_lower)) "benign" else "ambiguous"
      },
      "polyphen2" = {
        if (!is.na(score) && score > 0.908) "damaging" else if (!is.na(score) && score > 0.446) "ambiguous" else if (!is.na(score)) "benign" else if (grepl("^d|^p", pred_lower)) "damaging" else if (grepl("^b", pred_lower)) "benign" else "ambiguous"
      },
      "lrt" = {
        if (grepl("^d", pred_lower)) "damaging" else if (grepl("^n|^u", pred_lower)) "benign" else "ambiguous"
      },
      "fathmm" = {
        if (grepl("^d", pred_lower)) "damaging" else if (grepl("^t", pred_lower)) "benign" else if (!is.na(score) && score < -1.5) "damaging" else "ambiguous"
      },
      "provean" = {
        if (!is.na(score) && score < -2.5) "damaging" else if (!is.na(score)) "benign" else if (grepl("^d", pred_lower)) "damaging" else "ambiguous"
      },
      "mutationtaster" = {
        if (grepl("^a|^d", pred_lower)) "damaging" else if (grepl("^n|^p", pred_lower)) "benign" else "ambiguous"
      },
      "revel" = {
        if (!is.na(score) && score >= 0.75) "damaging" else if (!is.na(score) && score >= 0.5) "ambiguous" else if (!is.na(score)) "benign" else "ambiguous"
      },
      "metasvm" = {
        if (grepl("^d", pred_lower)) "damaging" else if (grepl("^t", pred_lower)) "benign" else "ambiguous"
      },
      "metalr" = {
        if (grepl("^d", pred_lower)) "damaging" else if (grepl("^t", pred_lower)) "benign" else "ambiguous"
      },
      "metarnn" = {
        if (grepl("^d", pred_lower)) "damaging" else if (grepl("^t", pred_lower)) "benign" else "ambiguous"
      },
      "cadd" = {
        if (!is.na(score) && score >= 30) "damaging" else if (!is.na(score) && score >= 20) "ambiguous" else if (!is.na(score)) "benign" else "ambiguous"
      },
      "dann" = {
        if (!is.na(score) && score >= 0.96) "damaging" else if (!is.na(score) && score >= 0.9) "ambiguous" else if (!is.na(score)) "benign" else "ambiguous"
      },
      "alphamissense" = {
        if (grepl("pathogenic", pred_lower)) "damaging" else if (grepl("benign", pred_lower)) "benign" else if (!is.na(score) && score >= 0.564) "damaging" else if (!is.na(score) && score < 0.34) "benign" else "ambiguous"
      },
      "ambiguous"
    )
  }, error = function(e) "ambiguous")
  
  result
}

# ══════════════════════════════════════════════════════════════════════════════
# ACMG Evidence Tag Computation
# Based on Richards et al. (2015) Standards and Guidelines
# Tags: PS/PM/PP (pathogenic support) and BS/BP (benign support)
# ══════════════════════════════════════════════════════════════════════════════
# NOTE: Full ACMG tag computation (Richards 2015 + Tavtigian 2018 + Pejaver 2022) is done
# inline inside build_variant_table() where CCRS, domain, ClinVar match type, ConSurf,
# UniProt repeats, and all dbNSFP scores are all in scope together.
# Criteria covered: PM1/PM1_strong, PM2, PM4, PM5, PP2, PP3/PP3_moderate/PP3_strong,
#                   PP5, PS1, BP3, BP4, BP6, BP7, BS1, BS2



gnomad_freqplot <- function(gene_name, af_cutoff, highlight = data.frame(), prot_length = NULL) {
  json_parsed <- query_gnomad_api(gene_name)
  gnomad_data <- parse_gnomad_variants(json_parsed)
  prot_length <- if (!is.null(prot_length) && prot_length > 0) prot_length else max(gnomad_data$aa_pos, na.rm = TRUE)
  ymaxlim <- af_cutoff + af_cutoff * 0.1
  message("[gnomAD freqplot] af_cutoff=", signif(af_cutoff,6), 
          " | ymaxlim=", signif(ymaxlim,6),
          " | total variants=", nrow(gnomad_data),
          " | max AF in data=", signif(max(gnomad_data$exome_af, na.rm=TRUE), 6))
  
  gnomad_data <- gnomad_data %>%
    filter(!is.na(exome_af) & exome_af > 0) %>%   # drop zero/NA-AF rows — outside scale
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
  
  message("[gnomAD freqplot] AboveCutoff (triangles): ", sum(gnomad_data$freq_status == "AboveCutoff"),
          " | BelowCutoff (circles): ", sum(gnomad_data$freq_status == "BelowCutoff"))
  
  p <- ggplot(gnomad_data, aes(x = aa_pos, y = plot_freq)) +
    geom_point(aes(shape = factor(shape), color = color), size = 2, na.rm = TRUE) +
    scale_shape_manual(values = c(`16` = 16, `17` = 17)) +
    scale_color_identity() +
    scale_y_reverse(limits = c(ymaxlim, -0.0000001), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-prot_length * 0.01, prot_length + prot_length * 0.01), expand = c(0,0)) +
    geom_hline(yintercept = af_cutoff, linetype = "dashed", color = "red", linewidth = 0.4, alpha = 0.6) +
    ggplot2::annotate("text", x = prot_length * 0.98, y = af_cutoff, 
             label = paste0("AF cutoff = ", signif(af_cutoff, 3)), 
             hjust = 1, vjust = -0.5, size = 2.5, color = "red") +
    labs(y = "gnomAD Freq") +
    theme_bw(base_size = vv_medium) +
    theme_vv() +
    theme(axis.title.x = element_blank())
  
  # Color-coded variant lines: blue if variant AF is below cutoff, red if above or not found
  if (!is.null(highlight) && nrow(highlight) > 0 && "prot_pos" %in% colnames(highlight)) {
    h <- highlight[highlight$prot_pos >= 0 & highlight$prot_pos <= prot_length, , drop = FALSE]
    if (nrow(h) > 0) {
      h$vline_color <- sapply(h$prot_pos, function(pos) {
        matched <- gnomad_data[gnomad_data$aa_pos == pos, , drop = FALSE]
        if (nrow(matched) == 0) {
          message("[gnomAD vline] pos=", pos, " -> NOT FOUND in gnomAD -> RED")
          return("red")
        }
        max_af <- max(matched$exome_af, na.rm = TRUE)
        if (max_af > af_cutoff) {
          message("[gnomAD vline] pos=", pos, " max_af=", signif(max_af,3), " > cutoff=", signif(af_cutoff,3), " -> BLUE (common)")
          return("blue")
        }
        message("[gnomAD vline] pos=", pos, " max_af=", signif(max_af,3), " <= cutoff=", signif(af_cutoff,3), " -> RED (rare)")
        return("red")
      })
      h_red  <- h[h$vline_color == "red",  , drop = FALSE]
      h_blue <- h[h$vline_color == "blue", , drop = FALSE]
      if (nrow(h_red) > 0) {
        p <- p + geom_vline(data = h_red, aes(xintercept = prot_pos),
                            linetype = "dashed", color = "red", linewidth = 0.5, alpha = 0.8)
      }
      if (nrow(h_blue) > 0) {
        p <- p + geom_vline(data = h_blue, aes(xintercept = prot_pos),
                            linetype = "dashed", color = "blue", linewidth = 0.5, alpha = 0.8)
      }
    }
  }
  
  return(p)
}

densityplot <- function(gnomad_data,clinvar_data,prot_length,allele_count,highlight,gene_ptm_data,gene_ccrs_data,user_path_variants) 
{
  kde.adjust<-0.05
  p <- ggplot2::ggplot()
  p <- p + ggplot2::geom_density(data=gnomad_data[gnomad_data$gnomad_allele_count>as.numeric(allele_count),], aes(x=prot_pos),adjust=kde.adjust,fill='blue',alpha = 0.5, color="blue") 
  p <- p + scale_x_continuous(limits = c(-prot_length * 0.01, prot_length + prot_length * 0.01), expand = c(0,0))
  p <- p + ggplot2::labs(y = "Mutation Density") 
  p <- p + ggplot2::labs(x = "") 
  p <- p + theme_vv()
  if(!is.null(clinvar_data) && is.data.frame(clinvar_data) && nrow(clinvar_data) > 2 ){
    p = p + ggplot2::geom_density(aes(x=unique(clinvar_data$prot_pos)),adjust=kde.adjust, fill='red', alpha = 0.5,color="red") 
  }
  if(!is.null(user_path_variants) && is.data.frame(user_path_variants) && nrow(user_path_variants) > 2){
    p = p + ggplot2::geom_density(aes(x=unique(user_path_variants$V1)),adjust=kde.adjust, fill='green', alpha = 0.5,color="green") 
  }
  
  if(nrow(highlight) > 0 ){
    p <- add_variant_vlines(p, highlight, xmax = prot_length)
  }
  return(p)
}

clinvar_ccrsplot <- function(pfam_data,uniprot_data,gene_clinvar_data,gene_ptm_data,gene_ccrs_data,highlight=data.frame()) 
{
  begin = end = NULL
  L <- pfam_data$sequence$length   # protein length shorthand
  label_x <- L * 0.02              # inside-plot row labels: just right of x=0
  p <- ggplot2::ggplot()
  p <- p + ggplot2::ylim(0, 2.1)
  p <- p + scale_x_continuous(limits = c(0, L + L * 0.01), expand = c(0,0))
  p <- p + ggplot2::labs(y = "ClinVar/PTMs/CCRs")
  p <- p + ggplot2::theme(
    axis.title.y = element_text(size = vv_medium, face = "bold"),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y  = element_blank(),
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x  = element_blank(),
    plot.margin = unit(c(0.15, 0.3, 0.15, 0.3), "cm"),
    legend.position="none")
  # Row bands — labels sit BELOW each band (vjust=1 anchors text top at y=band_bottom-gap)
  lbl_sz  <- vv_small * 0.27
  lbl_col <- "#555555"
  p <- p + ggplot2::annotate(geom = "rect", xmin=0, xmax=L, ymin=1.7, ymax=1.85, fill="#C0C0C0", alpha=0.2)
  p <- p + ggplot2::annotate(geom = "text", label = "Mis", x = label_x, y = 1.69, size = lbl_sz, colour = lbl_col, hjust = 0, vjust = 1)
  p <- p + ggplot2::annotate(geom = "rect", xmin=0, xmax=L, ymin=1.2, ymax=1.35, fill="#C0C0C0", alpha=0.2)
  p <- p + ggplot2::annotate(geom = "text", label = "LOF", x = label_x, y = 1.19, size = lbl_sz, colour = lbl_col, hjust = 0, vjust = 1)
  p <- p + ggplot2::annotate(geom = "rect", xmin=0, xmax=L, ymin=0.7, ymax=0.85, fill="#C0C0C0", alpha=0.2)
  p <- p + ggplot2::annotate(geom = "text", label = "PTM", x = label_x, y = 0.69, size = lbl_sz, colour = lbl_col, hjust = 0, vjust = 1)
  p <- p + ggplot2::annotate(geom = "rect", xmin=0, xmax=L, ymin=0.1, ymax=0.25, fill="#C0C0C0", alpha=0.2)
  p <- p + ggplot2::annotate(geom = "text", label = "RS",  x = label_x, y = 0.09, size = lbl_sz, colour = lbl_col, hjust = 0, vjust = 1)
  # Color scale: ClinVar gold stars + PTM types
  ptm_colors <- c(
    "Phosphorylation"              = "#E69F00",
    "Phosphorylation/Modification" = "#D55E00",
    "Acetylation"                  = "#0072B2",
    "Methylation"                  = "#009E73",
    "Ubiquitination"               = "#CC79A7",
    "SUMOylation"                  = "#56B4E9",
    "Glycosylation"                = "#F0E442",
    "Lipidation"                   = "#999999",
    "Cross-link"                   = "#882255",
    "Disulfide bond"               = "#AA4499",
    "Modified residue"             = "#D55E00"
  )
  clinvar_colors <- c('0' = '#808080', '1' = '#FF8C00', '2' = '#FF00FF', '3' = '#0000FF', '4' = '#00FF00')
  all_colors <- c(clinvar_colors, ptm_colors)
  p <- p + ggplot2::scale_color_manual(values = all_colors, na.value = "#888888")
  
  if(!is.null(gene_clinvar_data) && is.data.frame(gene_clinvar_data) && nrow(gene_clinvar_data) > 0 && "type" %in% colnames(gene_clinvar_data)) {
    if( nrow(gene_clinvar_data[gene_clinvar_data$type=="missense_variant",]) > 0){
      mis_data <- gene_clinvar_data[gene_clinvar_data$type=="missense_variant",]
      mis_data$clinvar_goldstar <- as.character(mis_data$clinvar_goldstar)
      p <- p + ggplot2::geom_segment(data = mis_data,aes(x=prot_pos,xend=prot_pos,y=1.7,yend=1.85,color=clinvar_goldstar),linewidth=0.5) 
    }
    if( nrow(gene_clinvar_data[gene_clinvar_data$type!="missense_variant",]) > 0){
      lof_data <- gene_clinvar_data[gene_clinvar_data$type!="missense_variant",]
      lof_data$clinvar_goldstar <- as.character(lof_data$clinvar_goldstar)
      p <- p + ggplot2::geom_segment(data = lof_data,aes(x=prot_pos,xend=prot_pos,y=1.2,yend=1.35,color=clinvar_goldstar),linewidth=0.5) 
    }
  }
  if(!is.null(gene_ptm_data) && is.data.frame(gene_ptm_data) && nrow(gene_ptm_data) > 0 ){
    gene_ptm_data$final_ptm_group <- as.character(gene_ptm_data$final_ptm_group)
    p <- p + ggplot2::geom_segment(data=gene_ptm_data,aes(x=location,xend=location,y=0.7,yend=0.85,color=final_ptm_group),linewidth=0.5) 
  }
  if(!is.null(gene_ccrs_data) && is.data.frame(gene_ccrs_data) && nrow(gene_ccrs_data) > 0 && "prot_pos" %in% colnames(gene_ccrs_data)){
    p <- p + ggplot2::geom_segment(data=gene_ccrs_data,aes(x=prot_pos,xend=prot_pos,y=0.1,yend=0.25),colour="#ADFF2F",linewidth=0.5) 
  }
  # Variant dashed lines — spans full y range so all rows are marked
  if (!is.null(highlight) && is.data.frame(highlight) && nrow(highlight) > 0 && "prot_pos" %in% colnames(highlight)) {
    p <- p + ggplot2::geom_vline(data = highlight,
                                  ggplot2::aes(xintercept = prot_pos),
                                  linetype = "dashed", color = "black",
                                  linewidth = 0.45, alpha = 0.7)
  }
  return(p)
}

# ─────────────────────────────────────────────────────────────────────────────
# fetch_consurf_best_res: auto-fetch ConSurf evolutionary conservation scores
# from ConSurf-DB (consurfdb.tau.ac.il) using the best available PDB structure.
# Sorted by crystallographic resolution; first successful download is returned.
# Falls back to NULL if no pre-calculated ConSurf entry exists for any structure.
# ─────────────────────────────────────────────────────────────────────────────
fetch_consurf_best_res <- function(uniprot_id) {
  cached <- cache_get("consurf", uniprot_id)
  if (!is.null(cached)) {
    if (identical(cached, NA)) return(NULL)   # cached miss — don't retry
    message("[ConSurf] Cache hit for ", uniprot_id)
    return(cached)
  }

  # 1. Get PDB cross-references from UniProt
  uniprot_url  <- paste0("https://rest.uniprot.org/uniprotkb/", uniprot_id, ".json")
  uniprot_resp <- tryCatch(httr::GET(uniprot_url, httr::timeout(15)),
                           error = function(e) NULL)
  if (is.null(uniprot_resp) || httr::status_code(uniprot_resp) != 200) {
    message("[ConSurf] UniProt fetch failed for ", uniprot_id)
    return(NULL)
  }

  data        <- jsonlite::fromJSON(httr::content(uniprot_resp, as = "text", encoding = "UTF-8"))
  xrefs       <- data$uniProtKBCrossReferences
  pdb_entries <- xrefs[xrefs$database == "PDB", ]
  if (nrow(pdb_entries) == 0) {
    message("[ConSurf] No PDB structures for ", uniprot_id)
    cache_set("consurf", uniprot_id, NA)
    return(NULL)
  }
  message("[ConSurf] ", nrow(pdb_entries), " PDB structures found for ", uniprot_id)

  # 2. Sort by resolution — vapply guarantees numeric(1) per row
  #    Non-numeric entries (NMR, EM, N/A) get 99.0 so X-ray structures sort first
  resolutions <- vapply(pdb_entries$properties, function(p) {
    val <- p$value[p$key == "Resolution"]
    if (length(val) == 0) return(99.0)
    num <- suppressWarnings(as.numeric(val))
    if (is.na(num)) 99.0 else num
  }, numeric(1))
  pdb_entries <- pdb_entries[order(resolutions), ]

  # 3. Helper: extract the first single-letter chain from "A=1-300" or "A/B=1-300"
  parse_chain <- function(props) {
    raw <- props$value[props$key == "Chains"]
    if (length(raw) == 0 || nchar(raw[1]) == 0) return(NA_character_)
    before_eq <- sub("=.*", "", raw[1])          # strip residue range
    chain      <- strsplit(before_eq, "/")[[1]][1]  # take first of multi-chain
    toupper(trimws(chain))
  }

  # 4. Iterate sorted structures — cap at 3 attempts, stop after 2 consecutive timeouts
  max_attempts      <- 3L
  attempts_made     <- 0L
  consecutive_fails <- 0L
  for (i in seq_len(nrow(pdb_entries))) {
    if (attempts_made >= max_attempts) {
      message("[ConSurf] Reached ", max_attempts, " attempt limit — stopping")
      break
    }
    if (consecutive_fails >= 2L) {
      message("[ConSurf] 2 consecutive timeouts — server unreachable, stopping")
      break
    }

    pdb_id    <- toupper(pdb_entries$id[i])
    props     <- pdb_entries$properties[[i]]
    chain_id  <- parse_chain(props)
    res_val   <- props$value[props$key == "Resolution"]
    res_label <- if (length(res_val) > 0) res_val[1] else "N/A"

    if (is.na(chain_id) || nchar(chain_id) != 1) {
      message("[ConSurf] Skipping ", pdb_id, " — unparseable chain")
      next  # skipped entries don't count toward attempt limit
    }
    attempts_made <- attempts_made + 1L
    message("[ConSurf] Trying ", pdb_id, " chain ", chain_id,
            " (resolution: ", res_label, ", attempt ", attempts_made, "/", max_attempts, ")")

    consurf_url  <- paste0("https://consurfdb.tau.ac.il/DB/", pdb_id, "/",
                           chain_id, "/consurf_summary.txt")
    consurf_resp <- tryCatch(httr::GET(consurf_url, httr::timeout(4)),
                             error = function(e) {
                               message("[ConSurf] Timeout/error: ", e$message); NULL
                             })
    if (is.null(consurf_resp) || httr::status_code(consurf_resp) != 200) {
      consecutive_fails <- consecutive_fails + 1L
      next
    }
    consecutive_fails <- 0L

    raw_text <- httr::content(consurf_resp, as = "text", encoding = "UTF-8")
    lines    <- strsplit(raw_text, "\\n")[[1]]

    # 5. Locate header dynamically — more robust than fixed skip=N
    header_idx <- grep("^\\s*(POS|#POS|# POS)", lines)[1]
    if (is.na(header_idx)) {
      message("[ConSurf] Header not found in ", pdb_id, " — skipping")
      next
    }

    parsed <- tryCatch({
      df_raw <- read.table(
        text         = paste(lines[header_idx:length(lines)], collapse = "\\n"),
        header       = TRUE,
        fill         = TRUE,
        comment.char = "#",
        stringsAsFactors = FALSE
      )
      # Standardise to the columns the rest of VarViz expects
      needed <- c("POS", "SEQ", "SCORE", "COLOR")
      if (!all(needed %in% colnames(df_raw))) stop("Missing expected columns")

      df_raw$POS   <- as.integer(gsub("[^0-9]",    "", trimws(as.character(df_raw$POS))))
      df_raw$SCORE <- as.numeric(gsub("[^0-9.\\-]", "", trimws(as.character(df_raw$SCORE))))
      df_raw$COLOR <- as.numeric(gsub("[^0-9.]",    "", trimws(as.character(df_raw$COLOR))))
      df_raw$SEQ   <- gsub("[^A-Za-z]", "", trimws(as.character(df_raw$SEQ)))
      if ("B.E" %in% colnames(df_raw)) {
        df_raw$BE <- tolower(gsub("[^A-Za-z]", "", trimws(as.character(df_raw$B.E))))
      } else {
        df_raw$BE <- NA_character_
      }
      df_raw <- df_raw[!is.na(df_raw$POS) & !is.na(df_raw$COLOR), ]
      df_raw
    }, error = function(e) {
      message("[ConSurf] Parse error for ", pdb_id, ": ", e$message); NULL
    })

    if (!is.null(parsed) && nrow(parsed) > 0) {
      message("[ConSurf] Success — ", pdb_id, " chain ", chain_id,
              " (", nrow(parsed), " residues, res=", res_label, ")")
      cache_set("consurf", uniprot_id, parsed)
      return(parsed)
    }
  }

  if (consecutive_fails >= 2L) {
    message("[ConSurf] Server unreachable — not caching so retry is possible")
  } else {
    message("[ConSurf] No pre-calculated data found for ", uniprot_id)
    cache_set("consurf", uniprot_id, NA)
  }
  return(NULL)
}


conservplot <- function(consurf_score, prot_length, highlight = data.frame()) {
  
  consurf_colors <- c(
    "9" = "#a02560",
    "8" = "#f07dab",
    "7" = "#fac9de",
    "6" = "#fcedf4",
    "5" = "#ffffff",
    "4" = "#eaffff",
    "3" = "#d7ffff",
    "2" = "#8cffff",
    "1" = "#10c8d1"
  )
  
  # Prepare fill factor
  consurf_score$COLOR <- factor(as.integer(consurf_score$COLOR), levels = 1:9)
  
  # Precompute contour line coordinates (ggplotly-safe)
  consurf_score$x_start <- consurf_score$POS - 0.5
  consurf_score$x_end   <- consurf_score$POS + 0.5
  
  p <- ggplot2::ggplot(consurf_score, aes(x = POS, y = as.integer(as.character(COLOR)), fill = COLOR)) +
    geom_col(width = 1) +
    geom_segment(aes(x = x_start, xend = x_end, y = as.integer(as.character(COLOR)), yend = as.integer(as.character(COLOR))),
                 color = "black", size = 0.4) +
    scale_y_reverse(limits = c(9.5, 0), breaks = 0:9, expand = c(0, 0)) +
    scale_fill_manual(values = consurf_colors, na.value = "#ffff96") +
    scale_x_continuous(
      limits = c(-prot_length * 0.01, prot_length + prot_length * 0.01),
      expand = c(0, 0)
    ) +
    labs(y = "Conservation", x = "Position") +
    theme_minimal(base_size = vv_medium) +
    theme_vv()
  
  p <- add_variant_vlines(p, highlight, xmax = prot_length)
  
  return(p)
}


pfamplot <- function(pfam_data,uniprot_data,gene_clinvar_data,highlight,label,for_plotly=FALSE,revel_scores=NULL)
{
  begin = end = NULL
  prot_len <- pfam_data$sequence$length

  # ── Layout (tight, matches old working version) ───────────────────────────
  # ylim 0–4.0; backbone y=1.5; domain bar ymin=1.0 ymax=2.0
  # Below bar: topo y=0.50, topo labels y=0.38, PTM y=0.20
  # Above bar: disulfide from 2.0 up; lollipop stems 2.05->2.70; dot 2.70
  # Protein name: caption (theme) — not in data space
  # Labels added via plotly annotations (angle not supported in ggplotly)

  p <- ggplot2::ggplot()
  p <- p + ggplot2::ylim(-0.60, 4.5)
  p <- p + scale_x_continuous(limits = c(-prot_len * 0.01, prot_len + prot_len * 0.01), expand = c(0,0))
  p <- p + ggplot2::labs(
    y = "Uniprot Domains",
    x = "Amino acid number"
  )
  p <- p + ggplot2::annotate("text",
    x = pfam_data$sequence$length / 2, y = -0.50,
    label = paste(pfam_data$uniProtkbId, "-", pfam_data$primaryAccession, "-",
                  pfam_data$proteinDescription$recommendedName$fullName$value,
                  "(", pfam_data$sequence$length, "aa)", sep=" "),
    color = "black", size = 2.4, hjust = 0.5, inherit.aes = FALSE)

  # ── Backbone ─────────────────────────────────────────────────────────────
  p <- p + ggplot2::geom_segment(
    ggplot2::aes(x=0, xend=pfam_data$sequence$length, y=1.5, yend=1.5),
    linewidth=2, color="grey80", inherit.aes=FALSE)

  # ── Domain bars ──────────────────────────────────────────────────────────
  darken_hex <- function(hex) {
    m <- grDevices::col2rgb(hex) / 255
    grDevices::rgb(m[1]*0.65, m[2]*0.65, m[3]*0.65)
  }
  if (length(pfam_data$Domain) > 0) {
    dom_df <- as.data.frame(pfam_data$Domain, stringsAsFactors=FALSE)
    dom_df$border_col <- sapply(dom_df$color, darken_hex)
    p <- p + ggplot2::geom_rect(
      data=dom_df,
      ggplot2::aes(xmin=start, xmax=end, ymin=1.0, ymax=2.0,
                   fill=I(color), colour=I(border_col)),
      linewidth=0.4, inherit.aes=FALSE)
    dom_wide <- dom_df[(dom_df$end - dom_df$start) > prot_len * 0.07, , drop=FALSE]
    if (nrow(dom_wide) > 0) {
      dom_wide$bar_pct <- (dom_wide$end - dom_wide$start) / prot_len * 100
      dom_wide$max_chr <- pmax(4L, as.integer(dom_wide$bar_pct * 0.9))
      dom_wide$lbl     <- mapply(function(d, n) {
        if (nchar(d) <= n) d else paste0(substr(d, 1, n-1), "\u2026")
      }, dom_wide$description, dom_wide$max_chr)
      p <- p + ggplot2::geom_text(
        data=dom_wide,
        ggplot2::aes(x=(start+end)/2, y=1.5, label=lbl),
        size=2.0, fontface="bold", colour="#1e293b", inherit.aes=FALSE)
    }
  }

  # ── Region bars ──────────────────────────────────────────────────────────
  if (length(pfam_data$Region) > 0) {
    reg_df <- as.data.frame(pfam_data$Region, stringsAsFactors=FALSE)
    reg_df$border_col <- sapply(reg_df$color, darken_hex)
    p <- p + ggplot2::geom_rect(
      data=reg_df,
      ggplot2::aes(xmin=start, xmax=end, ymin=1.0, ymax=2.0,
                   fill=I(color), colour=I(border_col)),
      linewidth=0.4, inherit.aes=FALSE)
    reg_wide <- reg_df[(reg_df$end - reg_df$start) > prot_len * 0.07, , drop=FALSE]
    if (nrow(reg_wide) > 0) {
      reg_wide$bar_pct <- (reg_wide$end - reg_wide$start) / prot_len * 100
      reg_wide$max_chr <- pmax(4L, as.integer(reg_wide$bar_pct * 0.9))
      reg_wide$lbl     <- mapply(function(d, n) {
        if (nchar(d) <= n) d else paste0(substr(d, 1, n-1), "\u2026")
      }, reg_wide$description, reg_wide$max_chr)
      p <- p + ggplot2::geom_text(
        data=reg_wide,
        ggplot2::aes(x=(start+end)/2, y=1.5, label=lbl),
        size=2.0, fontface="bold", colour="#1e293b", inherit.aes=FALSE)
    }
  }

  # ── TM helices ───────────────────────────────────────────────────────────
  if ("transmem" %in% uniprot_data$type &&
      nrow(uniprot_data[uniprot_data$type=="transmem",]) > 0) {
    tm_data <- uniprot_data[uniprot_data$type=="transmem", , drop=FALSE]
    p <- p + ggplot2::geom_rect(data=tm_data,
      ggplot2::aes(xmin=start, xmax=end, ymin=1.05, ymax=1.95),
      fill="#e07b39", color="#b8621b", alpha=0.85, linewidth=0.3, inherit.aes=FALSE)
    tm_wide <- tm_data[(tm_data$end - tm_data$start) > prot_len * 0.03, , drop=FALSE]
    if (nrow(tm_wide) > 0) {
      p <- p + ggplot2::geom_text(data=tm_wide,
        ggplot2::aes(x=(start+end)/2, y=1.5,
                     label=ifelse(!is.na(description)&nchar(description)>0, description, "TM")),
        size=1.8, color="white", fontface="bold", inherit.aes=FALSE)
    }
  }

  # ── Signal peptide ────────────────────────────────────────────────────────
  if (nrow(uniprot_data[uniprot_data$type=="signal",]) > 0) {
    p <- p + ggplot2::geom_segment(
      data=uniprot_data[uniprot_data$type=="signal",],
      ggplot2::aes(x=start, xend=end, y=1.5, yend=1.5),
      linewidth=3, color="red", inherit.aes=FALSE)
  }

  # ── Topological domains ───────────────────────────────────────────────────
  if ("topo_dom" %in% uniprot_data$type &&
      nrow(uniprot_data[uniprot_data$type=="topo_dom",]) > 0) {
    td_data <- uniprot_data[uniprot_data$type=="topo_dom", , drop=FALSE]
    td_data$td_color <- sapply(td_data$description, function(d) {
      if (is.na(d)) return("#888888")
      dl <- tolower(d)
      if (grepl("extracellular|exoplasmic", dl)) return("#4a90d9")
      if (grepl("cytoplasmic|intracellular", dl)) return("#50b050")
      if (grepl("lumenal", dl)) return("#9b59b6")
      return("#888888")
    }, USE.NAMES=FALSE)
    p <- p + ggplot2::geom_segment(data=td_data,
      ggplot2::aes(x=start, xend=end, y=0.50, yend=0.50),
      color=td_data$td_color, linewidth=3.5, alpha=0.8, inherit.aes=FALSE)
    td_data$width <- td_data$end - td_data$start
    td_wide <- td_data[td_data$width > prot_len * 0.05, , drop=FALSE]
    if (nrow(td_wide) > 0) {
      p <- p + ggplot2::geom_text(data=td_wide,
        ggplot2::aes(x=(start+end)/2, y=0.35,
                     label=ifelse(!is.na(description), description, "")),
        size=1.6, color="#333333", inherit.aes=FALSE)
    }
  }

  # ── Disulfide bridges (above bar, base=2.0, taller brackets) ────────────
  if (nrow(uniprot_data[uniprot_data$type=="disulfid",]) > 0) {
    values <- uniprot_data[uniprot_data$type=="disulfid" & !is.na(uniprot_data$end),]
    values <- values[order(values$start),]
    values$y_end <- 2.40
    k <- 0
    for (i in seq_len(nrow(values))) {
      if (i %% 4 == 0) k <- 0
      values$y_end[i] <- 2.40 + k
      k <- k + 0.12
    }
    p <- p + ggplot2::geom_segment(data=values,
      ggplot2::aes(x=start, y=2.0, xend=start, yend=y_end),
      color="#3b82f6", linewidth=0.5, inherit.aes=FALSE)
    p <- p + ggplot2::geom_segment(data=values,
      ggplot2::aes(x=end, y=2.0, xend=end, yend=y_end),
      color="#3b82f6", linewidth=0.5, inherit.aes=FALSE)
    p <- p + ggplot2::geom_segment(data=values,
      ggplot2::aes(x=start, y=y_end, xend=end, yend=y_end),
      color="#3b82f6", linewidth=0.5, inherit.aes=FALSE)
  }

  # ── PTM / active sites ────────────────────────────────────────────────────
  if (nrow(uniprot_data[uniprot_data$type=="mod_res",]) > 0) {
    p <- p + ggplot2::geom_segment(
      data=uniprot_data[uniprot_data$type=="mod_res",],
      ggplot2::aes(x=start, xend=start, y=0.25, yend=0.55),
      color="black", linewidth=0.3, inherit.aes=FALSE)
    p <- p + ggplot2::geom_point(
      data=uniprot_data[uniprot_data$type=="mod_res",],
      ggplot2::aes(x=start, y=0.25, fill=mod_res_group),
      shape=23, size=1.6, colour="black", stroke=0.4, inherit.aes=FALSE)
  }
  if (nrow(uniprot_data[uniprot_data$type=="act_site",]) > 0) {
    p <- p + ggplot2::geom_segment(
      data=uniprot_data[uniprot_data$type=="act_site",],
      ggplot2::aes(x=start, xend=start, y=0.25, yend=0.55),
      color="black", linewidth=0.3, inherit.aes=FALSE)
    p <- p + ggplot2::geom_point(
      data=uniprot_data[uniprot_data$type=="act_site",],
      ggplot2::aes(x=start, y=0.25),
      shape=22, size=1.6, fill="pink", colour="black", stroke=0.4, inherit.aes=FALSE)
  }
  if (nrow(uniprot_data[uniprot_data$type=="mutagen",]) > 0) {
    p <- p + ggplot2::geom_point(
      data=uniprot_data[uniprot_data$type=="mutagen",],
      ggplot2::aes(x=start, y=0.25),
      shape=10, size=1.6, inherit.aes=FALSE)
  }

  p <- p + ggplot2::theme(
    axis.text.y   = element_blank(),
    axis.ticks.y  = element_blank(),
    axis.line.y   = element_blank(),
    axis.title.y  = element_text(size=vv_medium, face="bold"),
    axis.title.x  = element_text(size=vv_medium, face="bold"),
    axis.text.x   = element_text(size=vv_small),
    plot.margin   = unit(c(0.15, 0.3, 0.15, 0.3), "cm"),
    legend.position = "none")

  # ── Lollipops (ggplot geoms — drawn for BOTH ggplot and plotly paths) ─────
  # Labels are NOT drawn here for plotly path — added via layout(annotations=)
  # so angle=-65 works. For static path, ggplot geom_text handles labels.
  if (nrow(highlight) > 0) {
    prot_len_l <- pfam_data$sequence$length
    h <- highlight[!is.na(highlight$prot_pos) &
                   highlight$prot_pos >= 0 &
                   highlight$prot_pos <= prot_len_l, , drop=FALSE]
    if (nrow(h) > 0) {
      aa_col <- c(D="#e04040",E="#e04040",H="#4472c4",K="#4472c4",R="#4472c4",
                  S="#2ca02c",T="#2ca02c",N="#2ca02c",Q="#2ca02c",C="#2ca02c",
                  F="#e67e22",W="#e67e22",Y="#e67e22",
                  A="#8e44ad",V="#8e44ad",I="#8e44ad",L="#8e44ad",M="#8e44ad",G="#8e44ad",
                  P="#8B4513")
      bare        <- sub("^p\\.", "", as.character(h$Mutation))
      ref_aa      <- toupper(substr(bare, 1, 1))
      var_aa      <- toupper(substr(bare, nchar(bare), nchar(bare)))
      h$ref_color <- ifelse(ref_aa %in% names(aa_col), aa_col[ref_aa], "#555555")
      h$var_color <- ifelse(var_aa %in% names(aa_col), aa_col[var_aa], "#555555")
      # Baked y-coords: stem 2.85->3.60 (clears disulfide ceiling ~2.76)
      h$stem_base <- 1.5
      h$stem_dot  <- 3.60
      p <- p + ggplot2::geom_segment(
        data=h, ggplot2::aes(x=prot_pos, y=stem_base, xend=prot_pos, yend=stem_dot),
        color=h$ref_color, linewidth=0.5, lineend="round", inherit.aes=FALSE)
      p <- p + ggplot2::geom_point(
        data=h, ggplot2::aes(x=prot_pos, y=stem_dot),
        color=h$var_color, fill=h$var_color, size=2, shape=21, stroke=0.5,
        inherit.aes=FALSE)

      # Static-only labels (plotly path gets them via layout(annotations))
      if (label == "Yes" && !for_plotly) {
        h$ldr_top <- 3.78
        p <- p + ggplot2::geom_segment(
          data=h, ggplot2::aes(x=prot_pos, y=stem_dot, xend=prot_pos, yend=ldr_top),
          color="grey60", linewidth=0.3, inherit.aes=FALSE)
        prot_len_local <- pfam_data$sequence$length
        h_sorted <- h[order(h$prot_pos), , drop=FALSE]
        stagger_levels <- c(3.85, 4.05, 4.22, 4.36, 4.42, 4.46)
        min_gap <- prot_len_local * 0.05
        h_sorted$label_y <- stagger_levels[1]
        h_sorted$stagger_idx <- 1L
        if (nrow(h_sorted) > 1) {
          for (j in 2:nrow(h_sorted)) {
            if (h_sorted$prot_pos[j] - h_sorted$prot_pos[j-1] < min_gap) {
              h_sorted$stagger_idx[j] <- (h_sorted$stagger_idx[j-1] %% length(stagger_levels)) + 1L
            } else {
              h_sorted$stagger_idx[j] <- 1L
            }
            h_sorted$label_y[j] <- stagger_levels[h_sorted$stagger_idx[j]]
          }
        }
        p <- p + ggplot2::geom_text(
          data=h_sorted,
          ggplot2::aes(x=prot_pos, y=label_y, label=Mutation),
          size=2.1, hjust=0, angle=65,
          color=h_sorted$var_color,
          family="Times", fontface="bold", inherit.aes=FALSE)
      }
    }
  }

  return(p)
}







# Define server logic for slider examples
# ============================================================
# ACMG Classification Engine
# Implements BOTH:
#   (a) Rule-based ACMG/AMP 2015 combination logic (Richards et al.)
#   (b) Tavtigian 2020 Bayesian point thresholds
# Returns: list(classification, rule, pts)
#   classification: "Pathogenic"|"Likely Pathogenic"|"VUS"|"Likely Benign"|"Benign"
#   rule: human-readable rule string, e.g. "1 PS + 2 PM"
#   pts: integer Bayesian score
# ============================================================
classify_acmg <- function(tags_vec) {
  # ── Point weights (Tavtigian 2020) ──────────────────────────────────────
  tag_pts_map <- c(
    PVS1=8,
    PS1=4, PS1_moderate=2, PS1_supporting=1, PS2=4, PS3=4, PS3_supporting=1, PS4=4,
    PM1_strong=4, PP3_strong=4, PP1_strong=4,
    PM1=2, PM2=2, PM3=1, PM3_moderate=2, PM3_strong=4, PM4=2, PM5=2, PM6=2, PP3_moderate=2, PP1_moderate=2,
    PP1=1, PP2=1, PP3=1, PP4=1, PP5=1,
    BA1=-8,
    BS1=-4, BS2=-4, BS3=-4, BS4=-4, BP6=-4,
    BP1=-1, BP2=-1, BP3=-1, BP4=-1, BP5=-1, BP7=-1
  )
  tags <- trimws(tags_vec)
  pts  <- sum(tag_pts_map[intersect(tags, names(tag_pts_map))], na.rm = TRUE)

  # ── Tag counts by tier ───────────────────────────────────────────────────
  n_pvs <- sum(grepl("^PVS",        tags))
  n_ps  <- sum(grepl("^PS[0-9]",    tags))
  # PM1_strong counts as PS-equivalent (4 pts)
  n_ps_eq <- n_ps + sum(tags %in% c("PM1_strong","PP3_strong"))
  n_pm  <- sum(grepl("^PM",         tags) & !tags %in% c("PM1_strong"))
  n_pp  <- sum(grepl("^PP",         tags) & !tags %in% c("PP3_strong","PP3_moderate"))
  n_pp_eq <- n_pp + sum(tags == "PP3_moderate")   # PP3_moderate counts as 2 PP
  n_bs  <- sum(grepl("^BS",         tags))
  n_bp  <- sum(grepl("^BP",         tags))
  has_ba1 <- "BA1" %in% tags

  # ── Rule-based classifier (Richards 2015) ────────────────────────────────
  # Returns list(class, rule) or NULL if no rule matched
  rule_classify <- function() {
    # ── Pathogenic rules ──
    if (n_pvs >= 1 && n_ps_eq >= 1)
      return(list(class="Pathogenic", rule=paste0("PVS1 + ", n_ps_eq, " PS")))
    if (n_pvs >= 1 && n_pm >= 2)
      return(list(class="Pathogenic", rule=paste0("PVS1 + ", n_pm, " PM")))
    if (n_pvs >= 1 && n_pm == 1 && n_pp_eq >= 2)
      return(list(class="Pathogenic", rule=paste0("PVS1 + 1 PM + ", n_pp_eq, " PP")))
    if (n_pvs >= 1 && n_pp_eq >= 4)
      return(list(class="Pathogenic", rule=paste0("PVS1 + ", n_pp_eq, " PP")))
    if (n_ps_eq >= 2)
      return(list(class="Pathogenic", rule=paste0(n_ps_eq, " PS")))
    if (n_ps_eq == 1 && n_pm >= 3)
      return(list(class="Pathogenic", rule=paste0("1 PS + ", n_pm, " PM")))
    if (n_ps_eq == 1 && n_pm == 2 && n_pp_eq >= 2)
      return(list(class="Pathogenic", rule=paste0("1 PS + 2 PM + ", n_pp_eq, " PP")))
    if (n_ps_eq == 1 && n_pm == 1 && n_pp_eq >= 4)
      return(list(class="Pathogenic", rule=paste0("1 PS + 1 PM + ", n_pp_eq, " PP")))
    # ── Likely Pathogenic rules ──
    if (n_pvs >= 1 && n_pm == 1)
      return(list(class="Likely Pathogenic", rule="PVS1 + 1 PM"))
    if (n_ps_eq == 1 && n_pm >= 1 && n_pm <= 2)
      return(list(class="Likely Pathogenic", rule=paste0("1 PS + ", n_pm, " PM")))
    if (n_ps_eq == 1 && n_pp_eq >= 2)
      return(list(class="Likely Pathogenic", rule=paste0("1 PS + ", n_pp_eq, " PP")))
    if (n_pm >= 3)
      return(list(class="Likely Pathogenic", rule=paste0(n_pm, " PM")))
    if (n_pm == 2 && n_pp_eq >= 2)
      return(list(class="Likely Pathogenic", rule=paste0("2 PM + ", n_pp_eq, " PP")))
    if (n_pm == 1 && n_pp_eq >= 4)
      return(list(class="Likely Pathogenic", rule=paste0("1 PM + ", n_pp_eq, " PP")))
    # ── Benign rules ──
    if (has_ba1)
      return(list(class="Benign", rule="BA1 (standalone)"))
    if (n_bs >= 2)
      return(list(class="Benign", rule=paste0(n_bs, " BS")))
    if (n_bs >= 1 && n_bp >= 2)
      return(list(class="Benign", rule=paste0("1 BS + ", n_bp, " BP")))
    # ── Likely Benign rules ──
    if (n_bs >= 1 && n_bp >= 1)
      return(list(class="Likely Benign", rule=paste0("1 BS + 1 BP")))
    if (n_bp >= 2)
      return(list(class="Likely Benign", rule=paste0(n_bp, " BP")))
    return(NULL)
  }

  rb <- rule_classify()

  # ── Tavtigian Bayesian thresholds (fallback / confirmation) ──────────────
  # VUS sub-tiers (SVC v4.0 direction, Tavtigian 2020 points 0–5):
  #   VUS-High  (4–5 pts): leans pathogenic
  #   VUS-Mid   (2–3 pts): truly uncertain
  #   VUS-Low   (0–1 pts): leans benign
  bay_class <- if (has_ba1)         "Benign" else
               if (pts >= 10)       "Pathogenic" else
               if (pts >= 6)        "Likely Pathogenic" else
               if (pts <= -7)       "Benign" else
               if (pts <= -4)       "Likely Benign" else
               if (pts >= 4)        "VUS-High" else
               if (pts >= 2)        "VUS-Mid" else
                                    "VUS-Low"

  # Use rule-based when it fires; fall back to Bayesian with pts note
  if (!is.null(rb)) {
    classification <- rb$class
    rule_str       <- rb$rule
  } else {
    classification <- bay_class
    rule_str       <- paste0("score ", pts, " pts")
  }

  # Format tag summary string, e.g. "PM1↑ PM2 PP2 PP3"
  tag_summary <- if (length(tags) > 0) {
    paste(sub("_strong$","\u2b06", sub("_moderate$","\u2191", tags)), collapse=" ")
  } else "—"

  # Build per-tag breakdown string: e.g. "PM1_strong(+4) + PM2(+2) + PP2(+1) + PP3(+1)"
  # Separate pathogenic (+) from benign (-) contributions
  path_parts   <- c()
  benign_parts <- c()
  for (t in tags) {
    v <- tag_pts_map[t]
    if (!is.na(v) && v != 0) {
      lbl <- sub("_strong$","⬆", sub("_moderate$","↑", t))
      if (v > 0) path_parts  <- c(path_parts,  paste0(lbl, "(+", v, ")")) else       benign_parts <- c(benign_parts, paste0(lbl, "(", v, ")"))
    }
  }
  breakdown_parts <- c(path_parts, benign_parts)
  # Numeric-only string for compact display: e.g. "+4+2+1+1" or "+4+2-1"
  nums_str <- if (length(breakdown_parts) > 0) {
    nums <- sapply(c(path_parts, benign_parts), function(x) {
      m <- regmatches(x, regexpr("[+-][0-9]+", x)); if (length(m)) m else ""
    })
    paste0("(", paste(nums, collapse=""), ")")
  } else ""

  list(
    classification = classification,
    rule           = rule_str,
    pts            = pts,
    tag_summary    = tag_summary,
    breakdown      = if (length(breakdown_parts)>0) paste(breakdown_parts, collapse=" + ") else "",
    pts_str        = nums_str
  )
}

# ============================================================
# ACMG Clinical Comment Generator
# Produces natural-language variant interpretation text from ACMG tags + data
# ============================================================
generate_acmg_comment <- function(acmg_tags_str, gene, mut, gnomad_af, gnomad_ac,
                                   gnomad_nhomalt,
                                   clinvar_sig, clinvar_name, clinvar_vcv, clinvar_vcv_pos,
                                   clinvar_trait, clinvar_match_type, domain, pp2_applies,
                                   revel_score, am_score, cadd_score, metasvm_v,
                                   phylop_val, phastcons_val, gerp_val, consurf_grade,
                                   variant_pos,
                                   clingen_class = "", clingen_disease = "", clingen_moi = "",
                                   ps3_proxy = FALSE, bs3_proxy = FALSE,
                                   gevir_gene_pct = NA_real_) {

  safe1 <- function(x, default = "") {
    if (is.null(x) || length(x) == 0) return(default)
    x <- x[1]; if (is.na(x)) return(default); as.character(x)
  }
  safe_num <- function(x) {
    if (is.null(x) || length(x) == 0) return(NA_real_)
    suppressWarnings(as.numeric(x[1]))
  }

  acmg_tags_str <- safe1(acmg_tags_str)
  gene          <- safe1(gene, "gene")
  mut           <- safe1(mut)
  gnomad_af     <- safe_num(gnomad_af)
  gnomad_nhomalt <- safe_num(gnomad_nhomalt)
  clinvar_name  <- safe1(clinvar_name)
  clinvar_vcv   <- safe1(clinvar_vcv)
  clinvar_vcv_pos <- safe1(clinvar_vcv_pos)
  clinvar_trait <- safe1(clinvar_trait)
  domain        <- safe1(domain)
  revel_score   <- safe_num(revel_score)
  am_score      <- safe_num(am_score)
  cadd_score    <- safe_num(cadd_score)
  metasvm_v     <- safe1(metasvm_v)
  consurf_grade <- safe1(consurf_grade)
  variant_pos   <- safe_num(variant_pos)

  tags <- trimws(strsplit(acmg_tags_str, ",")[[1]])
  tags <- tags[nchar(tags) > 0]
  if (length(tags) == 0) return("")
  sentences <- character(0)
  s <- function(x) sentences <<- c(sentences, x)

  # Clean domain: strip UniProt interaction boilerplate, keep short type names
  clean_domain <- function(dom) {
    if (nchar(dom) == 0) return("")
    parts <- trimws(strsplit(dom, ";")[[1]])
    non_structural_pattern <- paste0(
      "Interaction with|Required for|Sufficient for|Involved in|Necessary for|",
      "Binds to|Mediates|Responsible for|Important for|Essential for|",
      "Needed for|Critical for|Associates with|",
      "Disordered|disorder|",
      "Pro residues|Gly residues|Ala residues|Ser residues|compositional|",
      "basic residues|acidic residues|polar residues"
    )
    keep <- parts[
      nchar(parts) < 60 &
      !grepl(non_structural_pattern, parts, ignore.case = TRUE)
    ]
    # If nothing structural remains, return "" — do NOT fall back to interaction text
    if (length(keep) == 0) return("")
    paste(keep, collapse = "; ")
  }
  domain_clean <- clean_domain(domain)

  # Clean trait: reject transcript/accession IDs, gene symbols, too-short strings
  clean_trait <- function(trait, gene_sym) {
    if (nchar(trait) == 0) return("")
    if (grepl("^(NM_|NP_|NG_|NC_|ENST|ENSG|LRG_|chr)", trait)) return("")
    if (tolower(trimws(trait)) == tolower(gene_sym)) return("")
    if (nchar(trimws(trait)) < 5) return("")
    trait
  }
  trait_clean <- clean_trait(clinvar_trait, gene)
  disease_str <- if (nchar(trait_clean) > 0) trait_clean else paste0(gene, "-related disorder")

  vcv_cite <- function(vcv) {
    vcv <- safe1(vcv)
    if (nchar(vcv) == 0) return("")
    ids <- trimws(strsplit(vcv, ",")[[1]])
    ids <- ids[nchar(ids) > 0 & !grepl("^uid:", ids)]
    if (length(ids) == 0) return("")
    paste0("ClinVar ID: ", paste(ids, collapse = ", "))
  }

  fmt_af <- function(af) {
    af <- safe_num(af)
    if (is.na(af)) return(NULL)
    pct <- af * 100
    if (pct == 0)      return("0%")
    # Enough decimal places to show 2–3 significant figures at each scale
    if (pct < 0.0001)  return(paste0(formatC(pct, format = "e", digits = 2), "%"))  # e.g. 2.50e-05%
    if (pct < 0.01)    return(sprintf("%.5f%%", pct))   # e.g. 0.00171%  (covers 1e-06 to 1e-04 AF)
    if (pct < 0.1)     return(sprintf("%.3f%%", pct))   # e.g. 0.040%
    return(sprintf("%.2f%%", pct))                       # e.g. 5.00%
  }

  # Population frequency
  if ("PM2" %in% tags) {
    af_txt <- fmt_af(gnomad_af)
    if (!is.null(af_txt) && af_txt != "0%")
      s(paste0("The variant is observed at an extremely low frequency in the gnomAD v4 dataset (total allele frequency: ", af_txt, ").")) else
      s("The variant is absent from the gnomAD v4 population database.")
  } else if ("BA1" %in% tags) {
    af_txt <- fmt_af(gnomad_af)
    s(paste0("The variant is extremely common in the general population (gnomAD AF: ", if (!is.null(af_txt)) af_txt else ">5%", "), meeting the BA1 threshold for standalone Benign classification (Richards 2015)."))
  } else if ("BS2" %in% tags) {
    af_txt <- fmt_af(gnomad_af)
    nhom   <- safe_num(gnomad_nhomalt)
    if (!is.na(nhom) && nhom > 0)
      s(paste0("The variant has been observed in homozygous state in ", nhom, " unaffected individuals in gnomAD (gnomAD AF: ", if (!is.null(af_txt)) af_txt else "present", "), providing strong evidence against pathogenicity (BS2).")) else
      s(paste0("The variant is common in the general population (gnomAD AF: ", if (!is.null(af_txt)) af_txt else "high", "), suggesting it is likely benign (BS2)."))
  } else if ("BS1" %in% tags) {
    af_txt <- fmt_af(gnomad_af)
    s(paste0("The variant frequency exceeds what is expected for a pathogenic variant (gnomAD AF: ", if (!is.null(af_txt)) af_txt else "elevated", ") (BS1)."))
  }

  # Functional domain
  if ("PM1_strong" %in% tags) {
    if (nchar(domain_clean) > 0)
      s(paste0("The variant is located in a critical and highly conserved functional domain (", domain_clean, "), supporting PM1 at strong evidence strength.")) else
      s("The variant is located in a highly constrained region with strong conservation evidence, supporting PM1 at strong evidence strength.")
  } else if ("PM1" %in% tags) {
    if (nchar(domain_clean) > 0)
      s(paste0("The variant is located in a known functional domain (", domain_clean, ") (PM1).")) else
      s("The variant falls within a mutational hotspot or constrained genomic region (PM1).")
  }

  # Missense mechanism
  if ("PP2" %in% tags) {
    gevir_pp2_str <- if (!is.na(gevir_gene_pct) && gevir_gene_pct < 25)
      paste0(" GeVIR percentile = ", round(gevir_gene_pct, 1),
             " indicates the gene is highly intolerant to missense variation in gnomAD.")
      else ""
    if (nchar(clingen_class) > 0 && clingen_class %in% c("Definitive", "Strong", "Moderate")) {
      dis_txt <- if (nchar(clingen_disease) > 0) paste0(" for ", clingen_disease) else ""
      moi_txt <- if (nchar(clingen_moi) > 0 && !grepl("^HP:", clingen_moi))
                   paste0(" (", clingen_moi, ")") else ""
      s(paste0("Missense variants are a common disease-causing mechanism for ", gene,
               dis_txt, moi_txt,
               ", supported by ClinGen Gene-Disease Validity classification: ", clingen_class,
               " (PP2).", gevir_pp2_str))
    } else {
      s(paste0("Missense changes are a common disease-causing mechanism for this gene (PP2).",
               gevir_pp2_str))
    }
  }

  # In silico predictions
  has_pp3 <- isTRUE(any(grepl("^PP3", tags)))
  has_bp4 <- isTRUE("BP4" %in% tags)
  if (has_pp3 || has_bp4) {
    scores <- character(0)
    if (!is.na(revel_score)) scores <- c(scores, paste0("REVEL: ", round(revel_score, 2)))
    if (!is.na(am_score))    scores <- c(scores, paste0("AlphaMissense: ", round(am_score, 3)))
    if (!is.na(cadd_score))  scores <- c(scores, paste0("CADD: ", round(cadd_score, 1)))
    if (nchar(metasvm_v) > 0 && metasvm_v != "NA")
      scores <- c(scores, paste0("MetaSVM: ", metasvm_v))
    score_str <- if (length(scores) > 0) paste0(" (", paste(scores, collapse = "; "), ")") else ""
    strength  <- if ("PP3_strong" %in% tags)   " at strong evidence strength" else if ("PP3_moderate" %in% tags) " at moderate evidence strength" else ""
    if (has_pp3)
      s(paste0("In silico tool predictions suggest damaging effect of the variant on the gene or gene product", strength, score_str, ".")) else
      s(paste0("In silico tool predictions suggest the variant is likely tolerated", score_str, " (BP4)."))
    # PS3 proxy note — when AM >= 0.90 + REVEL >= 0.773, note that PP3_strong is used as proxy
    if (isTRUE(ps3_proxy))
      s(paste0("Convergent structural (AlphaMissense) and ensemble (REVEL) evidence supports damaging effect; PP3_strong applied as proxy for PS3_supporting pending experimental functional data."))
    if (isTRUE(bs3_proxy))
      s(paste0("Convergent benign structural (AlphaMissense < 0.10) and ensemble (REVEL <= 0.29) evidence suggests tolerated effect; consistent with BS3_supporting pending DMS functional data."))
  }

  # Conservation standalone
  cs_gr <- suppressWarnings(as.integer(consurf_grade))
  if (!has_pp3 && length(cs_gr) == 1L && !is.na(cs_gr) && cs_gr >= 8L)
    s(paste0("The affected residue is highly conserved across species (ConSurf grade ", cs_gr, "/9)."))

  # PS1: same AA change
  if ("PS1" %in% tags) {
    cite     <- vcv_cite(clinvar_vcv)
    cite_str <- if (nchar(cite) > 0) paste0(" (", cite, ")") else ""
    s(paste0("Same nucleotide change resulting in same amino acid change has been previously reported to be associated with ", disease_str, cite_str, "."))
  }

  # PP5: ClinVar pathogenic
  if ("PP5" %in% tags && !"PS1" %in% tags) {
    cite     <- vcv_cite(clinvar_vcv)
    cite_str <- if (nchar(cite) > 0) paste0(" (", cite, ")") else ""
    s(paste0("This variant has been reported as pathogenic in ClinVar", cite_str, ", associated with ", disease_str, "."))
  }

  # PM5: different AA at same codon — MUST match variant's own position
  if ("PM5" %in% tags) {
    alt_changes <- character(0)
    if (nchar(clinvar_name) > 0 && !is.na(variant_pos)) {
      matches <- regmatches(clinvar_name,
                            gregexpr("p\\.([A-Z][a-z]{0,2}\\d+[A-Z][a-z]{0,2})", clinvar_name))[[1]]
      # Only keep entries whose numeric position matches this variant's position
      pos_str <- as.character(as.integer(variant_pos))
      matches <- matches[grepl(paste0("[A-Za-z]", pos_str, "[A-Za-z]"), matches)]
      mut_bare    <- sub("^p\\.", "", mut)
      alt_changes <- unique(matches[!grepl(mut_bare, matches, fixed = TRUE)])
    }
    # Use position-filtered VCVs (pre-filtered at data extraction time to match this codon only)
    cite     <- vcv_cite(clinvar_vcv_pos)
    alt_str  <- if (length(alt_changes) > 0) paste0(" (", paste(alt_changes, collapse = ", "), ")") else ""
    cite_str <- if (nchar(cite) > 0) paste0(" (", cite, ")") else ""
    s(paste0("Different missense changes at the same codon", alt_str,
             " have been reported to be associated with ", disease_str, cite_str, "."))
  }

  # PM4
  if ("PM4" %in% tags)
    s("The variant causes a protein length change (in-frame indel or stop-loss), which may disrupt protein function (PM4).")

  # PS2: confirmed de novo
  if ("PS2" %in% tags)
    s("The variant has been confirmed as de novo (maternity and paternity both verified) in an affected individual with no relevant family history (PS2).")

  # PM6: assumed de novo
  if ("PM6" %in% tags)
    s("The variant is assumed to be de novo in an affected individual; parental testing has not been performed to confirm absence in both parents (PM6).")

  # PS3_supporting: variant disrupts annotated PTM site
  if ("PS3_supporting" %in% tags)
    s("The variant position coincides with a UniProt-annotated post-translational modification site (phosphorylation, ubiquitination, disulfide bond, or equivalent). Disruption of such sites provides supporting functional evidence for pathogenicity (PS3_supporting).")

  # BP1
  if ("BP1" %in% tags) {
    gevir_bp1_str <- if (!is.na(gevir_gene_pct))
      paste0(" GeVIR percentile = ", round(gevir_gene_pct, 1),
             " (>75th percentile), indicating the gene accumulates missense variants freely in the gnomAD population.")
      else ""
    s(paste0("This is a missense variant in a gene where the primary disease mechanism is not missense variation, reducing the prior probability of pathogenicity (BP1).",
             gevir_bp1_str))
  }

  # BP3
  if ("BP3" %in% tags)
    s("The in-frame indel is located in a repetitive region with no known functional importance (BP3).")

  # BP6
  if ("BP6" %in% tags) {
    cite     <- vcv_cite(clinvar_vcv)
    cite_str <- if (nchar(cite) > 0) paste0(" (", cite, ")") else ""
    s(paste0("This variant has been classified as benign or likely benign by a reputable source", cite_str, " (BP6)."))
  }

  # BP7
  if ("BP7" %in% tags)
    s("The variant is synonymous and is not predicted to impact splicing (BP7).")

  # Classification — use shared classify_acmg() engine
  acmg_result    <- classify_acmg(tags)
  classification <- tolower(acmg_result$classification)
  rule_str       <- acmg_result$rule
  pts            <- acmg_result$pts

  s(paste0("Therefore, this variant is classified as ",
           classification,
           " (rule: ", rule_str, "; score: ", pts, " pts)",
           " according to the recommendation of the ACMG/AMP guideline."))

  paste(sentences, collapse = " ")
}
# ============================================================
# Build Variant Intersection Table
# For each user-input variant, look up data from all tracks
# ============================================================
build_variant_table <- function(highlight_df, af_data, mean_data, afs_data, gnomad_data, 
                                 clinvar_data, pfam_data, uniprot_data, ccrs_data,
                                 af_cutoff = NULL, ac_cutoff = NULL,
                                 clinvar_missense = NULL, consurf_data = NULL,
                                 denovo_status = "not_denovo",
                                 inh_param = "monoallelic",
                                 segregation = "none",
                                 # ── Analysis session parameters ──────────────────
                                 cutoff_method = "calc_af",
                                 prevalence_1_in_n = 2000,
                                 allelic_het = 0.5,
                                 genetic_het = 1.0,
                                 penetrance = 1.0,
                                 pop_size = 125748,
                                 conf_interval = 0.95,
                                 clingen_disease_param = "",
                                 clingen_moi_param = "",
                                 consurf_file_name = "") {
  if (is.null(highlight_df) || nrow(highlight_df) == 0) return(NULL)
  
  # Pre-compute density curves for gnomAD (blue) and ClinVar missense (red)
  # so we can look up values at each variant position
  kde_adjust <- 0.05
  prot_length <- 0
  gnomad_density_fn <- NULL
  clinvar_density_fn <- NULL
  
  # gnomAD density (blue): variants above AC cutoff
  if (!is.null(gnomad_data) && nrow(gnomad_data) > 0) {
    prot_length <- max(gnomad_data$prot_pos, na.rm = TRUE)
    if (!is.null(ac_cutoff)) {
      gd_above <- gnomad_data$prot_pos[gnomad_data$gnomad_allele_count > as.numeric(ac_cutoff)]
    } else {
      gd_above <- gnomad_data$prot_pos
    }
    gd_above <- gd_above[!is.na(gd_above)]
    if (length(gd_above) > 2) {
      gnomad_density_fn <- tryCatch(
        density(gd_above, adjust = kde_adjust, from = 0, to = prot_length),
        error = function(e) NULL
      )
    }
  }
  
  # ClinVar missense density (red)
  if (!is.null(clinvar_missense) && nrow(clinvar_missense) > 0) {
    cv_pos <- unique(clinvar_missense$prot_pos[!is.na(clinvar_missense$prot_pos)])
    if (length(cv_pos) > 2 && prot_length > 0) {
      clinvar_density_fn <- tryCatch(
        density(cv_pos, adjust = kde_adjust, from = 0, to = prot_length),
        error = function(e) NULL
      )
    }
  }
  
  # Helper: lookup density value at a position from a density object
  density_at_pos <- function(dens, pos) {
    if (is.null(dens)) return(NA_real_)
    idx <- which.min(abs(dens$x - pos))
    if (length(idx) > 0) return(dens$y[idx[1]])
    NA_real_
  }
  
  # Gene name for API calls — derived once at gene level from highlight_df
  gene_name_for_api <- if ("gene" %in% colnames(highlight_df) && nrow(highlight_df) > 0)
                         as.character(highlight_df$gene[1]) else ""

  # PP2: Gene-level flag — missense is a common disease mechanism if ≥3 P/LP in ClinVar
  # Computed once per gene, applied to all variants
  pp2_applies <- FALSE
  if (!is.null(clinvar_data) && nrow(clinvar_data) > 0 && "ClinicalSignificance" %in% colnames(clinvar_data)) {
    n_path_lp <- sum(grepl("pathogenic", clinvar_data$ClinicalSignificance, ignore.case = TRUE) &
                     !grepl("conflicting|uncertain|benign", clinvar_data$ClinicalSignificance, ignore.case = TRUE),
                     na.rm = TRUE)
    pp2_applies <- n_path_lp >= 3
    message("[PP2] Gene has ", n_path_lp, " P/LP ClinVar variants — PP2 applies: ", pp2_applies)
  }

  # ClinGen Gene-Disease Validity — fetched once per gene (cached)
  # Used to: (1) strengthen PP2 if ClinGen Definitive/Strong even when ClinVar count < 3
  #          (2) provide authoritative disease name + MOI in ACMG comment
  clingen_validity <- tryCatch(
    fetch_clingen_validity(gene_name_for_api),
    error = function(e) list(classification = NA_character_, moi = "", disease = "", source = "ClinGen LDH")
  )
  clingen_class <- if (!is.null(clingen_validity$classification) &&
                        !is.na(clingen_validity$classification))
                     clingen_validity$classification else ""

  # If ClinGen says Definitive or Strong, PP2 applies regardless of ClinVar count
  if (clingen_class %in% c("Definitive", "Strong") && !pp2_applies) {
    pp2_applies <- TRUE
    message("[PP2] ClinGen validity (", clingen_class, ") -> PP2 enabled for ", gene_name_for_api)
  }

  # GeVIR gene-level percentiles — looked up once per gene from pre-merged gene_data
  # Low GeVIR_pct = missense intolerant (strengthens PP2)
  # High GeVIR_pct = missense tolerant  (supports BP1 for missense variants)
  gevir_row <- gene_data[gene_data$gene_name == gene_name_for_api, , drop = FALSE]
  gevir_pct <- if (nrow(gevir_row) > 0 && !is.na(gevir_row$GeVIR_pct[1])) gevir_row$GeVIR_pct[1] else NA_real_
  message("[GeVIR] ", gene_name_for_api, " — GeVIR_pct=", if (!is.na(gevir_pct)) round(gevir_pct,1) else "NA")

  # GeVIR_pct < 25 — gene is highly intolerant to missense variation in gnomAD
  # This is population-genetics evidence independent of ClinVar, so enables PP2
  # even when ClinVar has fewer than 3 P/LP entries
  if (!is.na(gevir_pct) && gevir_pct < 25 && !pp2_applies) {
    pp2_applies <- TRUE
    message("[PP2] GeVIR_pct=", round(gevir_pct, 1), " < 25 -> PP2 enabled for ", gene_name_for_api)
  }

  # BP1 gene-level flag — fires for missense variants in genes tolerant to missense
  # GeVIR_pct > 75 means the gene accumulates missense variants freely in gnomAD,
  # consistent with missense not being the primary disease mechanism
  # Suppressed if ClinGen validity is Definitive/Strong (contradicts BP1 for missense genes)
  bp1_gene_applies <- !is.na(gevir_pct) &&
                      gevir_pct > 75 &&
                      !clingen_class %in% c("Definitive", "Strong")
  message("[BP1] GeVIR_pct=", if (!is.na(gevir_pct)) round(gevir_pct, 1) else "NA",
          " -> BP1 gene flag: ", bp1_gene_applies)

  rows <- lapply(seq_len(nrow(highlight_df)), function(i) {
    mut <- as.character(highlight_df$Mutation[i])
    pos <- highlight_df$prot_pos[i]
    if (is.na(pos)) return(NULL)
    
    # --- pLDDT ---
    plddt_val <- NA_real_
    plddt_cat <- ""
    if (!is.null(af_data) && nrow(af_data) > 0) {
      idx <- which(af_data$residueNumber == pos)
      if (length(idx) > 0) {
        plddt_val <- af_data$confidenceScore[idx[1]]
        plddt_cat <- as.character(af_data$confidenceCategory[idx[1]])
      }
    }
    
    # --- AF Mean Pathogenicity (position-level mean across all substitutions) ---
    af_mean_path <- NA_real_
    af_path_cat <- ""
    if (!is.null(mean_data) && nrow(mean_data) > 0) {
      idx <- which(mean_data$Position == pos)
      if (length(idx) > 0) {
        af_mean_path <- round(mean_data$mean_pathogenicity[idx[1]], 3)
        af_path_cat <- as.character(mean_data$category[idx[1]])
      }
    }

    # --- AlphaMissense exact variant score from AlphaFold substitution CSV ---
    # CSV protein_variant column is 1-letter bare: e.g. "P72R"
    # mut is HGVSp 3-letter: "p.Pro72Arg" — convert with aa3to1() then strip "p."
    am_exact_score <- NA_real_
    am_exact_class <- NA_character_
    if (!is.null(afs_data) && is.data.frame(afs_data) && nrow(afs_data) > 0 &&
        all(c("protein_variant", "am_pathogenicity") %in% colnames(afs_data))) {
      mut_1l <- tryCatch(sub("^p\\.", "", aa3to1(mut)), error = function(e) "")
      if (nchar(mut_1l) > 0) {
        am_idx <- which(afs_data$protein_variant == mut_1l)
        if (length(am_idx) > 0) {
          am_exact_score <- round(as.numeric(afs_data$am_pathogenicity[am_idx[1]]), 3)
          am_exact_class <- if (!is.na(am_exact_score)) {
            if      (am_exact_score >= 0.564) "pathogenic" else if (am_exact_score >= 0.340) "ambiguous" else                              "benign"
          } else NA_character_
        }
      }
    }
    
    # --- gnomAD: exact ref+alt match first, position fallback, absent = AF 0 ---
    gnomad_af      <- NA_real_
    gnomad_ac      <- NA_integer_
    gnomad_nhomalt <- NA_integer_   # homozygote count — used for BS2
    gnomad_match_type <- "absent"   # "exact", "position", or "absent"
    if (!is.null(gnomad_data) && nrow(gnomad_data) > 0) {
      # Parse ref and alt from user variant string (3-letter or 1-letter)
      mut_bare <- sub("^p\\.", "", as.character(mut))
      m3 <- regmatches(mut_bare,
                       regexec("^([A-Za-z]{3})([0-9]+)([A-Za-z]{3}|Ter|\\*)", mut_bare))[[1]]
      m1 <- regmatches(mut_bare,
                       regexec("^([A-Z])([0-9]+)([A-Z\\*])", mut_bare))[[1]]
      three_to_one_lk <- c(Ala="A",Arg="R",Asn="N",Asp="D",Cys="C",
                           Gln="Q",Glu="E",Gly="G",His="H",Ile="I",
                           Leu="L",Lys="K",Met="M",Phe="F",Pro="P",
                           Ser="S",Thr="T",Trp="W",Tyr="Y",Val="V",
                           Ter="*",Sec="U")
      if (length(m3) == 4) {
        user_ref <- three_to_one_lk[m3[2]]; if (is.na(user_ref)) user_ref <- m3[2]
        user_alt <- if (m3[4] %in% c("Ter","*")) "*" else three_to_one_lk[m3[4]]
        if (is.na(user_alt)) user_alt <- m3[4]
      } else if (length(m1) == 4) {
        user_ref <- m1[2]; user_alt <- m1[4]
      } else {
        user_ref <- NA_character_; user_alt <- NA_character_
      }

      # Step 1: exact match — same position AND same ref AND same alt
      idx_pos <- which(gnomad_data$prot_pos == pos)
      idx_exact <- integer(0)
      if (length(idx_pos) > 0 &&
          !is.na(user_ref) && !is.na(user_alt) &&
          "aa_ref" %in% colnames(gnomad_data) &&
          "aa_alt" %in% colnames(gnomad_data)) {
        idx_exact <- idx_pos[
          !is.na(gnomad_data$aa_ref[idx_pos]) &
          !is.na(gnomad_data$aa_alt[idx_pos]) &
          gnomad_data$aa_ref[idx_pos] == user_ref &
          gnomad_data$aa_alt[idx_pos] == user_alt
        ]
      }

      if (length(idx_exact) > 0) {
        # Exact match found — use it
        best <- idx_exact[which.max(gnomad_data$gnomad_allele_freq[idx_exact])]
        gnomad_af         <- signif(gnomad_data$gnomad_allele_freq[best], 4)
        gnomad_ac         <- gnomad_data$gnomad_allele_count[best]
        gnomad_match_type <- "exact"
        if ("exome_ac_hom" %in% colnames(gnomad_data))
          gnomad_nhomalt <- as.integer(gnomad_data$exome_ac_hom[best])
        message("[gnomAD] Exact match for ", mut, " at pos ", pos,
                ": AF=", gnomad_af, " AC=", gnomad_ac)

      } else if (length(idx_pos) > 0 && (is.na(user_ref) || is.na(user_alt) ||
                 !("aa_ref" %in% colnames(gnomad_data)))) {
        # No ref/alt info available — fall back to max-AF at position
        best <- idx_pos[which.max(gnomad_data$gnomad_allele_freq[idx_pos])]
        gnomad_af         <- signif(gnomad_data$gnomad_allele_freq[best], 4)
        gnomad_ac         <- gnomad_data$gnomad_allele_count[best]
        gnomad_match_type <- "position"
        if ("exome_ac_hom" %in% colnames(gnomad_data))
          gnomad_nhomalt <- as.integer(gnomad_data$exome_ac_hom[best])
        message("[gnomAD] Position-only match for ", mut, " at pos ", pos,
                ": AF=", gnomad_af, " AC=", gnomad_ac)

      } else {
        # No match at this position — variant is absent from gnomAD
        gnomad_af         <- 0
        gnomad_ac         <- 0L
        gnomad_nhomalt    <- 0L
        gnomad_match_type <- "absent"
        message("[gnomAD] No match for ", mut, " at pos ", pos,
                " — treating as absent (AF=0)")
      }
    }
    
    # --- ClinVar ---
    clinvar_sig <- ""
    clinvar_stars <- ""
    clinvar_name <- ""
    clinvar_vcv     <- ""   # all VCV accessions at this position
    clinvar_vcv_pos <- ""   # position-filtered VCVs (only those matching exact codon)
    clinvar_trait <- ""  # disease / trait name
    clinvar_match_type <- ""  # "exact", "position", or ""
    if (!is.null(clinvar_data) && nrow(clinvar_data) > 0) {
      # First: try exact amino acid match (same position AND same type AND same protein_change)
      # The user variant format is e.g. p.Cys29Trp — we need to match against ClinVar's protein_change or name
      exact_idx <- integer(0)
      if ("name" %in% colnames(clinvar_data)) {
        # Try matching the 3-letter form of the user's variant in the ClinVar name
        mut_1 <- aa3to1(mut)  # e.g. p.C29W
        mut_3 <- aa1to3(mut)  # e.g. p.Cys29Trp
        # Check if any ClinVar name contains the exact protein change
        for (form in c(mut, mut_3, mut_1)) {
          form_clean <- gsub("^p\\.", "", form)  # remove p. prefix for partial match
          name_match <- which(
            clinvar_data$prot_pos == pos & 
            (grepl(form, clinvar_data$name, fixed = TRUE) | 
             grepl(form_clean, clinvar_data$name, fixed = TRUE))
          )
          if (length(name_match) > 0) {
            exact_idx <- name_match
            break
          }
        }
      }
      
      if (length(exact_idx) > 0) {
        # Exact match found — user's variant IS in ClinVar
        clinvar_sig <- paste(unique(clinvar_data$ClinicalSignificance[exact_idx]), collapse = "; ")
        clinvar_stars <- paste(unique(clinvar_data$clinvar_goldstar[exact_idx]), collapse = ", ")
        clinvar_name <- paste(unique(clinvar_data$name[exact_idx]), collapse = "; ")
        clinvar_vcv  <- if ("vcv_id" %in% colnames(clinvar_data))
          paste(unique(clinvar_data$vcv_id[exact_idx]), collapse = ", ") else ""
        clinvar_vcv_pos <- clinvar_vcv   # exact match VCVs are always position-correct
        clinvar_trait <- if ("trait_name" %in% colnames(clinvar_data)) {
          tr <- unique(clinvar_data$trait_name[exact_idx])
          tr <- tr[nchar(tr) > 0]
          if (length(tr) > 0) tr[1] else ""
        } else ""
        clinvar_match_type <- "exact"
      } else {
        # No exact match — fall back to same-position variants (could be different AA change)
        pos_idx <- which(clinvar_data$prot_pos == pos)
        if (length(pos_idx) > 0) {
          # PM5 requires a structurally comparable variant at the same codon.
          # Frameshift, nonsense, splice, in-frame indel at the same position do NOT
          # constitute evidence that a missense at that codon is pathogenic — their
          # mechanism is entirely different (loss-of-function vs amino acid substitution).
          # Only missense variants are eligible for PM5 position-level matching.
          # Two-layer filter:
          #   (1) type field must be "missense_variant" — excludes frameshift, nonsense, splice
          #   (2) ClinVar name must contain a missense-style protein change pattern p.XxxNNNXxx
          #       This catches entries where type may be mislabelled or absent, and explicitly
          #       excludes names like "NM_xxx:c.1220delA (p.Pro407fs)" from triggering PM5.
          missense_by_type <- if ("type" %in% colnames(clinvar_data))
            pos_idx[clinvar_data$type[pos_idx] == "missense_variant"]
          else
            pos_idx  # if no type column, fall through to name-pattern filter only

          # Name-level filter: require p.Xxx123Xxx pattern (3-letter AA + position + 3-letter AA)
          # This catches frameshifts/LOF even when stored as "missense" by mistake
          missense_pattern <- "[Pp]\\.[A-Z][a-z]{2}\\d+(?!(?:Ter|Ext|fs|\\*|=))[A-Z][a-z]{2}(?![A-Za-z*\\d])"
          missense_by_name <- if ("name" %in% colnames(clinvar_data))
            missense_by_type[grepl(missense_pattern, clinvar_data$name[missense_by_type], perl = TRUE)]
          else
            missense_by_type

          missense_idx <- missense_by_name

          # If no clean missense ClinVar entry at this codon, skip position match entirely
          if (length(missense_idx) == 0) {
            pos_idx <- integer(0)   # zero out — prevents clinvar_match_type <- "position" below
          }
          use_idx <- missense_idx

          # Only proceed if there are actually missense variants at this codon
          if (length(use_idx) > 0) {
            clinvar_sig <- paste(unique(clinvar_data$ClinicalSignificance[use_idx]), collapse = "; ")
            clinvar_stars <- paste(unique(clinvar_data$clinvar_goldstar[use_idx]), collapse = ", ")
            clinvar_name <- paste(unique(clinvar_data$name[use_idx]), collapse = "; ")
            clinvar_vcv  <- if ("vcv_id" %in% colnames(clinvar_data))
              paste(unique(clinvar_data$vcv_id[use_idx]), collapse = ", ") else ""
            # Position-filtered VCV: only keep VCVs whose name contains this exact position number
            # (e.g. "150" matches p.Thr150Ala but not p.Arg282Trp)
            clinvar_vcv_pos <- if ("vcv_id" %in% colnames(clinvar_data) && "name" %in% colnames(clinvar_data)) {
              pos_str_filter <- as.character(pos)
              pos_filter_idx <- use_idx[grepl(paste0("[A-Za-z]", pos_str_filter, "[A-Za-z]"),
                                              clinvar_data$name[use_idx])]
              if (length(pos_filter_idx) > 0)
                paste(unique(clinvar_data$vcv_id[pos_filter_idx]), collapse = ", ") else ""
            } else ""
            clinvar_trait <- if ("trait_name" %in% colnames(clinvar_data)) {
              tr <- unique(clinvar_data$trait_name[use_idx])
              tr <- tr[nchar(tr) > 0]
              if (length(tr) > 0) tr[1] else ""
            } else ""
            clinvar_match_type <- "position"
          }
          # If use_idx is empty (no missense at this codon), all clinvar fields remain ""
          # and clinvar_match_type stays "" — PM5 will not fire
        }
      }
    }
    
    # --- Protein Domain ---
    domain <- ""
    # pfam_data is a complex UniProt JSON list; processed domains are stored
    # as pfam_data$Domain, pfam_data$Region, etc. (dataframes with start, end, description)
    if (!is.null(pfam_data)) {
      # Check all processed feature types stored in pfam_data
      feature_types <- c("Domain", "Region", "Repeat", "Zinc finger", "Motif", 
                          "Coiled coil", "Compositional bias")
      for (ft in feature_types) {
        ft_df <- pfam_data[[ft]]
        if (!is.null(ft_df) && is.data.frame(ft_df) && nrow(ft_df) > 0 &&
            "start" %in% colnames(ft_df) && "end" %in% colnames(ft_df)) {
          hits <- ft_df[!is.na(ft_df$start) & !is.na(ft_df$end) & 
                         ft_df$start <= pos & ft_df$end >= pos, , drop = FALSE]
          if (nrow(hits) > 0 && "description" %in% colnames(hits)) {
            descs <- unique(hits$description[!is.na(hits$description) & nchar(hits$description) > 0])
            if (length(descs) > 0) {
              domain <- paste(c(domain[nchar(domain) > 0], descs), collapse = "; ")
            }
          }
        }
      }
    }
    # Fallback: check uniprot_data for domain info
    if (nchar(domain) == 0 && !is.null(uniprot_data) && is.data.frame(uniprot_data) && nrow(uniprot_data) > 0) {
      dom_rows <- uniprot_data[uniprot_data$type == "domain" & 
                                 !is.na(uniprot_data$start) & !is.na(uniprot_data$end) &
                                 uniprot_data$start <= pos & uniprot_data$end >= pos, , drop = FALSE]
      if (nrow(dom_rows) > 0) {
        domain <- paste(unique(dom_rows$description), collapse = "; ")
      }
    }
    
    # Also check transmembrane and topological domains from uniprot_data
    if (!is.null(uniprot_data) && is.data.frame(uniprot_data) && nrow(uniprot_data) > 0 &&
        "type" %in% colnames(uniprot_data)) {
      # Transmembrane helix
      tm_rows <- uniprot_data[uniprot_data$type == "transmem" & 
                                !is.na(uniprot_data$start) & !is.na(uniprot_data$end) &
                                uniprot_data$start <= pos & uniprot_data$end >= pos, , drop = FALSE]
      if (nrow(tm_rows) > 0) {
        tm_desc <- if ("description" %in% colnames(tm_rows) && any(!is.na(tm_rows$description) & nchar(tm_rows$description) > 0)) {
          paste("TM:", unique(tm_rows$description[!is.na(tm_rows$description) & nchar(tm_rows$description) > 0]))
        } else "Transmembrane helix"
        domain <- paste(c(domain[nchar(domain) > 0], tm_desc), collapse = "; ")
      }
      # Topological domain (extracellular, cytoplasmic, etc.)
      td_rows <- uniprot_data[uniprot_data$type == "topo_dom" & 
                                !is.na(uniprot_data$start) & !is.na(uniprot_data$end) &
                                uniprot_data$start <= pos & uniprot_data$end >= pos, , drop = FALSE]
      if (nrow(td_rows) > 0) {
        td_desc <- if ("description" %in% colnames(td_rows) && any(!is.na(td_rows$description) & nchar(td_rows$description) > 0)) {
          paste(unique(td_rows$description[!is.na(td_rows$description) & nchar(td_rows$description) > 0]))
        } else "Topological domain"
        domain <- paste(c(domain[nchar(domain) > 0], td_desc), collapse = "; ")
      }
    }
    
    # --- CCRS (is position in a constrained region?) ---
    ccrs_info <- lookup_ccrs_at_position(highlight_df$gene[1], pos)
    in_ccrs <- if (ccrs_info$in_ccrs) {
      paste0("Yes (", ccrs_info$percentile, "th)")
    } else if (!is.na(ccrs_info$percentile)) {
      paste0("No (", ccrs_info$percentile, "th)")
    } else {
      "N/A"
    }
    # ScoreCons inter-species conservation from CCRStoAAC (bonus column)
    scorecons_val <- if (!is.na(ccrs_info$conservation)) round(ccrs_info$conservation, 3) else ""
    
    # --- PTM at position ---
    ptm_info     <- ""
    ptm_acmg     <- ""
    ptm_strength <- ""
    if (!is.null(uniprot_data) && is.data.frame(uniprot_data) && nrow(uniprot_data) > 0 &&
        "type" %in% colnames(uniprot_data) && "start" %in% colnames(uniprot_data)) {
      ptm_types <- c("mod_res", "lipid", "carbohyd", "crosslnk", "disulfid")
      ptm_rows <- uniprot_data[uniprot_data$type %in% ptm_types &
                                !is.na(uniprot_data$start) &
                                uniprot_data$start == pos, , drop = FALSE]
      if (nrow(ptm_rows) > 0) {
        if ("description" %in% colnames(ptm_rows)) {
          ptm_info <- paste(unique(na.omit(ptm_rows$description)), collapse = "; ")
        }
        if (nchar(ptm_info) == 0) {
          ptm_info <- paste(unique(ptm_rows$type), collapse = "; ")
        }
        # Classify PTM functional impact strength
        # Strong: phosphorylation, ubiquitination, acetylation, disulfide bond, cross-link
        # Moderate: lipidation, glycosylation, methylation, sumoylation
        desc_lc <- tolower(ptm_info)
        types_present <- unique(ptm_rows$type)
        is_strong_ptm <- any(types_present == "disulfid") ||
                         any(types_present == "crosslnk") ||
                         grepl("phospho|ubiquit|acetyl|sumoyl", desc_lc)
        is_mod_ptm    <- any(types_present %in% c("lipid","carbohyd","mod_res")) &&
                         grepl("methyl|glycosyl|lipid|palmitoyl|myristoyl|GPI", desc_lc)
        if (is_strong_ptm) {
          ptm_acmg     <- "PS3_supporting"
          ptm_strength <- "Strong functional site"
        } else if (nrow(ptm_rows) > 0) {
          ptm_acmg     <- "PP_PTM"
          ptm_strength <- "Moderate functional site"
        }
        message("[PTM] Position ", pos, ": ", ptm_info, " | ACMG proxy: ", ptm_acmg)
      }
    }
    
    # --- Density values at this position ---
    gnomad_dens_val <- density_at_pos(gnomad_density_fn, pos)
    clinvar_dens_val <- density_at_pos(clinvar_density_fn, pos)
    # Determine which is dominant (for coloring)
    dens_dominant <- "none"
    if (!is.na(gnomad_dens_val) && !is.na(clinvar_dens_val)) {
      dens_dominant <- if (clinvar_dens_val > gnomad_dens_val) "red" else "blue"
    } else if (!is.na(clinvar_dens_val)) {
      dens_dominant <- "red"
    } else if (!is.na(gnomad_dens_val)) {
      dens_dominant <- "blue"
    }
    
    # --- gnomAD filter status: is this position above or below the AC/AF cutoff? ---
    gnomad_filter <- "absent"
    if (!is.na(gnomad_af) && gnomad_match_type != "absent") {
      # Variant found in gnomAD — classify by AF/AC
      if (!is.null(af_cutoff) && !is.na(af_cutoff) && gnomad_af > as.numeric(af_cutoff)) {
        gnomad_filter <- "common"
      } else if (!is.null(ac_cutoff) && !is.na(gnomad_ac) && gnomad_ac > as.numeric(ac_cutoff)) {
        gnomad_filter <- "common"
      } else {
        gnomad_filter <- "rare"
      }
    } else if (gnomad_match_type == "absent") {
      # Variant genuinely absent from gnomAD — AF=0, AC=0
      gnomad_af <- 0
      gnomad_ac <- 0L
      gnomad_filter <- "absent"
    }
    
    # --- ConSurf conservation ---
    consurf_score_val <- ""
    consurf_grade <- ""
    consurf_burial <- ""
    if (!is.null(consurf_data) && is.data.frame(consurf_data) && nrow(consurf_data) > 0 &&
        "POS" %in% colnames(consurf_data)) {
      cs_idx <- which(consurf_data$POS == pos)
      if (length(cs_idx) > 0) {
        cs_row <- consurf_data[cs_idx[1], , drop = FALSE]
        if ("SCORE" %in% colnames(cs_row) && !is.na(cs_row$SCORE)) {
          consurf_score_val <- round(cs_row$SCORE, 3)
        }
        if ("COLOR" %in% colnames(cs_row) && !is.na(cs_row$COLOR)) {
          grade <- as.integer(cs_row$COLOR)
          consurf_grade <- grade
        }
        if ("BE" %in% colnames(cs_row) && !is.na(cs_row$BE) && nchar(cs_row$BE) > 0) {
          consurf_burial <- ifelse(cs_row$BE == "b", "Buried", 
                            ifelse(cs_row$BE == "e", "Exposed", cs_row$BE))
        }
      }
    }
    
    # --- dbNSFP via MyVariant.info API ---
    gene_name_for_api <- if ("gene" %in% colnames(highlight_df)) highlight_df$gene[1] else ""
    dbnsfp_raw <- tryCatch(
      fetch_dbnsfp(gene_name_for_api, mut),
      error = function(e) { message("[dbNSFP] Error for ", mut, ": ", e$message); NULL }
    )
    dbnsfp <- parse_dbnsfp_scores(dbnsfp_raw)

    # Extract rsID from raw MyVariant.info result for LitVar2 links
    rsid_val <- tryCatch({
      id_field <- dbnsfp_raw$`_id`
      if (!is.null(id_field) && grepl("^rs[0-9]+$", id_field)) {
        id_field
      } else {
        rs2 <- dbnsfp_raw$dbsnp$rsid
        if (!is.null(rs2) && grepl("^rs[0-9]+$", rs2)) rs2 else ""
      }
    }, error = function(e) "")

    # Ancestry AFs: fetched via rsID since dbNSFP object has no population fields
    if (dbnsfp$has_data) {
      pop_data <- tryCatch(fetch_pop_by_rsid(dbnsfp_raw), error = function(e) NULL)
      if (!is.null(pop_data)) {
        if (!is.na(pop_data$gnomAD_AFR)) dbnsfp$population$gnomAD_AFR <- pop_data$gnomAD_AFR
        if (!is.na(pop_data$gnomAD_NFE)) dbnsfp$population$gnomAD_NFE <- pop_data$gnomAD_NFE
        if (!is.na(pop_data$gnomAD_EAS)) dbnsfp$population$gnomAD_EAS <- pop_data$gnomAD_EAS
        if (!is.na(pop_data$gnomAD_SAS)) dbnsfp$population$gnomAD_SAS <- pop_data$gnomAD_SAS
        if (!is.na(pop_data$OneKG_AF))   dbnsfp$population$OneKG_AF   <- pop_data$OneKG_AF
      }
    }

    # Helper to extract score + verdict from a parsed section item
    sv <- function(section, tool) {
      if (!dbnsfp$has_data) return(list(s = "", v = ""))
      item <- dbnsfp[[section]][[tool]]
      if (is.null(item)) return(list(s = "", v = ""))
      sc <- item$score
      vd <- item$verdict
      list(
        s = if (!is.null(sc) && !is.na(sc)) as.character(round(as.numeric(sc), 3)) else "",
        v = if (!is.null(vd) && !is.na(vd)) as.character(vd) else ""
      )
    }
    # Helper for conservation (no verdict, just numeric)
    cv_num <- function(name) {
      if (!dbnsfp$has_data) return("")
      val <- dbnsfp$conservation[[name]]
      if (is.null(val) || is.na(val) || identical(val, "null")) return("")
      as.character(round(as.numeric(val), 3))
    }
    # Helper for population freq
    pf <- function(name) {
      if (!dbnsfp$has_data) return("")
      val <- dbnsfp$population[[name]]
      if (is.null(val) || is.na(val) || identical(val, "null")) return("")
      v <- as.numeric(val)
      if (is.na(v)) return("")
      if (v == 0) "0" else if (v < 0.0001) formatC(v, format = "e", digits = 2) else as.character(signif(v, 3))
    }
    
    # ── Group 1: Pathogenicity Predictors ──
    sift       <- sv("pathogenicity", "SIFT")
    pp2_hdiv   <- sv("pathogenicity", "PolyPhen2_HDIV")
    pp2_hvar   <- sv("pathogenicity", "PolyPhen2_HVAR")
    lrt        <- sv("pathogenicity", "LRT")
    fathmm     <- sv("pathogenicity", "FATHMM")
    provean    <- sv("pathogenicity", "PROVEAN")
    mt         <- sv("pathogenicity", "MutationTaster")
    
    # ── Group 2: Ensemble / Meta Scores ──
    revel      <- sv("ensemble", "REVEL")
    metasvm    <- sv("ensemble", "MetaSVM")
    metalr     <- sv("ensemble", "MetaLR")
    metarnn    <- sv("ensemble", "MetaRNN")
    cadd       <- sv("ensemble", "CADD")
    dann       <- sv("ensemble", "DANN")
    am         <- sv("ensemble", "AlphaMissense")
    dbscsnv_ada_sv <- sv("ensemble", "dbscSNV_ADA")
    dbscsnv_rf_sv  <- sv("ensemble", "dbscSNV_RF")
    
    # ── Group 3: Conservation ──
    gerp_rs    <- cv_num("GERP_RS")
    phylop100  <- cv_num("PhyloP_100V")
    phylop30   <- cv_num("PhyloP_470M")
    phastcons  <- cv_num("PhastCons_100V")
    
    # ── Group 4: Population Frequencies ──
    # Primary source: gnomad_data (already fetched for rainfall plot).
    # Ancestry AFs computed from populations{ac,an} in the GraphQL response.
    get_gnomad_pop <- function(col) {
      if (is.null(gnomad_data) || nrow(gnomad_data) == 0 ||
          !col %in% colnames(gnomad_data)) return(NA_real_)
      # Use the same exact-match index determined above
      if (gnomad_match_type == "absent") return(0)
      idx_use <- if (gnomad_match_type == "exact") idx_exact else idx_pos
      if (length(idx_use) == 0) return(NA_real_)
      best <- idx_use[which.max(gnomad_data$gnomad_allele_freq[idx_use])]
      v <- gnomad_data[[col]][best]
      # NA  = population data not available (old gnomAD data or missing field)
      # 0   = variant genotyped in this ancestry cohort but AC=0 (absent)
      # >0  = variant observed in this ancestry
      if (is.na(v)) NA_real_ else signif(v, 3)   # preserves 0 exactly
    }
    pop_gnomad <- if (!is.na(gnomad_af) && gnomad_af > 0) as.character(signif(gnomad_af, 3)) else if (!is.na(gnomad_af) && gnomad_af == 0) "0" else ""
    pop_afr    <- get_gnomad_pop("af_afr")
    pop_nfe    <- get_gnomad_pop("af_nfe")
    pop_eas    <- get_gnomad_pop("af_eas")
    pop_sas    <- get_gnomad_pop("af_sas")
    pop_fin    <- get_gnomad_pop("af_fin")
    pop_1kg    <- pf("OneKG_AF")
    pop_exac   <- pf("ExAC_AF")
    
    # ACMG tags — compute from BOTH VarViz's own API data AND dbNSFP scores
    # VarViz's own data is primary (already fetched from ClinVar, gnomAD, CCRS)
    # dbNSFP computational predictions supplement PP3/BP4
    acmg_tags <- character(0)
    
    # ── ACMG Tag Computation ──────────────────────────────────────────────────
    # Based on Richards et al. (2015) + Tavtigian et al. (2018) Bayesian framework
    # + Pejaver et al. (2022) calibrated in silico thresholds

    # PM1: Variant in mutational hotspot / critical functional domain
    # Upgrade to PM1_strong when conservation evidence also supports critical residue:
    #   - ConSurf grade 8-9 (highly/fully conserved), OR
    #   - PhyloP 100V > 4.0 (strong vertebrate conservation), OR
    #   - PhastCons 100V > 0.9 (high phastCons), OR
    #   - GERP RS > 4.4 (strong purifying selection)
    #
    # Only fire PM1 for STRUCTURAL domains — NOT for:
    #   - Intrinsically disordered / compositional bias regions
    #   - Interaction annotations ("Interaction with X")
    #   - Binding sites / active sites (handled separately)
    #   - Pro/Gly/Ala compositional bias regions
    domain_is_structural <- function(dom) {
      # PM1 should only fire for functionally specific, well-characterised domains.
      # Broad topological annotations (cytoplasmic loop, extracellular segment, TM region,
      # coiled-coil, signal peptide, propeptide, transit peptide, disordered regions) are
      # NOT sufficient — they span large stretches without per-residue functional constraint.
      # Whitelist approach: domain MUST match a specific functional term AND not be a
      # broad structural/topological label.
      if (nchar(dom) == 0) return(FALSE)
      parts <- trimws(strsplit(dom, ";")[[1]])

      functional_whitelist <- paste0(
        "active site|catalytic|ATP.bind|GTP.bind|nucleotide.bind|substrate.bind|",
        "ligand.bind|cofactor.bind|metal.bind|zinc.finger|iron-sulfur|heme.bind|",
        "FAD.bind|NAD.bind|SAM.bind|coenzyme|prosthetic|",
        "WD repeat|HEAT repeat|ARM repeat|TPR repeat|ANK repeat|kelch|",
        "RING|BRCT|SH2|SH3|PH domain|pleckstrin|RRM|KH domain|",
        "Tudor|chromo|bromo|PHD finger|PWWP|MBD|SET domain|",
        "kinase|phosphatase.active|ATPase|GTPase|protease.active|lipase.active|",
        "Walker A|Walker B|P-loop|activation loop|DFG motif|",
        "switch.I|switch.II|effector.loop|",
        "hotspot|mutational.cluster|critical.residue|invariant residue"
      )

      broad_blocklist <- paste0(
        "cytoplasmic|extracellular|transmembrane|signal peptide|propeptide|",
        "transit peptide|coiled.coil|disordered|low.complexity|",
        "compositional bias|Pro residues|Gly residues|Ala residues|Ser residues|",
        "basic residues|acidic residues|polar residues|",
        "^Region|loop|segment|stretch|linker|spacer|junction"
      )

      qualifies <- sapply(parts, function(p) {
        hits_wl <- grepl(functional_whitelist, p, ignore.case = TRUE, perl = TRUE)
        hits_bl <- grepl(broad_blocklist,      p, ignore.case = TRUE, perl = TRUE)
        hits_wl && !hits_bl
      })
      any(qualifies)
    }

    # Pre-compute benign AF flags — used by PM1, PP2, PM5, PS1, PP5, PP3 below
    # Thresholds are inheritance-aware (Roberts 2024 / Ware 2018 framework):
    #   Monoallelic (dominant):  BA1 > 5%, BS1 > 1%,  PM2 < 0.0001
    #   Biallelic  (recessive):  BA1 > 5%, BS1 > 5%,  PM2 < 0.01
    #                            (BS1 relaxed; high carrier freq expected for recessive)
    # De novo overrides to monoallelic regardless of inh_param
    inh_mode <- if (!is.null(denovo_status) && denovo_status %in% c("denovo_confirmed","denovo_assumed")) {
      "monoallelic"
    } else if (!is.null(inh_param) && inh_param == "biallelic") {
      "biallelic"
    } else {
      "monoallelic"
    }

    # Set thresholds based on inheritance mode
    # PM2 threshold: use user's prevalence-computed af_cutoff when available (from data() reactive)
    # This respects the user's prevalence/penetrance/hetA/hetG settings.
    # For biallelic: always use 0.01 (recessive carrier frequency standard).
    # For monoallelic: use af_cutoff from calc_af method if set, else fall back to 0.0001 (ACMG default).
    ba1_thresh <- 0.05
    bs1_thresh <- if (inh_mode == "biallelic") 0.05 else 0.01
    pm2_thresh <- if (inh_mode == "biallelic") {
      0.01
    } else if (!is.null(af_cutoff) && !is.na(af_cutoff) && as.numeric(af_cutoff) > 0) {
      as.numeric(af_cutoff)   # user's prevalence-computed threshold (e.g. 4e-05)
    } else {
      0.0001                  # ACMG default fallback
    }

    ba1_fires <- !is.na(gnomad_af) && gnomad_af > ba1_thresh
    bs1_fires <- !is.na(gnomad_af) && gnomad_af > bs1_thresh
    bs2_fires <- ba1_fires ||
                 (!is.na(gnomad_nhomalt) && gnomad_nhomalt > 0)

    # PM1: Mutational hotspot / functional constraint
    # Three-pathway approach addressing both over- and under-firing:
    #
    # Path 1 — CCRS >= 90th percentile (Samocha 2017): empirically depleted missense region.
    #   Sufficient alone for PM1; no specific domain annotation required.
    #   PM1_strong upgrade if also cons_strong (convergent empirical + phylogenetic evidence).
    #
    # Path 2 — CCRS 85–89th + cons_strong: borderline depletion zone. Requires conservation
    #   to corroborate, preventing over-calling at permissive percentiles.
    #   Catches residues like His133 (88th), Thr145, Met147 that have real constraint but
    #   fall below the hard 90th cutoff. PM1 only (not _strong) at this tier.
    #
    # Path 3 — Specific functional domain (whitelist) without CCRS support.
    #   Requires cons_strong to fire. PM1_strong if domain + cons_strong, else PM1.
    #
    # Note: broad topological labels (cytoplasmic loop, extracellular segment) are excluded
    # from Path 3 by domain_is_structural() whitelist. Paths 1 and 2 bypass domain entirely
    # since CCRS is annotation-agnostic.
    if (!bs1_fires) {
      cs_grade_num  <- { cg <- suppressWarnings(as.integer(consurf_grade));  if (length(cg)==0||is.na(cg[1])) NA_integer_ else cg[1] }
      phylop_v      <- { pv <- if (dbnsfp$has_data) dbnsfp$conservation$PhyloP_100V    else NA_real_; if (is.null(pv)||length(pv)==0) NA_real_ else suppressWarnings(as.numeric(pv[1])) }
      phastcons_v   <- { pv <- if (dbnsfp$has_data) dbnsfp$conservation$PhastCons_100V else NA_real_; if (is.null(pv)||length(pv)==0) NA_real_ else suppressWarnings(as.numeric(pv[1])) }
      gerp_v        <- { pv <- if (dbnsfp$has_data) dbnsfp$conservation$GERP_RS        else NA_real_; if (is.null(pv)||length(pv)==0) NA_real_ else suppressWarnings(as.numeric(pv[1])) }
      scorecons_v   <- { sv <- suppressWarnings(as.numeric(scorecons_val));  if (length(sv)==0||is.na(sv[1])) NA_real_ else sv[1] }

      cons_strong <- isTRUE(!is.na(cs_grade_num)  && cs_grade_num  >= 8)   ||
                     isTRUE(!is.na(phylop_v)       && phylop_v       > 4.0)  ||
                     isTRUE(!is.na(phastcons_v)    && phastcons_v    > 0.9)  ||
                     isTRUE(!is.na(gerp_v)         && gerp_v         > 4.4)  ||
                     isTRUE(!is.na(scorecons_v)    && scorecons_v    > 0.8)

      has_domain  <- domain_is_structural(domain)
      ccrs_pct    <- if (!is.null(ccrs_info$percentile) && !is.na(ccrs_info$percentile) &&
                         ccrs_info$percentile > 0) ccrs_info$percentile else 0

      cons_used_for_pm1 <- FALSE  # tracks whether conservation already spent on PM1_strong
      if (ccrs_pct >= 90) {
        # Path 1: high confidence CCRS — PM1 unconditionally, PM1_strong with conservation
        if (cons_strong) {
          acmg_tags <- c(acmg_tags, "PM1_strong")
          cons_used_for_pm1 <- TRUE
        } else {
          acmg_tags <- c(acmg_tags, "PM1")
        }
      } else if (ccrs_pct >= 85 && cons_strong) {
        # Path 2: borderline CCRS — conservation required as second line of evidence
        acmg_tags <- c(acmg_tags, "PM1")
      } else if (has_domain && cons_strong) {
        # Path 3a: specific domain, conserved — PM1_strong
        # Flag: conservation was used here; suppress duplicate PP3 conservation tier
        acmg_tags <- c(acmg_tags, "PM1_strong")
        cons_used_for_pm1 <- TRUE
      } else if (has_domain) {
        # Path 3b: specific domain, not strongly conserved — PM1 supporting
        acmg_tags <- c(acmg_tags, "PM1")
      }
    }

    # PM2: Absent or extremely low frequency in population databases
    # Threshold: monoallelic < 0.0001, biallelic < 0.01 (Roberts 2024 / Ware 2018)
    # PM1 15aa neighborhood: independent of PM2 frequency gate.
    # A variant can be in a ClinVar hotspot regardless of population frequency.
    if (!bs1_fires && !is.null(clinvar_data) && nrow(clinvar_data) > 0 &&
        "prot_pos" %in% colnames(clinvar_data) && !is.na(pos)) {
      win_idx <- which(abs(clinvar_data$prot_pos - pos) <= 15)
      if (length(win_idx) > 0) {
        win_sig <- clinvar_data$ClinicalSignificance[win_idx]
        n_path_win <- sum(grepl("pathogenic", win_sig, ignore.case=TRUE) &
                         !grepl("conflicting|uncertain|benign", win_sig, ignore.case=TRUE), na.rm=TRUE)
        n_benign_win <- sum(grepl("benign", win_sig, ignore.case=TRUE) &
                           !grepl("pathogenic|conflicting", win_sig, ignore.case=TRUE), na.rm=TRUE)
        ratio_ok <- n_benign_win == 0 || (n_path_win / (n_path_win + n_benign_win)) >= 0.75
        if (n_path_win >= 3 && ratio_ok) {
          if ("PM1_strong" %in% acmg_tags) {
            # already strong from CCRS/domain+cons — neighborhood corroborates, no change
          } else if ("PM1" %in% acmg_tags) {
            # domain + neighborhood = upgrade to PM1_strong
            acmg_tags <- acmg_tags[acmg_tags != "PM1"]
            acmg_tags <- c(acmg_tags, "PM1_strong")
            message("[PM1] Neighborhood \u00b115aa upgrade to PM1_strong: ", n_path_win, " P/LP, ", n_benign_win, " B/LB")
          } else {
            # no domain/CCRS hit but neighborhood is strong — fire PM1
            acmg_tags <- c(acmg_tags, "PM1")
            message("[PM1] Neighborhood \u00b115aa fires PM1: ", n_path_win, " P/LP, ", n_benign_win, " B/LB")
          }
        }
      }
    }

    # PM2: Absent or extremely low frequency
    if (is.na(gnomad_af) || gnomad_af < pm2_thresh) {
      acmg_tags <- c(acmg_tags, "PM2")
    }

    # PS2 / PM6: De novo variant evidence
    # PS2: confirmed de novo (paternity verified) — Strong Pathogenic (+4 pts)
    # PM6: assumed de novo (paternity not confirmed) — Moderate Pathogenic (+2 pts)
    # Both fire only when variant is absent/ultra-rare (not BA1/BS1)
    if (!ba1_fires && !bs1_fires) {
      if (!is.null(denovo_status) && denovo_status == "denovo_confirmed") {
        acmg_tags <- c(acmg_tags, "PS2")
      } else if (!is.null(denovo_status) && denovo_status == "denovo_assumed") {
        acmg_tags <- c(acmg_tags, "PM6")
      }
    }

    # PP1: Cosegregation with disease in affected family members
    if (!is.null(segregation) && segregation != "none" && !ba1_fires && !bs1_fires) {
      pp1_tag <- switch(segregation,
        pp1          = "PP1",
        pp1_moderate = "PP1_moderate",
        pp1_strong   = "PP1_strong",
        NULL
      )
      if (!is.null(pp1_tag)) acmg_tags <- c(acmg_tags, pp1_tag)
    }

    # BS1/BS2 pre-check: common variants (AF > 1%) cannot also carry pathogenic moderate evidence
    # Compute benign AF flags early so downstream criteria can use them
    # PM4: Protein length change (in-frame indel or stop-loss)
    # Detect from variant name: contains "del", "ins", "dup", "fs", or "ext"
    mut_str <- tolower(as.character(mut))
    if (grepl("del|ins|dup|ext|\\*[0-9]", mut_str) && !grepl("fs|frameshift", mut_str)) {
      acmg_tags <- c(acmg_tags, "PM4")
    }

    # PP2: Missense in gene where missense is common disease mechanism (gene-level, pre-computed)
    # Suppressed when variant is common (BS1/BS2) — a 71% AF variant cannot use PP2
    if (pp2_applies && !bs1_fires) {
      acmg_tags <- c(acmg_tags, "PP2")
    }

    # ── PS1 vs PM5 disambiguation ────────────────────────────────────────────
    # PS1: SAME amino acid change as an established P/LP variant in ClinVar.
    #   - Requires clinvar_match_type == "exact" (same protein change matched)
    #   - Requires "Pathogenic" (not just LP) for full PS1; LP ClinVar → downgrade to PM5
    #   - Cannot co-exist with PM5 (mutually exclusive by definition)
    #   - Note: nucleotide-level distinction not available from protein data;
    #     clinicians should confirm different nucleotide context if relevant.
    #
    # PM5: DIFFERENT amino acid change at the SAME codon as an established P variant.
    #   - Requires clinvar_match_type == "position" (different AA at same codon)
    #   - Or: exact match but ClinVar is only LP (not full Pathogenic)
    #   - Cannot fire if PS1 fires

    ps1_fired <- FALSE
    pm5_fired <- FALSE

    if (!bs1_fires && nchar(clinvar_sig) > 0 &&
        grepl("pathogenic", clinvar_sig, ignore.case = TRUE) &&
        !grepl("conflicting|uncertain", clinvar_sig, ignore.case = TRUE)) {

      is_path_only <- grepl("^[Pp]athogenic", clinvar_sig) &&
                      !grepl("likely", clinvar_sig, ignore.case = TRUE) &&
                      !grepl("benign", clinvar_sig, ignore.case = TRUE)
      is_lp        <- grepl("likely pathogenic", clinvar_sig, ignore.case = TRUE) &&
                      !grepl("benign", clinvar_sig, ignore.case = TRUE)

      if (clinvar_match_type == "exact") {
        if (is_path_only) {
          # Same AA change, established Pathogenic → PS1 star-weighted
          stars_n <- suppressWarnings(as.integer(clinvar_stars))
          ps1_tag <- if (!is.na(stars_n) && stars_n >= 2) "PS1" else
                     if (!is.na(stars_n) && stars_n == 1) "PS1_moderate" else "PS1_supporting"
          acmg_tags <- c(acmg_tags, ps1_tag)
          ps1_fired <- TRUE
        } else if (is_lp) {
          # Same AA change but only Likely Pathogenic in ClinVar → downgrade to PM5
          # (LP is not "established" per Richards 2015 PS1 criteria)
          acmg_tags <- c(acmg_tags, "PM5")
          pm5_fired <- TRUE
        }
      } else if (clinvar_match_type == "position" && !ps1_fired) {
        # Different AA change at same codon, established P or LP → PM5
        # (PS1 requires identical AA change, so this cannot be PS1)
        if (is_path_only || is_lp) {
          acmg_tags <- c(acmg_tags, "PM5")
          pm5_fired <- TRUE
        }
      }
    }

    # PP5: Reputable source (ClinVar) reports pathogenic — exact match
    # Suppressed when PS1 fires: both use the same ClinVar assertion (double-dipping).
    # PP5 fires only when PS1 has NOT fired — i.e. when ClinVar has a pathogenic
    # assertion but it is not the exact same AA change (PS1 requires exact match).
    # In practice PP5 fires alongside PM5 (position match) or for multi-star LP variants.
    if (!bs1_fires && !ps1_fired &&
        clinvar_match_type == "exact" &&
        nchar(clinvar_sig) > 0 && grepl("pathogenic", clinvar_sig, ignore.case = TRUE) &&
        !grepl("conflicting|uncertain|benign", clinvar_sig, ignore.case = TRUE)) {
      acmg_tags <- c(acmg_tags, "PP5")
    }

    # PP3 / BP4: Computational evidence — Pejaver 2022 calibrated thresholds
    # Priority order: REVEL (best calibrated) → MetaSVM/MetaLR/MetaRNN → CADD → DANN → vote
    # Conservation scores (GERP, PhyloP, PhastCons, ConSurf, ScoreCons) add independent PP3

    # safe_sc/safe_pd: collapse any length-0/NULL/list to a single NA scalar
    # This is the root fix for "argument is of length zero" from character(0) pred fields
    safe_sc <- function(x) {
      if (is.null(x) || length(x) == 0) return(NA_real_)
      x <- suppressWarnings(as.numeric(x[[1]]))
      if (length(x) == 0 || is.nan(x)) NA_real_ else x
    }
    safe_pd <- function(x) {
      if (is.null(x) || length(x) == 0) return(NA_character_)
      x <- as.character(x[[1]])
      if (length(x) == 0 || is.na(x)) NA_character_ else x
    }

    revel_sc   <- safe_sc(if (dbnsfp$has_data) dbnsfp$ensemble$REVEL$score             else NULL)
    metasvm_pd <- safe_pd(if (dbnsfp$has_data) dbnsfp$ensemble$MetaSVM$pred            else NULL)
    metalr_pd  <- safe_pd(if (dbnsfp$has_data) dbnsfp$ensemble$MetaLR$pred             else NULL)
    metarnn_pd <- safe_pd(if (dbnsfp$has_data) dbnsfp$ensemble$MetaRNN$pred            else NULL)
    cadd_sc    <- safe_sc(if (dbnsfp$has_data) dbnsfp$ensemble$CADD$score              else NULL)
    dann_sc    <- safe_sc(if (dbnsfp$has_data) dbnsfp$ensemble$DANN$score              else NULL)
    # AlphaMissense score: prefer MyVariant (dbnsfp), fall back to AlphaFold CSV exact match
    am_sc <- {
      raw <- safe_sc(if (dbnsfp$has_data) dbnsfp$pathogenicity$AlphaMissense$score else NULL)
      if (!is.na(raw)) raw else am_exact_score   # am_exact_score initialised to NA_real_ above
    }

    # Helper: pick highest PP3 level seen so far
    pp3_level <- function(tags) {
      if (isTRUE(any(tags == "PP3_strong")))   return(3L)
      if (isTRUE(any(tags == "PP3_moderate"))) return(2L)
      if (isTRUE(any(tags == "PP3")))          return(1L)
      return(0L)
    }
    add_pp3 <- function(level) {
      tag <- c("1" = "PP3", "2" = "PP3_moderate", "3" = "PP3_strong")[level]
      acmg_tags <<- unique(c(acmg_tags[!grepl("^PP3", acmg_tags)], tag))
    }

    # 1. REVEL — best single calibrated predictor (Pejaver 2022 Table 2)
    if (!is.na(revel_sc)) {
      if      (revel_sc >= 0.932) add_pp3("3") else if (revel_sc >= 0.773) add_pp3("2") else if (revel_sc >= 0.644) add_pp3("1")
    }

    # 2. MetaSVM / MetaLR / MetaRNN — consensus meta-predictors
    # Per Pejaver 2022: 2/3 meta-predictors = moderate ONLY when a calibrated
    # predictor (REVEL ≥ 0.644 or CADD ≥ 28.1) corroborates.
    # Without calibrated corroboration, cap at PP3 supporting regardless of meta count.
    meta_dam <- sum(c(
      isTRUE(!is.na(metasvm_pd) && grepl("^[Dd]", metasvm_pd)),
      isTRUE(!is.na(metalr_pd)  && grepl("^[Dd]", metalr_pd)),
      isTRUE(!is.na(metarnn_pd) && grepl("^[Dd]", metarnn_pd))
    ))
    calibrated_sc <- isTRUE((!is.na(revel_sc) && revel_sc >= 0.644) ||
                             (!is.na(cadd_sc)  && cadd_sc  >= 28.1))
    if      (meta_dam >= 2 && calibrated_sc  && pp3_level(acmg_tags) < 2L) add_pp3("2") else if (meta_dam >= 1                   && pp3_level(acmg_tags) < 1L) add_pp3("1")

    # 3. CADD (Pejaver 2022: ≥ 28.1 = supporting; ≥ 35 = moderate)
    if (!is.na(cadd_sc)) {
      if      (cadd_sc >= 35   && pp3_level(acmg_tags) < 2L) add_pp3("2") else if (cadd_sc >= 28.1 && pp3_level(acmg_tags) < 1L) add_pp3("1")
    }

    # 4. DANN (≥ 0.96 = supporting)
    if (!is.na(dann_sc) && dann_sc >= 0.96 && pp3_level(acmg_tags) < 1L) add_pp3("1")

    # 5. AlphaMissense (≥ 0.564 = supporting)
    if (!is.na(am_sc) && am_sc >= 0.564 && pp3_level(acmg_tags) < 1L) add_pp3("1")

    # ── PS3 / BS3: Functional evidence ───────────────────────────────────────
    # PS3: well-established functional studies show damaging effect (+4 pts, Strong Pathogenic)
    # BS3: well-established functional studies show no damaging effect (-4 pts, Strong Benign)
    #
    # Ideal: direct DMS (deep mutational scanning) functional scores for the gene under study.
    # For PPP2R5C and related B56 subunit genes, saturation mutagenesis data would provide At near-saturation DMS coverage, each
    # variant has an experimentally measured functional score that maps directly to PS3/BS3.
    # VarViz does not yet ingest DMS data tables, so these criteria cannot fire from
    # experimental data. This is the core limitation identified for PPP2R5C analysis.
    #
    # Proxy tier (interim, lower confidence): AlphaMissense structural predictions.
    # AlphaMissense pathogenicity > 0.90 with REVEL >= 0.773 provides convergent structural +
    # evolutionary evidence approaching functional-study-level confidence for some variants,
    # but this is NOT equivalent to experimentally measured DMS scores.
    # PS3_supporting (not full PS3) is awarded; add a note in the comment.
    # BS3_supporting: AM < 0.10 with REVEL <= 0.290 — benign structural + computational.
    #
    # When DMS data are integrated, replace this block with direct Pd-score thresholds:
    #   PS3       : Pd >= 0.90  (damaging, high confidence)
    #   PS3_moderate: Pd >= 0.75
    #   BS3       : Pd <= 0.10  (functional, high confidence)
    #   BS3_supporting: Pd <= 0.25
    ps3_proxy_fired <- FALSE
    bs3_proxy_fired <- FALSE
    if (!bs1_fires && !ba1_fires) {
      am_high <- isTRUE(!is.na(am_sc) && am_sc >= 0.90)
      am_low  <- isTRUE(!is.na(am_sc) && am_sc <= 0.10)
      revel_dam  <- isTRUE(!is.na(revel_sc) && revel_sc >= 0.773)
      revel_ben  <- isTRUE(!is.na(revel_sc) && revel_sc <= 0.290)
      if (am_high && revel_dam) {
        # Convergent structural + ensemble evidence — PS3_supporting proxy
        # Do not add a new tag; instead upgrade PP3_moderate -> PP3_strong as proxy for PS3_supporting
        # (avoids tag proliferation while reflecting the stronger signal)
        if (pp3_level(acmg_tags) < 3L) add_pp3("3")
        ps3_proxy_fired <- TRUE
      } else if (am_low && revel_ben) {
        # Convergent benign structural + ensemble evidence — BS3_supporting proxy
        # Flag for comment generation; does not add a new ACMG tag at this time
        bs3_proxy_fired <- TRUE
      }
    }

    # ── Pre-compute BP4 votes BEFORE conservation gate ───────────────────────
    # BP4 votes must be tallied here so the conservation gate can check whether
    # benign predictors are contradicting pathogenic conservation evidence.
    # Actual BP4 tag assignment happens after conservation (see below).
    bp4_votes <- sum(c(
      isTRUE(!is.na(revel_sc)   && revel_sc   <= 0.290),
      isTRUE(!is.na(metasvm_pd) && grepl("^[Tt]", metasvm_pd)),
      isTRUE(!is.na(metalr_pd)  && grepl("^[Tt]", metalr_pd)),
      isTRUE(!is.na(metarnn_pd) && grepl("^[Tt]", metarnn_pd)),
      isTRUE(!is.na(cadd_sc)    && cadd_sc    < 15),
      isTRUE(!is.na(dann_sc)    && dann_sc    < 0.5),
      isTRUE(!is.na(am_sc)      && am_sc      <= 0.34)
    ))
    benign_predictors_active <- bp4_votes >= 2   # would BP4 fire?

    # 6. Conservation scores — independent evidence stream
    # Per Pejaver 2022 + ClinGen SVI guidance:
    #
    # Upgrade PP3→PP3_moderate:
    #   Requires ≥2 conservation tools AND calibrated seq predictor (REVEL ≥ 0.644
    #   or CADD ≥ 28.1) already fired PP3. Conservation corroborates, does NOT lead.
    #
    # Add PP3 supporting:
    #   Requires ≥1 conservation hit AND benign predictors are NOT contradicting
    #   (i.e. BP4 would not fire). If predictors call benign, conservation is
    #   overridden — conflicting computational signals cancel out.
    #
    # PP3_strong via conservation: NOT allowed (requires REVEL ≥ 0.932 alone).
    cs_grade_cons <- { cg <- suppressWarnings(as.integer(consurf_grade)); if (length(cg)==0) NA_integer_ else cg[1] }
    phylop_pp3    <- safe_sc(if (dbnsfp$has_data) dbnsfp$conservation$PhyloP_100V    else NA_real_)
    phastcons_pp3 <- safe_sc(if (dbnsfp$has_data) dbnsfp$conservation$PhastCons_100V else NA_real_)
    gerp_pp3      <- safe_sc(if (dbnsfp$has_data) dbnsfp$conservation$GERP_RS        else NA_real_)
    scorecons_pp3 <- { sv <- suppressWarnings(as.numeric(scorecons_val)); if (length(sv)==0) NA_real_ else sv[1] }

    n_cons_hits <- sum(c(
      isTRUE(!is.na(cs_grade_cons)  && cs_grade_cons  >= 8),
      isTRUE(!is.na(phylop_pp3)     && phylop_pp3      > 4.0),
      isTRUE(!is.na(phastcons_pp3)  && phastcons_pp3   > 0.9),
      isTRUE(!is.na(gerp_pp3)       && gerp_pp3        > 4.4),
      isTRUE(!is.na(scorecons_pp3)  && scorecons_pp3   > 0.8)
    ))

    calibrated_seq_fired <- calibrated_sc
    if (n_cons_hits >= 2 && calibrated_seq_fired && pp3_level(acmg_tags) < 2L &&
        !cons_used_for_pm1) {
      # ≥2 conservation + calibrated predictor → upgrade to moderate
      # Suppressed if conservation already spent on PM1_strong (double-dipping)
      add_pp3("2")
    } else if (n_cons_hits >= 1 && !benign_predictors_active && pp3_level(acmg_tags) < 1L &&
               !cons_used_for_pm1) {
      # ≥1 conservation, benign predictors silent → add supporting PP3
      # Suppressed if same conservation signal already contributed to PM1_strong
      add_pp3("1")
    }
    # If n_cons_hits >= 1 BUT benign_predictors_active: conflicting computational
    # signals — conservation says conserved, predictors say benign. Neither PP3
    # nor BP4 fires from conservation in this case. The predictor-based BP4 will
    # be handled below (PP3 vs BP4 mutual exclusion step).

    # BP1: Missense in gene where ONLY truncating variants cause disease
    # Cannot auto-assess without gene-level LOF mechanism knowledge — not implemented

    # BP3: In-frame indel in repeat region (no known function)
    if (grepl("del|ins|dup", mut_str) && !grepl("fs|frameshift", mut_str)) {
      in_repeat <- FALSE
      if (!is.null(uniprot_data) && is.data.frame(uniprot_data) && nrow(uniprot_data) > 0 &&
          "type" %in% colnames(uniprot_data)) {
        rep_rows <- uniprot_data[uniprot_data$type == "repeat" &
                                   !is.na(uniprot_data$start) & !is.na(uniprot_data$end) &
                                   uniprot_data$start <= pos & uniprot_data$end >= pos, , drop = FALSE]
        in_repeat <- nrow(rep_rows) > 0
      }
      if (in_repeat) acmg_tags <- c(acmg_tags, "BP3")
    }

    # ── PP3 vs BP4 mutual exclusion (InterVar conflict resolution) ────────────
    # PP3 and BP4 represent opposing computational verdicts and CANNOT co-exist.
    # Precedence rule (strongest evidence wins):
    #
    #   PP3_strong  fired → PP3 wins unconditionally; BP4 suppressed
    #   PP3_moderate fired → PP3 wins; BP4 suppressed
    #   PP3         fired → PP3 wins; BP4 suppressed
    #   PP3 not fired AND bp4_votes ≥ 2 → BP4 added
    #   PP3 not fired AND bp4_votes < 2 → neither fires
    #
    # This prevents the misleading tag list PP3+BP4 with net 0 pts.
    pp3_active <- isTRUE(any(grepl("^PP3", acmg_tags)))
    if (!pp3_active && bp4_votes >= 2) {
      acmg_tags <- c(acmg_tags, "BP4")
    }
    # pp3_active + bp4_votes ≥ 2 → PP3 wins, BP4 silently discarded

    # BA1: gnomAD AF > 5% — standalone Benign (Richards 2015 ACMG rule 2a)
    # This criterion alone is sufficient for Benign classification.
    # Remove ALL pathogenic tags; BA1 takes precedence.
    if (ba1_fires) {
      acmg_tags <- c(acmg_tags, "BA1")
      acmg_tags <- acmg_tags[!acmg_tags %in% c(
        "PM1", "PM1_strong", "PM2", "PM4", "PM5",
        "PP2", "PP3", "PP3_moderate", "PP3_strong", "PS1", "PS1_moderate", "PS1_supporting", "PP5"
      )]
    }

    # BS1: Allele frequency greater than expected for disorder (gnomAD AF 1–5%)
    # Note: variants with AF > 5% are already BA1; BS1 is the 1–5% range.
    # Also remove pathogenic moderate/supporting tags — a common variant cannot be PM/PP.
    if (bs1_fires && !ba1_fires) {
      acmg_tags <- c(acmg_tags, "BS1")
      acmg_tags <- acmg_tags[!acmg_tags %in% c("PM1", "PM1_strong", "PP2", "PP3", "PP3_moderate", "PP3_strong", "PM5", "PS1", "PS1_moderate", "PS1_supporting", "PP5")]
    }

    # BS2: Observed in unaffected adults (gnomAD homozygotes > 0, or AF > 5% as proxy)
    # Independent of BS1/BA1 — evaluates carrier observation in healthy population.
    if (bs2_fires) {
      acmg_tags <- c(acmg_tags, "BS2")
    }

    # BP6: Reputable source (ClinVar) reports benign — exact match
    # Mutual exclusion with PS1/PP5: if ClinVar has conflicting submissions,
    # neither PS1/PP5 nor BP6 should fire (already blocked upstream by
    # grepl("conflicting") guards, but make explicit here).
    # BP6 also cannot co-exist with PS1 — if PS1 fired, ClinVar said Pathogenic.
    if (clinvar_match_type == "exact" &&
        nchar(clinvar_sig) > 0 &&
        grepl("benign", clinvar_sig, ignore.case = TRUE) &&
        !grepl("pathogenic|conflicting", clinvar_sig, ignore.case = TRUE) &&
        !ps1_fired) {   # PS1 and BP6 are mutually exclusive
      acmg_tags <- c(acmg_tags, "BP6")
    }

    # BP7: Synonymous variant with no predicted splice impact
    # Detect from variant name: same ref and alt AA (e.g. p.Pro72Pro)
    syn_match <- regmatches(mut, regexec("^p\\.([A-Z][a-z]{0,2})(\\d+)([A-Z][a-z]{0,2})$", mut))[[1]]
    if (length(syn_match) == 4) {
      ref_aa_bp7 <- toupper(substr(syn_match[2], 1, 1))
      alt_aa_bp7 <- toupper(substr(syn_match[4], 1, 1))
      if (ref_aa_bp7 == alt_aa_bp7) acmg_tags <- c(acmg_tags, "BP7")
    }

    # BP1: Missense variant in gene where missense is NOT a common disease mechanism
    # Triggered by GeVIR_pct > 75 — gene accumulates missense variants freely in gnomAD
    # population, indicating missense variation is generally tolerated.
    # Suppressed when:
    #   - BA1/BS1 already fire (frequency evidence supersedes mechanism evidence)
    #   - PS1 fired (ClinVar confirms pathogenic at this exact position)
    #   - PP5 fired (ClinVar reports pathogenic — contradicts BP1)
    #   - ClinGen validity is Definitive/Strong (missense IS established mechanism)
    if (bp1_gene_applies &&
        !ba1_fires && !bs1_fires &&
        !any(c("PS1","PS1_moderate","PS1_supporting") %in% acmg_tags) &&
        !"PP5" %in% acmg_tags) {
      acmg_tags <- c(acmg_tags, "BP1")
      message("[BP1] Fired for ", mut, " — GeVIR_pct=", round(gevir_pct, 1))
    }

    acmg_tags <- unique(acmg_tags)
    acmg_str <- if (length(acmg_tags) > 0) paste(acmg_tags, collapse = ", ") else ""
    
    # Serialize the full dbNSFP scores as JSON for the detail card rendering
    dbnsfp_json <- if (dbnsfp$has_data) {
      tryCatch(jsonlite::toJSON(dbnsfp, auto_unbox = TRUE, null = "null", na = "null"),
               error = function(e) "")
    } else ""
    
    # ── Generate ACMG clinical comment ────────────────────────────────────────
    acmg_comment <- generate_acmg_comment(
      acmg_tags_str   = acmg_str,
      gene             = highlight_df$gene[1],
      mut              = mut,
      gnomad_af        = gnomad_af,
      gnomad_ac        = gnomad_ac,
      gnomad_nhomalt   = gnomad_nhomalt,
      clinvar_sig      = clinvar_sig,
      clinvar_name     = clinvar_name,
      clinvar_vcv      = clinvar_vcv,
      clinvar_vcv_pos  = clinvar_vcv_pos,
      clinvar_trait    = clinvar_trait,
      clinvar_match_type = clinvar_match_type,
      domain           = domain,
      pp2_applies      = pp2_applies,
      revel_score      = if (dbnsfp$has_data) dbnsfp$ensemble$REVEL$score else NA,
      am_score         = if (dbnsfp$has_data) dbnsfp$pathogenicity$AlphaMissense$score else NA,
      cadd_score       = if (dbnsfp$has_data) dbnsfp$ensemble$CADD$score else NA,
      metasvm_v        = if (dbnsfp$has_data) dbnsfp$ensemble$MetaSVM$pred else "",
      phylop_val       = if (dbnsfp$has_data) dbnsfp$conservation$PhyloP_100V else NA,
      phastcons_val    = if (dbnsfp$has_data) dbnsfp$conservation$PhastCons_100V else NA,
      gerp_val         = if (dbnsfp$has_data) dbnsfp$conservation$GERP_RS else NA,
      consurf_grade    = consurf_grade,
      variant_pos      = pos,
      clingen_class    = clingen_class,
      clingen_disease  = if (!is.null(clingen_validity$disease))  clingen_validity$disease  else "",
      clingen_moi      = if (!is.null(clingen_validity$moi))      clingen_validity$moi      else "",
      ps3_proxy        = ps3_proxy_fired,
      bs3_proxy        = bs3_proxy_fired,
      gevir_gene_pct   = gevir_pct
    )

    data.frame(
      Variant = mut,
      Position = pos,
      GeVIR_Gene_Pct = if (!is.na(gevir_pct)) round(gevir_pct, 1) else "",
      pLDDT = ifelse(is.na(plddt_val), "", round(plddt_val, 1)),
      pLDDT_Category = plddt_cat,
      AF_Mean_Pathogenicity = ifelse(is.na(af_mean_path), "", af_mean_path),
      AF_Path_Category = af_path_cat,
      gnomAD_AF = ifelse(is.na(gnomad_af), "", gnomad_af),
      gnomAD_AC = ifelse(is.na(gnomad_ac), "", gnomad_ac),
      gnomAD_Match = if (exists("gnomad_match_type")) gnomad_match_type else "",
      gnomAD_Filter = gnomad_filter,
      gnomAD_Nhomalt = ifelse(is.na(gnomad_nhomalt), "", as.integer(gnomad_nhomalt)),
      Density_gnomAD = ifelse(is.na(gnomad_dens_val), "", signif(gnomad_dens_val, 4)),
      Density_ClinVar = ifelse(is.na(clinvar_dens_val), "", signif(clinvar_dens_val, 4)),
      Density_Dominant = dens_dominant,
      ClinVar = clinvar_sig,
      ClinVar_Stars = clinvar_stars,
      ClinVar_Name = clinvar_name,
      ClinVar_VCV = clinvar_vcv,
      ClinVar_VCV_Pos = clinvar_vcv_pos,
      RSID = rsid_val,
      ClinVar_Trait = clinvar_trait,
      ClinVar_Match = clinvar_match_type,
      Domain = domain,
      CCRS = in_ccrs,
      ScoreCons = scorecons_val,
      PTM = ptm_info,
      PTM_ACMG = ptm_acmg,
      PTM_Strength = ptm_strength,
      ConSurf_Score = consurf_score_val,
      ConSurf_Grade = consurf_grade,
      ConSurf_Burial = consurf_burial,
      # ── Group 1: Pathogenicity ──
      SIFT_Score = sift$s, SIFT_V = sift$v,
      PP2_HDIV_Score = pp2_hdiv$s, PP2_HDIV_V = pp2_hdiv$v,
      PP2_HVAR_Score = pp2_hvar$s, PP2_HVAR_V = pp2_hvar$v,
      LRT_Score = lrt$s, LRT_V = lrt$v,
      FATHMM_Score = fathmm$s, FATHMM_V = fathmm$v,
      PROVEAN_Score = provean$s, PROVEAN_V = provean$v,
      MutTaster_Score = mt$s, MutTaster_V = mt$v,
      # ── Group 2: Ensemble ──
      REVEL_Score = revel$s, REVEL_V = revel$v,
      MetaSVM_Score = metasvm$s, MetaSVM_V = metasvm$v,
      MetaLR_Score = metalr$s, MetaLR_V = metalr$v,
      MetaRNN_Score = metarnn$s, MetaRNN_V = metarnn$v,
      CADD_Score = cadd$s, CADD_V = cadd$v,
      DANN_Score = dann$s, DANN_V = dann$v,
      AM_Score  = if (!is.na(am_exact_score)) am_exact_score else am$s,
      AM_V      = if (!is.na(am_exact_class)) am_exact_class  else am$v,
      AM_Source = if (!is.na(am_exact_score)) "AlphaFold" else if (!is.na(am$s)) "MyVariant" else NA_character_,
      dbscSNV_ADA = dbscsnv_ada_sv$s,
      dbscSNV_RF  = dbscsnv_rf_sv$s,
      # ── Group 3: Conservation ──
      GERP_RS = gerp_rs,
      PhyloP_100V = phylop100,
      PhyloP_470M = phylop30,
      PhastCons = phastcons,
      # ── Group 4: Population ──
      Pop_gnomAD = pop_gnomad,
      Pop_AFR = pop_afr,
      Pop_NFE = pop_nfe,
      Pop_EAS = pop_eas,
      Pop_SAS = pop_sas,
      Pop_FIN = pop_fin,
      Pop_1KG = pop_1kg,
      Pop_ExAC = pop_exac,
      # ── ACMG + Comment + JSON ──
      ACMG_Tags = acmg_str,
      Comment = acmg_comment,
      dbNSFP_JSON = as.character(dbnsfp_json),
      # ── Session analysis parameters (same for all variants in run) ──
      Analysis_Inheritance = inh_param,
      Analysis_Cutoff_Method = cutoff_method,
      Analysis_PM2_Threshold = signif(pm2_thresh, 4),
      Analysis_BS1_Threshold = bs1_thresh,
      Analysis_AF_Cutoff = ifelse(!is.null(af_cutoff) && !is.na(af_cutoff), signif(as.numeric(af_cutoff), 4), ""),
      Analysis_AC_Cutoff = ifelse(!is.null(ac_cutoff) && !is.na(ac_cutoff), as.integer(ac_cutoff), ""),
      Analysis_Prevalence = paste0("1 in ", prevalence_1_in_n),
      Analysis_Allelic_Het = allelic_het,
      Analysis_Genetic_Het = genetic_het,
      Analysis_Penetrance = penetrance,
      Analysis_Pop_Size = pop_size,
      Analysis_CI = conf_interval,
      ClinGen_Disease = clingen_disease_param,
      ClinGen_MOI = clingen_moi_param,
      ClinGen_Class = clingen_class,
      ConSurf_File = consurf_file_name,
      stringsAsFactors = FALSE
    )
  })
  
  rows <- rows[!sapply(rows, is.null)]
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}


shinyServer(function(input, output, session) {
  
  # Serve static files (help.html) from app directory
  addResourcePath("static", getwd())
  session$onFlushed(function() {
    # Load choices server-side (fast — just the gene name list).
    # Do NOT pre-select a default gene: selecting "TP53" here immediately triggers
    # all downstream eventReactives (gnomAD, ClinVar, AlphaFold CSV download, etc.)
    # which makes every session cold-start slow regardless of what the user wants.
    # The user picks a gene, then presses Go — nothing heavy runs before that.
    updateSelectizeInput(session, "gene_name", choices = gene_list,
                         selected = "CASR", server = TRUE,
                         options = list(placeholder = "Type or select a gene..."))
  }, once = TRUE)
  
  # ── Variant file upload -> populate text area ──────────────────────────
  observeEvent(input$variant_file, {
    req(input$variant_file)
    tryCatch({
      lines <- readLines(input$variant_file$datapath)
      lines <- trimws(lines[nchar(trimws(lines)) > 0])
      all_vars <- unlist(strsplit(paste(lines, collapse = ","), "[,\n\r]+"))
      all_vars <- trimws(all_vars[nchar(trimws(all_vars)) > 0])
      all_vars <- ifelse(grepl("^p\\.", all_vars, ignore.case = TRUE),
                         all_vars, paste0("p.", all_vars))
      updateTextAreaInput(session, "variants", value = paste(all_vars, collapse = ", "))
      message("[VariantFile] Loaded ", length(all_vars), " variants from: ", input$variant_file$name)
    }, error = function(e) {
      message("[VariantFile] Error: ", e$message)
    })
  })
  
  observeEvent(input$launchApp,  { updateTabsetPanel(session, "nav", selected = "Protein View") })
  observeEvent(input$launchApp2, { updateTabsetPanel(session, "nav", selected = "Protein View") })
  
  # Advanced section now integrated into Calculate AF / AC Filter panels
  
    # consurf_state: tracks fetch status for UI status panel
    # Values: "idle" | "loading" | list(status="ok", source, pdb_id, chain, resolution)
    #       | list(status="absent", uid) | list(status="file_ok", rows)
    #       | list(status="file_error", msg)
    consurf_state <- reactiveVal("idle")

    # parse_consurf_file: shared parser for both file upload and auto-fetch
    # Handles real ConSurf grades files which use TAB separators with a blank
    # column between SCORE and COLOR, and variable column positions for B/E.
    # Strategy: find the header row, map column NAMES to indices, then extract
    # by name rather than fixed position so blank/extra columns don't break it.
  parse_consurf_file <- function(lines) {
    # Find header line — matches " POS", "#POS", "# POS", "\tPOS" etc.
    # POSIX [:space:] works in base R grep without perl=TRUE
    header_index <- grep("^[[:space:]#]*POS\\b", lines)[1]
    if (is.na(header_index)) stop("Header line with POS column not found.")

    # Strip leading whitespace/# from header, split on TAB
    header_raw  <- gsub("^[#[:space:]]+", "", lines[header_index])
    hdr_cols    <- trimws(strsplit(header_raw, "\\t")[[1]])

    # Locate the columns we need by name
    pos_col   <- which(hdr_cols == "POS")[1]
    seq_col   <- which(hdr_cols == "SEQ")[1]
    score_col <- which(hdr_cols == "SCORE")[1]
    color_col <- which(hdr_cols == "COLOR")[1]
    be_col    <- which(hdr_cols %in% c("B/E", "BE"))[1]
    if (any(is.na(c(pos_col, seq_col, score_col, color_col))))
      stop("Required columns (POS, SEQ, SCORE, COLOR) not all found in header.")

    # Data lines: skip the units row (contains "normalized") and blank/comment lines
    data_lines <- lines[(header_index + 1):length(lines)]
    data_lines <- data_lines[!grepl("normalized|^[[:space:]]*#|^[[:space:]]*$", data_lines)]
    if (length(data_lines) == 0) stop("No data lines found after header.")

    # Parse each data row by TAB — extract by column index
    safe_col <- function(row, idx) {
      if (is.na(idx) || length(row) < idx) NA_character_ else row[idx]
    }
    rows <- lapply(data_lines, function(line) {
      row <- strsplit(line, "\\t")[[1]]
      c(safe_col(row, pos_col),
        safe_col(row, seq_col),
        safe_col(row, score_col),
        safe_col(row, color_col),
        safe_col(row, be_col))
    })
    rows <- rows[sapply(rows, function(r) !is.na(r[1]) && nchar(trimws(r[1])) > 0)]
    if (length(rows) == 0) stop("No valid data rows extracted.")

    df <- as.data.frame(do.call(rbind, rows), stringsAsFactors = FALSE)
    colnames(df) <- c("POS", "SEQ", "SCORE", "COLOR", "BE")
    df$POS   <- suppressWarnings(as.integer(gsub("[^0-9]",     "", trimws(df$POS))))
    df$SEQ   <- gsub("[^A-Za-z]",                                  "", trimws(df$SEQ))
    df$SCORE <- suppressWarnings(as.numeric(gsub("[^0-9.\\-]", "", trimws(df$SCORE))))
    df$COLOR <- suppressWarnings(as.numeric(gsub("[^0-9.]",       "", trimws(df$COLOR))))
    df$BE    <- tolower(gsub("[^A-Za-z]",                          "", trimws(df$BE)))
    df[!is.na(df$POS) & !is.na(df$COLOR), ]
  }

  consurf_score <- reactive({
    # Depend on both Go button and file upload so uploading a ConSurf file
    # after Go is clicked immediately updates the plot without re-clicking Go
    if (input$goButton == 0) return(data.frame())
    input$consurf_file  # reactive dependency — invalidates when file changes
    consurf_state("idle")

    # Priority 1: user-uploaded file (explicit user choice always wins)
    if (!is.null(input$consurf_file)) {
      message("[ConSurf] Reading uploaded file: ", input$consurf_file$name)
      lines <- tryCatch(readLines(input$consurf_file$datapath),
                        error = function(e) { message("[ConSurf] File read error: ", e$message); NULL })
      if (!is.null(lines)) {
        df <- tryCatch(parse_consurf_file(lines),
                       error = function(e) { message("[ConSurf] Parse error: ", e$message); NULL })
        if (!is.null(df) && nrow(df) > 0) {
          message("[ConSurf] Loaded ", nrow(df), " residues from uploaded file")
          consurf_state(list(status = "file_ok", rows = nrow(df),
                             fname = input$consurf_file$name))
          return(df)
        }
      }
      consurf_state(list(status = "file_error",
                         msg = paste0("Could not parse '", input$consurf_file$name, "'")))
      return(data.frame())
    }

    # Priority 2: auto-fetch disabled — ConSurf track shows only when file is uploaded
    return(data.frame())
  })
  
  # ── Multi-Conservation (UCSC) status panel ──
  output$multicons_status <- renderUI({
    state <- multicons_state()
    if (identical(state, "idle")) return(NULL)

    box_style <- function(bg, border)
      paste0("margin:4px 0 8px 0; padding:8px 10px; background:", bg,
             "; border-left:3px solid ", border, "; border-radius:4px; font-size:12px;")

    if (identical(state, "loading")) {
      return(tags$div(style = box_style("#f0f9ff", "#3b82f6"),
        tags$span(style = "color:#2563eb;",
          HTML('<svg xmlns="http://www.w3.org/2000/svg" width="13" height="13"
                     viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"
                     style="animation:spin 1s linear infinite; vertical-align:middle;
                     margin-right:5px;"><path d="M21 12a9 9 0 1 1-6.219-8.56"/></svg>'),
          "Fetching PhyloP / PhastCons / GERP++ from UCSC..."
        ),
        tags$style("@keyframes spin{from{transform:rotate(0deg)}to{transform:rotate(360deg)}}")
      ))
    }

    if (is.list(state) && identical(state$status, "ok")) {
      pct    <- round(100 * state$covered / max(state$total, 1))

      return(tags$div(style = box_style("#f0fdf4", "#22c55e"),
        tags$span(style = "color:#15803d; font-weight:600;",
                  "✓ UCSC conservation loaded"),
        tags$br(),
        tags$span(style = "color:#4b5563;",
          paste0(state$covered, "/", state$total, " residues (", pct,
                 "%) • PhyloP 100V • PhyloP 30M • PhastCons"))
      ))
    }

    if (is.list(state) && identical(state$status, "absent")) {
      return(tags$div(style = box_style("#fef2f2", "#ef4444"),
        tags$span(style = "color:#b91c1c; font-weight:600;",
                  "⚠ No conservation data returned"),
        tags$br(),
        tags$span(style = "color:#4b5563;",
          paste0("UCSC BigWig tracks returned no scores for ", state$gene, ".")),
        tags$br(),
        tags$span(style = "color:#6b7280; font-size:11px;",
          "Gene may be novel, on an untracked chromosome, or Ensembl coordinates unavailable.")
      ))
    }

    if (is.list(state) && identical(state$status, "error")) {
      return(tags$div(style = box_style("#fef2f2", "#ef4444"),
        tags$span(style = "color:#b91c1c; font-weight:600;", "✘ Fetch failed"),
        tags$br(),
        tags$span(style = "color:#4b5563;", state$msg)
      ))
    }

    NULL
  })


  # ── ConSurf status panel: shown below the file input in the sidebar ──
  output$consurf_status <- renderUI({
    state <- consurf_state()

    # Not yet run or track not selected — show nothing
    if (identical(state, "idle")) return(NULL)

    # Loading spinner while fetch is in progress
    if (identical(state, "loading")) {
      return(tags$div(
        style = "margin: 4px 0 8px 0; padding: 8px 10px; background:#f0f9ff;
                 border-left: 3px solid #3b82f6; border-radius: 4px; font-size:12px;",
        tags$span(style = "color:#2563eb;",
          HTML('<svg xmlns="http://www.w3.org/2000/svg" width="13" height="13"
                     viewBox="0 0 24 24" fill="none" stroke="currentColor"
                     stroke-width="2" class="spin" style="animation:spin 1s linear infinite;
                     vertical-align:middle; margin-right:5px;">
                 <path d="M21 12a9 9 0 1 1-6.219-8.56"/>
               </svg>'),
          "Searching ConSurf-DB..."
        ),
        tags$style("@keyframes spin { from{transform:rotate(0deg)} to{transform:rotate(360deg)} }")
      ))
    }

    # Auto-fetch success
    if (is.list(state) && identical(state$status, "ok")) {
      return(tags$div(
        style = "margin: 4px 0 8px 0; padding: 8px 10px; background:#f0fdf4;
                 border-left: 3px solid #22c55e; border-radius: 4px; font-size:12px;",
        tags$span(style = "color:#15803d; font-weight:600;", "✓ ConSurf-DB loaded"),
        tags$br(),
        tags$span(style = "color:#4b5563;",
          paste0(state$rows, " residues • ", state$uid))
      ))
    }

    # File upload success
    if (is.list(state) && identical(state$status, "file_ok")) {
      return(tags$div(
        style = "margin: 4px 0 8px 0; padding: 8px 10px; background:#f0fdf4;
                 border-left: 3px solid #22c55e; border-radius: 4px; font-size:12px;",
        tags$span(style = "color:#15803d; font-weight:600;", "✓ File loaded"),
        tags$br(),
        tags$span(style = "color:#4b5563;",
          paste0(state$rows, " residues • ", state$fname))
      ))
    }

    # Auto-fetch absent — no PDB / not in ConSurf-DB
    if (is.list(state) && identical(state$status, "absent")) {
      return(tags$div(
        style = "margin: 4px 0 8px 0; padding: 8px 10px; background:#fef2f2;
                 border-left: 3px solid #ef4444; border-radius: 4px; font-size:12px;",
        tags$span(style = "color:#b91c1c; font-weight:600;",
          paste0("⚠ Not in ConSurf-DB")),
        tags$br(),
        tags$span(style = "color:#4b5563;",
          paste0(state$gene, " (", state$uid, ") has no PDB structure or ",
                 "ConSurf-DB entry.")),
        tags$br(),
        tags$a(
          href = paste0("https://consurf.tau.ac.il/?pdb_id=ALPHAFOLD&uniprot_id=",
                        state$uid),
          target = "_blank",
          style = "color:#2563eb; font-size:11px;",
          "Run ConSurf on AlphaFold model ↗"
        ),
        tags$span(style = "color:#6b7280; font-size:11px;",
          " then upload grades file above.")
      ))
    }

    # File parse error
    if (is.list(state) && identical(state$status, "file_error")) {
      return(tags$div(
        style = "margin: 4px 0 8px 0; padding: 8px 10px; background:#fef2f2;
                 border-left: 3px solid #ef4444; border-radius: 4px; font-size:12px;",
        tags$span(style = "color:#b91c1c; font-weight:600;", "✘ Parse failed"),
        tags$br(),
        tags$span(style = "color:#4b5563;", state$msg),
        tags$br(),
        tags$span(style = "color:#6b7280; font-size:11px;",
          "Expected: ConSurf grades file (.txt) in tab-separated format.")
      ))
    }

    NULL
  })


  user_path_variants <- eventReactive(input$goButton, {
    # Read custom pathogenic positions file.
    # Accepts two formats:
    #   (a) p.notation variant names (p.Arg175His, one per line or TSV)
    #       -> numeric position extracted via extract_protein_position()
    #   (b) raw numeric positions (one per line or TSV V1 column)
    # Returns data.frame with column V1 containing numeric positions.
    if (!is.null(input$user_file)) {
      tryCatch({
        raw_lines <- readLines(input$user_file$datapath)
        raw_lines <- trimws(raw_lines[nchar(trimws(raw_lines)) > 0])
        # Detect p.notation format
        if (any(grepl("^p\\.", raw_lines, ignore.case = TRUE))) {
          # Extract numeric positions from p.notation strings
          positions <- sapply(raw_lines, extract_protein_position)
          positions <- as.numeric(positions[!is.na(positions)])
          message("[UserFile] Parsed ", length(positions), " p.notation variants -> positions")
        } else {
          # Try as TSV with numeric positions in first column
          df <- tryCatch(
            data.table::fread(input$user_file$datapath, sep = "\t", header = FALSE, quote = ""),
            error = function(e) NULL
          )
          positions <- if (!is.null(df) && nrow(df) > 0)
            suppressWarnings(as.numeric(df$V1))
          else as.numeric(raw_lines)
          positions <- positions[!is.na(positions)]
          message("[UserFile] Parsed ", length(positions), " numeric positions from TSV")
        }
        if (length(positions) == 0) return(data.frame())
        data.frame(V1 = positions)
      }, error = function(e) {
        message("[UserFile] Parse error: ", e$message)
        data.frame()
      })
    } else {
      data.frame()
    }
  })

  # Reset uploaded files when gene changes — ignoreInit=TRUE so startup
  # auto-selection of CASR does not wipe files uploaded before clicking Go
  observeEvent(input$gene_name, {
    req(!is.null(input$gene_name) && nchar(trimws(input$gene_name)) > 0)
    shinyjs::reset("consurf_file")
  }, ignoreInit = TRUE)
  observeEvent(input$gene_name, {
    req(!is.null(input$gene_name) && nchar(trimws(input$gene_name)) > 0)
    shinyjs::reset("user_file")
  }, ignoreInit = TRUE)
  
  variants <- eventReactive(input$goButton,{
    raw <- strsplit(input$variants, "[,]")[[1]]
    raw <- trimws(raw)
    raw <- raw[nchar(raw) > 0]
    # Normalize: add p. prefix if missing
    raw <- ifelse(grepl("^p\\.", raw, ignore.case = TRUE), raw,
                  paste0("p.", raw))
    as.data.frame(raw, stringsAsFactors = FALSE) |>
      setNames("x")
  })
  highlight = eventReactive(input$goButton,{ 
    highlight=data.frame(variants(), sapply(variants(), extract_protein_position), input$gene_name) 
    colnames(highlight) <- c("Mutation","prot_pos","gene")
    highlight
  })
  
  # ============================================================
  # METHOD 1: Calculate AF (by prevalence / genetic architecture)
  # ============================================================
  maxAF_calc <- reactive({
    myPrev = 1/input$prev
    if(input$inh=="monoallelic"){
      myMaxAF = (1/2) * myPrev * input$hetA * input$hetG * (1/input$pen)
    }
    if(input$inh=="biallelic"){
      myMaxAF = sqrt(myPrev) * input$hetA * sqrt(input$hetG) * (1/sqrt(input$pen))
    }
    myMaxAF
  })
  
  maxAC_from_AF <- reactive({
    qpois(p = as.numeric(input$CI_af), lambda = input$popSize_af * maxAF_calc())
  })
  
  # ============================================================
  # METHOD 2: AC Filter (direct max population AF input)
  # ============================================================
  maxAC_calc <- reactive({
    qpois(p = as.numeric(input$CI), lambda = input$popSize * input$maxPopAF)
  })
  
  # ============================================================
  # Compact summary (always visible below dropdown)
  # ============================================================
  output$cutoff_summary <- renderText({
    if (input$cutoff_method == "calc_af") {
      af <- maxAF_calc()
      ac <- maxAC_from_AF()
      paste0(
        '<div style="background:#fff; border:1px solid #e2e8f0; border-radius:8px; padding:10px; margin:6px 0; display:flex; justify-content:space-around; text-align:center;">',
        '<div><span style="font-size:11px; color:#475569;">Max AF</span><br>',
        '<span style="font-size:18px; font-weight:700; color:#ef4444;">', signif(af, 3), '</span></div>',
        '<div><span style="font-size:11px; color:#475569;">Max AC</span><br>',
        '<span style="font-size:18px; font-weight:700; color:#ef4444;">', ac, '</span></div>',
        '</div>')
    } else {
      ac <- maxAC_calc()
      paste0(
        '<div style="background:#fff; border:1px solid #e2e8f0; border-radius:8px; padding:10px; margin:6px 0; text-align:center;">',
        '<span style="font-size:11px; color:#475569;">Max tolerated AC</span><br>',
        '<span style="font-size:18px; font-weight:700; color:#ef4444;">', ac, '</span>',
        '</div>')
    }
  })
  
  # Toggle cutoff details panel
  shinyjs::onclick("toggleCutoff", {
    shinyjs::toggle(id = "cutoff_details", anim = TRUE)
    shinyjs::runjs('
      var lnk = document.getElementById("toggleCutoff");
      var spans = lnk.querySelectorAll("span:not(.fa)");
      var txt = lnk.textContent.trim();
      if (txt.indexOf("Hide") !== -1) {
        lnk.childNodes[lnk.childNodes.length-1].textContent = " Edit parameters";
      } else {
        lnk.childNodes[lnk.childNodes.length-1].textContent = " Hide parameters";
      }
    ')
  })
  
  # ============================================================
  # UNIFIED: data() reactive routes to the selected method for plotting
  # ============================================================
  data <- reactive({
    if (input$cutoff_method == "calc_af") {
      list(maxAF_calc(), maxAC_from_AF())
    } else {
      list(input$maxPopAF, maxAC_calc())
    }
  }) 
  
  gene_attrib <- eventReactive(input$goButton,{ gene_data[gene_data$gene_name==input$gene_name, ] })
  uniprotID <- eventReactive(input$goButton,{ as.character(gene_attrib()$uniprot_id)  })
  uniprot_data <- eventReactive(input$goButton,{ extract_uniprot_feature_data(uniprotID())  })
  
  # AlphaFold data extraction - FIXED
  af <- eventReactive(input$goButton, { 
    extract_alphafold_plddt(uniprotID()) 
  })
  
  # AlphaFold Mean Pathogenicity - FIXED
  mean_data <- eventReactive(input$goButton, {
    get_mean_pathogenicity(uniprotID())
  })
  # afs_data: raw per-substitution AlphaMissense CSV rows, loaded separately
  # (same file that get_mean_pathogenicity downloads — read from disk cache)
  afs_data <- eventReactive(input$goButton, {
    uid <- uniprotID()
    dest <- file.path(getwd(), paste0(uid, "-F1-aa-substitutions.csv"))
    if (!file.exists(dest)) {
      # Try to download if not present (may already be there from mean_data fetch)
      urls <- c(
        paste0("https://alphafold.ebi.ac.uk/files/AF-", uid, "-F1-aa-substitutions.csv"),
        paste0("https://alphafold.ebi.ac.uk/files/AF-", uid, "-F1-aa-substitutions-v1.csv")
      )
      for (url in urls) {
        ok <- tryCatch({
          download.file(url, dest, mode = "wb", quiet = TRUE)
          file.exists(dest) && file.size(dest) > 0
        }, error = function(e) FALSE, warning = function(w) FALSE)
        if (ok) break
      }
    }
    if (file.exists(dest) && file.size(dest) > 0) {
      tryCatch(data.table::fread(dest), error = function(e) NULL)
    } else NULL
  })
  
  pfam_data <- eventReactive(input$goButton,{ extract_pfam(uniprotID()) })
  
  # Gene info from UniProt (for GeneInfo tab)
  gene_info_uniprot <- eventReactive(input$goButton, {
    extract_gene_info_uniprot(uniprotID(), input$gene_name)
  })
  
  gene_clinvar_data <- eventReactive(input$goButton,{ extract_clinvar(input$gene_name) })

  # ClinGen gene validity — fetched once per gene at top level so Gene Info card can show badge
  clingen_validity_reactive <- eventReactive(input$goButton, {
    # Get HGNC ID from UniProt gi (most reliable — fetched from UniProt cross-refs)
    # Falls back to NULL if gene_info_uniprot not yet available
    hgnc_for_clingen <- tryCatch({
      gi <- gene_info_uniprot()
      if (!is.null(gi$hgnc_id) && nchar(gi$hgnc_id) > 0) gi$hgnc_id else NULL
    }, error = function(e) NULL)
    tryCatch(
      fetch_clingen_validity(input$gene_name, hgnc_id = hgnc_for_clingen),
      error = function(e) list(classification = NA_character_, moi = "", disease = "", source = "ClinGen LDH")
    )
  })
  gene_gnomad_data <- eventReactive(input$goButton,{ extract_gnomad(input$gene_name) })
  gene_ccrs_data <- eventReactive(input$goButton,{ extract_ccrs(input$gene_name, pfam_data()$primaryAccession) })

  # Multi-conservation track: UCSC REST API for true per-base PhyloP/PhastCons/GERP++ scores
  # Conservation scores: reactive on EITHER Go button OR multiconservation checkbox toggle.
  # Using reactiveVal + observeEvent so checking the box after Go also triggers fetch.
  conservation_scores_val <- reactiveVal(NULL)
  multicons_state         <- reactiveVal("idle")  # idle|loading|ok|absent|error


  observe({
    # Re-evaluate whenever Go is clicked OR multiconservation checkbox changes
    input$goButton
    input$plotselection

    isolate({
      if (!"multiconservation" %in% input$plotselection) {
        multicons_state("idle")
        return()
      }
      if (input$goButton == 0) return()
      gene <- trimws(input$gene_name)
      if (nchar(gene) == 0) return()

      # Cache hit — instant
      cached <- cache_get("ucsc_cons", gene)
      if (!is.null(cached)) {
        conservation_scores_val(cached)
        n_ok <- sum(!is.na(cached$phylop100))
        multicons_state(list(status = "ok", gene = gene, covered = n_ok,
                             total = nrow(cached)))
        return()
      }

      # Get protein length
      prot_len <- tryCatch(pfam_data()$sequence$length, error = function(e) NULL)
      if (is.null(prot_len) || prot_len == 0) {
        prot_len <- tryCatch({
          gd <- gene_gnomad_data()
          if (!is.null(gd) && nrow(gd) > 0) max(gd$prot_pos, na.rm = TRUE) else NULL
        }, error = function(e) NULL)
      }
      if (is.null(prot_len) || prot_len == 0 || !is.finite(prot_len)) {
        message("[UCSC Cons] No protein length available, skipping conservation fetch")
        multicons_state(list(status = "error", msg = "Protein length unavailable"))
        return()
      }

      multicons_state("loading")
      withProgress(message = "Fetching conservation scores from UCSC…", value = 0.5, {
        result <- tryCatch(
          fetch_conservation_scores(gene, as.integer(prot_len)),
          error = function(e) { message("[UCSC Cons] Fetch error: ", e$message); NULL }
        )
        conservation_scores_val(result)
        if (!is.null(result) && nrow(result) > 0) {
        n_ok   <- sum(!is.na(result$phylop100))
        n_gerp <- sum(!is.na(result$gerp_raw))
        multicons_state(list(status = "ok", gene = gene, covered = n_ok,
                             total = nrow(result)))
        } else {
          multicons_state(list(status = "absent", gene = gene))
        }
      })
    })
  })

  conservation_scores <- reactive({ conservation_scores_val() })

  # ============================================================
  # BACKGROUND PREFETCH: On Go click, start all API calls and
  # show GeneInfo tab immediately. Data is cached for plotting.
  # ============================================================
  prefetch_state <- reactiveValues(
    done = FALSE,
    variant_table = NULL
  )
  
  observeEvent(input$goButton, {
    # Guard: require both gene and at least one variant before doing any work
    missing_gene <- is.null(input$gene_name) || nchar(trimws(input$gene_name)) == 0
    missing_vars <- is.null(input$variants)  || nchar(trimws(input$variants))  == 0

    if (missing_gene || missing_vars) {
      missing_items <- c(
        if (missing_gene) "<li><strong>Gene symbol</strong> — select a gene from the dropdown</li>",
        if (missing_vars) "<li><strong>Variant(s)</strong> — enter at least one pMut (e.g. p.Arg175His)</li>"
      )
      showModal(modalDialog(
        title = tags$span(
          tags$span(style = "color:#dc2626; margin-right:8px;", HTML("&#9888;")),
          "Required fields missing"
        ),
        HTML(paste0(
          '<div style="font-size:14px; line-height:1.7; color:#374151;">',
          '<p style="margin:0 0 10px;">Please fill in the following before running the analysis:</p>',
          '<ul style="margin:0; padding-left:20px; color:#1e3a5f;">',
          paste(missing_items, collapse = ""),
          '</ul>',
          '</div>'
        )),
        footer = modalButton("Got it"),
        easyClose = TRUE,
        size = "s"
      ))
      return()
    }

    # Reset state
    prefetch_state$done <- FALSE
    prefetch_state$variant_table <- NULL
    
    # Stay on GeneInfo tab so user sees gene info while data loads
    updateTabsetPanel(session, "main_tabs", selected = "GeneInfo")
    
    # Show persistent notification toast (bottom-right) like the plot loading bar
    nid <- showNotification(
      HTML('<div style="font-size:13px;"><strong>Loading gene data...</strong> Fetching UniProt info...</div>'),
      duration = NULL, closeButton = TRUE, type = "message"
    )
    
    # Helper to update the notification toast
    update_note <- function(detail) {
      showNotification(
        HTML(paste0('<div style="font-size:13px;"><strong>Loading gene data...</strong> ', detail, '</div>')),
        id = nid, duration = NULL, closeButton = TRUE, type = "message"
      )
    }
    
    # Kick off all API calls sequentially (they cache internally)
    # This warms the cache so the Plot tab renders instantly
    tryCatch({
      gi <- gene_info_uniprot()    # Fast - shows GeneInfo immediately
    }, error = function(e) { message("[Prefetch] UniProt error: ", e$message) })
    
    update_note("Fetching protein domains...")
    tryCatch({
      pfam_data()
    }, error = function(e) { message("[Prefetch] Pfam error: ", e$message) })
    
    update_note("Fetching gnomAD variants...")
    tryCatch({
      gd <- gene_gnomad_data()
      if (is.null(gd) || nrow(gd) == 0) {
        showNotification(
          HTML('<div style="font-size:13px;">&#9888; <strong>gnomAD</strong> returned no variants for this gene.</div>'),
          duration = 8, closeButton = TRUE, type = "warning"
        )
      }
    }, error = function(e) { message("[Prefetch] gnomAD error: ", e$message) })
    
    update_note("Fetching ClinVar variants...")
    tryCatch({
      cv <- gene_clinvar_data()
      if (is.null(cv) || nrow(cv) == 0) {
        showNotification(
          HTML('<div style="font-size:13px;">&#9888; <strong>ClinVar</strong> returned no pathogenic variants for this gene.</div>'),
          duration = 8, closeButton = TRUE, type = "warning"
        )
      }
    }, error = function(e) { message("[Prefetch] ClinVar error: ", e$message) })
    
    update_note("Fetching AlphaFold pLDDT...")
    af_ok <- tryCatch({
      af_result <- af()
      !is.null(af_result) && is.data.frame(af_result) && nrow(af_result) > 0
    }, error = function(e) { message("[Prefetch] AlphaFold pLDDT error: ", e$message); FALSE })
    if (!af_ok) {
      showNotification(
        HTML('<div style="font-size:13px;">&#9888; <strong>AlphaFold pLDDT</strong> data not available for this protein. Track will show placeholder.</div>'),
        duration = 8, closeButton = TRUE, type = "warning"
      )
    }
    
    update_note("Fetching AF Pathogenicity...")
    afm_ok <- tryCatch({
      md_result <- mean_data()
      !is.null(md_result) && is.data.frame(md_result) && nrow(md_result) > 0
    }, error = function(e) { message("[Prefetch] AF Pathogenicity error: ", e$message); FALSE })
    if (!afm_ok) {
      showNotification(
        HTML('<div style="font-size:13px;">&#9888; <strong>AlphaFold Pathogenicity</strong> data not available for this protein. Track will show placeholder.</div>'),
        duration = 8, closeButton = TRUE, type = "warning"
      )
    }
    
    update_note("Fetching constrained regions (CCRS)...")
    tryCatch({
      gene_ccrs_data()
    }, error = function(e) { message("[Prefetch] CCRS error: ", e$message) })
    
    # Prefetch dbNSFP scores for each variant via MyVariant.info
    update_note("Fetching dbNSFP pathogenicity scores...")
    tryCatch({
      h <- highlight()
      if (nrow(h) > 0) {
        for (vi in seq_len(nrow(h))) {
          tryCatch({
            fetch_dbnsfp(input$gene_name, as.character(h$Mutation[vi]))
          }, error = function(e) message("[Prefetch] dbNSFP error for ", h$Mutation[vi], ": ", e$message))
        }
      }
    }, error = function(e) { message("[Prefetch] dbNSFP batch error: ", e$message) })
    
    # Now build the variant intersection table
    update_note("Building variant annotation table + ACMG tags...")
    tryCatch({
      h <- highlight()
      message("[VariantTable] highlight has ", nrow(h), " variants")
      if (nrow(h) > 0) {
        # Get cutoff values for gnomAD filter color-coding
        cutoff_vals <- tryCatch(data(), error = function(e) list(NA, NA))
        current_af_cutoff <- cutoff_vals[[1]]
        current_ac_cutoff <- cutoff_vals[[2]]
        
        # Get ClinVar missense subset for red density curve
        cv_data <- tryCatch(gene_clinvar_data(), error = function(e) NULL)
        cv_missense <- NULL
        if (!is.null(cv_data) && is.data.frame(cv_data) && nrow(cv_data) > 0 && "type" %in% colnames(cv_data)) {
          cv_missense <- cv_data[cv_data$type == "missense_variant", , drop = FALSE]
        }
        
        cs_data <- tryCatch(consurf_score(), error = function(e) { 
          message("[VariantTable] ConSurf error: ", e$message); NULL 
        })
        if (!is.null(cs_data) && is.data.frame(cs_data) && nrow(cs_data) > 0) {
          message("[VariantTable] ConSurf data available: ", nrow(cs_data), " positions, cols: ", 
                  paste(colnames(cs_data), collapse=", "))
        } else {
          message("[VariantTable] No ConSurf data available")
        }
        
        # Capture all session analysis parameters for provenance tracking
        cg_res     <- tryCatch(clingen_validity(), error = function(e) NULL)
        cg_disease <- if (!is.null(cg_res) && length(cg_res$disease) > 0) cg_res$disease[1] else ""
        cg_moi     <- if (!is.null(cg_res) && length(cg_res$moi)     > 0) cg_res$moi[1]     else ""
        cs_fname   <- if (!is.null(input$consurf_file)) input$consurf_file$name else ""
        vtbl <- build_variant_table(
          h, af(), mean_data(), afs_data(), gene_gnomad_data(), cv_data,
          pfam_data(), uniprot_data(), gene_ccrs_data(),
          af_cutoff             = current_af_cutoff,
          ac_cutoff             = current_ac_cutoff,
          clinvar_missense       = cv_missense,
          consurf_data           = cs_data,
          denovo_status          = "not_denovo",  # handled per-variant in card dropdowns
          inh_param              = isolate(input$inh),
          cutoff_method          = isolate(input$cutoff_method),
          prevalence_1_in_n      = isolate(input$prev),
          allelic_het            = isolate(input$hetA),
          genetic_het            = isolate(input$hetG),
          penetrance             = isolate(input$pen),
          pop_size               = isolate(if (input$cutoff_method == "calc_af") input$popSize_af else input$popSize),
          conf_interval          = isolate(if (input$cutoff_method == "calc_af") input$CI_af else input$CI),
          clingen_disease_param  = cg_disease,
          clingen_moi_param      = cg_moi,
          consurf_file_name      = cs_fname
        )
        message("[VariantTable] Built table with ", 
                if(!is.null(vtbl)) nrow(vtbl) else 0, " rows")
        prefetch_state$variant_table <- vtbl
      }
    }, error = function(e) {
      message("[VariantTable] Error building table: ", e$message)
      prefetch_state$variant_table <- paste("Error:", e$message)
    })
    
    prefetch_state$done <- TRUE
    
    # Replace loading toast with completion toast (auto-dismiss after 4s)
    showNotification(
      HTML('<div style="font-size:13px;"><strong>&#10003; All data loaded</strong> — Plot tab ready!</div>'),
      id = nid, duration = 4, closeButton = TRUE, type = "message"
    )
  })


  # --- Variant Intersection Table (below GeneInfo) ---
  output$variant_table_section <- renderUI({
    if (input$goButton == 0) return(NULL)
    if (!prefetch_state$done) return(NULL)
    
    vtbl <- prefetch_state$variant_table
    
    # Handle error case (string message stored)
    if (is.character(vtbl)) {
      return(HTML(paste0(
        '<div style="background:#fef2f2; border:1px solid #fecaca; border-radius:12px; padding:18px; margin-bottom:18px;">',
        '<h4 style="color:#ef4444; margin:0 0 8px;">Variant Table Error</h4>',
        '<p style="color:#991b1b; font-size:13px;">', htmltools::htmlEscape(vtbl), '</p>',
        '</div>'
      )))
    }
    
    if (is.null(vtbl) || !is.data.frame(vtbl) || nrow(vtbl) == 0) {
      return(HTML(paste0(
        '<div style="background:#f8fafc; border:1px solid #e2e8f0; border-radius:12px; padding:18px; margin-bottom:18px;">',
        '<p style="color:#94a3b8; font-style:italic;">No variant annotation data to display. ',
        'Ensure variants are entered in pMut format (e.g., p.I554N) before clicking Go.</p>',
        '</div>'
      )))
    }
    
    # Build a styled HTML table
    esc <- function(x) {
      x <- as.character(x)
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      x <- gsub(">", "&gt;", x, fixed = TRUE)
      x
    }
    
    # Color coding helpers
    plddt_color <- function(val) {
      v <- suppressWarnings(as.numeric(val))
      if (is.na(v)) return("#64748b")
      if (v >= 90) return("#0053d6")
      if (v >= 70) return("#65cbf3")
      if (v >= 50) return("#ffdb13")
      return("#ff7d45")
    }
    
    table_html <- ""   # Removed — all data shown in Variant Summary cards below

    # ── Variant Summary Cards (one per variant, grouped field sections) ─────
    comment_html <- ""
    if (nrow(vtbl) > 0) {

      # ── helper: render one label+value cell ──────────────────────────────
      # Score/numeric cells — nowrap, compact (predictor scores, population AFs)
      sum_cell <- function(label, value, color = "#374151", note = "", bold = FALSE) {
        bw <- if (bold) "font-weight:700;" else "font-weight:400;"
        note_html <- if (nchar(note) > 0)
          paste0('<br><span style="font-size:9px;color:#94a3b8;white-space:nowrap;">', note, '</span>') else ""
        paste0(
          '<td style="padding:3px 6px; vertical-align:top; border:1px solid #e5e7eb; ',
          'background:#fff; min-width:55px; max-width:110px;">',
          '<div style="font-size:8px;color:#94a3b8;text-transform:uppercase;',
          'letter-spacing:0.3px;margin-bottom:1px;white-space:nowrap;overflow:hidden;text-overflow:ellipsis;">', label, '</div>',
          '<div style="font-size:11px;color:', color, ';', bw, 'white-space:nowrap;">', value, note_html, '</div>',
          '</td>'
        )
      }
      # Info/text cells — wrapping allowed, used for Variant Info row fields
      info_cell <- function(label, value, color = "#374151", note = "", bold = FALSE) {
        bw <- if (bold) "font-weight:700;" else "font-weight:400;"
        note_html <- if (nchar(note) > 0)
          paste0('<br><span style="font-size:9px;color:#94a3b8;">', note, '</span>') else ""
        paste0(
          '<td style="padding:3px 7px; vertical-align:top; border:1px solid #e5e7eb; ',
          'background:#fff; min-width:70px; max-width:150px; word-break:break-word;">',
          '<div style="font-size:8px;color:#94a3b8;text-transform:uppercase;',
          'letter-spacing:0.3px;margin-bottom:1px;white-space:nowrap;">', label, '</div>',
          '<div style="font-size:11px;color:', color, ';', bw, 'line-height:1.3;">', value, note_html, '</div>',
          '</td>'
        )
      }

      # ── helper: score verdict color ──────────────────────────────────────
      vd_col <- function(v) switch(v, damaging="#ef4444", benign="#10b981", ambiguous="#f59e0b", "#94a3b8")

      # ── helper: format score+verdict pair ────────────────────────────────
      score_val <- function(s, v) {
        if (is.na(s) || nchar(as.character(s)) == 0) return('<span style="color:#cbd5e1;">—</span>')
        vc <- vd_col(as.character(v))
        paste0('<span style="color:', vc, ';font-weight:600;">', s,
               '</span> <span style="font-size:9px;color:', vc, ';">', substr(v,1,3), '</span>')
      }

      # ── helper: ACMG tag badge ────────────────────────────────────────────
      acmg_badge <- function(tag) {
        tag <- trimws(tag)
        tc <- if (grepl("_strong$", tag) && grepl("^P", tag)) "#9b1c1c" else if (grepl("_moderate$", tag) && grepl("^P", tag)) "#dc2626" else if (grepl("^PS", tag)) "#dc2626" else if (grepl("^PM", tag)) "#ef4444" else if (grepl("^PP", tag)) "#f97316" else if (grepl("^BS", tag)) "#065f46" else if (grepl("^BP", tag)) "#059669" else if (tag == "BA1") "#1d4ed8" else "#64748b"
        disp <- sub("_strong$", "\u2b06", sub("_moderate$", "\u2191", tag))
        paste0('<span title="', tag, '" style="display:inline-block;background:', tc, '22;color:', tc,
               ';padding:2px 7px;border-radius:4px;font-weight:700;font-size:11px;',
               'margin:2px 2px;white-space:nowrap;border:1px solid ', tc, '44;">', disp, '</span>')
      }

      # ── VarSome-style criterion grid ─────────────────────────────────────
      tag_pts_map <- c(
        PVS1=8, PS1=4, PS1_moderate=2, PS1_supporting=1, PS2=4, PS3=4, PS3_supporting=1, PS4=4,
        PM1_strong=4, PP3_strong=4, PP1_strong=4,
        PM1=2, PM2=2, PM3=1, PM3_moderate=2, PM3_strong=4, PM4=2, PM5=2, PM6=2, PP3_moderate=2, PP1_moderate=2,
        PP1=1, PP2=1, PP3=1, PP4=1, PP5=1,
        BA1=-8, BS1=-4, BS2=-4, BS3=-4, BS4=-4, BP6=-4,
        BP1=-1, BP2=-1, BP3=-1, BP4=-1, BP5=-1, BP7=-1
      )
      strength_label <- function(tag) {
        if (grepl("_strong$", tag))    return("Strong")
        if (grepl("_moderate$", tag))  return("Moderate")
        if (tag %in% c("PVS1"))        return("Very Strong")
        if (tag %in% c("PS1","PS2","PS3","PS4","BS1","BS2","BS3","BS4","BP6")) return("Strong")
        if (tag %in% c("PM1","PM2","PM4","PM5","PM6")) return("Moderate")
        if (tag == "PM3_strong")   return("Strong")
        if (tag == "PM3_moderate") return("Moderate")
        if (tag == "PM3")          return("Supporting")
        if (tag == "PS3_supporting")  return("Supporting")
        if (tag == "PS1_moderate")   return("Moderate")
        if (tag == "PS1_supporting")  return("Supporting")
        if (tag %in% c("PP1","PP2","PP3","PP4","PP5","BP1","BP2","BP3","BP4","BP5","BP7")) return("Supporting")
        if (tag == "BA1") return("Stand Alone")
        return("Supporting")
      }
      tag_display <- function(tag) sub("_strong$","",sub("_moderate$","",tag))
      tag_color <- function(tag, is_path) {
        if (!is_path) return(list(bg="#dcfce7", border="#16a34a", text="#14532d", badge_bg="#16a34a"))
        if (grepl("^PVS",tag)||grepl("_strong$",tag)) return(list(bg="#fee2e2",border="#dc2626",text="#7f1d1d",badge_bg="#dc2626"))
        if (grepl("^PS",tag))  return(list(bg="#fee2e2",border="#dc2626",text="#7f1d1d",badge_bg="#dc2626"))
        if (grepl("^PM",tag))  return(list(bg="#fff7ed",border="#ea580c",text="#7c2d12",badge_bg="#ea580c"))
        if (grepl("^PP",tag))  return(list(bg="#fef9c3",border="#ca8a04",text="#713f12",badge_bg="#ca8a04"))
        return(list(bg="#f1f5f9",border="#64748b",text="#1e293b",badge_bg="#64748b"))
      }

      make_criterion_card <- function(tag, gevir_pct_val = NA) {
        tag  <- trimws(tag)
        pts  <- if (tag %in% names(tag_pts_map)) tag_pts_map[[tag]] else 0
        is_p <- pts > 0
        lbl  <- strength_label(tag)
        disp <- tag_display(tag)
        col  <- tag_color(tag, is_p)
        pts_str <- if (pts > 0) paste0("+", pts, "pts") else paste0(pts, "pts")
        # Tooltip: generic criterion description, with GeVIR detail for BP1/PP2
        tip <- switch(tag,
          BP1 = if (!is.na(gevir_pct_val))
                  paste0("Missense variant in gene tolerant to missense variation. ",
                         "GeVIR percentile = ", round(gevir_pct_val, 1),
                         " (>75 = tolerant → BP1)")
                else "Missense variant in gene where only truncating variants cause disease",
          PP2 = if (!is.na(gevir_pct_val))
                  paste0("Missense variant in gene intolerant to missense variation. ",
                         "GeVIR percentile = ", round(gevir_pct_val, 1),
                         " (<25 = intolerant → PP2)")
                else "Missense variant in gene with high rate of pathogenic missense variants",
          PS1 = "Same amino acid change as established pathogenic variant (ClinVar)",
          PM1 = "Variant in functional domain (UniProt); conservation evidence may upgrade strength",
          PM2 = "Absent or rare in gnomAD at disease-prevalence-adjusted threshold",
          PM5 = "Novel missense at codon with different established pathogenic change (ClinVar)",
          PP3 = "Computational evidence of deleteriousness (Pejaver 2022 calibrated thresholds)",
          PP5 = "Reported as pathogenic in ClinVar with at least 1-star review",
          BA1 = "Allele frequency >5% in gnomAD — standalone Benign",
          BS1 = "Allele frequency above disease-prevalence-adjusted threshold",
          BS2 = "Observed homozygous in gnomAD in healthy individuals",
          BP3 = "In-frame indel in repetitive region without known function",
          BP4 = "Multiple computational tools predict benign effect",
          BP6 = "Reported as benign in ClinVar",
          BP7 = "Synonymous variant with no predicted splice impact",
          tag  # fallback: show tag name
        )
        paste0(
          '<div title="', gsub('"', "&quot;", tip), '" ',
          'style="display:inline-flex;flex-direction:column;align-items:center;',
          'background:', col$bg, ';border:2px solid ', col$border, ';',
          'border-radius:10px;padding:5px 10px;margin:3px;min-width:62px;',
          'box-shadow:0 1px 3px rgba(0,0,0,0.08);vertical-align:top;cursor:help;">',
          # Tag name
          '<div style="font-weight:800;font-size:13px;color:', col$text, ';',
          'letter-spacing:0.3px;margin-bottom:2px;">', disp, '</div>',
          # Strength label
          '<div style="font-size:9px;color:', col$text, ';opacity:0.8;',
          'margin-bottom:3px;white-space:nowrap;">', lbl, '</div>',
          # Points badge
          '<div style="background:', col$badge_bg, ';color:#fff;',
          'font-size:9px;font-weight:700;padding:1px 6px;border-radius:10px;',
          'white-space:nowrap;">', pts_str, '</div>',
          '</div>'
        )
      }

      # ── section header row ────────────────────────────────────────────────
      sec_hdr <- function(label, bg, icon = "") {
        paste0('<tr><td colspan="99" style="padding:5px 10px;background:', bg,
               ';color:#fff;font-size:10px;font-weight:700;text-transform:uppercase;',
               'letter-spacing:0.6px;border:none;white-space:nowrap;">', icon, ' ', label, '</td></tr>')
      }

      # ── classification badge ──────────────────────────────────────────────
      cls_colors <- list(
        "Pathogenic"        = list(bg="#fee2e2", col="#9b1c1c"),
        "Likely Pathogenic" = list(bg="#fef3c7", col="#b45309"),
        "VUS-High"          = list(bg="#fff7ed", col="#c2410c"),   # warm orange — leans pathogenic
        "VUS-Mid"           = list(bg="#f1f5f9", col="#475569"),   # neutral grey — truly uncertain
        "VUS-Low"           = list(bg="#f0fdf4", col="#166534"),   # cool green — leans benign
        "VUS"               = list(bg="#f1f5f9", col="#475569"),   # fallback legacy
        "Likely Benign"     = list(bg="#d1fae5", col="#047857"),
        "Benign"            = list(bg="#dcfce7", col="#065f46")
      )
      class_badge <- function(acmg_res) {
        lbl <- acmg_res$classification
        pal <- if (lbl %in% names(cls_colors)) cls_colors[[lbl]] else list(bg="#f1f5f9", col="#475569")
        paste0('<span style="background:', pal$bg, ';color:', pal$col,
               ';padding:3px 10px;border-radius:6px;font-weight:800;font-size:12px;">',
               lbl, '</span>',
               ' <span style="font-size:10px;color:', pal$col,
               ';background:', pal$bg,
               ';padding:2px 7px;border-radius:4px;margin-left:4px;',
               'border:1px solid ', pal$col, '44;">',
               esc(acmg_res$rule), '</span>'
        )
      }

      summary_blocks <- lapply(seq_len(nrow(vtbl)), function(i) {
        r <- vtbl[i, ]

        # ── Per-card cosegregation input (unique ID per variant) ──────────
        seg_input_id <- paste0("seg_", i)
        card_seg     <- if (!is.null(input[[seg_input_id]])) input[[seg_input_id]] else "none"
        pm3_input_id <- paste0("pm3_", i)
        card_pm3     <- if (!is.null(input[[pm3_input_id]])) input[[pm3_input_id]] else "none"
        dn_input_id  <- paste0("dn_", i)
        card_dn      <- if (!is.null(input[[dn_input_id]])) input[[dn_input_id]] else "not_denovo"

        # ── classify using shared engine + per-card PP1 ───────────────────
        tags_vec <- if (nchar(as.character(r$ACMG_Tags)) > 0)
                      trimws(strsplit(as.character(r$ACMG_Tags), ",")[[1]]) else character(0)
        # Add PP1 tag dynamically based on per-card segregation input
        pp1_tag <- switch(card_seg,
          pp1          = "PP1",
          pp1_moderate = "PP1_moderate",
          pp1_strong   = "PP1_strong",
          NULL
        )
        if (!is.null(pp1_tag)) tags_vec <- c(tags_vec, pp1_tag)
        pm3_tag <- switch(card_pm3,
          pm3          = "PM3",
          pm3_moderate = "PM3_moderate",
          pm3_strong   = "PM3_strong",
          NULL
        )
        if (!is.null(pm3_tag)) tags_vec <- c(tags_vec, pm3_tag)
        # PTM hit — add PS3_supporting to tags if strong PTM site
        if (nchar(as.character(r$PTM_ACMG)) > 0 && r$PTM_ACMG == "PS3_supporting")
          tags_vec <- c(tags_vec, "PS3_supporting")
        # Per-card de novo — overrides global denovo_status for this variant
        dn_tag <- switch(card_dn,
          denovo_confirmed = "PS2",
          denovo_assumed   = "PM6",
          NULL
        )
        if (!is.null(dn_tag)) tags_vec <- c(tags_vec, dn_tag)
        acmg_res <- classify_acmg(tags_vec)
        pts      <- acmg_res$pts

        # ── gnomAD AF color — uses same inheritance-aware thresholds as ACMG engine ─
        # Re-derive thresholds from current UI inputs for display consistency
        card_inh_mode <- if (isTRUE(card_dn %in% c("denovo_confirmed","denovo_assumed"))) "monoallelic" else
                         if (isTRUE(input$inh == "biallelic")) "biallelic" else "monoallelic"
        card_bs1_thresh <- if (card_inh_mode == "biallelic") 0.05 else 0.01
        # Pull data() ONCE per card — used for both colour coding and pill display
        cutoff_card   <- tryCatch(data(), error = function(e) list(NA_real_, NA_real_))
        cutoff_live_af <- tryCatch(as.numeric(cutoff_card[[1]]), error = function(e) NA_real_)
        cutoff_live_ac <- tryCatch(as.numeric(cutoff_card[[2]]), error = function(e) NA_real_)
        card_pm2_thresh <- if (card_inh_mode == "biallelic") 0.01 else
                           if (!is.na(cutoff_live_af) && cutoff_live_af > 0) cutoff_live_af else 0.0001

        gaf_v   <- suppressWarnings(as.numeric(r$gnomAD_AF))
        gaf_col <- if (!is.na(gaf_v) && gaf_v > 0.05)             "#2563eb" else
                   if (!is.na(gaf_v) && gaf_v > card_bs1_thresh)  "#f59e0b" else
                   if (!is.na(gaf_v) && gaf_v > 0)                "#ef4444" else "#94a3b8"
        gaf_disp <- if (!is.na(gaf_v) && gaf_v > 0) formatC(gaf_v, format="g", digits=3) else if (nchar(as.character(r$gnomAD_AF)) > 0) as.character(r$gnomAD_AF) else "Absent"

        # ── pLDDT color ───────────────────────────────────────────────────
        plddt_v   <- suppressWarnings(as.numeric(r$pLDDT))
        plddt_col <- if (!is.na(plddt_v) && plddt_v >= 90) "#1d4ed8" else if (!is.na(plddt_v) && plddt_v >= 70) "#0d9488" else if (!is.na(plddt_v) && plddt_v >= 50) "#f59e0b" else "#ef4444"

        # ── ClinVar color ─────────────────────────────────────────────────
        cv_sig  <- as.character(r$ClinVar)
        cv_col  <- if (grepl("pathogenic", tolower(cv_sig)) && !grepl("benign", tolower(cv_sig))) "#ef4444" else if (grepl("benign", tolower(cv_sig)) && !grepl("pathogenic", tolower(cv_sig))) "#10b981" else if (grepl("uncertain|conflicting", tolower(cv_sig))) "#f59e0b" else "#64748b"

        # ── Section rows ──────────────────────────────────────────────────

        # — Row 0: title bar —
        rsid_link <- if (nchar(as.character(r$RSID)) > 0 && !is.na(r$RSID))
          paste0(' &nbsp;<a href="https://www.ncbi.nlm.nih.gov/snp/', esc(r$RSID),
                 '" target="_blank" style="font-size:11px;color:#93c5fd;">',
                 esc(r$RSID), '</a>') else ""

        # ── PubMed MeSH variant search link ───────────────────────────────
        # Query: humans[MeSH Terms] AND GENE AND (variation OR variant OR mutation) AND Ref###Alt
        # Handles both 1-letter (p.R175H) and 3-letter (p.Arg175His) input forms.
        pubmed_link <- tryCatch({
          gene_name_card <- isolate(input$gene_name)
          var_bare <- sub("^p\\.", "", as.character(r$Variant))
          pos_n    <- gsub("[^0-9]", "", var_bare)
          aa1to3_m <- c(A="Ala",R="Arg",N="Asn",D="Asp",C="Cys",E="Glu",Q="Gln",G="Gly",
                        H="His",I="Ile",L="Leu",K="Lys",M="Met",F="Phe",P="Pro",S="Ser",
                        T="Thr",W="Trp",Y="Tyr",V="Val")
          if (grepl("^[A-Z][a-z]{2}[0-9]", var_bare)) {
            m3 <- regmatches(var_bare,
                    regexec("^([A-Z][a-z]{2,})([0-9]+)([A-Z][a-z]{2,}|Ter|del|dup|ins|fs|ext)",
                            var_bare))[[1]]
            ref3 <- if (length(m3) == 4) m3[2] else substr(var_bare, 1, 3)
            alt3 <- if (length(m3) == 4) m3[4] else substr(var_bare, nchar(pos_n) + 4, nchar(var_bare))
          } else {
            ref1 <- substr(var_bare, 1, 1)
            alt1 <- substr(var_bare, nchar(var_bare), nchar(var_bare))
            ref3 <- if (!is.na(aa1to3_m[ref1])) aa1to3_m[[ref1]] else ref1
            alt3 <- if (!is.na(aa1to3_m[alt1])) aa1to3_m[[alt1]] else alt1
          }
          three_letter <- paste0(ref3, pos_n, alt3)
          mesh_query <- paste0("humans[MeSH Terms] AND ", gene_name_card,
                               " AND (variation OR variant OR mutation) AND ", three_letter)
          pubmed_url <- paste0("https://pubmed.ncbi.nlm.nih.gov/?term=",
                               URLencode(mesh_query, reserved = TRUE))
          paste0('&nbsp;<a href="', pubmed_url, '" target="_blank" rel="noopener" ',
                 'style="font-size:10px;color:#bfdbfe;text-decoration:none;',
                 'background:rgba(255,255,255,0.12);border-radius:4px;',
                 'padding:1px 7px;white-space:nowrap;border:1px solid rgba(255,255,255,0.2);">',
                 '&#128269; PubMed</a>')
        }, error = function(e) "")

        title_bar <- paste0(
          '<tr><td colspan="99" style="padding:10px 14px;background:#003f5c;border:none;">',
          '<span style="font-size:16px;font-weight:800;color:#fff;">',
          esc(as.character(r$Variant)), '</span>',
          rsid_link,
          pubmed_link,
          ' &nbsp;&nbsp;', class_badge(acmg_res),
          '<span style="float:right;font-size:11px;color:#7dd3fc;margin-top:2px;text-align:right;">',
          'Score: <strong style="font-size:13px;color:#fff;">', acmg_res$pts, ' pts</strong>',
          if (nchar(acmg_res$pts_str) > 0)
            paste0('<br><span style="font-size:10px;color:#93c5fd;letter-spacing:0.2px;">',
                   esc(acmg_res$pts_str), '</span>') else "",
          '</span>',
          '</td></tr>'
        )

        # — Row 1: Variant Info (Variant → PTM) —
        row1 <- paste0(
          sec_hdr("Variant Info", "#1e3a5f", "&#9432;"),
          '<tr>',
          info_cell("Variant",   esc(as.character(r$Variant)), "#003f5c", bold=TRUE),
          info_cell("Position",  esc(as.character(r$Position))),
          info_cell("rsID",      if (nchar(as.character(r$RSID))>0 && !is.na(r$RSID))
                                  paste0('<a href="https://www.ncbi.nlm.nih.gov/snp/',esc(r$RSID),
                                         '" target="_blank" style="color:#0369a1;">',esc(r$RSID),'</a>') else '<span style="color:#cbd5e1;">—</span>'),
          info_cell("pLDDT",     if (nchar(as.character(r$pLDDT))>0)
                                  paste0('<span style="color:',plddt_col,';font-weight:600;">',r$pLDDT,'</span> <span style="font-size:9px;color:',plddt_col,';">',r$pLDDT_Category,'</span>') else '<span style="color:#cbd5e1;">—</span>'),
          info_cell("AF Path",   if (nchar(as.character(r$AF_Mean_Pathogenicity))>0)
                                  paste0(r$AF_Mean_Pathogenicity,' <span style="font-size:9px;color:#94a3b8;">',r$AF_Path_Category,'</span>') else '<span style="color:#cbd5e1;">—</span>'),
          info_cell("Density",   if (nchar(as.character(r$Density_Dominant))>0)
                                  paste0('<span style="color:', if(r$Density_Dominant=="red")"#ef4444" else "#3b82f6", ';font-weight:600;">',
                                         if(r$Density_Dominant=="red") r$Density_ClinVar else r$Density_gnomAD,
                                         '</span> <span style="font-size:9px;color:#94a3b8;">',r$Density_Dominant,'</span>') else '<span style="color:#cbd5e1;">—</span>'),
          info_cell("Domain",    if (nchar(as.character(r$Domain))>0) esc(r$Domain) else '<span style="color:#cbd5e1;">—</span>'),
          info_cell("CCRS",      esc(as.character(r$CCRS))),
          info_cell("ScoreCons", if (nchar(as.character(r$ScoreCons))>0) {
                                  sc <- suppressWarnings(as.numeric(r$ScoreCons))
                                  sc_col <- if (!is.na(sc) && sc>=0.7)"#a02560" else if (!is.na(sc) && sc>=0.4)"#f59e0b" else "#10b981"
                                  paste0('<span style="color:',sc_col,';font-weight:600;">',r$ScoreCons,'</span>')
                                } else '<span style="color:#cbd5e1;">—</span>'),
          info_cell("PTM",       if (nchar(as.character(r$PTM))>0) esc(r$PTM) else '<span style="color:#cbd5e1;">—</span>'),
          info_cell("ClinVar",   paste0('<span style="color:',cv_col,';font-weight:600;">',esc(cv_sig),'</span>'),
                                note=paste0(
                                  {
                                    n_stars <- suppressWarnings(as.integer(r$ClinVar_Stars))
                                    if (!is.na(n_stars) && n_stars >= 0) {
                                      star_labels <- c("0"="No assertion criteria",
                                                        "1"="Criteria provided, single submitter",
                                                        "2"="Criteria provided, multiple submitters",
                                                        "3"="Reviewed by expert panel",
                                                        "4"="Practice guideline")
                                      star_tip <- if (as.character(n_stars) %in% names(star_labels))
                                                    star_labels[[as.character(n_stars)]] else ""
                                      filled   <- paste(rep("\u2b50", n_stars),      collapse="")
                                      empty    <- paste(rep("\u2606", 4L - n_stars), collapse="")
                                      paste0('<span title="', star_tip, '" style="font-size:11px;letter-spacing:1px;">',
                                             filled, empty, '</span> ')
                                    } else ""
                                  },
                                  if(nchar(as.character(r$ClinVar_Trait))>0) paste0(" · ", r$ClinVar_Trait) else ""
                                )),
          '</tr>'
        )

        # — Row 2: Pathogenicity —
        row2 <- paste0(
          sec_hdr("&#9881; Pathogenicity", "#b91c1c"),
          '<tr>',
          sum_cell("SIFT",      score_val(r$SIFT_Score, r$SIFT_V)),
          sum_cell("PP2 HDIV",  score_val(r$PP2_HDIV_Score, r$PP2_HDIV_V)),
          sum_cell("PP2 HVAR",  score_val(r$PP2_HVAR_Score, r$PP2_HVAR_V)),
          sum_cell("LRT",       score_val(r$LRT_Score, r$LRT_V)),
          sum_cell("FATHMM",    score_val(r$FATHMM_Score, r$FATHMM_V)),
          sum_cell("PROVEAN",   score_val(r$PROVEAN_Score, r$PROVEAN_V)),
          sum_cell("MutTaster", score_val(r$MutTaster_Score, r$MutTaster_V)),
          '</tr>'
        )

        # — Row 3: Ensemble —
        row3 <- paste0(
          sec_hdr("&#9733; Ensemble", "#6d28d9"),
          '<tr>',
          sum_cell("REVEL",     score_val(r$REVEL_Score, r$REVEL_V)),
          sum_cell("MetaSVM",   score_val(r$MetaSVM_Score, r$MetaSVM_V)),
          sum_cell("MetaLR",    score_val(r$MetaLR_Score, r$MetaLR_V)),
          sum_cell("MetaRNN",   score_val(r$MetaRNN_Score, r$MetaRNN_V)),
          sum_cell("CADD",      score_val(r$CADD_Score, r$CADD_V)),
          sum_cell("DANN",      score_val(r$DANN_Score, r$DANN_V)),
          sum_cell("AlphaMiss", score_val(r$AM_Score, r$AM_V),
                   note = if (!is.na(r$AM_Source)) as.character(r$AM_Source) else ""),
          {
            ada_v <- suppressWarnings(as.numeric(r$dbscSNV_ADA))
            rf_v  <- suppressWarnings(as.numeric(r$dbscSNV_RF))
            ada_col <- if (!is.na(ada_v) && ada_v>0.6)"#ef4444" else if (!is.na(ada_v))"#10b981" else "#94a3b8"
            rf_col  <- if (!is.na(rf_v)  && rf_v >0.6)"#ef4444" else if (!is.na(rf_v)) "#10b981" else "#94a3b8"
            paste0(
              sum_cell("dbscSNV ADA", if (!is.na(ada_v)) paste0('<span style="color:',ada_col,';font-weight:600;">',round(ada_v,3),'</span>') else '<span style="color:#cbd5e1;">—</span>',
                       note=if (!is.na(ada_v) && ada_v>0.6)"splice" else ""),
              sum_cell("dbscSNV RF",  if (!is.na(rf_v))  paste0('<span style="color:',rf_col, ';font-weight:600;">',round(rf_v,3),'</span>') else '<span style="color:#cbd5e1;">—</span>',
                       note=if (!is.na(rf_v)  && rf_v >0.6)"splice" else "")
            )
          },
          '</tr>'
        )

        # — Row 4: Conservation —
        gerp_v <- suppressWarnings(as.numeric(r$GERP_RS))
        gerp_col <- if (!is.na(gerp_v) && gerp_v>4.4)"#059669" else if (!is.na(gerp_v) && gerp_v>2)"#0d9488" else "#64748b"
        row4 <- paste0(
          sec_hdr("&#127795; Conservation", "#047857"),
          '<tr>',
          sum_cell("PhyloP 100V", {
            v <- suppressWarnings(as.numeric(r$PhyloP_100V))
            if (!is.na(v)) paste0('<span style="color:',if(abs(v)>4)"#059669" else if(abs(v)>2)"#0d9488" else "#64748b",';font-weight:600;">',round(v,3),'</span>') else '<span style="color:#cbd5e1;">—</span>'
          }),
          sum_cell("PhyloP 470M", {
            v <- suppressWarnings(as.numeric(r$PhyloP_470M))
            if (!is.na(v)) paste0('<span style="color:',if(abs(v)>4)"#059669" else if(abs(v)>2)"#0d9488" else "#64748b",';font-weight:600;">',round(v,3),'</span>') else '<span style="color:#cbd5e1;">—</span>'
          }),
          sum_cell("PhastCons", {
            v <- suppressWarnings(as.numeric(r$PhastCons))
            if (!is.na(v)) paste0('<span style="color:',if(v>0.9)"#059669" else if(v>0.5)"#0d9488" else "#64748b",';font-weight:600;">',round(v,3),'</span>') else '<span style="color:#cbd5e1;">—</span>'
          }),
          sum_cell("GERP++",    if (!is.na(gerp_v)) paste0('<span style="color:',gerp_col,';font-weight:600;">',round(gerp_v,2),'</span>') else '<span style="color:#cbd5e1;">—</span>',
                                note=if (!is.na(gerp_v) && gerp_v>4.4)"conserved (PP3)" else ""),
          sum_cell("ScoreCons", if (nchar(as.character(r$ScoreCons))>0) {
                                  sc <- suppressWarnings(as.numeric(r$ScoreCons))
                                  sc_col <- if (!is.na(sc)&&sc>=0.7)"#a02560" else if (!is.na(sc)&&sc>=0.4)"#f59e0b" else "#10b981"
                                  paste0('<span style="color:',sc_col,';font-weight:600;">',r$ScoreCons,'</span>')
                                } else '<span style="color:#cbd5e1;">—</span>',
                                note="inter-species"),
          if (nchar(as.character(r$ConSurf_Grade))>0) {
            cs_colors <- c("1"="#10c8d1","2"="#8cffff","3"="#d7ffff","4"="#eaffff",
                           "5"="#ffffff","6"="#fcedf4","7"="#fac9de","8"="#f07dab","9"="#a02560")
            cs_fg     <- c("1"="#0a7a80","2"="#2a8a8a","3"="#4a8a8a","4"="#6a8a8a",
                           "5"="#888","6"="#8a4a6a","7"="#9a2a5a","8"="#8b0057","9"="#fff")
            csg <- as.character(r$ConSurf_Grade)
            bg_c <- if (csg %in% names(cs_colors)) cs_colors[[csg]] else "#e2e8f0"
            fg_c <- if (csg %in% names(cs_fg))     cs_fg[[csg]]     else "#333"
            paste0(
              sum_cell("ConSurf",
                paste0('<span style="background:',bg_c,';color:',fg_c,
                       ';padding:1px 8px;border-radius:4px;font-weight:700;border:1px solid #ccc;">',
                       csg,'</span>'),
                note=paste0("score: ",r$ConSurf_Score," · ",r$ConSurf_Burial))
            )
          } else "",
          '</tr>'
        )

        # — Row 5: Population —
        fmt_pop <- function(v) {
          n <- suppressWarnings(as.numeric(v))
          if (is.na(n) || length(n)==0) return('<span style="color:#cbd5e1;">—</span>')
          if (n == 0) return('<span style="color:#6b7280;font-style:italic;font-size:10px;">Absent</span>')
          col <- if (n < 0.0001)"#ef4444" else if (n < 0.01)"#f59e0b" else "#3b82f6"
          paste0('<span style="color:',col,';font-weight:600;">',formatC(n,format="e",digits=2),'</span>')
        }
        nhom_v <- suppressWarnings(as.integer(r$gnomAD_Nhomalt))
        row5 <- paste0(
          sec_hdr("&#127758; Population", "#1d4ed8"),
          '<tr>',
          {
            gv <- suppressWarnings(as.numeric(r$GeVIR_Gene_Pct))
            gv_col  <- if (!is.na(gv) && gv < 25)  "#ef4444" else
                       if (!is.na(gv) && gv > 75)  "#3b82f6" else "#64748b"
            gv_note <- if (!is.na(gv) && gv < 25)  "intolerant → PP2" else
                       if (!is.na(gv) && gv > 75)  "tolerant → BP1"  else ""
            gv_disp <- if (!is.na(gv)) paste0('<span style="color:',gv_col,
                         ';font-weight:600;">',round(gv,1),'th pct</span>') else
                         '<span style="color:#cbd5e1;">—</span>'
            sum_cell("GeVIR", gv_disp, note=gv_note)
          },
          sum_cell("gnomAD AF",  paste0('<span style="color:',gaf_col,';font-weight:600;">',gaf_disp,'</span>'),
                                 note=if (!is.na(gaf_v) && gaf_v>0.05)"BA1 (common)" else if (!is.na(gaf_v) && gaf_v>card_bs1_thresh)"BS1" else ""),
          sum_cell("gnomAD AC",  if (nchar(as.character(r$gnomAD_AC))>0) esc(as.character(r$gnomAD_AC)) else '<span style="color:#cbd5e1;">—</span>'),
          sum_cell("Nhomalt",    if (!is.na(nhom_v) && nchar(as.character(r$gnomAD_Nhomalt))>0)
                                   paste0('<span style="color:',if(nhom_v>0)"#2563eb" else "#94a3b8",';font-weight:600;">',nhom_v,'</span>') else '<span style="color:#cbd5e1;">—</span>',
                                 note=if (!is.na(nhom_v) && nhom_v>0) "BS2" else ""),
          sum_cell("AFR",  fmt_pop(r$Pop_AFR)),
          sum_cell("NFE",  fmt_pop(r$Pop_NFE)),
          sum_cell("EAS",  fmt_pop(r$Pop_EAS)),
          sum_cell("SAS",  fmt_pop(r$Pop_SAS)),
          sum_cell("FIN",  fmt_pop(r$Pop_FIN)),
          sum_cell("1000G",fmt_pop(r$Pop_1KG)),
          sum_cell("ExAC", fmt_pop(r$Pop_ExAC)),
          '</tr>'
        )

        # — Row 6: ACMG — VarSome-style criterion grid
        path_tags  <- tags_vec[tags_vec %in% names(tag_pts_map) & sapply(tags_vec, function(t) if(t %in% names(tag_pts_map)) tag_pts_map[[t]] > 0 else FALSE)]
        benign_tags <- tags_vec[tags_vec %in% names(tag_pts_map) & sapply(tags_vec, function(t) if(t %in% names(tag_pts_map)) tag_pts_map[[t]] < 0 else FALSE)]

        path_grid  <- if (length(path_tags) > 0)
          paste0('<div style="margin-bottom:4px;">',
                 '<div style="font-size:9px;font-weight:700;color:#dc2626;text-transform:uppercase;',
                 'letter-spacing:0.5px;margin-bottom:3px;">&#9650; Pathogenic</div>',
                 paste(sapply(path_tags,  function(t) make_criterion_card(t, gevir_pct_val=r$GeVIR_Gene_Pct)), collapse=""),
                 '</div>') else ""

        benign_grid <- if (length(benign_tags) > 0)
          paste0('<div>',
                 '<div style="font-size:9px;font-weight:700;color:#16a34a;text-transform:uppercase;',
                 'letter-spacing:0.5px;margin-bottom:3px;">&#9660; Benign</div>',
                 paste(sapply(benign_tags, function(t) make_criterion_card(t, gevir_pct_val=r$GeVIR_Gene_Pct)), collapse=""),
                 '</div>') else ""

        badges_html <- if (nchar(path_grid) > 0 || nchar(benign_grid) > 0)
          paste0(path_grid, benign_grid)
        else '<span style="color:#cbd5e1;font-size:12px;">No ACMG criteria fired</span>'

        # De novo context label
        denovo_note <- if (!is.null(card_dn) && card_dn != "not_denovo") {
          dn_label <- if (card_dn == "denovo_confirmed") "De Novo (confirmed) — PS2 applied" else "De Novo (assumed) — PM6 applied"
          dn_color <- if (card_dn == "denovo_confirmed") "#1e40af" else "#6d28d9"
          paste0(
            '<span style="display:inline-flex;align-items:center;gap:4px;',
            'background:', dn_color, '15;color:', dn_color,
            ';padding:2px 9px;border-radius:12px;font-size:11px;font-weight:700;',
            'border:1px solid ', dn_color, '44;margin-bottom:4px;margin-right:6px;">',
            '&#9698; ', dn_label, '</span>'
          )
        } else ""

        seg_note <- if (card_seg != "none") {
          seg_label <- switch(card_seg,
            pp1          = "1\u20132 affected relatives \u2014 PP1 (Supporting)",
            pp1_moderate = "3\u20134 affected relatives \u2014 PP1 Moderate",
            pp1_strong   = "\u22655 affected relatives \u2014 PP1 Strong",
            ""
          )
          paste0(
            '<span style="display:inline-flex;align-items:center;gap:4px;',
            'background:#f0fdf415;color:#15803d;',
            'padding:2px 9px;border-radius:12px;font-size:11px;font-weight:700;',
            'border:1px solid #16a34a44;margin-bottom:4px;margin-right:6px;">',
            '&#9651; Cosegregation: ', seg_label, '</span>'
          )
        } else ""

        # Inheritance mode context label — show actual computed AF/AC cutoffs from data() reactive
        live_af     <- cutoff_live_af   # already fetched above
        live_ac     <- cutoff_live_ac
        live_af_str <- if (!is.na(live_af) && live_af > 0)
                         paste0("PM2 &lt; ", formatC(live_af, format="e", digits=1))
                       else if (card_inh_mode == "biallelic") "PM2 &lt; 1%" else "PM2 &lt; 0.01%"
        live_bs1_str <- if (card_inh_mode == "biallelic") "BS1 &gt; 5%" else "BS1 &gt; 1%"
        live_ac_str  <- if (!is.na(live_ac) && live_ac > 0)
                          paste0(" | Max AC: ", as.integer(live_ac)) else ""

        inh_note <- paste0(
          '<span style="display:inline-flex;align-items:center;gap:4px;',
          'background:#f1f5f9;color:#475569;',
          'padding:2px 9px;border-radius:12px;font-size:11px;',
          'border:1px solid #cbd5e1;margin-bottom:4px;">',
          if (card_inh_mode == "biallelic") "&#9898; Biallelic" else "&#9898; Monoallelic",
          ' | ', live_bs1_str, ' | ', live_af_str, live_ac_str,
          '</span>'
        )

        # ── Regenerate comment live from current tags (includes PP1/PS2/PM3) ──
        # Derive clingen_class from PP2 tag presence + ClinVar data
        # (clingen_class is not stored in vtbl, infer from tags and GeVIR)
        live_clingen_class <- if ("PP2" %in% tags_vec) {
          gv <- suppressWarnings(as.numeric(r$GeVIR_Gene_Pct))
          if (!is.na(gv) && gv < 25) "Definitive" else ""
        } else ""
        cmt <- tryCatch(
          generate_acmg_comment(
            acmg_tags_str      = paste(tags_vec, collapse = ", "),
            gene               = as.character(r$Variant),
            mut                = as.character(r$Variant),
            gnomad_af          = suppressWarnings(as.numeric(r$gnomAD_AF)),
            gnomad_ac          = suppressWarnings(as.integer(r$gnomAD_AC)),
            gnomad_nhomalt     = suppressWarnings(as.integer(r$gnomAD_Nhomalt)),
            clinvar_sig        = as.character(r$ClinVar),
            clinvar_name       = as.character(r$ClinVar_Name),
            clinvar_vcv        = as.character(r$ClinVar_VCV),
            clinvar_vcv_pos    = as.character(r$ClinVar_VCV_Pos),
            clinvar_trait      = as.character(r$ClinVar_Trait),
            clinvar_match_type = as.character(r$ClinVar_Match),
            domain             = as.character(r$Domain),
            pp2_applies        = isTRUE("PP2" %in% tags_vec),
            revel_score        = suppressWarnings(as.numeric(r$REVEL_Score)),
            am_score           = suppressWarnings(as.numeric(r$AM_Score)),
            cadd_score         = suppressWarnings(as.numeric(r$CADD_Score)),
            metasvm_v          = as.character(r$MetaSVM_V),
            phylop_val         = suppressWarnings(as.numeric(r$PhyloP_100V)),
            phastcons_val      = suppressWarnings(as.numeric(r$PhastCons)),
            gerp_val           = suppressWarnings(as.numeric(r$GERP_RS)),
            consurf_grade      = as.character(r$ConSurf_Grade),
            variant_pos        = suppressWarnings(as.integer(r$Position)),
            clingen_class      = if ("ClinGen_Class" %in% colnames(vtbl)) as.character(r$ClinGen_Class) else live_clingen_class,
            clingen_disease    = as.character(r$ClinGen_Disease),
            clingen_moi        = as.character(r$ClinGen_MOI),
            ps3_proxy          = isTRUE("PS3_supporting" %in% tags_vec),
            gevir_gene_pct     = suppressWarnings(as.numeric(r$GeVIR_Gene_Pct))
          ),
          error = function(e) {
            message("[Comment] Live regen error: ", e$message)
            if ("Comment" %in% colnames(vtbl)) as.character(r$Comment) else ""
          }
        )

        # VUS reclassification note (fires for any VUS tier)
        vus_note <- if (grepl("^VUS", acmg_res$classification)) {
          subtier_msg <- switch(acmg_res$classification,
            "VUS-High" = " This VUS leans pathogenic — prioritize for functional follow-up.",
            "VUS-Mid"  = " Evidence is balanced — monitor literature every 12–18 months.",
            "VUS-Low"  = " This VUS leans benign — unlikely to require immediate action.",
            ""
          )
          paste0(
            '<div style="margin-top:6px;padding:5px 9px;background:#fef9c3;border-left:3px solid #ca8a04;',
            'border-radius:0 5px 5px 0;font-size:11px;color:#78350f;line-height:1.5;">',
            '<strong>&#9203; VUS Reclassification:</strong> Studies show 20–30% of VUS are reclassified ',
            'within 2 years (PMID: 32803765). Re-evaluate every 12–18 months or when new evidence emerges.',
            subtier_msg,
            '</div>'
          )
        } else ""
        # ── Analysis Parameters card row ──────────────────────────────
        inh_label  <- switch(as.character(r$Analysis_Inheritance),
          monoallelic = "Monoallelic (AD / XL)",
          biallelic   = "Biallelic (AR)",
          as.character(r$Analysis_Inheritance))
        pm2_disp   <- if (nchar(as.character(r$Analysis_PM2_Threshold)) > 0)
                        paste0("PM2 AF < ", r$Analysis_PM2_Threshold) else ""
        prev_disp  <- if (nchar(as.character(r$Analysis_Prevalence)) > 0)
                        paste0("Prevalence ", r$Analysis_Prevalence) else ""
        param_pills <- paste0(
          '<span style="background:#f0f9ff;border:1px solid #bae6fd;border-radius:4px;',
          'padding:2px 7px;font-size:10px;color:#0369a1;margin-right:4px;">', inh_label, '</span>',
          if (nchar(pm2_disp) > 0) paste0(
            '<span style="background:#fef9c3;border:1px solid #fde68a;border-radius:4px;',
            'padding:2px 7px;font-size:10px;color:#854d0e;margin-right:4px;">', pm2_disp, '</span>') else "",
          if (nchar(prev_disp) > 0) paste0(
            '<span style="background:#f0fdf4;border:1px solid #bbf7d0;border-radius:4px;',
            'padding:2px 7px;font-size:10px;color:#166534;margin-right:4px;">', prev_disp, '</span>') else "",
          if (nchar(as.character(r$ClinGen_Disease)) > 0) paste0(
            '<span style="background:#fdf4ff;border:1px solid #e9d5ff;border-radius:4px;',
            'padding:2px 7px;font-size:10px;color:#7e22ce;margin-right:4px;">',
            r$ClinGen_Disease,
            if (nchar(as.character(r$ClinGen_MOI)) > 0) paste0(" (", r$ClinGen_MOI, ")") else "",
            '</span>') else "",
          if (nchar(as.character(r$ConSurf_File)) > 0) paste0(
            '<span style="background:#fff7ed;border:1px solid #fed7aa;border-radius:4px;',
            'padding:2px 7px;font-size:10px;color:#9a3412;margin-right:4px;">ConSurf: ',
            r$ConSurf_File, '</span>') else ""
        )
        row5b <- paste0(
          '<tr><td colspan="99" style="padding:4px 10px 4px;border:1px solid #e5e7eb;',
          'background:#f8fafc;font-size:11px;color:#64748b;">',
          '<strong style="color:#374151;">&#9881; Analysis parameters:</strong> ',
          param_pills,
          '</td></tr>'
        )

        row6 <- paste0(
          sec_hdr("&#9670; ACMG Classification", "#b45309"),
          '<tr>',
          '<td colspan="99" style="padding:8px 10px;border:1px solid #e5e7eb;background:#fff;">',
          if (nchar(denovo_note) > 0 || nchar(seg_note) > 0 || TRUE)
            paste0('<div style="margin-bottom:5px;">', denovo_note, seg_note, inh_note, '</div>') else "",
          # Per-card de novo selector
          '<div style="margin-bottom:8px;display:flex;align-items:center;gap:8px;">',
          '<span style="font-size:11px;font-weight:600;color:#374151;white-space:nowrap;">',
          '&#9800; De novo:</span>',
          as.character(selectInput(
            inputId  = paste0("dn_", i),
            label    = NULL,
            choices  = list(
              "Not assessed"                                  = "not_denovo",
              "Confirmed de novo — both parents tested (PS2 +4)"  = "denovo_confirmed",
              "Assumed de novo — parents not tested (PM6 +2)"     = "denovo_assumed"
            ),
            selected = card_dn,
            width    = "320px"
          )),
          '</div>',
          # Per-card cosegregation selector
          '<div style="margin-bottom:8px;display:flex;align-items:center;gap:8px;">',
          '<span style="font-size:11px;font-weight:600;color:#374151;white-space:nowrap;">',
          '&#9651; Cosegregation:</span>',
          as.character(selectInput(
            inputId  = paste0("seg_", i),
            label    = NULL,
            choices  = list(
              "Not assessed"                   = "none",
              "1\u20132 affected relatives (PP1 Supporting)" = "pp1",
              "3\u20134 relatives (PP1 Moderate)"           = "pp1_moderate",
              "\u22655 relatives (PP1 Strong)"              = "pp1_strong"
            ),
            selected = card_seg,
            width    = "280px"
          )),
          '</div>',
          '<div style="margin-bottom:8px;display:flex;align-items:center;gap:8px;">',
          '<span style="font-size:11px;font-weight:600;color:#374151;white-space:nowrap;">',
          '&#9650; In trans (PM3 — AR genes):</span>',
          as.character(selectInput(
            inputId  = paste0("pm3_", i),
            label    = NULL,
            choices  = list(
              "Not assessed (PM3 N/A)"                          = "none",
              "1 compound het observation (PM3 Supporting +1)"  = "pm3",
              "2 compound het observations (PM3 Moderate +2)"   = "pm3_moderate",
              "\u22653 compound het observations (PM3 Strong +4)" = "pm3_strong"
            ),
            selected = card_pm3,
            width    = "320px"
          )),
          '</div>',
          '<div style="margin-bottom:6px;">', badges_html, '</div>',
          if (nchar(trimws(cmt)) > 0)
            paste0('<div style="font-size:12px;color:#374151;line-height:1.7;',
                   'border-top:1px solid #f1f5f9;padding-top:6px;margin-top:4px;">',
                   cmt, '</div>') else "",
          vus_note,
          '</td>',
          '</tr>'
        )

        # — Assemble card —
        paste0(
          '<div style="background:#fff;border:1px solid #e5e7eb;border-radius:14px;',
          'margin-bottom:16px;box-shadow:0 1px 4px rgba(0,0,0,0.06);overflow:hidden;">',
          # title bar always full-width (no scroll)
          '<table style="width:100%;border-collapse:collapse;">', title_bar, '</table>',
          # scrollable body for data rows
          '<div style="overflow-x:auto;-webkit-overflow-scrolling:touch;">',
          '<table style="border-collapse:collapse;table-layout:auto;width:auto;">',
          row1, row2, row3, row4, row5, row5b, row6,
          '</table>',
          '</div>',
          '</div>'
        )
      })

      summary_blocks <- summary_blocks[nchar(summary_blocks) > 0]
      if (length(summary_blocks) > 0) {

        # ── Collapsible cards: show first 2 expanded, rest hidden behind toggle ──
        n_blocks   <- length(summary_blocks)
        first_two  <- paste(summary_blocks[seq_len(min(2, n_blocks))], collapse = "\n")
        rest_cards <- if (n_blocks > 2) paste(summary_blocks[3:n_blocks], collapse = "\n") else ""
        n_hidden   <- max(0, n_blocks - 2)

        collapsed_section <- if (n_hidden > 0) {
          paste0(
            # Toggle button row
            '<div id="vv-more-toggle" ',
            'onclick="',
              'var s=document.getElementById(\'vv-more-cards\');',
              'var b=document.getElementById(\'vv-more-btn\');',
              'var open=s.style.display!==\'none\';',
              'if(open){s.style.display=\'none\';}else{s.style.display=\'block\';}',
              'b.innerHTML=open',
              '?\'&#9660; Show ', n_hidden, ' more variant', if(n_hidden!=1)"s" else "", '\'',
              ':\'&#9650; Hide ', n_hidden, ' variant', if(n_hidden!=1)"s" else "", '\';',
            '" ',
            'style="cursor:pointer;margin:4px 0 12px;padding:10px 18px;',
            'background:#f8fafc;border:1px dashed #cbd5e1;border-radius:10px;',
            'display:flex;align-items:center;justify-content:space-between;',
            'user-select:none;">',
              '<span style="font-size:13px;font-weight:600;color:#475569;">',
                '<span id="vv-more-btn">&#9660; Show ', n_hidden, ' more variant', if(n_hidden!=1)"s" else "", '</span>',
              '</span>',
              '<span style="font-size:11px;color:#94a3b8;">click to expand</span>',
            '</div>',
            # Hidden cards container
            '<div id="vv-more-cards" style="display:none;">',
            rest_cards,
            '</div>'
          )
        } else ""

        comment_html <- paste0(
          '<div style="background:#f8fafc;border:1px solid #e5e7eb;border-radius:16px;',
          'padding:22px;margin-top:14px;margin-bottom:18px;">',
          '<div style="display:flex;align-items:center;margin-bottom:14px;">',
          '<div style="height:6px;width:72px;border-radius:6px;background:#bc5090;margin-right:14px;"></div>',
          '<h4 style="color:#003f5c;margin:0;font-size:17px;">Variant Summary</h4>',
          '</div>',
          '<p style="font-size:11px;color:#64748b;margin:0 0 12px;">',
          'Cross-referencing input variants with all tracks + dbNSFP scores (MyVariant.info). ',
          'ACMG evidence tags per ',
          '<strong>Richards et al. (2015)</strong> + <strong>Tavtigian et al. (2018)</strong> + <strong>Pejaver et al. (2022)</strong>. ',
          'PP3_strong/PM1_strong shown with ⬆/↑. ',
          'Score verdicts: ',
          '<span style="color:#ef4444;font-weight:700;">dam</span>=damaging &nbsp;',
          '<span style="color:#f59e0b;font-weight:700;">amb</span>=ambiguous &nbsp;',
          '<span style="color:#10b981;font-weight:700;">ben</span>=benign. &nbsp;',
          'Sections: ',
          '<span style="background:#b91c1c;color:#fff;padding:1px 5px;border-radius:3px;font-size:10px;">Pathogenicity</span> ',
          '<span style="background:#6d28d9;color:#fff;padding:1px 5px;border-radius:3px;font-size:10px;">Ensemble</span> ',
          '<span style="background:#047857;color:#fff;padding:1px 5px;border-radius:3px;font-size:10px;">Conservation</span> ',
          '<span style="background:#1d4ed8;color:#fff;padding:1px 5px;border-radius:3px;font-size:10px;">Population</span>. ',
          'Clinical review required.',
          '</p>',
          first_two,
          collapsed_section,
          '</div>'
        )
      }
    }
    HTML(paste0(table_html, comment_html))
  })
  
  # --- Variant Table: Download as TSV ---
  output$variant_table_ready <- reactive({
    vtbl <- prefetch_state$variant_table
    !is.null(vtbl) && is.data.frame(vtbl) && nrow(vtbl) > 0
  })
  outputOptions(output, "variant_table_ready", suspendWhenHidden = FALSE)
  
  output$download_variant_table <- downloadHandler(
    filename = function() {
      gene <- tryCatch(input$gene_name, error = function(e) "variants")
      paste0("VarViz_", gene, "_variant_summary_", format(Sys.Date(), "%Y%m%d"), ".tsv")
    },
    content = function(file) {
      vtbl <- prefetch_state$variant_table
      if (is.null(vtbl) || !is.data.frame(vtbl) || nrow(vtbl) == 0) {
        writeLines("No variant data available.", file)
        return()
      }

      # ── Enrich with per-card cosegregation selections ─────────────────────
      seg_labels <- c(
        none         = "Not assessed",
        pp1          = "1-2 affected relatives (PP1 Supporting)",
        pp1_moderate = "3-4 affected relatives (PP1 Moderate)",
        pp1_strong   = ">=5 affected relatives (PP1 Strong)"
      )
      vtbl$Cosegregation <- sapply(seq_len(nrow(vtbl)), function(i) {
        val <- tryCatch(input[[paste0("seg_", i)]], error = function(e) "none")
        if (is.null(val) || !val %in% names(seg_labels)) val <- "none"
        seg_labels[[val]]
      })
      vtbl$PP1_Applied <- sapply(seq_len(nrow(vtbl)), function(i) {
        val <- tryCatch(input[[paste0("seg_", i)]], error = function(e) "none")
        switch(val, pp1="PP1", pp1_moderate="PP1_moderate", pp1_strong="PP1_strong", "")
      })
      dn_labels <- c(
        not_denovo       = "Not assessed",
        denovo_confirmed = "Confirmed de novo, both parents tested (PS2)",
        denovo_assumed   = "Assumed de novo, parents not tested (PM6)"
      )
      vtbl$DeNovo_Evidence <- sapply(seq_len(nrow(vtbl)), function(i) {
        val <- tryCatch(input[[paste0("dn_", i)]], error = function(e) "not_denovo")
        if (is.null(val) || !val %in% names(dn_labels)) val <- "not_denovo"
        dn_labels[[val]]
      })
      vtbl$DeNovo_Applied <- sapply(seq_len(nrow(vtbl)), function(i) {
        val <- tryCatch(input[[paste0("dn_", i)]], error = function(e) "not_denovo")
        switch(val, denovo_confirmed="PS2", denovo_assumed="PM6", "")
      })
      pm3_labels <- c(
        none         = "Not assessed",
        pm3          = "1 compound het observation (PM3 Supporting)",
        pm3_moderate = "2 compound het observations (PM3 Moderate)",
        pm3_strong   = ">=3 compound het observations (PM3 Strong)"
      )
      vtbl$PM3_Evidence <- sapply(seq_len(nrow(vtbl)), function(i) {
        val <- tryCatch(input[[paste0("pm3_", i)]], error = function(e) "none")
        if (is.null(val) || !val %in% names(pm3_labels)) val <- "none"
        pm3_labels[[val]]
      })
      vtbl$PM3_Applied <- sapply(seq_len(nrow(vtbl)), function(i) {
        val <- tryCatch(input[[paste0("pm3_", i)]], error = function(e) "none")
        switch(val, pm3="PM3", pm3_moderate="PM3_moderate", pm3_strong="PM3_strong", "")
      })

      # Recompute Final_ACMG_Tags and Final_Classification with PP1
      vtbl$Final_ACMG_Tags <- sapply(seq_len(nrow(vtbl)), function(i) {
        base_tags <- if (nchar(as.character(vtbl$ACMG_Tags[i])) > 0)
                       trimws(strsplit(as.character(vtbl$ACMG_Tags[i]), ",")[[1]])
                     else character(0)
        pp1 <- vtbl$PP1_Applied[i]
        if (nchar(pp1) > 0) base_tags <- c(base_tags, pp1)
        pm3 <- vtbl$PM3_Applied[i]
        if (nchar(pm3) > 0) base_tags <- c(base_tags, pm3)
        dn <- vtbl$DeNovo_Applied[i]
        if (nchar(dn) > 0) base_tags <- c(base_tags, dn)
        ptm_tag_export <- if ("PTM_ACMG" %in% colnames(vtbl) &&
                              nchar(as.character(vtbl$PTM_ACMG[i])) > 0 &&
                              vtbl$PTM_ACMG[i] == "PS3_supporting")
                            "PS3_supporting" else ""
        if (nchar(ptm_tag_export) > 0) base_tags <- c(base_tags, ptm_tag_export)
        paste(base_tags, collapse=", ")
      })
      vtbl$Final_Classification <- sapply(seq_len(nrow(vtbl)), function(i) {
        tags <- if (nchar(vtbl$Final_ACMG_Tags[i]) > 0)
                  trimws(strsplit(vtbl$Final_ACMG_Tags[i], ",")[[1]])
                else character(0)
        tryCatch(classify_acmg(tags)$classification, error = function(e) "")
      })
      vtbl$Final_Points <- sapply(seq_len(nrow(vtbl)), function(i) {
        tags <- if (nchar(vtbl$Final_ACMG_Tags[i]) > 0)
                  trimws(strsplit(vtbl$Final_ACMG_Tags[i], ",")[[1]])
                else character(0)
        tryCatch(classify_acmg(tags)$pts, error = function(e) NA_integer_)
      })

      # Export all columns except internal ones
      exclude_cols <- c("dbNSFP_JSON", "SIFT_V", "PP2_HDIV_V", "PP2_HVAR_V", "LRT_V",
                        "FATHMM_V", "PROVEAN_V", "MutTaster_V", "REVEL_V", "MetaSVM_V",
                        "MetaLR_V", "MetaRNN_V", "CADD_V", "DANN_V", "AM_V")
      export_cols <- setdiff(colnames(vtbl), exclude_cols)
      write.table(vtbl[, export_cols, drop = FALSE], file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
    }
  )

  # --- Gene Info: Summary card (left column) ---
  output$geneinfo_summary <- renderText({
    gi <- tryCatch(gene_info_uniprot(), error = function(e) NULL)
    if (is.null(gi)) return("")
    
    # Helper to escape HTML
    esc <- function(x) {
      x <- as.character(x)
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      x <- gsub(">", "&gt;", x, fixed = TRUE)
      x <- gsub('"', "&quot;", x, fixed = TRUE)
      x
    }

    # ClinGen validity badge — uses richer all_assertions from upgraded fetch_clingen_validity()
    clingen_val <- tryCatch(clingen_validity_reactive(), error = function(e) NULL)
    clingen_cls <- if (!is.null(clingen_val) && !is.na(clingen_val$classification) &&
                       nchar(clingen_val$classification) > 0) clingen_val$classification else "Unknown"
    clingen_disease_txt <- if (!is.null(clingen_val) && !is.na(clingen_val$disease) &&
                                nchar(clingen_val$disease) > 0) clingen_val$disease else ""
    clingen_moi_txt     <- if (!is.null(clingen_val) && !is.na(clingen_val$moi) &&
                                nchar(clingen_val$moi) > 0) clingen_val$moi else ""
    clingen_source      <- if (!is.null(clingen_val$source)) clingen_val$source else "ClinGen"

    # Aggregate stats from all_assertions
    all_a <- clingen_val$all_assertions
    n_assertions  <- if (!is.null(all_a) && nrow(all_a) > 0) nrow(all_a) else 0
    n_submitters  <- if (!is.null(all_a) && nrow(all_a) > 0) length(unique(all_a$submitter)) else 0
    n_diseases    <- if (!is.null(all_a) && nrow(all_a) > 0) length(unique(all_a$disease[nchar(all_a$disease)>0])) else 0

    # Best-tier disease list (show up to 3 diseases at the top tier)
    top_diseases <- if (!is.null(all_a) && nrow(all_a) > 0) {
      tier_map <- c("Definitive"=6,"Strong"=5,"Moderate"=4,"Limited"=3,"Disputed"=2,"Refuted"=1)
      all_a$tn <- tier_map[all_a$classification]
      top_tier <- max(all_a$tn, na.rm = TRUE)
      top_dis  <- unique(all_a$disease[!is.na(all_a$tn) & all_a$tn == top_tier & nchar(all_a$disease)>0])
      if (length(top_dis) > 3) paste0(paste(head(top_dis,3), collapse=", "), " +", length(top_dis)-3, " more")
      else if (length(top_dis) > 0) paste(top_dis, collapse=", ")
      else clingen_disease_txt
    } else clingen_disease_txt

    # Badge color by tier
    clingen_badge_style <- switch(clingen_cls,
      "Definitive"  = "background:#d1fae5;color:#065f46;border:1.5px solid #6ee7b7;",
      "Strong"      = "background:#dbeafe;color:#1e40af;border:1.5px solid #93c5fd;",
      "Moderate"    = "background:#fef9c3;color:#713f12;border:1.5px solid #fde68a;",
      "Limited"     = "background:#fef3c7;color:#92400e;border:1.5px solid #fcd34d;",
      "Disputed"    = "background:#fee2e2;color:#991b1b;border:1.5px solid #fca5a5;",
      "Refuted"     = "background:#f3f4f6;color:#374151;border:1.5px solid #d1d5db;",
                      "background:#f1f5f9;color:#64748b;border:1.5px solid #cbd5e1;"
    )
    clingen_icon <- switch(clingen_cls,
      "Definitive" = "&#9989;", "Strong" = "&#128308;", "Moderate" = "&#128992;",
      "Limited" = "&#128993;", "Disputed" = "&#9888;", "Refuted" = "&#10060;", "&#9679;"
    )

    # Build compact multi-disease assertion table (for tooltip / hover detail)
    tooltip_parts <- c(paste0("ClinGen Gene-Disease Validity: ", clingen_cls))
    if (nchar(clingen_moi_txt) > 0) tooltip_parts <- c(tooltip_parts, paste0("MOI: ", clingen_moi_txt))
    if (n_submitters > 0) tooltip_parts <- c(tooltip_parts, paste0(n_assertions, " assertions from ", n_submitters, " submitter(s)"))
    clingen_tooltip <- paste(tooltip_parts, collapse = " | ")

    # Submitter/disease stats pills
    stat_pills <- if (n_assertions > 0) {
      paste0(
        '<span style="font-size:11px;color:#64748b;margin-left:8px;">',
        n_assertions, ' assertion', if(n_assertions!=1)"s" else "", ' &middot; ',
        n_diseases, ' disease', if(n_diseases!=1)"s" else "",
        ' &middot; <span style="font-size:10px;color:#94a3b8;">', esc(clingen_source), '</span>',
        '</span>'
      )
    } else ""

    # Disease name line
    disease_line <- if (nchar(top_diseases) > 0) {
      paste0('<div style="font-size:11px;color:#64748b;margin-top:2px;font-style:italic;">',
             esc(top_diseases),
             if (nchar(clingen_moi_txt) > 0) paste0(' &middot; ', esc(clingen_moi_txt)) else "",
             '</div>')
    } else ""

    # ClinGen search link — uses HGNC ID for direct gene page (e.g. HGNC:950 → BAP1)
    clingen_search_url <- if (!is.null(gi$hgnc_id) && nchar(gi$hgnc_id) > 0) {
      paste0("https://search.clinicalgenome.org/kb/genes/", URLencode(gi$hgnc_id, reserved = TRUE))
    } else if (!is.null(gi$gene_name) && nchar(gi$gene_name) > 0) {
      paste0("https://search.clinicalgenome.org/kb/genes?search=", URLencode(gi$gene_name, reserved = TRUE))
    } else NULL

    vcep_link <- if (!is.null(clingen_search_url) && n_assertions > 0) {
      paste0(' <a href="', clingen_search_url, '" target="_blank" ',
             'style="font-size:10px;color:#0369a1;text-decoration:none;white-space:nowrap;">',
             'ClinGen &#8599;</a>')
    } else ""

    clingen_badge_html <- paste0(
      '<div>',
      '<span title="', esc(clingen_tooltip), '" style="',
      clingen_badge_style,
      'display:inline-flex;align-items:center;gap:5px;',
      'padding:4px 12px;border-radius:20px;font-size:12px;font-weight:700;',
      'cursor:help;white-space:nowrap;">',
      clingen_icon, ' ClinGen: ', esc(clingen_cls),
      '</span>',
      stat_pills,
      vcep_link,
      disease_line,
      '</div>'
    )
    
    paste0(
      '<div style="background:#fff; border:1px solid #e5e7eb; border-radius:16px; padding:22px; margin-bottom:18px;">',
        '<div style="display:flex; align-items:center; margin-bottom:14px;">',
          '<div style="height:6px; width:72px; border-radius:6px; background:#ff6361; margin-right:14px;"></div>',
          '<h3 style="color:#003f5c; margin:0; font-size:22px;">Gene</h3>',
        '</div>',
        '<table style="font-size:15px; margin-bottom:4px; line-height:1.8;">',
          '<tr><td style="padding:4px 16px 4px 0; font-weight:700; color:#58508d; white-space:nowrap;">Gene:</td>',
          '<td class="txt-big" style="padding:4px 0; font-weight:700;">', esc(gi$gene_name), '</td></tr>',
          '<tr><td style="padding:4px 16px 4px 0; font-weight:700; color:#58508d; white-space:nowrap;">UniProt ID:</td>',
          '<td class="txt-small" style="padding:4px 0;"><a href="https://www.uniprot.org/uniprotkb/', esc(gi$uniprot_id), '/entry" target="_blank" style="color:#bc5090; font-weight:600;">', esc(gi$uniprot_id), '</a></td></tr>',
          '<tr><td style="padding:4px 16px 4px 0; font-weight:700; color:#58508d; white-space:nowrap;">Synonyms:</td>',
          '<td class="txt-small" style="padding:4px 0;">', esc(gi$synonyms), '</td></tr>',
          '<tr><td style="padding:4px 16px 4px 0; font-weight:700; color:#58508d; white-space:nowrap;">Protein:</td>',
          '<td class="txt-medium" style="padding:4px 0;">', esc(gi$protein_name), '</td></tr>',
          '<tr><td style="padding:4px 16px 4px 0; font-weight:700; color:#58508d; white-space:nowrap;">Length:</td>',
          '<td class="txt-small" style="padding:4px 0;">', esc(as.character(gi$length)), ' aa</td></tr>',
          if (clingen_cls != "Unknown") paste0(
          '<tr><td style="padding:4px 16px 4px 0; font-weight:700; color:#58508d; white-space:nowrap;">ClinGen:</td>',
          '<td style="padding:4px 0;">', clingen_badge_html, '</td></tr>'
          ) else "",
          {
            gr <- gene_data[gene_data$gene_name == gi$gene_name, , drop=FALSE]
            gv <- if (nrow(gr)>0 && !is.na(gr$GeVIR_pct[1])) gr$GeVIR_pct[1] else NA_real_
            if (!is.na(gv)) {
              gv_col  <- if (gv < 25) "#ef4444" else if (gv > 75) "#3b82f6" else "#64748b"
              gv_bg   <- if (gv < 25) "#fee2e2" else if (gv > 75) "#dbeafe" else "#f1f5f9"
              gv_lbl  <- if (gv < 25) paste0(round(gv,1), "th pct • missense intolerant")
                         else if (gv > 75) paste0(round(gv,1), "th pct • missense tolerant")
                         else paste0(round(gv,1), "th pct")
              gv_tip  <- paste0("GeVIR gene percentile for missense variation intolerance. ",
                                "Low = intolerant (supports PP2). High = tolerant (supports BP1). ",
                                "Source: gevirank.org")
              paste0(
                '<tr><td style="padding:4px 16px 4px 0; font-weight:700; color:#58508d; white-space:nowrap;">GeVIR:</td>',
                '<td style="padding:4px 0;">',
                '<span title="', gv_tip, '" style="background:', gv_bg, ';color:', gv_col, ';',
                'padding:3px 10px;border-radius:12px;font-size:12px;font-weight:700;cursor:help;">',
                gv_lbl, '</span>',
                ' <a href="https://www.gevirank.org/gevir/" target="_blank" ',
                'style="font-size:10px;color:#0369a1;text-decoration:none;">GeVIR &#8599;</a>',
                '</td></tr>'
              )
            } else ""
          },
        '</table>',
      '</div>'
    )
  })

  # --- Gene Info: Details cards (function, disease, links) ---
  output$geneinfo_details <- renderText({
    gi <- tryCatch(gene_info_uniprot(), error = function(e) NULL)
    if (is.null(gi)) return("")
    
    esc <- function(x) {
      x <- as.character(x)
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      x <- gsub(">", "&gt;", x, fixed = TRUE)
      x <- gsub('"', "&quot;", x, fixed = TRUE)
      x
    }
    
    # --- Build External Links ---
    link_html <- paste(sapply(gi$links, function(lk) {
      paste0('<a href="', esc(lk$url), '" target="_blank" ',
             'style="display:inline-block; padding:6px 14px; margin:4px 6px 4px 0;',
             'background:linear-gradient(135deg,#003f5c,#58508d); color:#fff;',
             'border-radius:8px; font-size:13px; font-weight:600;',
             'text-decoration:none;">', esc(lk$name), ' &#8599;</a>')
    }), collapse = "\n")
    
    # --- Build Disease Table (UniProt source) ---
    if (!is.null(gi$diseases) && nrow(gi$diseases) > 0) {
      disease_rows_html <- paste(sapply(seq_len(nrow(gi$diseases)), function(i) {
        row <- gi$diseases[i, ]
        mim_cell <- if (!is.na(row$MIM) && row$MIM != "—" && row$MIM != "") {
          paste0('<a href="https://www.omim.org/entry/', esc(row$MIM), 
                 '" target="_blank" style="color:#bc5090; font-weight:600;">',
                 'MIM:', esc(row$MIM), ' &#8599;</a>')
        } else { "—" }
        
        paste0('<tr>',
               '<td style="padding:8px 10px; border-bottom:1px solid #e5e7eb; font-weight:600; min-width:140px;">', esc(row$Disease), '</td>',
               '<td style="padding:8px 10px; border-bottom:1px solid #e5e7eb; font-size:13px; max-width:280px;">', esc(row$Description), '</td>',
               '<td style="padding:8px 10px; border-bottom:1px solid #e5e7eb; white-space:nowrap;">', mim_cell, '</td>',
               '<td style="padding:8px 10px; border-bottom:1px solid #e5e7eb; font-size:12px; color:#64748b; max-width:300px;">', esc(row$Notes), '</td>',
               '</tr>')
      }), collapse = "\n")
      
      disease_html <- paste0(
        '<p style="font-size:11px;color:#94a3b8;margin:0 0 8px;">Source: UniProt</p>',
        '<div style="overflow-x:auto;">',
        '<table style="width:100%; border-collapse:collapse; font-size:14px;">',
        '<thead><tr style="background:#f1f5f9;">',
        '<th style="padding:10px; text-align:left; color:#003f5c; white-space:nowrap;">Disease (UniProt ID)</th>',
        '<th style="padding:10px; text-align:left; color:#003f5c;">Description</th>',
        '<th style="padding:10px; text-align:left; color:#003f5c; white-space:nowrap;">OMIM</th>',
        '<th style="padding:10px; text-align:left; color:#003f5c;">Clinical Notes</th>',
        '</tr></thead><tbody>', disease_rows_html, '</tbody></table>',
        '</div>')
    } else {
      disease_html <- '<p style="color:#64748b; font-style:italic; font-size:13px;">No disease associations found in UniProt for this gene.</p>'
    }

    # --- ClinGen Gene-Disease Validity section (from all_assertions) ---
    clingen_val2 <- tryCatch(clingen_validity_reactive(), error = function(e) NULL)
    all_a2 <- if (!is.null(clingen_val2)) clingen_val2$all_assertions else NULL

    clingen_assertions_html <- if (!is.null(all_a2) && nrow(all_a2) > 0) {
      tier_order2 <- c("Definitive"=6,"Strong"=5,"Moderate"=4,"Limited"=3,"Disputed"=2,"Refuted"=1)
      all_a2$tn <- tier_order2[all_a2$classification]
      all_a2 <- all_a2[order(-all_a2$tn, all_a2$disease), ]

      badge_cls <- function(cls) {
        sty <- switch(cls,
          "Definitive" = "background:#d1fae5;color:#065f46;border:1px solid #6ee7b7;",
          "Strong"     = "background:#dbeafe;color:#1e40af;border:1px solid #93c5fd;",
          "Moderate"   = "background:#fef9c3;color:#713f12;border:1px solid #fde68a;",
          "Limited"    = "background:#fef3c7;color:#92400e;border:1px solid #fcd34d;",
          "Disputed"   = "background:#fee2e2;color:#991b1b;border:1px solid #fca5a5;",
          "Refuted"    = "background:#f3f4f6;color:#374151;border:1px solid #d1d5db;",
                         "background:#f1f5f9;color:#64748b;border:1px solid #cbd5e1;"
        )
        paste0('<span style="', sty, 'padding:2px 8px;border-radius:10px;font-size:11px;font-weight:700;white-space:nowrap;">', esc(cls), '</span>')
      }

      rows_html <- paste(sapply(seq_len(nrow(all_a2)), function(i) {
        row <- all_a2[i,]
        paste0('<tr>',
               '<td style="padding:7px 10px;border-bottom:1px solid #f1f5f9;font-size:13px;font-weight:600;">', esc(row$disease), '</td>',
               '<td style="padding:7px 10px;border-bottom:1px solid #f1f5f9;">', badge_cls(row$classification), '</td>',
               '<td style="padding:7px 10px;border-bottom:1px solid #f1f5f9;font-size:12px;color:#64748b;">', esc(row$moi), '</td>',
               '<td style="padding:7px 10px;border-bottom:1px solid #f1f5f9;font-size:12px;color:#94a3b8;">', esc(row$submitter), '</td>',
               '</tr>')
      }), collapse="\n")

      paste0(
        '<p style="font-size:11px;color:#94a3b8;margin:0 0 8px;">Source: ', esc(clingen_val2$source), ' · ',
        nrow(all_a2), ' assertions</p>',
        '<div style="overflow-x:auto;">',
        '<table style="width:100%;border-collapse:collapse;font-size:14px;">',
        '<thead><tr style="background:#f8fafc;">',
        '<th style="padding:8px 10px;text-align:left;color:#003f5c;">Disease</th>',
        '<th style="padding:8px 10px;text-align:left;color:#003f5c;">Validity</th>',
        '<th style="padding:8px 10px;text-align:left;color:#003f5c;">MOI</th>',
        '<th style="padding:8px 10px;text-align:left;color:#003f5c;">Submitter</th>',
        '</tr></thead><tbody>', rows_html, '</tbody></table>',
        '</div>'
      )
    } else {
      '<p style="color:#94a3b8;font-size:12px;font-style:italic;">ClinGen validity data unavailable — check ClinGen directly.</p>'
    }
    
    paste0(
      '<div style="background:#fff; border:1px solid #e5e7eb; border-radius:16px; padding:22px; margin-bottom:18px;">',
        '<div style="display:flex; align-items:center; margin-bottom:14px;">',
          '<div style="height:6px; width:72px; border-radius:6px; background:#ffa600; margin-right:14px;"></div>',
          '<h4 style="color:#003f5c; margin:0; font-size:18px;">Function</h4>',
        '</div>',
        '<p style="font-size:14px; line-height:1.7; color:#1e293b;">', esc(gi$function_text), '</p>',
        '<p style="font-size:12px; color:#94a3b8; margin-top:8px; margin-bottom:0;">Source: UniProt</p>',
      '</div>',
      
      '<div style="background:#fff; border:1px solid #e5e7eb; border-radius:16px; padding:22px; margin-bottom:18px;">',
        '<div style="display:flex; align-items:center; margin-bottom:14px;">',
          '<div style="height:6px; width:72px; border-radius:6px; background:#003f5c; margin-right:14px;"></div>',
          '<h4 style="color:#003f5c; margin:0; font-size:18px;">Involvement in Disease</h4>',
        '</div>',
        disease_html,
      '</div>',

      '<div style="background:#fff; border:1px solid #e5e7eb; border-radius:16px; padding:22px; margin-bottom:18px;">',
        '<div style="display:flex; align-items:center; justify-content:space-between; margin-bottom:14px;">',
          '<div style="display:flex; align-items:center;">',
            '<div style="height:6px; width:72px; border-radius:6px; background:#2563eb; margin-right:14px;"></div>',
            '<h4 style="color:#003f5c; margin:0; font-size:18px;">ClinGen Gene-Disease Validity</h4>',
          '</div>',
          '<a href="',
          if (!is.null(gi$hgnc_id) && nchar(gi$hgnc_id) > 0)
            paste0('https://search.clinicalgenome.org/kb/genes/', esc(gi$hgnc_id))
          else
            paste0('https://search.clinicalgenome.org/kb/genes?search=', esc(isolate(input$gene_name))),
          '" target="_blank" ',
             'style="font-size:12px;color:#0369a1;text-decoration:none;padding:3px 10px;',
             'border:1px solid #bfdbfe;border-radius:6px;white-space:nowrap;">',
             esc(isolate(input$gene_name)), ' on ClinGen &#8599;</a>',
        '</div>',
        clingen_assertions_html,
      '</div>'
    )
  })

  output$external_links_section <- renderUI({
    if (input$goButton == 0) return(NULL)
    gi <- tryCatch(gene_info_uniprot(), error = function(e) NULL)
    if (is.null(gi) || is.null(gi$links) || length(gi$links) == 0) return(NULL)
    esc <- function(x) {
      x <- as.character(x)
      x <- gsub("&", "&amp;", x, fixed = TRUE)
      x <- gsub("<", "&lt;", x, fixed = TRUE)
      x <- gsub(">", "&gt;", x, fixed = TRUE)
      x <- gsub('"', "&quot;", x, fixed = TRUE)
      x
    }
    link_html <- paste(sapply(gi$links, function(lk) {
      paste0('<a href="', esc(lk$url), '" target="_blank" ',
             'style="display:inline-block; padding:6px 14px; margin:4px 6px 4px 0;',
             'background:linear-gradient(135deg,#003f5c,#58508d); color:#fff;',
             'border-radius:8px; font-size:13px; font-weight:600;',
             'text-decoration:none;">', esc(lk$name), ' &#8599;</a>')
    }), collapse = "\n")
    HTML(paste0(
      '<div style="background:#fff; border:1px solid #e5e7eb; border-radius:16px; padding:22px; margin-bottom:18px;">',
        '<div style="display:flex; align-items:center; margin-bottom:14px;">',
          '<div style="height:6px; width:72px; border-radius:6px; background:#58508d; margin-right:14px;"></div>',
          '<h4 style="color:#003f5c; margin:0; font-size:18px;">External Links</h4>',
        '</div>',
        '<div style="display:flex; flex-wrap:wrap; gap:4px;">', link_html, '</div>',
      '</div>'
    ))
  })

  # --- NCBI gene summary for wordcloud ---
  get_ncbi_gene_summary <- function(gene_symbol) {
    cached <- cache_get("ncbi_summary", gene_symbol)
    if (!is.null(cached)) { message("[NCBI Summary] Cache hit for: ", gene_symbol); return(cached) }
    tryCatch({
      search_term <- URLencode(paste0(gene_symbol, "[Gene Name] AND human[Organism]"))
      search_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gene&term=",
                           search_term, "&retmode=json")
      search_res <- jsonlite::fromJSON(httr::content(httr::GET(search_url), "text", encoding = "UTF-8"))
      gene_id <- search_res$esearchresult$idlist[1]
      if (is.null(gene_id)) return("")
      summary_url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gene&id=",
                            gene_id, "&retmode=json")
      summary_res <- jsonlite::fromJSON(httr::content(httr::GET(summary_url), "text", encoding = "UTF-8"))
      result <- summary_res$result[[gene_id]]$summary
      cache_set("ncbi_summary", gene_symbol, result)
      return(result)
    }, error = function(e) return(""))
  }

  # Wordcloud plot
  output$wordcloud_plot <- renderPlot({
    gi <- tryCatch(gene_info_uniprot(), error = function(e) NULL)
    if (is.null(gi)) {
      plot.new(); text(0.5, 0.5, "Click Go to generate", col = "#94a3b8", cex = 1.4)
      return()
    }

    # Combine UniProt text + NCBI summary
    uniprot_text <- paste(gi$function_text, collapse = " ")
    if (!is.null(gi$diseases) && nrow(gi$diseases) > 0) {
      uniprot_text <- paste(uniprot_text,
                            paste(gi$diseases$Disease, gi$diseases$Description, collapse = " "))
    }
    ncbi_text <- get_ncbi_gene_summary(gi$gene_name)
    combined_text <- paste(uniprot_text, ncbi_text)

    # Clean
    clean_text <- gsub("\\(PubMed:[^)]+\\)", "", combined_text)
    clean_text <- gsub("\\[[^]]*\\]", "", clean_text)
    clean_text <- gsub("[^[:alpha:][:space:]]", " ", clean_text)

    if (nchar(trimws(clean_text)) < 10) {
      plot.new(); text(0.5, 0.5, "Not enough text for word cloud", col = "#94a3b8", cex = 1.4)
      return()
    }

    docs <- tm::Corpus(tm::VectorSource(clean_text))
    docs <- tm::tm_map(docs, tm::content_transformer(tolower))
    docs <- tm::tm_map(docs, tm::removeNumbers)
    docs <- tm::tm_map(docs, tm::removeWords, tm::stopwords("english"))
    extra_stops <- c("gene", "protein", "encoded", "human", "results", "associated",
                     "refseq", "provided", "also", "may", "can", "two", "one",
                     "role", "involved", "required", "member", "family", "known")
    docs <- tm::tm_map(docs, tm::removeWords, extra_stops)
    docs <- tm::tm_map(docs, tm::removePunctuation)
    docs <- tm::tm_map(docs, tm::stripWhitespace)

    par(mar = c(0, 0, 0, 0))
    tryCatch({
      wordcloud::wordcloud(
        words = docs,
        min.freq = 1,
        max.words = 80,
        random.order = FALSE,
        rot.per = 0.35,
        colors = c("#003f5c", "#58508d", "#bc5090", "#ff6361", "#ffa600"),
        scale = c(3, 0.4)
      )
    }, error = function(e) {
      plot.new(); text(0.5, 0.5, "Could not generate word cloud", col = "#94a3b8", cex = 1.4)
    })
  })

  output$debug_cols <- renderPrint({
   colnames(consurf_score())
   head(consurf_score(), 5)
  })
   
   
  
  final_plot <- reactive({
    
    if (input$goButton == 0)
      return()
    
    withProgress(message = "Loading gene data...", value = 0, {
    
    incProgress(0.05, detail = "Fetching ClinVar variants...")
    gene_clinvar_data = gene_clinvar_data()
    if (!is.null(gene_clinvar_data) && nrow(gene_clinvar_data) > 0 && "clinvar_goldstar" %in% colnames(gene_clinvar_data)) {
      if (input$clinvar_filter == "1 Star or more"){
        gene_clinvar_data_filtered = gene_clinvar_data[as.numeric(gene_clinvar_data$clinvar_goldstar) > 0, , drop = FALSE]
      } else if (input$clinvar_filter == "2 Stars or more") {
        gene_clinvar_data_filtered = gene_clinvar_data[as.numeric(gene_clinvar_data$clinvar_goldstar) > 1, , drop = FALSE]
      } else {
        gene_clinvar_data_filtered = gene_clinvar_data
      }
    } else {
      gene_clinvar_data_filtered = gene_clinvar_data
    }
    incProgress(0.15, detail = paste0("ClinVar: ", nrow(gene_clinvar_data), " variants"))
    
    incProgress(0.05, detail = "Fetching gnomAD variants...")
    gene_gnomad_data=gene_gnomad_data()
    incProgress(0.15, detail = paste0("gnomAD: ", nrow(gene_gnomad_data), " variants"))
    
    incProgress(0.05, detail = "Fetching UniProt domains & features...")
    pfam_data=pfam_data()
    
    # Derive PTM data from uniprot_data (single source of truth)
    ud <- uniprot_data()
    ptm_types <- c("mod_res", "lipid", "carbohyd", "crosslnk", "disulfid")
    if (!is.null(ud) && nrow(ud) > 0) {
      ptm_rows <- ud[ud$type %in% ptm_types, , drop = FALSE]
      if (nrow(ptm_rows) > 0) {
        gene_ptm_data <- data.frame(
          uniprot_id = ptm_rows$uniprot_id,
          location = ptm_rows$start,
          final_ptm_group = ifelse(!is.na(ptm_rows$mod_res_group), ptm_rows$mod_res_group,
                                   dplyr::case_when(
                                     ptm_rows$type == "disulfid" ~ "Disulfide bond",
                                     ptm_rows$type == "lipid"    ~ "Lipidation",
                                     ptm_rows$type == "carbohyd" ~ "Glycosylation",
                                     ptm_rows$type == "crosslnk" ~ "Cross-link",
                                     TRUE ~ "Other modification"
                                   )),
          stringsAsFactors = FALSE
        )
        # For disulfide bonds, also add the end position
        ds_rows <- ptm_rows[ptm_rows$type == "disulfid" & !is.na(ptm_rows$end) & ptm_rows$end != ptm_rows$start, ]
        if (nrow(ds_rows) > 0) {
          ds_end <- data.frame(
            uniprot_id = ds_rows$uniprot_id,
            location = ds_rows$end,
            final_ptm_group = "Disulfide bond",
            stringsAsFactors = FALSE
          )
          gene_ptm_data <- rbind(gene_ptm_data, ds_end)
        }
      } else {
        gene_ptm_data <- data.frame(uniprot_id = character(), location = numeric(), final_ptm_group = character(), stringsAsFactors = FALSE)
      }
    } else {
      gene_ptm_data <- data.frame(uniprot_id = character(), location = numeric(), final_ptm_group = character(), stringsAsFactors = FALSE)
    }
    incProgress(0.1, detail = paste0("UniProt features: ", nrow(ud), " | PTMs: ", nrow(gene_ptm_data)))
    
    # Generate all plots including AlphaFold
    incProgress(0.05, detail = "Fetching AlphaFold data...")
    # Build named REVEL score vector: position -> score (for lollipop height)
    revel_vec <- tryCatch({
      vtbl <- prefetch_state$variant_table
      if (is.data.frame(vtbl) && "Position" %in% colnames(vtbl) && "REVEL_Score" %in% colnames(vtbl)) {
        rv <- suppressWarnings(as.numeric(vtbl$REVEL_Score))
        names(rv) <- as.character(vtbl$Position)
        rv[!is.na(rv)]
      } else NULL
    }, error = function(e) NULL)

    p1 <- pfamplot(pfam_data(),uniprot_data(),gene_clinvar_data,highlight(),input$label,
                   for_plotly=TRUE,  revel_scores=revel_vec)
    p1_static <- pfamplot(pfam_data(),uniprot_data(),gene_clinvar_data,highlight(),input$label,
                          for_plotly=FALSE, revel_scores=revel_vec)
    p2 <- plot_afmps(mean_data(), highlight(), prot_length = pfam_data()$sequence$length) 
    p6 <- plot_pLDDT(af(), highlight(), prot_length = pfam_data()$sequence$length)
    incProgress(0.1, detail = "Building ClinVar/PTM/CCRS track...")
    p7 <- clinvar_ccrsplot(pfam_data(),uniprot_data(),gene_clinvar_data,gene_ptm_data,gene_ccrs_data(),highlight=highlight())
    message("[Plot] clinvar_ccrsplot generated with ", nrow(gene_ptm_data), " PTM sites")
    
    incProgress(0.1, detail = "Building gnomAD & density plots...")
    if( (nrow(gene_gnomad_data)>0) ) {
      message("[gnomAD freqplot] Using AF cutoff = ", signif(data()[[1]], 6), " | AC cutoff = ", data()[[2]])
      p5 <- gnomad_freqplot(input$gene_name, data()[[1]], highlight(), prot_length = pfam_data()$sequence$length)
    } else {
      p5 <- ggplot() + ggplot2::annotate(geom = "text", x = 0.5, y = 0.5, 
              label = "gnomAD data not available — API may be temporarily overloaded, try again", 
              size = 4.5, color = "#94a3b8", fontface = "italic") +
            labs(y = "gnomAD Freq") +
            theme_bw(base_size = vv_medium) + theme_vv() +
            theme(axis.title.x = element_blank(), axis.text = element_blank(), 
                  axis.ticks = element_blank(), panel.grid = element_blank())
    }
    
    if( (nrow(gene_gnomad_data)>0) ) {
      # Safely subset ClinVar to missense only
      cv_missense <- if (!is.null(gene_clinvar_data_filtered) && is.data.frame(gene_clinvar_data_filtered) && 
                         nrow(gene_clinvar_data_filtered) > 0 && "type" %in% colnames(gene_clinvar_data_filtered)) {
        gene_clinvar_data_filtered[gene_clinvar_data_filtered$type=="missense_variant", , drop = FALSE]
      } else {
        data.frame(prot_pos = numeric(), stringsAsFactors = FALSE)
      }
      p4 <- densityplot(gene_gnomad_data, cv_missense, pfam_data$sequence$length, data()[[2]], highlight(), gene_ptm_data, gene_ccrs_data(), user_path_variants())
    } else {
      p4 <- ggplot() + ggplot2::annotate(geom = "text", x = 0.5, y = 0.5, 
              label = "Mutation density data not available — requires gnomAD data", 
              size = 4.5, color = "#94a3b8", fontface = "italic") +
            labs(y = "Mutation Density") +
            theme_bw(base_size = vv_medium) + theme_vv() +
            theme(axis.title.x = element_blank(), axis.text = element_blank(), 
                  axis.ticks = element_blank(), panel.grid = element_blank())
    }
    
    cs <- tryCatch(consurf_score(), error = function(e) data.frame())
    if (is.data.frame(cs) && nrow(cs) > 0) {
      p3 <- conservplot(cs, pfam_data$sequence$length, highlight())
    } else {
      cs_msg <- if (!is.null(input$consurf_file)) {
        "Could not parse uploaded ConSurf file.\nCheck format and try again."
      } else {
        uid <- tryCatch(as.character(gene_attrib()$uniprot_id[1]), error = function(e) "")
        paste0("No ConSurf-DB entry found for ", input$gene_name,
               if (nchar(uid) > 0) paste0(" (", uid, ")") else "",
               ".\nNo PDB structure available or gene not in ConSurf-DB.\n",
               "Upload a ConSurf grades file manually to view this track.")
      }
      p3 <- ggplot() +
            ggplot2::annotate(geom = "text", x = 0.5, y = 0.5,
              label = cs_msg, size = 3.8, color = "#ef4444",
              fontface = "italic", hjust = 0.5, vjust = 0.5) +
            labs(y = "Conservation") +
            theme_bw(base_size = vv_medium) + theme_vv() +
            theme(axis.title.x = element_blank(), axis.text = element_blank(),
                  axis.ticks = element_blank(), panel.grid = element_blank())
    }

    # Build p8: Multi-Conservation track (UCSC: PhyloP100V, PhyloP470M, PhastCons, GERP++)
    # Only built when checkbox is selected — avoids firing UCSC fetch otherwise
    # p8        = plotly widget  (for interactive HPlot via subplot)
    # p8_static = ggplot object  (for static SPlot via cowplot — htmlwidgets not allowed)
    p8 <- NULL
    p8_static <- NULL
    if ("multiconservation" %in% input$plotselection) {
      cons_data <- tryCatch(conservation_scores(), error = function(e) {
        message("[MultiCons] reactive error: ", e$message); NULL
      })
      prot_len_mc <- tryCatch(pfam_data()$sequence$length, error = function(e) 600)
      p8 <- plot_multiconservation(cons_data, prot_len_mc, highlight())

      # Static ggplot version for PDF/PNG download — plotly widgets cannot be
      # passed to cowplot::plot_grid(), so we build a matching ggplot version
      p8_static <- plot_multiconservation_static(cons_data, prot_len_mc, highlight())
    }

    # Build plot list - AlphaFold always at top, then selected plots
    plotlist <- list(p6)  # Start with AlphaFold pLDDT (always shown)
    
    # Add selected plots based on checkbox inputs
    if("freq" %in% input$plotselection) {
      plotlist <- c(plotlist, list(p5))  # gnomAD Frequency
    }
    
    if("density" %in% input$plotselection) {
      plotlist <- c(plotlist, list(p4))  # Mutation Density
    }
    
    if("clinvar" %in% input$plotselection) {
      plotlist <- c(plotlist, list(p2))  # AF Mean Pathogenicity
    }
    
    if("clinvar_ptm" %in% input$plotselection) {
      plotlist <- c(plotlist, list(p7))  # ClinVar / PTM / CCRS
    }
    
    if (is.data.frame(cs) && nrow(cs) > 0) {
      plotlist <- c(plotlist, list(p3))  # ConSurf (auto-shown when file uploaded)
    }

    # Multi-conservation track — inserted between CCRS and UniProt Domains
    if("multiconservation" %in% input$plotselection && !is.null(p8)) {
      plotlist <- c(plotlist, list(p8))  # PhyloP / PhastCons / GERP++ from UCSC
    }
    
    # Always include the Uniprot Domains plot at bottom
    plotlist <- c(plotlist, list(p1))
    
    # Remove any NULL plots
    plotlist <- plotlist[!sapply(plotlist, is.null)]
    
    # Create final plot — equal panel heights like original working version
    HPlot <- plotly::subplot(
      plotlist,
      nrows  = length(plotlist),
      shareX = TRUE,
      titleY = TRUE,
      margin = c(0.02, 0.02, 0.02, 0.02)
    )
    
    # Plotly variant labels via layout(annotations) — matches old working approach.
    # Stems/dots are drawn by ggplot (ggplotly converts them correctly).
    # Only labels need special treatment because ggplotly ignores angle=.
    if (input$label == "Yes" && nrow(highlight()) > 0) {
      local({
        h <- tryCatch(highlight(), error = function(e) data.frame())
        if (nrow(h) == 0) return(invisible(NULL))
        prot_len <- tryCatch(pfam_data()$sequence$length, error = function(e) 600L)
        h <- h[!is.na(h$prot_pos) & h$prot_pos >= 0 & h$prot_pos <= prot_len, , drop=FALSE]
        if (nrow(h) == 0) return(invisible(NULL))

        n_pan      <- length(plotlist)
        y_axis_ref <- if (n_pan == 1) "y" else paste0("y", n_pan)

        aa_col <- c(D="#e04040",E="#e04040",H="#4472c4",K="#4472c4",R="#4472c4",
                    S="#2ca02c",T="#2ca02c",N="#2ca02c",Q="#2ca02c",C="#2ca02c",
                    F="#e67e22",W="#e67e22",Y="#e67e22",
                    A="#8e44ad",V="#8e44ad",I="#8e44ad",L="#8e44ad",M="#8e44ad",G="#8e44ad",
                    P="#8B4513")
        bare       <- sub("^p\\.", "", h$Mutation)
        var_aa     <- toupper(substr(bare, nchar(bare), nchar(bare)))
        h$lbl_col  <- ifelse(var_aa %in% names(aa_col), aa_col[var_aa], "#555555")

        # Stagger label y-levels (same as old script)
        stagger_levels <- c(3.85, 4.05, 4.22, 4.36, 4.42, 4.46)
        min_gap <- prot_len * 0.05
        h <- h[order(h$prot_pos), , drop=FALSE]
        h$label_y    <- stagger_levels[1]
        h$stagger_idx <- 1L
        if (nrow(h) > 1) {
          for (j in 2:nrow(h)) {
            gap <- h$prot_pos[j] - h$prot_pos[j-1]
            if (gap < min_gap) {
              h$stagger_idx[j] <- (h$stagger_idx[j-1] %% length(stagger_levels)) + 1L
            } else {
              h$stagger_idx[j] <- 1L
            }
            h$label_y[j] <- stagger_levels[h$stagger_idx[j]]
          }
        }

        annotations <- lapply(seq_len(nrow(h)), function(i) {
          list(x=h$prot_pos[i], y=h$label_y[i],
               text=paste0("<b>", as.character(h$Mutation[i]), "</b>"),
               xref="x", yref=y_axis_ref,
               showarrow=FALSE, textangle=-65,
               font=list(size=8, family="Times New Roman", color=h$lbl_col[i]),
               xanchor="left", yanchor="bottom")
        })

        # Protein name (caption not rendered by ggplotly)
        HPlot <<- HPlot %>% plotly::layout(
          annotations = annotations
        )
      })
    }
    
    incProgress(0.15, detail = "Assembling final visualization...")
    
    # Also create static version for download (with ggplot labels instead of plotly annotations)
    # IMPORTANT: plotly htmlwidgets cannot go into cowplot::plot_grid().
    # Rebuild the static plotlist using only ggplot objects:
    #   - p8 (plotly) is replaced with p8_static (ggplot)
    #   - p1 (plotly domain bar) is replaced with p1_static (ggplot)
    static_plotlist <- list(p6)  # pLDDT (ggplot)
    if ("freq"       %in% input$plotselection) static_plotlist <- c(static_plotlist, list(p5))
    if ("density"    %in% input$plotselection) static_plotlist <- c(static_plotlist, list(p4))
    if ("clinvar"    %in% input$plotselection) static_plotlist <- c(static_plotlist, list(p2))
    if ("clinvar_ptm" %in% input$plotselection) static_plotlist <- c(static_plotlist, list(p7))
    if (is.data.frame(cs) && nrow(cs) > 0)
      static_plotlist <- c(static_plotlist, list(p3))
    if ("multiconservation" %in% input$plotselection && !is.null(p8_static))
      static_plotlist <- c(static_plotlist, list(p8_static))  # ggplot, not plotly
    static_plotlist <- c(static_plotlist, list(p1_static))   # domain bar static version
    static_plotlist <- static_plotlist[!sapply(static_plotlist, is.null)]
    aligned <- cowplot::align_plots(plotlist = static_plotlist, align = "v", axis = "lr")
    SPlot <- plot_grid(plotlist = aligned, ncol = 1)
    
    incProgress(0.05, detail = "Done!")
    
    list(HTMLversion = HPlot, GGversion = SPlot, n_panels = length(plotlist))
    
    })  # end withProgress
  })
  
  output$mplot  <- renderPlotly({
    fp <- final_plot()
    p  <- fp$HTMLversion
    p %>% plotly::config(
      displayModeBar = TRUE,
      modeBarButtonsToAdd = list("resetScale2d"),
      modeBarButtonsToRemove = list("lasso2d", "select2d"),
      displaylogo = FALSE
    )
  })

  # Dynamic container height: non-domain panels ~130px each, domain panel ~280px
  output$mplot_container <- renderUI({
    fp <- tryCatch(final_plot(), error = function(e) NULL)
    n_pan <- if (!is.null(fp)) fp$n_panels else 7L
    px_h  <- (n_pan - 1L) * 130L + 280L
    withSpinner(
      plotlyOutput("mplot", width="100%", height=paste0(px_h, "px")),
      type=6
    )
  })
  
  highlight_data <- eventReactive(input$goButton,{ 
    highlight()
  })
  
  output$highlight <- DT::renderDataTable( highlight_data(), options = list(scrollX = TRUE))
  
  # ── gnomAD raw data download ──────────────────────────────────────────────
  output$download_gnomad_raw <- downloadHandler(
    filename = function() {
      gene  <- tryCatch(input$gene_name, error = function(e) "gene")
      paste0("VarViz_", gene, "_gnomAD_raw_", format(Sys.Date(), "%Y%m%d"), ".tsv")
    },
    content = function(file) {
      gd <- tryCatch(gene_gnomad_data(), error = function(e) NULL)
      if (is.null(gd) || nrow(gd) == 0) {
        write.table(data.frame(message = "No gnomAD data available for this gene."),
                    file, sep = "\t", row.names = FALSE, quote = FALSE)
        return()
      }
      # Re-fetch raw API JSON to get variant_id and all fields
      gene_name  <- tryCatch(input$gene_name, error = function(e) "")
      json_raw   <- tryCatch(query_gnomad_api(gene_name), error = function(e) NULL)
      if (is.null(json_raw) || length(json_raw$data$gene$variants) == 0) {
        # Fall back to already-parsed gene_gnomad_data
        out <- gd
        out$Gene <- gene_name
        out$Source <- "VarViz_parsed"
        write.table(out, file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
        return()
      }
      # Parse all fields from raw JSON including variant_id and both exome+genome
      variants <- json_raw$data$gene$variants
      pop_val <- function(pops, pop_id, field) {
        if (is.null(pops) || length(pops) == 0) return(NA_real_)
        for (p in pops) {
          if (identical(p$id, pop_id)) {
            v <- p[[field]]
            return(if (!is.null(v)) as.numeric(v) else NA_real_)
          }
        }
        NA_real_
      }
      rows <- lapply(variants, function(v) {
        ex <- v$exome; ge <- v$genome
        ex_pops <- if (!is.null(ex$populations)) ex$populations else list()
        ge_pops <- if (!is.null(ge$populations)) ge$populations else list()
        data.frame(
          Gene            = gene_name,
          Variant_ID      = if (!is.null(v$variant_id)) v$variant_id else NA_character_,
          HGVSp           = if (!is.null(v$hgvsp))     v$hgvsp      else NA_character_,
          Exome_AC        = if (!is.null(ex$ac))  as.integer(ex$ac)  else NA_integer_,
          Exome_AN        = if (!is.null(ex$an))  as.integer(ex$an)  else NA_integer_,
          Exome_AF        = if (!is.null(ex$af))  as.numeric(ex$af)  else NA_real_,
          Exome_Nhomalt   = if (!is.null(ex$ac_hom)) as.integer(ex$ac_hom) else NA_integer_,
          Genome_AC       = if (!is.null(ge$ac))  as.integer(ge$ac)  else NA_integer_,
          Genome_AN       = if (!is.null(ge$an))  as.integer(ge$an)  else NA_integer_,
          Genome_AF       = if (!is.null(ge$af))  as.numeric(ge$af)  else NA_real_,
          Genome_Nhomalt  = if (!is.null(ge$ac_hom)) as.integer(ge$ac_hom) else NA_integer_,
          Exome_AF_AFR    = pop_val(ex_pops, "afr", "af"),
          Exome_AF_NFE    = pop_val(ex_pops, "nfe", "af"),
          Exome_AF_EAS    = pop_val(ex_pops, "eas", "af"),
          Exome_AF_SAS    = pop_val(ex_pops, "sas", "af"),
          Exome_AF_AMR    = pop_val(ex_pops, "amr", "af"),
          Exome_AF_FIN    = pop_val(ex_pops, "fin", "af"),
          Genome_AF_AFR   = pop_val(ge_pops, "afr", "af"),
          Genome_AF_NFE   = pop_val(ge_pops, "nfe", "af"),
          Genome_AF_EAS   = pop_val(ge_pops, "eas", "af"),
          Genome_AF_SAS   = pop_val(ge_pops, "sas", "af"),
          Genome_AF_AMR   = pop_val(ge_pops, "amr", "af"),
          Genome_AF_FIN   = pop_val(ge_pops, "fin", "af"),
          stringsAsFactors = FALSE
        )
      })
      out <- do.call(rbind, rows)
      write.table(out, file, sep = "\t", row.names = FALSE, quote = FALSE, na = "")
      message("[gnomAD download] Wrote ", nrow(out), " variants for ", gene_name)
    }
  )

  # downloadHandler 
  output$down <- downloadHandler(
    filename =  function() {
      paste("VarViz",input$gene_name, input$format, sep=".")
    },
    content = function(file) {
      # Count number of subplot panels for height scaling
      n_panels <- final_plot()$n_panels
      plot_height <- max(8, 2.5 * n_panels)
      plot_width  <- 12
      
      if(input$format == "png")
        png(file, width = plot_width, height = plot_height, units = "in", res = 300) else if(input$format == "jpeg")
        jpeg(file, width = plot_width, height = plot_height, units = "in", res = 300, quality = 95) else
        pdf(file, width = plot_width, height = plot_height)
      
      # Use cowplot::ggdraw to wrap the composite plot and add a border
      bordered <- cowplot::ggdraw() +
        cowplot::draw_plot(final_plot()$GGversion, x = 0.01, y = 0.01, width = 0.98, height = 0.98) +
        ggplot2::theme(
          plot.background = ggplot2::element_rect(color = "#333333", linewidth = 1.2, fill = NA),
          plot.margin = ggplot2::margin(4, 4, 4, 4, "mm")
        )
      
      print(bordered)
      dev.off()
    } 
  )
  
})
