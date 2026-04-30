# local_predictors.R — bulk REVEL + AlphaMissense lookups for the dual-pass
# classification harness (analyses/05_classify_harness.R, Option C Task 2).
#
# Replaces per-variant MyVariant API calls (~1.5 s each, ~30 min per gene) with
# in-memory data.table joins (~150 MB per chrom, microsecond lookups). The
# harness pre-loads one (REVEL chrom, AlphaMissense UID) pair per gene at the
# start of classify_gene(), then issues binary-search lookups during the
# variant loop.
#
# Background data, all gitignored under analyses/raw/:
#   REVEL v1.3 (May 2021):  chr,hg19_pos,grch38_pos,ref,alt,aaref,aaalt,REVEL,Ensembl_transcriptid
#                           split per-chrom at analyses/raw/revel/by_chrom/chr_<N>.csv
#   AlphaMissense (EBI):    protein_variant,am_pathogenicity,am_class
#                           one CSV per UniProt accession at analyses/raw/alphamissense/<UID>-F1-aa-substitutions.csv
#
# Reference: analyses/OPTION_C_HANDOFF.md §2.

suppressMessages(library(data.table))

# ---------------------------------------------------------------------------
# AlphaMissense
# ---------------------------------------------------------------------------

load_alphamissense_csv <- function(uniprot_id,
                                   am_dir = "analyses/raw/alphamissense") {
  path <- file.path(am_dir, paste0(uniprot_id, "-F1-aa-substitutions.csv"))
  if (!file.exists(path)) return(NULL)
  dt <- data.table::fread(path, showProgress = FALSE)
  data.table::setkey(dt, protein_variant)
  dt
}

lookup_alphamissense <- function(am_dt, p_notation) {
  na_out <- list(score = NA_real_, class = NA_character_)
  if (is.null(am_dt) || length(p_notation) == 0L || is.na(p_notation)) {
    return(na_out)
  }
  key <- sub("^p\\.", "", as.character(p_notation))
  hit <- am_dt[.(key), on = "protein_variant", nomatch = NA]
  if (nrow(hit) == 0L || is.na(hit$am_pathogenicity[1])) return(na_out)
  list(score = as.numeric(hit$am_pathogenicity[1]),
       class = as.character(hit$am_class[1]))
}

# ---------------------------------------------------------------------------
# REVEL
# ---------------------------------------------------------------------------

load_revel_for_chrom <- function(chrom,
                                 revel_dir = "analyses/raw/revel/by_chrom") {
  path <- file.path(revel_dir, paste0("chr_", chrom, ".csv"))
  if (!file.exists(path)) {
    stop(sprintf("REVEL chrom file not found: %s", path), call. = FALSE)
  }
  # REVEL marks variants that fail to lift over to GRCh38 with grch38_pos="."
  # — coerce to NA so fread keeps the column integer, then drop those rows.
  dt <- data.table::fread(path, showProgress = FALSE, na.strings = c("", "NA", "."))
  dt <- dt[!is.na(grch38_pos)]
  # The bulk file emits one row per Ensembl transcript context; REVEL itself
  # is a position-level score so collapsing on (pos, ref, alt) is safe.
  dt <- unique(dt, by = c("grch38_pos", "ref", "alt"))
  data.table::setkey(dt, grch38_pos, ref, alt)
  attr(dt, "revel_chrom") <- as.character(chrom)
  dt
}

lookup_revel <- function(revel_dt, chrom, pos, ref, alt) {
  dt_chrom <- attr(revel_dt, "revel_chrom")
  if (!is.null(dt_chrom) && as.character(chrom) != dt_chrom) {
    stop(sprintf("revel_dt is for chromosome %s but lookup requested %s",
                 dt_chrom, chrom), call. = FALSE)
  }
  hit <- revel_dt[.(as.integer(pos), as.character(ref), as.character(alt)),
                  on = c("grch38_pos", "ref", "alt"),
                  nomatch = NA]
  if (nrow(hit) == 0L || is.na(hit$REVEL[1])) return(NA_real_)
  as.numeric(hit$REVEL[1])
}
