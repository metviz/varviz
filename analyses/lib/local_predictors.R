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

# Find REVEL for a missense at amino-acid position whose codon spans
# `codon_positions` (3 ints — caller computes via aa_to_genomic + strand).
# Filters by the protein-level alt amino acid (REVEL row's `aaalt`), so the
# caller does not need to reverse-translate or know the genomic ref nucleotide.
# When several rows match (rare — different alt nucleotides encoding the same
# amino acid), returns the median REVEL.
lookup_revel_by_aa <- function(revel_dt, chrom, codon_positions, alt_aa) {
  dt_chrom <- attr(revel_dt, "revel_chrom")
  if (!is.null(dt_chrom) && as.character(chrom) != dt_chrom) {
    stop(sprintf("revel_dt is for chromosome %s but lookup requested %s",
                 dt_chrom, chrom), call. = FALSE)
  }
  pos_int <- suppressWarnings(as.integer(codon_positions))
  pos_int <- pos_int[!is.na(pos_int)]
  if (length(pos_int) == 0L) return(NA_real_)
  hits <- revel_dt[grch38_pos %in% pos_int & aaalt == as.character(alt_aa)]
  if (nrow(hits) == 0L) return(NA_real_)
  as.numeric(stats::median(hits$REVEL, na.rm = TRUE))
}

# ---------------------------------------------------------------------------
# dbNSFP 4.9a (single-source, hg38)
#
# Single-file build at analyses/tmp/dbNSFPv4.9a_hg38_custombuild.bgz with a
# tabix index (.tbi). 458 columns including SIFT, Polyphen2 HDIV/HVAR, LRT,
# MutationTaster, FATHMM, PROVEAN, MetaSVM/LR/RNN, REVEL, AlphaMissense,
# CADD, DANN, GERP++ RS/NR, PhyloP 100way/470way/17way primate, PhastCons
# 100way/17way primate, SiPhy_29way pi.
#
# Strategy: tabix once per gene region (load_dbnsfp_for_region), hash by
# "<pos>_<ref>_<alt>" for O(1) per-variant lookup (lookup_dbnsfp_local).
# Returns a list shaped like the raw MyVariant `hit` so server.R's
# parse_dbnsfp_scores(hit) consumes it directly — no synthesis layer needed.
#
# Column indices below are pinned to dbNSFP 4.9a (Aug 2024 release). If the
# upstream schema shifts, recompute via:
#   zcat <chr>.gz | head -1 | tr '\t' '\n' | nl | grep -E '<col>'
# ---------------------------------------------------------------------------

DBNSFP_DEFAULT_PATH <- "analyses/tmp/dbNSFPv4.9a_hg38_custombuild.bgz"

load_dbnsfp_for_region <- function(chrom, start, end,
                                   dbnsfp_path = DBNSFP_DEFAULT_PATH) {
  if (!file.exists(dbnsfp_path)) return(NULL)
  if (!file.exists(paste0(dbnsfp_path, ".tbi"))) {
    stop(sprintf("tabix index missing: %s.tbi (build with: tabix -s1 -b2 -e2 %s)",
                 dbnsfp_path, dbnsfp_path), call. = FALSE)
  }
  region <- sprintf("%s:%d-%d", as.character(chrom),
                    as.integer(start), as.integer(end))
  cmd <- sprintf("tabix %s %s", shQuote(dbnsfp_path), shQuote(region))
  rows <- tryCatch(system(cmd, intern = TRUE, ignore.stderr = TRUE),
                   error = function(e) character(0))
  if (length(rows) == 0L) return(NULL)
  env   <- new.env(hash = TRUE, parent = emptyenv(), size = length(rows))
  env_p <- new.env(hash = TRUE, parent = emptyenv(), size = length(rows))
  for (line in rows) {
    p <- strsplit(line, "\t", fixed = TRUE)[[1]]
    if (length(p) < 4L) next
    # Genomic key: (pos, ref, alt) — for the rare case the harness has hgvs_g
    g_key <- paste0(p[2], "_", p[3], "_", p[4])
    assign(g_key, p, envir = env)
    # Protein key: (aapos, aaalt) — primary lookup path for hgvsp queries.
    # aapos (col 12) can be ";"-separated for multi-transcript variants;
    # index every value so isoform-specific HGVSp lookups also resolve.
    if (length(p) >= 12L) {
      aa_alt <- if (length(p) >= 6L) p[6] else "."
      aa_pos_raw <- p[12]
      if (!is.na(aa_alt) && nzchar(aa_alt) && aa_alt != "." &&
          !is.na(aa_pos_raw) && nzchar(aa_pos_raw) && aa_pos_raw != ".") {
        aa_pos_parts <- unique(strsplit(aa_pos_raw, ";", fixed = TRUE)[[1]])
        aa_pos_parts <- aa_pos_parts[aa_pos_parts != "" & aa_pos_parts != "."]
        for (ap in aa_pos_parts) {
          p_key <- paste0(ap, "_", aa_alt)
          # First-write-wins; multiple genomic encodings of the same
          # missense produce identical predictor scores in practice.
          if (!exists(p_key, envir = env_p, inherits = FALSE)) {
            assign(p_key, p, envir = env_p)
          }
        }
      }
    }
  }
  attr(env, "dbnsfp_chrom")   <- as.character(chrom)
  attr(env, "dbnsfp_region")  <- c(start = as.integer(start), end = as.integer(end))
  attr(env, "protein_env")    <- env_p
  env
}

# Private: convert a single dbNSFP row (character vector of 458 cols) into
# the parse_dbnsfp_scores-shaped list. Tolerates "." (dbNSFP missing marker)
# and ";"-separated transcript-level values (takes first non-empty entry —
# matches MyVariant.info's behaviour of surfacing the first transcript
# prediction).
.dbnsfp_row_to_hit <- function(p) {
  num_first <- function(i) {
    if (i > length(p)) return(NA_real_)
    v <- p[i]
    if (is.na(v) || v == "." || v == "") return(NA_real_)
    parts <- strsplit(v, ";", fixed = TRUE)[[1]]
    parts <- parts[parts != "" & parts != "."]
    if (length(parts) == 0L) return(NA_real_)
    suppressWarnings(as.numeric(parts[1]))
  }
  chr_first <- function(i) {
    if (i > length(p)) return(NA_character_)
    v <- p[i]
    if (is.na(v) || v == "." || v == "") return(NA_character_)
    parts <- strsplit(v, ";", fixed = TRUE)[[1]]
    parts <- parts[parts != "" & parts != "."]
    if (length(parts) == 0L) return(NA_character_)
    parts[1]
  }
  list(
    dbnsfp = list(
      sift            = list(score = num_first(38),  pred = chr_first(40)),
      polyphen2       = list(
        hdiv          = list(score = num_first(44),  pred = chr_first(46)),
        hvar          = list(score = num_first(47),  pred = chr_first(49))
      ),
      lrt             = list(score = num_first(50),  pred = chr_first(52)),
      mutationtaster  = list(score = num_first(54),  pred = chr_first(56)),
      fathmm          = list(score = num_first(62),  pred = chr_first(64)),
      provean         = list(score = num_first(65),  pred = chr_first(67)),
      metasvm         = list(score = num_first(70),  pred = chr_first(72)),
      metalr          = list(score = num_first(73),  pred = chr_first(75)),
      metarnn         = list(score = num_first(77),  pred = chr_first(79)),
      revel           = list(score = num_first(83)),
      alphamissense   = list(am_pathogenicity = num_first(138),
                             am_class         = chr_first(140)),
      cadd            = list(raw_rankscore = num_first(154),
                             phred         = num_first(155)),
      dann            = list(score = num_first(159)),
      `gerp++`        = list(rs = num_first(191), nr = num_first(190)),
      phylop          = list(
        `100way_vertebrate` = list(score = num_first(195)),
        `470way_mammalian`  = list(score = num_first(197)),
        `17way_primate`     = list(score = num_first(199))
      ),
      phastcons       = list(
        `100way_vertebrate` = list(score = num_first(201)),
        `17way_primate`     = list(score = num_first(205))
      ),
      siphy_29way     = list(pi = list(a = num_first(207)))
      # dbscsnv (splice ada/rf) intentionally absent — not in dbNSFP 4.9a;
      # safe_extract_num returns NA_real_ for the missing path.
    )
  )
}

lookup_dbnsfp_local <- function(env, chrom, pos, ref, alt) {
  if (is.null(env)) return(NULL)
  dt_chrom <- attr(env, "dbnsfp_chrom")
  if (!is.null(dt_chrom) && as.character(chrom) != dt_chrom) return(NULL)
  key <- paste0(as.integer(pos), "_", as.character(ref), "_", as.character(alt))
  if (!exists(key, envir = env, inherits = FALSE)) return(NULL)
  .dbnsfp_row_to_hit(get(key, envir = env, inherits = FALSE))
}

# Protein-key path — the harness's primary lookup since the universe TSV
# only carries hgvs_g for ~0.8% of variants. aa_pos comes from the HGVS-p
# (parsed by extract_pos), aa_alt is the one-letter alt code.
lookup_dbnsfp_by_aa <- function(env, chrom, aa_pos, aa_alt) {
  if (is.null(env)) return(NULL)
  dt_chrom <- attr(env, "dbnsfp_chrom")
  if (!is.null(dt_chrom) && as.character(chrom) != dt_chrom) return(NULL)
  env_p <- attr(env, "protein_env")
  if (is.null(env_p)) return(NULL)
  key <- paste0(as.integer(aa_pos), "_", as.character(aa_alt))
  if (!exists(key, envir = env_p, inherits = FALSE)) return(NULL)
  .dbnsfp_row_to_hit(get(key, envir = env_p, inherits = FALSE))
}
