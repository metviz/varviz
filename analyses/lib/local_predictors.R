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
    # Protein key: (aapos, aaref, aaalt) — primary lookup path for hgvsp
    # queries. aapos (col 12) can be ";"-separated for multi-transcript
    # variants, so ONE row is indexed at several positions. Keying on position
    # and alt alone therefore collides across transcripts: SNCA's canonical
    # G51D and an alternate-transcript N51D both land on "51_D", and whichever
    # was read first wins. That returned N51D's scores (MetaSVM tolerated,
    # REVEL absent) for a G51D query and cost the variant a PP3 tier. The
    # reference residue disambiguates them. The legacy (aapos, aaalt) key is
    # still written, first-write-wins, so a caller that cannot supply aaref
    # keeps the old behaviour.
    if (length(p) >= 12L) {
      aa_ref <- if (length(p) >= 5L) p[5] else "."
      aa_alt <- if (length(p) >= 6L) p[6] else "."
      aa_pos_raw <- p[12]
      if (!is.na(aa_alt) && nzchar(aa_alt) && aa_alt != "." &&
          !is.na(aa_pos_raw) && nzchar(aa_pos_raw) && aa_pos_raw != ".") {
        aa_pos_parts <- unique(strsplit(aa_pos_raw, ";", fixed = TRUE)[[1]])
        aa_pos_parts <- aa_pos_parts[aa_pos_parts != "" & aa_pos_parts != "."]
        for (ap in aa_pos_parts) {
          if (!is.na(aa_ref) && nzchar(aa_ref) && aa_ref != ".") {
            ref_key <- paste0(ap, "_", aa_ref, "_", aa_alt)
            if (!exists(ref_key, envir = env_p, inherits = FALSE)) {
              assign(ref_key, p, envir = env_p)
            }
          }
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

# ---------------------------------------------------------------------------
# MANE transcript resolution
# ---------------------------------------------------------------------------
# dbNSFP's per-transcript columns are only usable once you know which element
# belongs to the MANE Select transcript, and MANE is the first element for only
# half the benchmark genes (see analyses/reports/dbnsfp_mane_index_survey.tsv).
# The map is read once from the same gene_mane.tsv the app's MANE join uses, then
# memoised, because the row mapper runs once per variant.
.MANE_TX_ENV <- new.env(parent = emptyenv())

.dbnsfp_mane_tx <- function(gene, mane_path = "analyses/derived/gene_mane.tsv") {
  if (length(gene) == 0L) return(NA_character_)
  gene <- as.character(gene)[1]
  if (is.na(gene) || !nzchar(gene)) return(NA_character_)
  if (is.null(.MANE_TX_ENV$map)) {
    # Dual-cwd, matching the tests: run from the project root or from
    # analyses/tests. A missing file is not fatal -- it degrades to the
    # first-element behaviour rather than dropping every predictor.
    path <- if (file.exists(mane_path)) mane_path else file.path("..", mane_path)
    if (!file.exists(path)) {
      .MANE_TX_ENV$map <- character(0)
      return(NA_character_)
    }
    tb <- data.table::fread(path, select = c("gene_name", "mane_ensembl"),
                            showProgress = FALSE, na.strings = c("", "NA"))
    tx <- sub("\\..*$", "", as.character(tb$mane_ensembl))   # drop the version suffix
    nm <- as.character(tb$gene_name)
    ok <- !is.na(nm) & nzchar(nm) & !is.na(tx) & nzchar(tx)
    m  <- tx[ok]
    names(m) <- nm[ok]
    # First row wins: gene_mane.tsv carries one MANE Select per symbol, so a
    # repeated symbol is an alias row, not a second transcript.
    .MANE_TX_ENV$map <- m[!duplicated(names(m))]
  }
  m <- .MANE_TX_ENV$map
  if (!length(m)) return(NA_character_)
  i <- match(gene, names(m))
  if (is.na(i)) NA_character_ else unname(m[[i]])
}

# NA asks the row mapper to resolve MANE from the row; NULL is the pre-change
# first-element read. The option exists so a historical run can be reproduced
# byte for byte, which is otherwise impossible once selection changes.
.dbnsfp_mane_mode <- function() {
  if (isTRUE(getOption("varviz.dbnsfp_mane", TRUE))) NA_character_ else NULL
}

# Private: convert a single dbNSFP row (character vector of 458 cols) into
# the parse_dbnsfp_scores-shaped list. Tolerates "." (dbNSFP missing marker)
# and ";"-separated transcript-level values.
#
# mane_tx controls which element of a ";"-separated column is read:
#   NULL  take the first non-"." element. The pre-change behaviour, kept so a
#         caller that has no gene context, and every existing test, is unchanged.
#   NA    resolve the MANE transcript from this row's own genename column.
#   <id>  use this Ensembl transcript id.
#
# Only columns verified to be aligned with Ensembl_transcriptid (col 15) are
# indexed. Measured over 5,256 multi-transcript rows across eight genes, list
# lengths matched the transcript count in every row for aapos, Ensembl_proteinid,
# HGVSp_VEP, SIFT, SIFT4G, Polyphen2 HDIV/HVAR, FATHMM, PROVEAN, MetaRNN, REVEL
# and AlphaMissense, and never matched it for LRT, MetaSVM, MetaLR, CADD, DANN,
# GERP, PhyloP, PhastCons, SiPhy, 1000Gp3_AF and ExAC_AF, which carry one value.
# MutationTaster is the trap: 357 rows aligned, 331 single and 4,568 neither,
# because its list is per prediction model. Indexing it by transcript position
# would pair a transcript with another axis's value, which is the same class of
# error this change removes, so it keeps the first element.
.dbnsfp_row_to_hit <- function(p, mane_tx = NULL) {
  mane_idx <- NA_integer_
  if (!is.null(mane_tx)) {
    tx_raw <- if (length(p) >= 15L) p[15] else "."
    if (!is.na(tx_raw) && tx_raw != "." && tx_raw != "") {
      tx   <- strsplit(tx_raw, ";", fixed = TRUE)[[1]]
      want <- as.character(mane_tx)[1]
      if (is.na(want)) {
        # A dbNSFP row can span overlapping genes, so try each symbol and keep
        # the first whose MANE transcript this row actually lists.
        gn_raw <- if (length(p) >= 13L) p[13] else "."
        if (!is.na(gn_raw) && gn_raw != "." && gn_raw != "") {
          for (g in unique(strsplit(gn_raw, ";", fixed = TRUE)[[1]])) {
            cand <- .dbnsfp_mane_tx(g)
            if (!is.na(cand) && cand %in% tx) { want <- cand; break }
          }
        }
      }
      if (!is.na(want)) {
        i <- which(tx == want)
        if (length(i)) mane_idx <- i[1]
      }
    }
  }
  # use_mane = TRUE reads the MANE element and falls back to the first scored one,
  # because dbNSFP writes "." for a transcript a tool could not score and losing
  # that row's evidence entirely would be worse than reporting another isoform's.
  pick <- function(i, use_mane) {
    if (i > length(p)) return(NA_character_)
    v <- p[i]
    if (is.na(v) || v == "." || v == "") return(NA_character_)
    parts <- strsplit(v, ";", fixed = TRUE)[[1]]
    if (use_mane && !is.na(mane_idx) && mane_idx <= length(parts)) {
      cand <- parts[mane_idx]
      if (!is.na(cand) && cand != "." && cand != "") return(cand)
    }
    keep <- parts[parts != "" & parts != "."]
    if (length(keep) == 0L) return(NA_character_)
    keep[1]
  }
  num_tx    <- function(i) suppressWarnings(as.numeric(pick(i, TRUE)))
  chr_tx    <- function(i) pick(i, TRUE)
  num_first <- function(i) suppressWarnings(as.numeric(pick(i, FALSE)))
  chr_first <- function(i) pick(i, FALSE)
  list(
    dbnsfp = list(
      sift            = list(score = num_tx(38),    pred = chr_tx(40)),
      polyphen2       = list(
        hdiv          = list(score = num_tx(44),    pred = chr_tx(46)),
        hvar          = list(score = num_tx(47),    pred = chr_tx(49))
      ),
      lrt             = list(score = num_first(50), pred = chr_first(52)),
      mutationtaster  = list(score = num_first(54), pred = chr_first(56)),
      fathmm          = list(score = num_tx(62),    pred = chr_tx(64)),
      provean         = list(score = num_tx(65),    pred = chr_tx(67)),
      metasvm         = list(score = num_first(70), pred = chr_first(72)),
      metalr          = list(score = num_first(73), pred = chr_first(75)),
      metarnn         = list(score = num_tx(77),    pred = chr_tx(79)),
      revel           = list(score = num_tx(83)),
      alphamissense   = list(am_pathogenicity = num_tx(138),
                             am_class         = chr_tx(140)),
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
      siphy_29way     = list(pi = list(a = num_first(207))),
      # Population frequencies. Verified against the dbNSFP 4.9a chr1 header:
      # column 213 is 1000Gp3_AF and column 235 is ExAC_AF. These were previously
      # left NA and a per-variant rsID fetch was attempted instead, through a
      # function that does not exist; the values were already in the row.
      onekg           = list(af = num_first(213)),
      exac            = list(af = num_first(235))
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
  .dbnsfp_row_to_hit(get(key, envir = env, inherits = FALSE), mane_tx = .dbnsfp_mane_mode())
}

# Protein-key path — the harness's primary lookup since the universe TSV
# only carries hgvs_g for ~0.8% of variants. aa_pos comes from the HGVS-p
# (parsed by extract_pos), aa_alt is the one-letter alt code.
lookup_dbnsfp_by_aa <- function(env, chrom, aa_pos, aa_alt, aa_ref = NA_character_) {
  if (is.null(env)) return(NULL)
  dt_chrom <- attr(env, "dbnsfp_chrom")
  if (!is.null(dt_chrom) && as.character(chrom) != dt_chrom) return(NULL)
  env_p <- attr(env, "protein_env")
  if (is.null(env_p)) return(NULL)
  # With a reference residue, demand an exact (pos, ref, alt) match: a row
  # indexed at this position for a DIFFERENT reference residue is a different
  # substitution, and its scores must never stand in for this one. Returning
  # NULL sends the caller to its documented fallback, which is recoverable;
  # returning the wrong variant's predictors is not.
  if (!is.na(aa_ref) && nzchar(as.character(aa_ref))) {
    key <- paste0(as.integer(aa_pos), "_", as.character(aa_ref), "_", as.character(aa_alt))
    if (!exists(key, envir = env_p, inherits = FALSE)) return(NULL)
    return(.dbnsfp_row_to_hit(get(key, envir = env_p, inherits = FALSE),
                              mane_tx = .dbnsfp_mane_mode()))
  }
  key <- paste0(as.integer(aa_pos), "_", as.character(aa_alt))
  if (!exists(key, envir = env_p, inherits = FALSE)) return(NULL)
  .dbnsfp_row_to_hit(get(key, envir = env_p, inherits = FALSE), mane_tx = .dbnsfp_mane_mode())
}
