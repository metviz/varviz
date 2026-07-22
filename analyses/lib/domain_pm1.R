# domain_pm1.R — conservation-based domain constraint for a PM1 pathway.
#
# Motivation: the existing four PM1 pathways each have a dependency we would
# rather not carry. CCRS is derived from gnomAD depletion, so it is not
# independent of the PM2/BA1/BS1 frequency gate. The ClinVar 15-residue
# hotspot is ClinVar-circular by construction. DOLPHIN is an external REST
# service (and went dark for us when its TLS certificate lapsed on
# 2026-07-13). This module derives the same kind of evidence offline from
# cross-species conservation within InterPro-annotated domains, which is
# independent of BOTH ClinVar and population frequency — the property the
# Pass-Full / Pass-Blind design needs and that no existing domain method has
# (MetaDome aggregates gnomAD + ClinVar; see docs/CITATIONS.md).
#
# The score answers "is this residue unusually conserved FOR ITS OWN DOMAIN",
# which is the PM1 question ("critical / well-established functional domain"),
# rather than "is this residue conserved in absolute terms" — a whole-domain
# conservation floor would fire PM1 across entire domains and be useless for
# discrimination.
#
# Pure functions only; the extraction driver is analyses/13_build_domain_pm1.R.

# dbNSFP packs per-transcript values into ";"-separated fields, and uses "."
# for missing. A field can be ".;.;.;." — which is NOT missing-as-a-string and
# is why a naive x != "." test overcounts annotated rows roughly two-fold.
dbnsfp_first_present <- function(x) {
  if (length(x) != 1L || is.na(x)) return(NA_character_)
  parts <- strsplit(x, ";", fixed = TRUE)[[1]]
  parts <- trimws(parts)
  parts <- parts[nzchar(parts) & parts != "."]
  if (length(parts) == 0L) return(NA_character_)
  parts[1]
}

dbnsfp_numeric <- function(x) {
  v <- dbnsfp_first_present(x)
  if (is.na(v)) return(NA_real_)
  suppressWarnings(as.numeric(v))
}

# Within-domain percentile rank of conservation, in [0, 100].
#
# min_positions guards the small-domain case: a 4-residue domain would hand
# its top residue the 100th percentile on no evidence at all. Below the gate
# we return NA rather than a number, so the caller can distinguish "not
# constrained" from "not enough data to say" — the same distinction the
# DOLPHIN outage taught us to keep explicit.
domain_position_percentile <- function(conservation, min_positions = 20L,
                                       threshold = 90, max_fire_frac = 0.25) {
  n <- sum(!is.na(conservation))
  if (n < min_positions) return(rep(NA_real_, length(conservation)))
  # ties.method="max", NOT "average". GERP++ saturates: on the CASR
  # ligand-binding domain 23% of residues share its ceiling value of 6.17, so
  # under "average" the entire top plateau lands at the 88th percentile and a
  # 90th-percentile threshold becomes unreachable — the first build of this
  # scorer fired PM1 exactly zero times for that reason. With "max", every
  # member of a tied top group is treated as top-ranked, which is also the
  # correct semantics: conservation genuinely cannot discriminate within a
  # plateau, so it should not pretend to.
  r <- rank(conservation, na.last = "keep", ties.method = "max")
  pct <- 100 * (r - 1) / (n - 1)

  # Discrimination guard. "max" ties fix the saturation problem but introduce
  # the mirror one: in a uniformly conserved domain every residue ties at the
  # top and would fire PM1 everywhere. If more than max_fire_frac of the domain
  # clears the threshold, conservation is not discriminating within this domain
  # at all, so we return NA — "cannot say" — rather than blanketing it with
  # evidence. This is the same not-scorable-vs-negative distinction the small
  # domain gate makes, and that the DOLPHIN outage showed matters.
  if (mean(pct >= threshold, na.rm = TRUE) > max_fire_frac)
    return(rep(NA_real_, length(conservation)))
  pct
}

# PM1 decision. Supporting-only by default: this is positional evidence about
# a residue's context, not about the specific substitution, and ACMG PM1 is
# Moderate at most. Returning NA for an unscorable position keeps "unknown"
# separable from "no".
domain_pm1_call <- function(percentile, threshold = 90) {
  if (length(percentile) == 0L) return(logical(0))
  ifelse(is.na(percentile), NA, percentile >= threshold)
}

# Collapse per-variant rows to one row per (gene, domain, aapos). Conservation
# is a property of the genomic position, so every variant at a residue carries
# the same value; taking max() rather than mean() avoids letting the number of
# alternate alleles at a codon influence the score.
collapse_domain_positions <- function(df) {
  stopifnot(all(c("gene", "domain", "aapos", "cons") %in% names(df)))
  agg <- stats::aggregate(cons ~ gene + domain + aapos, data = df,
                          FUN = function(z) max(z, na.rm = TRUE))
  agg$cons[!is.finite(agg$cons)] <- NA_real_
  agg[order(agg$gene, agg$domain, agg$aapos), ]
}

# Domain/family regions for one UniProt accession, straight from InterPro.
#
# We do NOT use dbNSFP's Interpro_domain column for boundaries. That column
# carries only entries typed "domain", so it silently omits family-typed
# entries — SNCA has 1,126 rows in dbNSFP and zero domain annotations, because
# its InterPro entry IPR001058 (Synuclein, 1-130) is a *family*. Since SNCA is
# the intrinsically disordered case the application is meant to handle, a
# source that cannot see it is not usable. Querying InterPro directly also
# gives exact residue ranges, which the dbNSFP column never provided.
#
# This runs at BUILD time only; the resulting table ships offline, so the
# runtime path keeps no network dependency (the whole point after DOLPHIN).
fetch_interpro_regions <- function(uniprot_acc,
                                   types = c("domain", "family",
                                             "homologous_superfamily"),
                                   timeout_s = 30) {
  url <- sprintf(
    "https://www.ebi.ac.uk/interpro/api/entry/InterPro/protein/UniProt/%s/?page_size=100",
    uniprot_acc)
  txt <- tryCatch(
    paste(readLines(url(url), warn = FALSE), collapse = ""),
    error = function(e) NULL)
  if (is.null(txt) || !nzchar(txt)) return(NULL)
  js <- tryCatch(jsonlite::fromJSON(txt, simplifyVector = FALSE),
                 error = function(e) NULL)
  if (is.null(js) || is.null(js$results)) return(NULL)

  out <- list()
  for (r in js$results) {
    md <- r$metadata
    if (!is.null(types) && !(md$type %in% types)) next
    for (p in r$proteins) {
      for (loc in p$entry_protein_locations) {
        for (fr in loc$fragments) {
          out[[length(out) + 1L]] <- data.frame(
            acc = md$accession, type = md$type, name = md$name,
            start = as.integer(fr$start), end = as.integer(fr$end),
            stringsAsFactors = FALSE)
        }
      }
    }
  }
  if (!length(out)) return(NULL)
  do.call(rbind, out)
}

# Map residue positions onto InterPro regions. A residue inside several nested
# entries (common: a family spanning a domain spanning a conserved site) is
# assigned to the SMALLEST enclosing region, because the tightest annotation is
# the most specific statement about that residue's role — and scoring against a
# 900-residue family would wash out exactly the local signal PM1 wants.
assign_interpro_region <- function(aapos, regions) {
  if (is.null(regions) || !nrow(regions)) return(rep(NA_character_, length(aapos)))
  width <- regions$end - regions$start
  ord <- order(width)
  regions <- regions[ord, , drop = FALSE]
  vapply(aapos, function(p) {
    hit <- which(regions$start <= p & regions$end >= p)
    if (!length(hit)) return(NA_character_)
    regions$acc[hit[1]]
  }, character(1))
}

# Full pipeline for one collapsed table: percentile within each (gene, domain)
# instance, then the PM1 call.
score_domain_table <- function(pos_df, min_positions = 20L, threshold = 90) {
  stopifnot(all(c("gene", "domain", "aapos", "cons") %in% names(pos_df)))
  key <- paste(pos_df$gene, pos_df$domain, sep = "")
  pos_df$pctile <- unsplit(
    lapply(split(pos_df$cons, key), domain_position_percentile,
           min_positions = min_positions),
    key)
  pos_df$pm1 <- domain_pm1_call(pos_df$pctile, threshold = threshold)
  pos_df
}
