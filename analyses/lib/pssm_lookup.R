# pssm_lookup.R — runtime lookup against the offline Pfam PSSM table.
#
# Replaces the DOLPHIN API call in PM1 Path 4 with a local table read. The table
# (analyses/derived/pfam_pssm_human.rds, built by 15p) holds, per Pfam family, an
# int8-quantised score matrix [alignment column x 20 AA], plus a map from each
# human residue to its (family, column). The runtime cost is a data-frame filter
# and one matrix index — no network, no rate limit.
#
# Coordinate note: the map is keyed on the UniProt ENTRY NAME (e.g. SYUA_HUMAN),
# not the accession. VarViz holds the accession (P37840); bridge with the entry
# name UniProt returns as `uniProtkbId` on the same fetch it already does, or
# precompute an accession->entry map. Residue numbering is the UniProt canonical
# isoform — the same basis as VarViz's other UniProt-anchored layers.

# Load once and keep in the app (it is ~50 MB / 600 MB in RAM).
pssm_table_load <- function(path = "analyses/derived/pfam_pssm_human.rds") {
  stopifnot(file.exists(path))
  t <- readRDS(path)
  # Shrink the in-memory footprint (~640 MB -> ~280 MB) so it fits a 1 GB cloud
  # instance without OOM. The rds stores each int8 score matrix as R `integer`
  # (4 bytes/cell); re-pack as `raw` (1 byte). Values are in [-127, 127], so
  # store q+127 in [0, 254]; .pssm_dequant reads raw as already-offset. The map's
  # repeated character keys (id ~52k entries, family ~6.5k) become factors.
  if (length(t$pssm) && is.integer(t$pssm[[1L]])) {
    t$pssm <- lapply(t$pssm, function(m) {
      r <- as.raw(as.integer(m) + 127L)
      dim(r) <- dim(m); dimnames(r) <- dimnames(m); r
    })
  }
  if (is.character(t$map$id))     t$map$id     <- as.factor(t$map$id)
  if (is.character(t$map$family)) t$map$family <- as.factor(t$map$family)
  t
}

# int8 -> log-odds. The quant bounds live in the table so this stays in sync.
# Accepts either the raw-packed value (already q+127, from pssm_table_load) or a
# plain integer/numeric q (needs the +127 offset), so tests can call it directly.
.pssm_dequant <- function(q, tbl) {
  base <- if (is.raw(q)) as.numeric(q) else as.numeric(q) + 127
  base / 254 * (tbl$quant[["max"]] - tbl$quant[["min"]]) + tbl$quant[["min"]]
}

# All (family, column) a residue maps to. A residue can sit in more than one
# Pfam domain, so this may return several rows.
pssm_sites <- function(tbl, entry_name, residue) {
  m <- tbl$map
  m[m$id == entry_name & m$residue == residue, c("family", "column"), drop = FALSE]
}

# delta = M(column, mut) - M(column, wt), most-disfavoured across the residue's
# families (DOLPHIN reports negative = substitution disfavoured at that column).
# Returns list(delta, family, column) or delta = NA if the residue is in no
# family or the amino acids are non-standard.
pssm_delta <- function(tbl, entry_name, residue, wt, mut) {
  wt <- toupper(wt); mut <- toupper(mut)
  aa <- tbl$aa
  if (!(wt %in% aa) || !(mut %in% aa)) return(list(delta = NA_real_, family = NA, column = NA))
  sites <- pssm_sites(tbl, entry_name, residue)
  if (!nrow(sites)) return(list(delta = NA_real_, family = NA, column = NA))
  best <- list(delta = Inf, family = NA, column = NA)
  for (i in seq_len(nrow(sites))) {
    fam <- as.character(sites$family[i]); col <- sites$column[i]
    M <- tbl$pssm[[fam]]
    if (is.null(M) || col > nrow(M)) next
    d <- .pssm_dequant(M[col, mut], tbl) - .pssm_dequant(M[col, wt], tbl)
    if (d < best$delta) best <- list(delta = d, family = fam, column = col)
  }
  if (is.infinite(best$delta)) best$delta <- NA_real_
  best
}

# PM1 decision. DOLPHIN's own cutoff is server-side and undocumented, so we set
# our own: fire when the substitution is clearly disfavoured at a conserved
# column. -4 cleanly separates the DOLPHIN-pathogenic set (SNCA G14R/E46K/G51D,
# CASR G143E: -5 to -9) from the tolerated A53T (+1.2); calibrate against the
# dual-pass benchmark before production.
pssm_fires_pm1 <- function(tbl, entry_name, residue, wt, mut, threshold = -4) {
  d <- pssm_delta(tbl, entry_name, residue, wt, mut)$delta
  isTRUE(!is.na(d) && d <= threshold)
}

# ── demo / self-check ──────────────────────────────────────────────────────
if (sys.nframe() == 0L) {
  tbl <- pssm_table_load()
  cat(sprintf("table: %s families, %s residues\n\n",
              format(length(tbl$pssm), big.mark=","), format(nrow(tbl$map), big.mark=",")))
  # entry_name, residue, wt, mut, DOLPHIN reference delta
  cases <- list(
    c("SYUA_HUMAN", 14, "G", "R", -9.161),
    c("SYUA_HUMAN", 46, "E", "K", -8.873),
    c("SYUA_HUMAN", 51, "G", "D", -5.415),
    c("SYUA_HUMAN", 53, "A", "T",  1.165),
    c("CASR_HUMAN",143, "G", "E", -5.959))
  cat(sprintf("%-14s %-6s %-8s %-9s %-6s\n", "protein", "var", "ours", "DOLPHIN", "PM1?"))
  for (c0 in cases) {
    en <- c0[1]; pos <- as.integer(c0[2]); wt <- c0[3]; mut <- c0[4]; ref <- as.numeric(c0[5])
    r <- pssm_delta(tbl, en, pos, wt, mut)
    fire <- pssm_fires_pm1(tbl, en, pos, wt, mut)
    cat(sprintf("%-14s %s%d%s%-3s %-6.3f %-9.3f %s\n", en, wt, pos, mut, "",
                r$delta, ref, if (fire) "PM1" else "-"))
    stopifnot(!is.na(r$delta), abs(r$delta - ref) < 0.3)   # within int8 + drift
  }
  cat("\npssm_lookup demo: all deltas within tolerance of DOLPHIN\n")
}
