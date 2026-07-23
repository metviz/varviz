#!/usr/bin/env Rscript
# 15_build_pfam_pssm_table.R — genome-wide offline PM1 domain evidence.
#
# Streams Pfam-A.full.gz one family at a time, builds the Corcuff et al. 2023
# score matrix for each, and keeps only what human variant interpretation
# needs. Output replaces the DOLPHIN API call in PM1 Path 4.
#
# Storage design. The naive layouts do not fit:
#   per human residue x 19 alt AA   ~76M rows   too large
#   per residue, 20 doubles         ~640 MB     too large to bundle
# Instead we mirror DOLPHIN's own data model and store the matrix ONCE per
# (family, column) — deduplicating across every protein carrying that domain —
# plus a small map from human residue to (family, column):
#   pssm : family x column x 20, int8-quantised      ~25 MB
#   map  : uniprot, residue -> family, column         ~16 MB
# Runtime lookup is map -> pssm -> delta = M[col, mut] - M[col, wt].
#
# Usage:
#   Rscript analyses/15_build_pfam_pssm_table.R [path/to/Pfam-A.full.gz]
#
# Pfam 33.1 is the release DOLPHIN used; 38.2 is current. The table records
# which was used so the provenance is reproducible from the artifact itself.

source("analyses/lib/pfam_pssm.R")

args    <- commandArgs(trailingOnly = TRUE)
PFAM_GZ <- if (length(args)) args[1] else "analyses/tmp/pfam/Pfam33.1-A.full.gz"
OUT     <- "analyses/derived/pfam_pssm_human.rds"
# Quantisation: log-odds live roughly in [-12, +6]. int8 over that range costs
# ~0.07 per step, far below the ~1.0 release-to-release drift we measured
# against DOLPHIN, so it changes no PM1 call while cutting storage 8x.
Q_MIN <- -12; Q_MAX <- 6
# Preserve the column x 20 shape + AA names: as.integer() alone would flatten
# the matrix and break the runtime M[column, mut] lookup. See 15p for the
# parallel build that supersedes this script.
quantise   <- function(x) {
  q <- as.integer(pmax(-127, pmin(127,
         round((pmin(pmax(x, Q_MIN), Q_MAX) - Q_MIN) / (Q_MAX - Q_MIN) * 254 - 127))))
  dim(q) <- dim(x); dimnames(q) <- dimnames(x); q
}
dequantise <- function(q) (as.numeric(q) + 127) / 254 * (Q_MAX - Q_MIN) + Q_MIN

if (!file.exists(PFAM_GZ)) stop("Pfam alignment not found: ", PFAM_GZ)
dir.create("analyses/derived", showWarnings = FALSE, recursive = TRUE)

message("[pfam-pssm] streaming ", PFAM_GZ)
con <- gzfile(PFAM_GZ, "rt")
on.exit(close(con), add = TRUE)

pssm_rows <- list(); map_rows <- list()
acc <- NULL; buf <- character(0); nfam <- 0L; nkept <- 0L
t0 <- Sys.time()

flush_family <- function(acc, buf) {
  # Only families containing at least one human (UniProt _HUMAN) row matter.
  hum <- grep("_HUMAN/", buf, value = TRUE)
  if (!length(hum)) return(NULL)
  parts <- strsplit(trimws(buf), "[[:space:]]+")
  ids  <- vapply(parts, `[`, character(1), 1L)
  seqs <- vapply(parts, function(p) paste(p[-1], collapse = ""), character(1))
  aln  <- tapply(seqs, ids, paste, collapse = "")
  if (length(aln) < 2L) return(NULL)
  M <- build_pssm(aln)

  hids <- grep("_HUMAN/", names(aln), value = TRUE)
  mp <- list()
  for (hid in hids) {
    rng   <- sub("^.*/", "", hid)
    start <- suppressWarnings(as.integer(sub("-.*$", "", rng)))
    end   <- suppressWarnings(as.integer(sub("^.*-", "", rng)))
    if (is.na(start) || is.na(end)) next
    chars  <- strsplit(aln[[hid]], "", fixed = TRUE)[[1]]
    nongap <- which(chars %in% c(AA20, tolower(AA20)))
    n <- min(length(nongap), end - start + 1L)
    if (n < 1L) next
    mp[[length(mp) + 1L]] <- data.frame(
      id = sub("/.*$", "", hid), residue = start + seq_len(n) - 1L,
      column = nongap[seq_len(n)], stringsAsFactors = FALSE)
  }
  if (!length(mp)) return(NULL)
  list(M = M, map = do.call(rbind, mp))
}

repeat {
  line <- readLines(con, n = 1L, warn = FALSE)
  if (!length(line)) break
  if (startsWith(line, "#=GF AC")) {
    acc <- sub("^#=GF AC\\s+", "", trimws(line)); next
  }
  if (startsWith(line, "//")) {
    nfam <- nfam + 1L
    res <- if (!is.null(acc) && length(buf)) flush_family(acc, buf) else NULL
    if (!is.null(res)) {
      nkept <- nkept + 1L
      fam <- sub("\\..*$", "", acc)
      pssm_rows[[fam]] <- quantise(res$M)          # column x 20, int8
      res$map$family <- fam
      map_rows[[length(map_rows) + 1L]] <- res$map
    }
    if (nfam %% 1000L == 0L)
      message(sprintf("  %6d families scanned, %5d human-bearing, %.1f min",
                      nfam, nkept, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    acc <- NULL; buf <- character(0); next
  }
  if (startsWith(line, "#") || !nzchar(trimws(line))) next
  buf <- c(buf, line)
}

map <- do.call(rbind, map_rows)
out <- list(
  pssm = pssm_rows, map = map,
  quant = c(min = Q_MIN, max = Q_MAX),
  source = basename(PFAM_GZ),
  built = format(Sys.time(), "%Y-%m-%d"),
  aa = AA20)
saveRDS(out, OUT, compress = "xz")

message("\n[pfam-pssm] summary")
message(sprintf("  families scanned      : %s", format(nfam, big.mark = ",")))
message(sprintf("  human-bearing kept    : %s", format(nkept, big.mark = ",")))
message(sprintf("  human residues mapped : %s", format(nrow(map), big.mark = ",")))
message(sprintf("  in-RAM                : %s", format(object.size(out), units = "MB")))
message(sprintf("  on disk (xz)          : %.1f MB", file.size(OUT) / 1e6))
message("  written to ", OUT)
