#!/usr/bin/env Rscript
# 15p_build_pfam_pssm_parallel.R — parallel build of the genome-wide PM1 PSSM
# table. Same output as 15_build_pfam_pssm_table.R, but built to use many cores.
#
# The serial version was bottlenecked twice over: readLines(n=1) crawls a 11 GB
# gz one line at a time (R's slowest I/O), and every per-family build_pssm runs
# on a single core. Both are fixed here:
#
#   Stage 1 (stream)  pigz -dc | awk splits the alignment into one file per
#                     HUMAN-BEARING family (families with no _HUMAN row can never
#                     produce a human residue map, so they are dropped on the fly).
#                     awk in C does in seconds what per-line readLines did in
#                     minutes.
#   Stage 2 (build)   mclapply builds the Corcuff PSSM + human residue map for
#                     each family file across cores — families are independent.
#   Stage 3 (assemble) quantise, bind the maps, saveRDS. Serial, cheap.
#
# Usage:  Rscript analyses/15p_build_pfam_pssm_parallel.R [path/to/Pfam-A.full.gz] [ncores]

suppressMessages(library(parallel))
source("analyses/lib/pfam_pssm.R")

args    <- commandArgs(trailingOnly = TRUE)
PFAM_GZ <- if (length(args) >= 1) args[1] else "analyses/tmp/pfam/Pfam33.1-A.full.gz"
NCORES  <- if (length(args) >= 2) as.integer(args[2]) else 24L
FAMDIR  <- "analyses/tmp/pfam/fam"
OUT     <- "analyses/derived/pfam_pssm_human.rds"
Q_MIN <- -12; Q_MAX <- 6
# Preserve the column x 20 shape + AA column names: as.integer() alone would
# flatten the matrix to a vector, breaking the runtime M[column, mut] lookup.
quantise <- function(x) {
  q <- as.integer(pmax(-127, pmin(127,
         round((pmin(pmax(x, Q_MIN), Q_MAX) - Q_MIN) / (Q_MAX - Q_MIN) * 254 - 127))))
  dim(q) <- dim(x); dimnames(q) <- dimnames(x); q
}
# Cap sequences per family so a mega-family (PF00005 ABC transporter: 676k x
# 2369) cannot OOM a worker. All _HUMAN rows are always kept (they drive the
# residue map); the rest are subsampled. Column frequencies are stable well
# below this cap, so PM1 calls are unaffected.
MAX_SEQ <- 200000L

if (!file.exists(PFAM_GZ)) stop("Pfam alignment not found: ", PFAM_GZ)
dir.create("analyses/derived", showWarnings = FALSE, recursive = TRUE)
unlink(FAMDIR, recursive = TRUE); dir.create(FAMDIR, recursive = TRUE)

t0 <- Sys.time()

# ── Stage 1: split into per-family files, human-bearing only ───────────────
# Buffer each family's sequence lines; on "//" write them to fam/<AC>.sto but
# only if a _HUMAN row was seen. The AC is stripped to its family base (PF#####).
# gsub in awk keeps the version off the filename so downstream keys are clean.
message("[pssm|| ] stage 1: splitting ", PFAM_GZ, " by family (human-bearing only)")
awk_prog <- '
  /^#=GF AC/ { ac=$3; sub(/\\..*/,"",ac); n=0; hum=0; next }
  /^\\/\\// {
    if (hum && ac!="") { f=dir"/"ac".sto"; for(i=0;i<n;i++) print buf[i] > f; close(f); k++ }
    ac=""; next
  }
  /^#/  { next }
  /^$/  { next }
  { buf[n++]=$0; if ($0 ~ /_HUMAN\\//) hum=1 }
  END { print k > "/dev/stderr" }
'
# pigz decompresses with helper threads; awk streams the split in one pass.
cmd <- sprintf("pigz -dc %s | awk -v dir=%s %s",
               shQuote(PFAM_GZ), shQuote(FAMDIR), shQuote(awk_prog))
kept <- system(cmd, intern = FALSE, ignore.stderr = FALSE)
fam_files <- list.files(FAMDIR, pattern = "\\.sto$", full.names = TRUE)
message(sprintf("[pssm|| ] stage 1 done: %s human-bearing families, %.1f min",
                format(length(fam_files), big.mark = ","),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
if (!length(fam_files)) stop("no human-bearing family files written")

# ── Stage 2: build PSSM + human residue map per family, in parallel ────────
# Returns the quantised score matrix (column x 20) and the human residue->column
# map for one family. Mirrors flush_family() in the serial script exactly.
build_one <- function(path) {
  fam <- sub("\\.sto$", "", basename(path))
  buf <- readLines(path, warn = FALSE)
  if (length(buf) < 2L) return(NULL)
  # Subsample huge families before the O(nseq*ncol) matrix in build_pssm, always
  # retaining human rows so the residue map is complete.
  if (length(buf) > MAX_SEQ) {
    is_hum <- grepl("_HUMAN/", buf)
    others <- which(!is_hum)
    keep <- sort(c(which(is_hum), others[seq_len(MAX_SEQ - sum(is_hum))]))
    buf <- buf[keep]
  }
  parts <- strsplit(trimws(buf), "[[:space:]]+")
  ids   <- vapply(parts, `[`, character(1), 1L)
  seqs  <- vapply(parts, function(p) paste(p[-1], collapse = ""), character(1))
  aln   <- tapply(seqs, ids, paste, collapse = "")
  if (length(aln) < 2L) return(NULL)
  M <- tryCatch(build_pssm(aln), error = function(e) NULL)
  if (is.null(M)) return(NULL)

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
      column = nongap[seq_len(n)], family = fam, stringsAsFactors = FALSE)
  }
  if (!length(mp)) return(NULL)
  list(fam = fam, pssm = quantise(M), map = do.call(rbind, mp))
}

message(sprintf("[pssm|| ] stage 2: building %s families on %d cores",
                format(length(fam_files), big.mark = ","), NCORES))
# preschedule=FALSE so the few very large families are load-balanced across
# workers rather than piling onto one chunk.
res <- mclapply(fam_files, build_one, mc.cores = NCORES, mc.preschedule = FALSE)
res <- Filter(Negate(is.null), res)
message(sprintf("[pssm|| ] stage 2 done: %s families produced a PSSM, %.1f min",
                format(length(res), big.mark = ","),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

# ── Stage 3: assemble + save ───────────────────────────────────────────────
pssm_rows <- lapply(res, `[[`, "pssm")
names(pssm_rows) <- vapply(res, `[[`, character(1), "fam")
map <- do.call(rbind, lapply(res, `[[`, "map"))
out <- list(pssm = pssm_rows, map = map,
            quant = c(min = Q_MIN, max = Q_MAX),
            source = basename(PFAM_GZ),
            built  = format(Sys.time(), "%Y-%m-%d"),
            aa = AA20)
saveRDS(out, OUT, compress = "xz")

# ── Text companions (grep-able, no R needed) ───────────────────────────────
# Full human residue -> (family, column) lookup, and a per-family summary.
MAP_TSV <- sub("\\.rds$", "_map.tsv", OUT)
FAM_TSV <- sub("\\.rds$", "_families.tsv", OUT)
write.table(map[, c("id", "residue", "family", "column")], MAP_TSV,
            sep = "\t", row.names = FALSE, quote = FALSE)
fam_summary <- data.frame(
  family      = names(pssm_rows),
  columns     = vapply(pssm_rows, nrow, integer(1)),
  human_resid = as.integer(table(factor(map$family, levels = names(pssm_rows)))),
  stringsAsFactors = FALSE)
fam_summary <- fam_summary[order(-fam_summary$human_resid), ]
write.table(fam_summary, FAM_TSV, sep = "\t", row.names = FALSE, quote = FALSE)

message("\n[pssm|| ] summary")
message(sprintf("  human-bearing families : %s", format(length(fam_files), big.mark = ",")))
message(sprintf("  families with a PSSM   : %s", format(length(res), big.mark = ",")))
message(sprintf("  human residues mapped  : %s", format(nrow(map), big.mark = ",")))
message(sprintf("  in-RAM                 : %s", format(object.size(out), units = "MB")))
message(sprintf("  on disk (xz)           : %.1f MB", file.size(OUT) / 1e6))
message(sprintf("  wall time              : %.1f min",
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
message("  written to ", OUT)
message("            ", MAP_TSV)
message("            ", FAM_TSV)
