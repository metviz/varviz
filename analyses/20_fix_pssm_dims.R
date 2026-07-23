#!/usr/bin/env Rscript
# 20_fix_pssm_dims.R — repair the PSSM rds from the first 15p run and add the
# one family that OOM'd, without a full 21-minute rebuild.
#
# Two defects in that build:
#   1. quantise() used as.integer(M), which dropped the column x 20 matrix shape,
#      leaving each PSSM a flat vector — the runtime M[column, mut] lookup and the
#      families TSV both need the matrix back. as.integer flattens column-major,
#      so matrix(vec, ncol = 20) reconstructs it exactly; AA column names are
#      restored from AA20.
#   2. PF00005 (ABC transporter, 676k x 2369) exceeded worker memory and produced
#      no result. Rebuild it alone here, subsampling non-human rows to MAX_SEQ
#      (human rows always kept) so it fits, then splice it in.
#
# Usage:  Rscript analyses/20_fix_pssm_dims.R
source("analyses/lib/pfam_pssm.R")

RDS     <- "analyses/derived/pfam_pssm_human.rds"
FAMDIR  <- "analyses/tmp/pfam/fam"
Q_MIN <- -12; Q_MAX <- 6; MAX_SEQ <- 200000L
quantise <- function(x) {
  q <- as.integer(pmax(-127, pmin(127,
         round((pmin(pmax(x, Q_MIN), Q_MAX) - Q_MIN) / (Q_MAX - Q_MIN) * 254 - 127))))
  dim(q) <- dim(x); dimnames(q) <- dimnames(x); q
}

o <- readRDS(RDS)

# ── 1. Reshape every flat PSSM back to column x 20 ─────────────────────────
reshaped <- 0L
o$pssm <- lapply(o$pssm, function(v) {
  if (!is.null(dim(v))) return(v)                 # already a matrix
  reshaped <<- reshaped + 1L
  m <- matrix(v, ncol = length(AA20))             # column-major -> exact recon
  dimnames(m) <- list(NULL, AA20); m
})
message(sprintf("[fix] reshaped %d flat PSSMs to column x 20", reshaped))

# ── 2. Rebuild the dropped family (subsampled) and splice it in ────────────
build_one <- function(path) {
  fam <- sub("\\.sto$", "", basename(path))
  buf <- readLines(path, warn = FALSE)
  if (length(buf) > MAX_SEQ) {
    is_hum <- grepl("_HUMAN/", buf); others <- which(!is_hum)
    buf <- buf[sort(c(which(is_hum), others[seq_len(MAX_SEQ - sum(is_hum))]))]
    message(sprintf("[fix] %s subsampled to %s sequences", fam, format(length(buf), big.mark=",")))
  }
  parts <- strsplit(trimws(buf), "[[:space:]]+")
  ids   <- vapply(parts, `[`, character(1), 1L)
  seqs  <- vapply(parts, function(p) paste(p[-1], collapse = ""), character(1))
  aln   <- tapply(seqs, ids, paste, collapse = "")
  M <- build_pssm(aln)
  hids <- grep("_HUMAN/", names(aln), value = TRUE); mp <- list()
  for (hid in hids) {
    rng <- sub("^.*/", "", hid)
    start <- suppressWarnings(as.integer(sub("-.*$", "", rng)))
    end   <- suppressWarnings(as.integer(sub("^.*-", "", rng)))
    if (is.na(start) || is.na(end)) next
    chars  <- strsplit(aln[[hid]], "", fixed = TRUE)[[1]]
    nongap <- which(chars %in% c(AA20, tolower(AA20)))
    n <- min(length(nongap), end - start + 1L); if (n < 1L) next
    mp[[length(mp)+1L]] <- data.frame(id = sub("/.*$","",hid),
      residue = start + seq_len(n) - 1L, column = nongap[seq_len(n)],
      family = fam, stringsAsFactors = FALSE)
  }
  list(fam = fam, pssm = quantise(M), map = if (length(mp)) do.call(rbind, mp) else NULL)
}

allf <- sub("\\.sto$", "", list.files(FAMDIR, pattern = "\\.sto$"))
miss <- setdiff(allf, names(o$pssm))
message("[fix] missing families: ", if (length(miss)) paste(miss, collapse=", ") else "none")
for (m in miss) {
  r <- tryCatch(build_one(file.path(FAMDIR, paste0(m, ".sto"))),
                error = function(e) { message("[fix] ", m, " failed: ", e$message); NULL })
  if (is.null(r)) next
  o$pssm[[r$fam]] <- r$pssm
  if (!is.null(r$map)) o$map <- rbind(o$map, r$map[, colnames(o$map)])
  message(sprintf("[fix] added %s: %d columns, %d human residues",
                  r$fam, nrow(r$pssm), if (is.null(r$map)) 0L else nrow(r$map)))
}

# ── Re-save + text companions ──────────────────────────────────────────────
o$built <- format(Sys.Date(), "%Y-%m-%d")
saveRDS(o, RDS, compress = "xz")
MAP_TSV <- sub("\\.rds$", "_map.tsv", RDS); FAM_TSV <- sub("\\.rds$", "_families.tsv", RDS)
write.table(o$map[, c("id","residue","family","column")], MAP_TSV,
            sep = "\t", row.names = FALSE, quote = FALSE)
fam <- data.frame(family = names(o$pssm),
  columns = vapply(o$pssm, nrow, integer(1)),
  human_resid = as.integer(table(factor(o$map$family, levels = names(o$pssm)))),
  stringsAsFactors = FALSE)
fam <- fam[order(-fam$human_resid), ]
write.table(fam, FAM_TSV, sep = "\t", row.names = FALSE, quote = FALSE)

message(sprintf("\n[fix] done: %s families, %s human residues, %.1f MB rds",
                format(length(o$pssm), big.mark=","), format(nrow(o$map), big.mark=","),
                file.size(RDS)/1e6))
message("       ", MAP_TSV, " / ", FAM_TSV)
