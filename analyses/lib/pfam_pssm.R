# pfam_pssm.R — position-specific scoring matrices from Pfam alignments.
#
# Reimplements the score construction of Corcuff et al. 2023 (DOLPHIN,
# Front. Bioinform. 3:1127341) so the PM1 domain evidence can be computed
# offline instead of fetched from a third-party API. Their five steps:
#
#   C(p,l)   count of amino acid l at alignment column p
#   F'(p,l)  = (C(p,l) + c * f_l) / (sum_l C(p,l) + c)      pseudo-count c = 1
#   F''(p,l) = F'(p,l) / f_l
#   M(p,l)   = ln F''(p,l)
#   delta    = M(p, mut) - M(p, wt)
#
# f_l is the expected amino-acid frequency, taken from UniProt/Swiss-Prot
# composition (as the paper specifies).
#
# Why domain-scoped rather than a whole-protein PSSM: PM1 is a claim about
# LOCATION ("mutational hot spot or critical functional domain"). A per-protein
# profile scores every residue regardless of domain membership, so it would not
# mean what the criterion says — and it would double-count the alignment signal
# VarViz already consumes as PP3 via PROVEAN.

# Swiss-Prot amino-acid composition (release statistics). Used as the expected
# background f_l; the paper uses the same source.
SWISSPROT_AA_FREQ <- c(
  A = 0.0825, R = 0.0553, N = 0.0406, D = 0.0545, C = 0.0137,
  Q = 0.0393, E = 0.0675, G = 0.0707, H = 0.0227, I = 0.0596,
  L = 0.0966, K = 0.0584, M = 0.0242, F = 0.0386, P = 0.0470,
  S = 0.0656, T = 0.0534, W = 0.0108, Y = 0.0292, V = 0.0687
)
AA20 <- names(SWISSPROT_AA_FREQ)

# Parse a gzipped Stockholm alignment into a named character vector of
# gapped sequences. Stockholm wraps long alignments across blocks, so the same
# sequence id can appear repeatedly and must be concatenated in order.
read_stockholm <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, "rt") else file(path, "rt")
  on.exit(close(con))
  lines <- readLines(con, warn = FALSE)
  lines <- lines[!grepl("^#", lines) & !grepl("^//", lines) & nzchar(trimws(lines))]
  if (!length(lines)) return(character(0))
  parts <- strsplit(trimws(lines), "[[:space:]]+")
  ids <- vapply(parts, `[`, character(1), 1L)
  seqs <- vapply(parts, function(p) paste(p[-1], collapse = ""), character(1))
  tapply(seqs, ids, paste, collapse = "")
}

# Count matrix over alignment columns. Gaps and non-standard residues are not
# counted; upper/lower case in Stockholm both denote residues (lower case marks
# insertions relative to the model) so we fold case.
pssm_counts <- function(aln) {
  chars <- do.call(rbind, strsplit(toupper(unname(aln)), "", fixed = TRUE))
  ncol_aln <- ncol(chars)
  C <- matrix(0, nrow = ncol_aln, ncol = length(AA20),
              dimnames = list(NULL, AA20))
  for (j in seq_len(ncol_aln)) {
    tb <- table(chars[, j])
    tb <- tb[names(tb) %in% AA20]
    if (length(tb)) C[j, names(tb)] <- as.numeric(tb)
  }
  C
}

# The log score matrix M(p,l), following the paper's formulae exactly.
build_pssm <- function(aln, pseudocount = 1) {
  C <- pssm_counts(aln)
  f <- SWISSPROT_AA_FREQ[AA20]
  # F'(p,l) = (C + c*f_l) / (rowSums(C) + c)
  Fprime <- sweep(C, 2, pseudocount * f, `+`)
  Fprime <- Fprime / (rowSums(C) + pseudocount)
  # F''(p,l) = F'/f_l ; M = ln F''
  log(sweep(Fprime, 2, f, `/`))
}

# Map a residue index of one aligned sequence onto its alignment column.
# Stockholm ids carry the residue range (e.g. SYUA_HUMAN/1-131), so residue n
# of the protein is the (n - start + 1)-th non-gap character of that row.
alignment_column_for_residue <- function(aln, id, residue, start = NULL) {
  # aln comes from tapply(), so it is a named array: [[ ]] on an absent name
  # errors rather than returning NULL. Match the name explicitly.
  if (!(id %in% names(aln))) return(NA_integer_)
  s <- aln[[id]]
  if (is.null(s) || is.na(s)) return(NA_integer_)
  if (is.null(start)) {
    rng <- sub("^.*/", "", id)
    start <- suppressWarnings(as.integer(sub("-.*$", "", rng)))
    if (is.na(start)) start <- 1L
  }
  idx <- residue - start + 1L
  if (idx < 1L) return(NA_integer_)
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  nongap <- which(chars %in% c(AA20, tolower(AA20)))
  if (idx > length(nongap)) return(NA_integer_)
  nongap[idx]
}

# delta = M(p, mut) - M(p, wt); negative means the substitution is disfavoured
# at that column, which is the direction DOLPHIN reports.
pssm_delta <- function(M, column, wt, mut) {
  if (is.na(column) || column < 1 || column > nrow(M)) return(NA_real_)
  wt <- toupper(wt); mut <- toupper(mut)
  if (!(wt %in% AA20) || !(mut %in% AA20)) return(NA_real_)
  unname(M[column, mut] - M[column, wt])
}
