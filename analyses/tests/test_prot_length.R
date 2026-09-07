# Self-check for the protein-length resolution shared by the three harnesses.
#
# fetch_conservation_scores(gene, prot_length) builds its PhyloP/PhastCons
# arrays to exactly prot_length residues, so an under-estimate silently
# truncates every position-indexed annotation past it.
#
# The old order was "AlphaMissense CSV max position, else the highest variant
# position in the universe". For a curated cohort the second branch is not the
# protein length at all: on the 226-variant RASopathy set, 15 of 16 genes had no
# AlphaMissense file and every one of them was truncated -- SOS1 to 1237 aa
# instead of 1333, MRAS to 71 instead of 208.
#
# UniProt's sequence length (pfam_d$length) is authoritative and is already
# fetched before this point, so it now leads. This test evaluates the exact
# block lifted from analyses/ps_final_harness.R against stub inputs.
#
# Run from project root: Rscript analyses/tests/test_prot_length.R

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

src   <- readLines("analyses/ps_final_harness.R", warn = FALSE)
start <- grep("^  prot_length_for_gene <- \\{", src)
stopifnot(length(start) == 1L)
# The block ends at the first line that is exactly two spaces then a closing brace.
end   <- start - 1L + which(src[start:length(src)] == "  }")[1]
block <- paste(src[start:end], collapse = "\n")

extract_pos <- function(p) as.integer(sub("^p\\.[A-Za-z]+([0-9]+).*$", "\\1", p))

resolve <- function(uniprot_len, am_positions, universe_positions) {
  env <- new.env(parent = globalenv())
  env$pfam_d <- if (is.null(uniprot_len)) NULL else list(length = uniprot_len)
  env$am_dt_local <- if (is.null(am_positions)) NULL else
    data.frame(protein_variant = paste0("M", am_positions, "A"), stringsAsFactors = FALSE)
  env$gene_name <- "TESTGENE"
  env$universe <- data.frame(gene = rep("TESTGENE", length(universe_positions)),
                             p_notation = paste0("p.Met", universe_positions, "Ala"),
                             stringsAsFactors = FALSE)
  env$extract_pos <- extract_pos
  out <- capture.output(v <- eval(parse(text = block), envir = env))
  list(value = v, log = paste(out, collapse = "\n"))
}

# 1. UniProt wins even when a shorter universe maximum is available.
#    This is the SOS1 case: 1333 aa, highest curated variant at 846.
r <- resolve(1333, NULL, c(108, 552, 846))
check(identical(as.integer(r$value), 1333L), "UniProt length wins over the universe maximum")
check(!grepl("WARNING", r$log), "no warning when UniProt supplies the length")

# 2. UniProt wins over AlphaMissense too; they agree in practice, and when they
#    disagree the authoritative sequence is the right answer.
r <- resolve(764, c(1, 700, 764), c(36, 744))
check(identical(as.integer(r$value), 764L), "UniProt length wins over the AlphaMissense maximum")

# 3. AlphaMissense is the fallback when the UniProt payload has no length.
r <- resolve(NULL, c(1, 500, 1333), c(108, 846))
check(identical(as.integer(r$value), 1333L), "AlphaMissense maximum used when UniProt is absent")

# 4. Both gone: the universe maximum is used, but never silently.
r <- resolve(NULL, NULL, c(108, 552, 846))
check(identical(as.integer(r$value), 846L), "universe maximum is the last resort")
check(grepl("WARNING", r$log) && grepl("truncates", r$log),
      "the last resort warns that annotations past it are truncated")

# 5. The old behaviour is what this guards against: had the universe maximum won
#    for SOS1, conservation would have been built to 846 rather than 1333.
r <- resolve(1333, NULL, c(846))
check(as.integer(r$value) > 846L, "a curated cohort no longer shortens the protein")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
