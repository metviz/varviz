# Guards the on-disk gnomAD payload cache.
#
# Why this exists. gnomAD was the only per-gene source with no disk cache at all.
# ClinVar is primed to .rds, UniProt is stored as its verbatim payload, Ensembl
# exons are cached, but gnomAD was fetched live from the GraphQL API on every run
# and kept only in an in-memory list for the life of the process. Allele
# frequency drives PM2, BA1, BS1 and BS2, and the population-genetic panel is
# built from gnomAD singletons, so no run was reproducible with respect to it:
# re-running the same code on the same day could read different frequencies than
# the run whose numbers were published.
#
# The payload is stored verbatim, as the response text, not as a re-serialised
# parse. A re-serialisation is a derived view: field order, number formatting and
# null handling all shift, so it cannot be diffed against a later fetch to prove
# the source did or did not change. The ClinVar cache is the cautionary case in
# this repository, reduced to ten curated columns, which makes every question
# outside those ten a re-fetch.
#
# Round-trip fidelity is the assertion that matters here. If the bytes that come
# back are not the bytes that went in, the cache cannot be used as evidence of
# what the source said.
#
# Run from project root: Rscript analyses/tests/test_gnomad_cache.R

resolve <- function(rel) if (file.exists(rel) || dir.exists(rel)) rel else file.path("..", rel)
src_path <- resolve("server.R")
src <- readLines(src_path, warn = FALSE)

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# Pull the three helpers out of server.R rather than sourcing the whole app.
# Same readLines/parse/eval idiom as test_acmg_posterior.R; the input is a fixed
# line range of this repository's own tracked source, and reading it from
# server.R is the point, so the assertions land on the shipped functions.
want <- c("gnomad_cache_dir", "gnomad_cache_path", "gnomad_cache_read", "gnomad_cache_write")
missing <- character(0)
for (fn in want) {
  beg <- grep(paste0("^", fn, " <- function"), src)
  if (length(beg) != 1L) { missing <- c(missing, fn); next }
  end <- beg - 1L + which(src[beg:length(src)] == "}")[1]
  eval(parse(text = paste(src[beg:end], collapse = "\n")), envir = globalenv())
}
if (length(missing)) {
  cat("FAIL: not defined exactly once in server.R: ", paste(missing, collapse = ", "), "\n", sep = "")
  quit(status = 1L)
}

tmp <- file.path(tempdir(), paste0("gnomad_cache_", as.integer(runif(1, 1e6, 9e6))))
on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
old <- getOption("varviz.gnomad_cache_dir")
options(varviz.gnomad_cache_dir = tmp)
on.exit(options(varviz.gnomad_cache_dir = old), add = TRUE)

cat("-- path resolution --\n")
check(identical(gnomad_cache_dir(), tmp),
      "the cache directory follows varviz.gnomad_cache_dir")
check(identical(gnomad_cache_path("CASR"), file.path(tmp, "CASR.json")),
      "the payload path is <dir>/<GENE>.json")

cat("-- a miss is NULL, not an error --\n")
check(is.null(gnomad_cache_read("NOSUCHGENE")),
      "reading an absent gene returns NULL")

cat("-- verbatim round trip --\n")
# Deliberately awkward: trailing spaces, a tab, CRLF, unicode, deep nesting and
# a number that any re-serialisation would reformat. A parse-and-rewrite cache
# fails this; a byte cache passes it.
payload <- paste0(
  '{"data":{"gene":{"variants":[',
  '{"variant_id":"3-122283880-G-A","hgvsp":"p.Val104Ile",',
  '"exome":{"ac":1,"an":1461886,"af":6.840000000000001e-07,"ac_hom":0},',
  '"genome":null,"note":"tab\\there, trailing space  ","sym":"\\u00e5\\u00df"}',
  ']}}}\r\n'
)
gnomad_cache_write("CASR", payload,
                   meta = list(url = "https://gnomad.broadinstitute.org/api",
                               dataset = "gnomad_r4", http_status = 200L,
                               n_variants = 1L))
back <- gnomad_cache_read("CASR")
check(identical(back, payload),
      "the payload reads back byte-identical, including CRLF, tab and unicode escapes")
check(!grepl("6.84e-07", back, fixed = TRUE),
      "the float is not reformatted (a re-serialising cache would shorten it)")

cat("-- the sidecar records provenance --\n")
mp <- file.path(tmp, "CASR.meta.json")
check(file.exists(mp), "a <GENE>.meta.json sidecar is written")
if (file.exists(mp)) {
  meta <- paste(readLines(mp, warn = FALSE), collapse = "")
  for (k in c("url", "dataset", "http_status", "n_variants", "sha256", "fetched_at", "bytes")) {
    check(grepl(paste0('"', k, '"'), meta, fixed = TRUE),
          sprintf("the sidecar records %s", k))
  }
  check(grepl("gnomad_r4", meta, fixed = TRUE), "the sidecar names the dataset queried")
}

cat("-- the checksum is over the stored bytes --\n")
if (requireNamespace("digest", quietly = TRUE)) {
  meta <- paste(readLines(mp, warn = FALSE), collapse = "")
  sha <- sub('^.*"sha256"\\s*:\\s*"([0-9a-f]+)".*$', "\\1", meta)
  check(identical(sha, digest::digest(payload, algo = "sha256", serialize = FALSE)),
        "sha256 in the sidecar matches the payload bytes")
} else {
  cat("  [skip] digest not installed; sha256 value not verified\n")
}

cat("-- the fetch path consults the cache before the network --\n")
beg <- grep("^query_gnomad_api <- function", src)
if (length(beg) == 1L) {
  end <- beg - 1L + which(grepl("^}", src[beg:length(src)]))[1]
  body <- src[beg:end]
  ipost <- which(grepl("POST\\(", body))[1]
  iread <- which(grepl("gnomad_cache_read\\(", body))[1]
  iwrite <- which(grepl("gnomad_cache_write\\(", body))[1]
  check(!is.na(iread), "query_gnomad_api calls gnomad_cache_read")
  check(!is.na(iwrite), "query_gnomad_api calls gnomad_cache_write")
  if (!is.na(iread) && !is.na(ipost)) {
    check(iread < ipost, "the disk read happens before the POST, so a primed gene never hits the network")
  }
} else {
  check(FALSE, "query_gnomad_api found exactly once")
}

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
