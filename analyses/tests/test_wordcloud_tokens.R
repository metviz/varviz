# Check: base-R tokeniser reproduces the tm pipeline it replaced (top-80 terms).
if (!requireNamespace("tm", quietly = TRUE)) {
  cat("SKIP: tm not installed (it is no longer a VarViz dependency)\n"); quit(save = "no")
}

sample_text <- paste(
  "Calcium-sensing receptor (CaSR) senses extracellular calcium (PubMed:12345678).",
  "The encoded protein is a G-protein coupled receptor [PMID 999]; it is required for",
  "parathyroid hormone secretion. Mutations in this gene cause familial hypocalciuric",
  "hypercalcemia 1 and neonatal severe hyperparathyroidism. Provided by RefSeq, Jul 2020.",
  "Don't confuse it with CASR-related disorders; results associated with 2 alleles.")

extra_stops <- c("gene", "protein", "encoded", "human", "results", "associated",
                 "refseq", "provided", "also", "may", "can", "two", "one",
                 "role", "involved", "required", "member", "family", "known")

clean <- function(x) {
  x <- gsub("\\(PubMed:[^)]+\\)", "", x)
  x <- gsub("\\[[^]]*\\]", "", x)
  gsub("[^[:alpha:][:space:]]", " ", x)
}
clean_text <- clean(sample_text)

# --- old path (tm) ---
docs <- tm::Corpus(tm::VectorSource(clean_text))
docs <- tm::tm_map(docs, tm::content_transformer(tolower))
docs <- tm::tm_map(docs, tm::removeNumbers)
docs <- tm::tm_map(docs, tm::removeWords, tm::stopwords("english"))
docs <- tm::tm_map(docs, tm::removeWords, extra_stops)
docs <- tm::tm_map(docs, tm::removePunctuation)
docs <- tm::tm_map(docs, tm::stripWhitespace)
tdm  <- as.matrix(tm::TermDocumentMatrix(docs))
old  <- sort(rowSums(tdm), decreasing = TRUE)

# --- new path (base R, as now in server.R) ---
stop_en <- setdiff(tm::stopwords("english"), grep("'", tm::stopwords("english"), value = TRUE))
words <- strsplit(tolower(clean_text), "[[:space:]]+")[[1]]
words <- words[nchar(words) >= 3 & !words %in% c(stop_en, extra_stops)]
new   <- sort(table(words), decreasing = TRUE)

# tm's TermDocumentMatrix drops 1-2 char terms by default; compare on >=3 chars.
o <- old
n <- new
o <- o[order(names(o))]; n <- n[order(names(n))]

stopifnot(identical(names(o), names(n)))
stopifnot(identical(as.numeric(o), as.numeric(n)))
cat("OK: base-R tokeniser matches tm on", length(o), "terms\n")
