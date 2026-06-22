# dolphin_bulk.R — per-gene DOLPHIN precompute.
#
# Sequentially queries DOLPHIN (Corcuff et al. 2023, Front. Bioinform.
# 3:1127341) for every missense variant in a gene, respecting the documented
# 60-req/min server rate limit, and writes a per-gene TSV cache to
# analyses/raw/dolphin/by_gene/{GENE}.tsv. Resumable across crashes (rereads
# the TSV, skips variants already done). Used by the dual-pass classification
# harness so that the 60-req/min API is hit at most once per gene-variant pair.
#
# Output schema (TSV, no quoting, tab-separated):
#   p_notation    : the input p-notation (whatever format the caller passed)
#   p_single      : normalised single-letter form (e.g. "S296N")
#   ensembl       : Ensembl transcript ID used for the query
#   pm1           : TRUE/FALSE — does DOLPHIN fire PM1?
#   acmg_raw      : pipe-delimited concatenation of acmg fields across results[]
#                   (preserves DOLPHIN's PM2 / BP8 tags for audit, even though
#                   the engine consumes only PM1)
#   pfam_id       : pipe-delimited Pfam IDs of all hit domains
#   score_delta   : pipe-delimited score_delta values across results[]
#   fetched_at    : ISO 8601 UTC timestamp of the API call
#
# Dependencies (sourced relative to project root):
#   analyses/lib/dolphin.R
#
# Usage from harness:
#   source("analyses/lib/dolphin_bulk.R")
#   df <- fetch_dolphin_gene("CASR", c("p.S296N", "p.G509R", ...))

if (!exists("fetch_dolphin", mode = "function")) {
  .dbulk_here <- function(rel) if (file.exists(rel)) rel else file.path("..", rel)
  source(.dbulk_here("analyses/lib/dolphin.R"), local = FALSE)
}

#' @return data.frame with one row per input variant, columns as above.
fetch_dolphin_gene <- function(gene,
                                p_notations,
                                ensembl = NULL,
                                cache_dir = "analyses/raw/dolphin/by_gene",
                                rate_sleep_s = 1.05,
                                flush_every = 50,
                                verbose = TRUE) {
  gene <- trimws(as.character(gene))
  if (nchar(gene) == 0) stop("fetch_dolphin_gene: empty gene")
  p_notations <- unique(trimws(as.character(p_notations)))
  p_notations <- p_notations[nchar(p_notations) > 0]
  if (length(p_notations) == 0) {
    return(data.frame(p_notation=character(), p_single=character(),
                      ensembl=character(), pm1=logical(),
                      acmg_raw=character(), pfam_id=character(),
                      score_delta=character(), fetched_at=character(),
                      stringsAsFactors = FALSE))
  }

  if (is.null(ensembl)) {
    ensembl <- dolphin_canonical_enst(gene)
    if (is.na(ensembl) || nchar(ensembl) == 0) {
      stop("fetch_dolphin_gene: could not resolve canonical ENST for ", gene)
    }
  }

  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(cache_dir, paste0(gene, ".tsv"))

  # Resume: read existing TSV if present.
  done <- if (file.exists(out_path)) {
    tryCatch(
      utils::read.table(out_path, sep = "\t", header = TRUE,
                        stringsAsFactors = FALSE, quote = "",
                        comment.char = ""),
      error = function(e) {
        message("[DOLPHIN bulk] cannot parse existing ", out_path,
                " (", e$message, "); starting fresh")
        NULL
      }
    )
  } else {
    NULL
  }

  empty_df <- data.frame(p_notation=character(), p_single=character(),
                         ensembl=character(), pm1=logical(),
                         acmg_raw=character(), pfam_id=character(),
                         score_delta=character(), fetched_at=character(),
                         stringsAsFactors = FALSE)
  if (is.null(done) || nrow(done) == 0) done <- empty_df

  todo <- setdiff(p_notations, done$p_notation)

  if (verbose) {
    message(sprintf("[DOLPHIN bulk] gene=%s enst=%s todo=%d already=%d",
                    gene, ensembl, length(todo), nrow(done)))
  }

  if (length(todo) == 0) return(done)

  flush_buf <- list()
  flush_to_disk <- function() {
    if (length(flush_buf) == 0) return(invisible(NULL))
    chunk <- do.call(rbind, flush_buf)
    done <<- rbind(done, chunk)
    utils::write.table(done, out_path, sep = "\t", row.names = FALSE,
                       quote = FALSE, na = "")
    flush_buf <<- list()
  }

  on.exit(flush_to_disk(), add = TRUE)

  for (i in seq_along(todo)) {
    p <- todo[i]
    p_single <- .dolphin_to_single_letter(p)
    resp <- tryCatch(
      fetch_dolphin(ensembl = ensembl, p_notation = p, timeout_s = 15),
      error = function(e) {
        message("[DOLPHIN bulk] ", gene, " ", p, " error: ", e$message)
        NULL
      }
    )

    pm1 <- isTRUE(dolphin_fires_pm1(resp))

    # Extract per-domain summary for audit (works for the array case;
    # absent / string-sentinel cases yield empty strings)
    acmg_raw <- ""
    pfam_id <- ""
    score_delta <- ""
    if (!is.null(resp) && is.list(resp$results) && !is.character(resp$results)) {
      acmg_raw <- paste(vapply(resp$results, function(r) {
        if (is.list(r) && !is.null(r$acmg)) {
          paste(unlist(r$acmg), collapse = ";")
        } else ""
      }, character(1)), collapse = "|")
      pfam_id <- paste(vapply(resp$results, function(r) {
        if (is.list(r) && !is.null(r$pfid)) as.character(r$pfid) else ""
      }, character(1)), collapse = "|")
      score_delta <- paste(vapply(resp$results, function(r) {
        if (is.list(r) && !is.null(r$score$score_delta)) {
          format(as.numeric(r$score$score_delta), nsmall = 3)
        } else ""
      }, character(1)), collapse = "|")
    }

    row <- data.frame(
      p_notation  = p,
      p_single    = p_single,
      ensembl     = ensembl,
      pm1         = pm1,
      acmg_raw    = acmg_raw,
      pfam_id     = pfam_id,
      score_delta = score_delta,
      fetched_at  = format(Sys.time(), tz = "UTC", "%Y-%m-%dT%H:%M:%SZ"),
      stringsAsFactors = FALSE
    )
    flush_buf[[length(flush_buf) + 1]] <- row

    if (length(flush_buf) >= flush_every || i == length(todo)) {
      flush_to_disk()
      if (verbose) {
        message(sprintf("[DOLPHIN bulk] %s: %d/%d done (last=%s pm1=%s)",
                        gene, nrow(done), length(p_notations), p, pm1))
      }
    }

    if (i < length(todo)) Sys.sleep(rate_sleep_s)
  }

  done
}
