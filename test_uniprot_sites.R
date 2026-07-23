# Self-check for the UniProt functional-site evidence path.
# Run: Rscript test_uniprot_sites.R
#
# ponytail: slices the real block out of server.R and evaluates it, rather than
# restating the rules here — so a drift in the regex or the type list fails this.
# The eval() parses fixed line ranges of our own in-repo server.R (a trusted
# build artifact, never user input), same approach as test_clinvar_star_weighting.R.
src <- readLines("server.R", warn = FALSE)

# ── The API request actually asks for the fields we score on ──────────────
fields_line <- src[grep("req_url_query\\(fields = \"ft_signal", src)]
stopifnot(
  length(fields_line) == 1,
  grepl("ft_binding", fields_line),
  grepl("ft_site",    fields_line),
  grepl("ft_act_site", fields_line),
  grepl("ft_mutagen",  fields_line)
)

# ── UniProt's own type names map to the short codes the scorer matches on ──
stopifnot(
  any(grepl('"Binding site"\\s*=\\s*"binding"', src)),
  any(grepl('"Site"\\s*=\\s*"site"',            src)),
  # A binding site carries its ligand outside `description`; without this the
  # label would be empty and the PM1 sentence would name nothing.
  any(grepl('f\\$ligand\\$name', src))
)

# ── Evaluate the real site / mutagenesis block ────────────────────────────
beg <- grep("^    uniprot_site_desc  <- \"\"$", src)
end <- grep("^    # --- Density values at this position ---$", src) - 1L
stopifnot(length(beg) == 1, length(end) == 1, end > beg)
block <- paste(src[beg:end], collapse = "\n")

# Drive the block with a synthetic feature table; returns the three PTM slots
# plus the site label, exactly as the surrounding code would see them.
run_block <- function(features, position) {
  e <- new.env(parent = globalenv())
  e$uniprot_data <- features
  e$pos          <- position
  e$ptm_info <- ""; e$ptm_acmg <- ""; e$ptm_strength <- ""
  eval(parse(text = block), envir = e)
  list(site = e$uniprot_site_desc, mutagen = e$uniprot_mutagen,
       ptm_acmg = e$ptm_acmg, ptm_info = e$ptm_info)
}
ft <- function(type, start, end = start, description = "")
  data.frame(type = type, start = start, end = end,
             description = description, stringsAsFactors = FALSE)

# Active site: the criterion's canonical example, exact residue.
r <- run_block(ft("act_site", 145, 145, "Proton acceptor"), 145)
stopifnot(r$site == "Proton acceptor", r$ptm_acmg == "")

# Binding sites may span a short range; a residue inside it still counts.
r <- run_block(ft("binding", 87, 89, "Binds Ca(2+)"), 88)
stopifnot(r$site == "Binds Ca(2+)")
stopifnot(run_block(ft("binding", 87, 89, "Binds Ca(2+)"), 90)$site == "")

# A PTM is not a functional-site hit — it travels the separate PS3 channel.
stopifnot(run_block(ft("mod_res", 129, 129, "Phosphoserine"), 129)$site == "")

# Mutagenesis: only a reported LOSS is evidence for pathogenicity.
damaging <- c("Loss of activity", "Abolishes ligand binding",
              "Strongly reduces receptor signaling", "Impairs trafficking")
for (d in damaging) {
  r <- run_block(ft("mutagen", 50, 50, d), 50)
  stopifnot(r$ptm_acmg == "PS3_supporting", grepl("Mutagenesis:", r$ptm_info))
}
benign_phrasing <- c("No effect on binding", "No loss of activity",
                     "No significant change in signaling", "Does not affect folding")
for (d in benign_phrasing) {
  r <- run_block(ft("mutagen", 50, 50, d), 50)
  stopifnot(r$mutagen == "", r$ptm_acmg == "")
}

# An existing PTM call must survive alongside a mutagenesis note.
e <- new.env(parent = globalenv())
e$uniprot_data <- ft("mutagen", 50, 50, "Loss of activity")
e$pos <- 50; e$ptm_info <- "Phosphoserine"
e$ptm_acmg <- "PS3_supporting"; e$ptm_strength <- "Strong functional site"
eval(parse(text = block), envir = e)
stopifnot(grepl("^Phosphoserine; Mutagenesis: ", e$ptm_info))

# No UniProt data at all must be inert, not an error.
stopifnot(run_block(ft("act_site", 1, 1)[0, ], 145)$site == "")
stopifnot(run_block(NULL, 145)$site == "")

# ── The site pathway is consulted before the region-level ones ────────────
# PM1 is a claim about location; a residue that IS the catalytic position is
# more specific than one merely inside a constrained region, so it wins.
site_branch <- grep("if \\(nchar\\(uniprot_site_desc\\) > 0\\) \\{", src)
ccrs_branch <- grep("\\} else if \\(ccrs_pct >= 90\\) \\{", src)
stopifnot(length(site_branch) == 1, length(ccrs_branch) == 1,
          site_branch < ccrs_branch,
          any(grepl('pm1_pathway_val <- "uniprot_site"', src)))

# The pathway label must reach the export, and the comment must be able to name
# the residue in both the first-pass and the live-regeneration call sites.
stopifnot(
  any(grepl("UniProt_Site = uniprot_site_desc", src)),
  length(grep("uniprot_site\\s*=", src)) >= 3,   # signature + 2 call sites
  length(grep("ptm_note\\s*=", src))     >= 3
)

# ── Real mutagenesis evidence stands the PS3 proxy down ───────────────────
# PP3_strong is applied as a stand-in for PS3_supporting when no functional
# data exists. Once UniProt supplies published mutagenesis, both would occupy
# the same evidence slot, so the proxy must be gated off.
stopifnot(any(grepl(
  "if \\(am_high && revel_dam && nchar\\(uniprot_mutagen\\) == 0\\)", src)))

cat("uniprot sites: all checks pass\n")
