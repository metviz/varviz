# Self-check for acmg_blind(), the application's Pass-Blind scorer.
#
# The dual-pass output was claimed in the manuscript before it existed in the
# application: strip_clinvar_tags() was written for the benchmark harness and
# appeared in server.R only inside comments. acmg_blind() wires it in, so the
# card and the TSV export now report a call with and without ClinVar-derived
# evidence. This checks the wiring and the pathway-dependent PM1 handling.
#
# Run from project root: Rscript analyses/tests/test_acmg_blind.R

src <- readLines("server.R", warn = FALSE)
grab <- function(start_pat, end_pat) {
  i <- grep(start_pat, src)[1]
  j <- i - 1 + grep(end_pat, src[i:length(src)])[1]
  paste(src[i:j], collapse = "\n")
}
source("analyses/lib/clinvar_blind.R")
eval(parse(text = grab("^ACMG_TAG_PTS <- c\\(", "^\\)")))
eval(parse(text = grab("^PP1_PP4_CAP <- ", "^PP1_PP4_CAP <- ")))
eval(parse(text = grab("^classify_acmg <- function\\(tags_vec\\) \\{", "^\\}")))
# acmg_blind() applies the Pathogenic-boundary cap to its own tag set, so the
# cap has to come across too. This extraction is by regex rather than by
# sourcing server.R, so a new dependency is invisible until it is named here.
eval(parse(text = grab("^cap_am_pathogenic <- function", "^\\}")))
eval(parse(text = grab("^acmg_blind <- function\\(tags_vec", "^\\}")))

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# 1. A ClinVar-derived criterion is withheld and the score drops by its weight.
b <- acmg_blind(c("PM1_strong", "PM2", "PP2", "PM5_supporting", "PP3"), "uniprot_domain")
check(!("PM5_supporting" %in% b$tags), "PM5_supporting is withheld")
check(identical(b$withheld, "PM5_supporting"), "the withheld criterion is reported back")
check(b$pts == 7L, "score drops from 8 to 7 when PM5_supporting is withheld")

# 2. PM1 from a non-circular pathway survives the blinding.
b <- acmg_blind(c("PM1", "PM2", "PP2", "PP3_moderate"), "ccrs")
check("PM1" %in% b$tags, "PM1 from CCRS survives blinding")
check(length(b$withheld) == 0, "nothing withheld when no criterion is ClinVar-derived")

# 3. PM1 that fired only from the ClinVar hotspot is withheld entirely.
b <- acmg_blind(c("PM1", "PM2", "PP2", "PP3_moderate"), "clinvar_hotspot")
check(!("PM1" %in% b$tags), "PM1 from the ClinVar hotspot alone is withheld")
check(identical(as.character(b$classification), "VUS-High"),
      "the CASR p.Thr676Arg case drops from Likely Pathogenic to VUS-High")

# 4. A hotspot upgrade is undone but the non-circular base strength is kept.
b <- acmg_blind(c("PM1_strong", "PM2", "PP2", "PP3_moderate"), "mds+hotspot_upgrade(PM1)")
check("PM1" %in% b$tags && !("PM1_strong" %in% b$tags),
      "PM1_strong demotes to PM1 when the hotspot supplied the upgrade")

# 5. A variant with no ClinVar evidence classifies identically under both passes.
tags <- c("PM1", "PM2", "PP2", "PP3_moderate")
f <- classify_acmg(tags); b <- acmg_blind(tags, "mds")
check(identical(f$classification, b$classification) && identical(f$pts, b$pts),
      "SNCA p.G14R: both passes agree when ClinVar holds nothing to withhold")

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
