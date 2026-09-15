# Self-check for the Pathogenic-boundary cap on AlphaMissense corroboration.
#
# AlphaMissense may raise a PP3 level another tool assigned, but it may not be
# the evidence that carries a variant across the Pathogenic threshold. The cap
# undoes only the crossing: a raise that leaves the variant in the same bin, or
# moves it between lower bins, stands.
#
# Run from project root: Rscript analyses/tests/test_am_cap_pathogenic.R
suppressMessages(source("server.R"))

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (ok) "ok" else "FAIL", msg))
  if (!ok) fails <<- fails + 1L
}

# The cap as server.R applies it, in isolation from the fetch machinery.
apply_cap <- function(tags, raised_from, on = TRUE) {
  if (!on || is.na(raised_from)) return(tags)
  if (!identical(classify_acmg(tags)$classification, "Pathogenic")) return(tags)
  rev <- c("1"="PP3","2"="PP3_moderate","3"="PP3_moderate_plus","4"="PP3_strong")[as.character(raised_from)]
  without <- unique(c(tags[!grepl("^PP3", tags)], unname(rev)))
  if (!identical(classify_acmg(without)$classification, "Pathogenic")) without else tags
}

# A variant at 8 points with PP3_moderate, raised to PP3_strong, reaches 10.
base <- c("PM1_strong","PM2","PP2")
check(classify_acmg(c(base,"PP3_moderate"))$pts == 8L, "PM1_strong+PM2+PP2+PP3_moderate = 8 pts")
check(classify_acmg(c(base,"PP3_strong"))$pts == 10L, "the same tags with PP3_strong = 10 pts")
check(identical(classify_acmg(c(base,"PP3_strong"))$classification, "Pathogenic"),
      "10 pts is Pathogenic, so the raise crosses the boundary")

capped <- apply_cap(c(base,"PP3_strong"), 2L)
check(identical(classify_acmg(capped)$classification, "Likely Pathogenic"),
      "cap reverts the raise, leaving Likely Pathogenic")
check("PP3_moderate" %in% capped && !("PP3_strong" %in% capped),
      "the reverted tag is the level the other tool assigned, not a removal")

# A raise that does not cross the boundary must survive untouched.
low <- c("PM2","PP2")
check(identical(apply_cap(c(low,"PP3_strong"), 1L), c(low,"PP3_strong")),
      "a raise that lands below Pathogenic is left alone")

# A variant already Pathogenic without any AlphaMissense help keeps its call.
strong <- c("PS1","PM1_strong","PM2","PP2")
check(identical(classify_acmg(c(strong,"PP3_strong"))$classification, "Pathogenic"),
      "PS1-backed variant is Pathogenic with PP3_strong")
kept <- apply_cap(c(strong,"PP3_strong"), 2L)
check(identical(classify_acmg(kept)$classification, "Pathogenic"),
      "cap does not fire when the variant is Pathogenic without the raise")

# The option must be switchable, so the pre-cap behaviour stays measurable.
check(identical(apply_cap(c(base,"PP3_strong"), 2L, on = FALSE), c(base,"PP3_strong")),
      "with the option off, the raise stands")
check(isTRUE(getOption("varviz.am_cap_pathogenic", TRUE)), "the cap defaults on")

cat(sprintf("\n%s: %d failure(s)\n", if (fails == 0L) "PASS" else "FAIL", fails))
quit(status = if (fails == 0L) 0L else 1L)
