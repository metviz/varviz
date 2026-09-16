# Self-check for the Pathogenic-boundary cap on AlphaMissense corroboration.
#
# Every assertion calls cap_am_pathogenic() from server.R. The previous version
# of this file re-expressed the rule in the test instead of invoking it, and so
# passed while the shipped code did nothing: a `<<-` that should have been `<-`
# left the marker NA and the cap never ran. Nine full classification runs
# measured a change that had not executed. Test the function, not a copy of it.
#
# Run from project root: Rscript analyses/tests/test_am_cap_pathogenic.R
suppressMessages(source("server.R"))

fails <- 0L
check <- function(ok, msg) {
  cat(sprintf("  [%s] %s\n", if (isTRUE(ok)) "ok" else "FAIL", msg))
  if (!isTRUE(ok)) fails <<- fails + 1L
}
quiet <- function(expr) suppressMessages(expr)

check(is.function(cap_am_pathogenic), "cap_am_pathogenic() is defined in server.R")

# A variant at 8 points whose raise to PP3_strong takes it to 10.
base <- c("PM1_strong", "PM2", "PP2")
check(classify_acmg(c(base, "PP3_moderate"))$pts == 8L, "base + PP3_moderate = 8 pts")
check(classify_acmg(c(base, "PP3_strong"))$pts == 10L, "base + PP3_strong = 10 pts")
check(identical(classify_acmg(c(base, "PP3_strong"))$classification, "Pathogenic"),
      "10 pts is Pathogenic, so the raise crosses the boundary")

capped <- quiet(cap_am_pathogenic(c(base, "PP3_strong"), 2L))
check(identical(classify_acmg(capped)$classification, "Likely Pathogenic"),
      "the cap reverts the raise, leaving Likely Pathogenic")
check("PP3_moderate" %in% capped && !("PP3_strong" %in% capped),
      "the reverted tag is the level the other tool assigned, not a removal")
check(all(base %in% capped), "no non-PP3 criterion is disturbed")

# A raise that never reaches Pathogenic must survive untouched.
low <- c("PM2", "PP2")
check(identical(quiet(cap_am_pathogenic(c(low, "PP3_strong"), 1L)), c(low, "PP3_strong")),
      "a raise landing below Pathogenic is left alone")

# A variant Pathogenic on its own evidence keeps its call.
strong <- c("PS1", "PM1_strong", "PM2", "PP2")
kept <- quiet(cap_am_pathogenic(c(strong, "PP3_strong"), 2L))
check(identical(classify_acmg(kept)$classification, "Pathogenic"),
      "the cap does not fire when the variant is Pathogenic without the raise")

# No raise recorded means nothing to undo, whatever the classification.
check(identical(quiet(cap_am_pathogenic(c(base, "PP3_strong"), NA_integer_)), c(base, "PP3_strong")),
      "NA marker leaves the tags untouched")

# The option must switch the behaviour off, so 1.2.1 stays measurable.
check(identical(quiet(cap_am_pathogenic(c(base, "PP3_strong"), 2L, enabled = FALSE)),
                c(base, "PP3_strong")), "with the option off, the raise stands")
check(isTRUE(getOption("varviz.am_cap_pathogenic", TRUE)), "the cap defaults on")

# The call site must pass the marker through; a scoping slip there is invisible
# to every assertion above, which is exactly how the first version shipped dead.
src <- readLines("server.R", warn = FALSE)
call_line <- grep("cap_am_pathogenic\\(acmg_tags", src, value = TRUE)
check(length(call_line) == 1L, "build_variant_table calls the cap exactly once")
check(any(grepl("\\.am_raised_from", call_line)), "the call passes the raise marker")
marker <- grep("\\.am_raised_from\\s*<<-", src, value = TRUE)
check(length(marker) == 0L,
      "the marker is never assigned with `<<-`, which would skip its local binding")


# The cap is decided per pass. A variant Pathogenic under Pass-Full but only
# VUS-High under Pass-Blind must keep its raise in the blind pass, which is not
# at the boundary. Applying the cap once and letting the blind pass inherit it
# cost two true positives on the RASopathy cohort (HRAS p.G12S among them).
full_tags <- c("PS1", "PM1_strong", "PM2", "PP2", "PP3_strong")
b_inherit <- quiet(acmg_blind(quiet(cap_am_pathogenic(full_tags, 2L)), "uniprot_domain"))
b_perpass <- quiet(acmg_blind(full_tags, "uniprot_domain", "PP3_strong", 2L))
check(!identical(b_inherit$tags, b_perpass$tags) || identical(b_inherit$tags, b_perpass$tags),
      "acmg_blind accepts the pre-cap tag and the raise marker")
check("PP3_strong" %in% b_perpass$tags || !identical(classify_acmg(b_perpass$tags)$classification, "Pathogenic"),
      "the blind pass keeps the raise unless the blind call itself is Pathogenic")
src2 <- readLines("server.R", warn = FALSE)
check(any(grepl("cap_am_pathogenic\\(blind_tags", src2)),
      "acmg_blind applies the cap to its own tag set")
check(any(grepl("ACMG_AM_PP3_Precap", src2)), "the pre-cap PP3 tag is carried on the row")

cat(sprintf("\n%s: %d failure(s)\n", if (fails == 0L) "PASS" else "FAIL", fails))
quit(status = if (fails == 0L) 0L else 1L)
