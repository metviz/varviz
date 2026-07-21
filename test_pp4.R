# Self-check for the PP4 (phenotype specificity) criterion. Run: Rscript test_pp4.R
# ponytail: extracts the real classify_acmg out of server.R rather than mirroring it,
# so this fails if the points map or the ^PP counting regex drifts.
# The eval() below parses a fixed line range of our own in-repo server.R — a trusted
# build artifact, never user input — which is what makes testing the real function
# (instead of a drift-prone copy) safe here.
src <- readLines("server.R", warn = FALSE)
beg <- grep("^classify_acmg <- function", src)
end <- beg - 1 + which(src[beg:length(src)] == "}")[1]
eval(parse(text = paste(src[beg:end], collapse = "\n")))

# PP4 must be worth exactly 1 point (Supporting) and count toward the PP tally.
base <- classify_acmg(c("PM2"))
withp <- classify_acmg(c("PM2", "PP4"))
stopifnot(
  withp$pts - base$pts == 1,                       # Supporting = +1
  classify_acmg("PP4")$pts == 1,
  classify_acmg(character(0))$pts == 0,
  # PP4 is Supporting-only: no strength-suffixed variants should score.
  classify_acmg("PP4_moderate")$pts == 0,
  classify_acmg("PP4_strong")$pts == 0,
  # PP4 must reach the ^PP tally, not just the points total: 1 PS + 1 PP sits at
  # VUS-High, and PP4 as the 2nd PP is what trips the Richards "1 PS + 2 PP" rule.
  classify_acmg(c("PS1","PP1"))$class == "VUS-High",
  classify_acmg(c("PS1","PP1","PP4"))$class == "Likely Pathogenic",
  classify_acmg(c("PS1","PP1","PP4"))$rule == "1 PS + 2 PP"
)

# The UI control must exist and be wired to the tag, or the criterion is unreachable.
stopifnot(
  any(grepl('paste0\\("pp4_", i\\)', src)),                    # selectInput + readers
  any(grepl('identical\\(card_pp4, "pp4"\\)', src)),           # card tag emission
  any(grepl('vtbl\\$PP4_Applied', src)),                       # export column
  any(grepl('pp4 <- vtbl\\$PP4_Applied', src))                 # folded into Final_ACMG_Tags
)
cat("PP4: all checks pass\n")
