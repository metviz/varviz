# Self-check for ClinVar star-weighted PS1/PM5. Run: Rscript test_clinvar_star_weighting.R
# ponytail: extracts the real classify_acmg + strength helpers out of server.R rather
# than mirroring them, so this fails if the points map or suffix handling drifts.
# The eval() below parses fixed line ranges of our own in-repo server.R — a trusted
# build artifact, never user input — which is what makes testing the real code safe.
src <- readLines("server.R", warn = FALSE)
eval(parse(text = src[grep("^PP1_PP4_CAP <- ", src)]))
beg <- grep("^classify_acmg <- function", src)
end <- beg - 1 + which(src[beg:length(src)] == "}")[1]
eval(parse(text = paste(src[beg:end], collapse = "\n")))

# ── Star-weighted tiers score at their declared strengths ─────────────────
stopifnot(
  classify_acmg("PS1")$pts            == 4,   # >=2 stars: Strong
  classify_acmg("PS1_moderate")$pts   == 2,   # 1 star
  classify_acmg("PS1_supporting")$pts == 1,   # 0 stars / unparseable
  classify_acmg("PM5")$pts            == 2,   # >=2 stars: Moderate
  classify_acmg("PM5_supporting")$pts == 1,   # below 2 stars: Supporting
  # A weak-review PM5 must be worth strictly less than a well-reviewed one,
  # which is the whole point of star-weighting it.
  classify_acmg("PM5_supporting")$pts < classify_acmg("PM5")$pts
)

# ── The emission logic is actually star-gated, not flat ───────────────────
stopifnot(
  any(grepl('pm5_tag <- if \\(!is\\.na\\(stars_n\\) && stars_n >= 2\\) "PM5" else "PM5_supporting"', src)),
  # Both PM5 emission sites must use the star-weighted tag, never a bare "PM5".
  length(grep('acmg_tags, pm5_tag', src)) == 2,
  length(grep('acmg_tags, "PM5"', src)) == 0,
  # stars_n must be hoisted above the exact/position branches so both can read it.
  grep('stars_n <- suppressWarnings', src)[1] < grep('if \\(clinvar_match_type == "exact"', src)[1]
)

# ── Strength suffixes resolve for display, tier and tooltip ───────────────
# Regression guard: tag_display previously stripped only _strong/_moderate, so
# _supporting tags rendered as raw names and missed the tooltip switch.
sl  <- grep("strength_label <- function", src)
tdl <- grep("tag_display <- function", src)
sl_end <- sl - 1 + which(trimws(src[sl:length(src)]) == "}")[1]
eval(parse(text = paste(src[sl:sl_end], collapse = "\n")))
eval(parse(text = src[tdl]))
for (base in c("PS1", "PS3", "PM5")) {
  stopifnot(tag_display(paste0(base, "_supporting")) == base)
}
stopifnot(
  tag_display("PP3_strong")   == "PP3",
  tag_display("PP1_moderate") == "PP1",
  tag_display("PM2")          == "PM2",          # untouched when unsuffixed
  strength_label("PM5_supporting") == "Supporting",
  strength_label("PS1_supporting") == "Supporting",
  strength_label("PS3_supporting") == "Supporting",
  strength_label("PM5")            == "Moderate",
  strength_label("PS1")            == "Strong"
)
cat("clinvar star weighting: all checks pass\n")
