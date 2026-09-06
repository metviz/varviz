# Self-check for classify_acmg(): classification comes from summed Tavtigian
# points only. No Richards rule ladder, so
#   * a reduced-strength tag counts at its assigned strength, not its prefix
#     (reviewer: PS1_supporting + PM2 used to reach Likely Pathogenic at 3 pts);
#   * conflicting evidence sums instead of being ordered away
#     (reviewer: PS1 + PS2 + BA1 used to return Pathogenic at 0 pts);
#   * PM2 is Supporting (+1) per ClinGen SVI 2020-09-04; PVS1 + PM2 = 9 = LP;
#   * PP5 / BP6 do not score (retired by ClinGen SVI 2018).
# Bands (Tavtigian 2020): P >= 10, LP 6-9, VUS 0-5, LB -1..-6, B <= -7;
# BA1 is stand-alone Benign.
#
# Extracts the real classify_acmg out of server.R (same trick as test_pp4.R).
# Run from project root: Rscript analyses/tests/test_acmg_points_only.R

src <- readLines("server.R")
beg <- grep("^classify_acmg <- function", src)
end <- beg + which(src[beg:length(src)] == "}")[1] - 1
eval(parse(text = src[beg:end]))
eval(parse(text = src[grep("^PP1_PP4_CAP <- ", src)]))
pb <- grep("^ACMG_TAG_PTS <- c\\($", src); pe <- pb - 1 + which(src[pb:length(src)] == ")")[1]
eval(parse(text = paste(src[pb:pe], collapse = "\n")))                # points table classify_acmg reads

fails <- 0L
check <- function(tags, class, pts, note = "") {
  r <- classify_acmg(tags)
  ok <- identical(r$classification, class) && r$pts == pts
  cat(sprintf("  [%s] %-38s -> %-17s %3d pts  %s\n", if (ok) "ok" else "FAIL",
              paste(tags, collapse = "+"), r$classification, r$pts,
              if (ok) note else sprintf("EXPECTED %s %d", class, pts)))
  if (!ok) fails <<- fails + 1L
}

cat("reviewer cases\n")
check(c("PS1_supporting", "PM2"),     "VUS-Mid",           2, "was LP at 3")
check(c("PS1_moderate", "PS2"),       "Likely Pathogenic", 6, "was P at 6")
check(c("PS1", "PS2", "BA1"),         "Benign",            0, "BA1 stand-alone; was P at 0")
check(c("PS1", "PM2", "BS1", "BS2"),  "Likely Benign",    -3, "was LP at -2")
check(c("PP3_strong", "PM2"),         "VUS-High",          5, "was LP at 6")

cat("ClinGen PM2_Supporting\n")
check("PM2",                          "VUS-Low",           1)
check(c("PVS1", "PM2"),               "Likely Pathogenic", 9, "SVI: very strong + supporting = LP")
check(c("PVS1", "PS1"),               "Pathogenic",       12)

cat("retired criteria do not score\n")
check("PP5",                          "VUS-Low",           0)
check("BP6",                          "VUS-Low",           0)
check(c("PM1", "PP3", "BP6"),         "VUS-Mid",           3, "BP6 no longer -4")

cat("band edges\n")
check(c("PS1", "PM1", "PP3"),         "Likely Pathogenic", 7)
check(c("PS1", "PS3", "PM1"),         "Pathogenic",       10)
check(c("PS1", "PP1"),                "VUS-High",          5)
check("BP4",                          "Likely Benign",    -1, "Tavtigian LB starts at -1")
check(c("BS1", "BP4", "BP1", "BP7"),  "Benign",           -7)
check(character(0),                   "VUS-Low",           0)

cat("PP1/PP4 locus cap still applies\n")
check(c("PP1_strong", "PP4"),         "VUS-High",          5, "4 + 1 = 5, at cap")
check(c("PP1_strong", "PP1_moderate", "PP4"), "VUS-High",  5, "7 capped to 5")

cat("rule string reports the score\n")
r <- classify_acmg(c("PS1", "PP1", "PP4"))
ok <- identical(r$rule, "score 6 pts"); if (!ok) fails <- fails + 1L
cat(sprintf("  [%s] rule = %s\n", if (ok) "ok" else "FAIL", r$rule))

cat(sprintf("\n%s: %d failure(s)\n", if (fails) "FAIL" else "PASS", fails))
quit(status = if (fails) 1L else 0L)
