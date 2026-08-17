# Self-check for the MDS (Missense Disfavour Score) PM1 Path 4 wiring in server.R.
# Run: Rscript test_mds_pm1.R
#
# Three parts: (1) server.R is wired to use MDS as Path 4 and tiers it by LR+
# (-4 Moderate, -8 Moderate-plus, -12 Strong); (2) the tier DECISION is exercised
# at its boundaries against the real thresholds extracted from server.R; (3) the
# lookup helper reproduces DOLPHIN on the shipped table and real reference deltas
# land in the expected tier. Slices server.R as a trusted in-repo build artifact,
# never user input.
src <- readLines("server.R", warn = FALSE)

# ── server.R actually wires MDS as Path 4, tiered ──────────────────────────
stopifnot(
  any(grepl('source\\("analyses/lib/pssm_lookup.R"', src)),
  any(grepl('MDS_TABLE <- tryCatch\\(pssm_table_load\\("data/pfam_pssm_human.rds"', src)),
  any(grepl('^MDS_PM1_THRESHOLD *<-', src)),
  any(grepl('^MDS_PM1_MODPLUS_THRESHOLD *<-', src)),
  any(grepl('^MDS_PM1_STRONG_THRESHOLD *<-', src)),
  # entry-name bridge from the UniProt JSON already fetched
  any(grepl('mds_entry_name <-', src)) && any(grepl('pfam_data\\$uniProtkbId', src)),
  # the Path 4 gate, the three tier branches, and the unavailable sentinel
  any(grepl('pm1_pathway_val <- "mds"', src)),
  any(grepl('pm1_pathway_val <- "mds_unavailable"', src)),
  any(grepl('mds_val <= MDS_PM1_THRESHOLD', src)),
  any(grepl('mds_val <= MDS_PM1_STRONG_THRESHOLD', src)),
  any(grepl('mds_val <= MDS_PM1_MODPLUS_THRESHOLD', src)),
  any(grepl('"PM1_strong"', src)),
  any(grepl('"PM1_moderate_plus"', src)),
  # tier -> points map (Strong +4, Moderate-plus +3, Moderate +2)
  any(grepl('PM1_strong=4', src)),
  any(grepl('PM1_moderate_plus=3', src)),
  any(grepl('PM1=2', src)),
  # score is exported
  any(grepl('MDS_Score = ', src)),
  # the retired DOLPHIN call is no longer the runtime Path 4
  !any(grepl('pm1_pathway_val <- "dolphin"', src))
)

# ── Tier decision exercised at its boundaries against server.R's real thresholds ──
num_after <- function(pat) {
  rhs <- sub(".*<-", "", grep(pat, src, value = TRUE)[1])
  as.numeric(regmatches(rhs, regexpr("-?[0-9]+\\.?[0-9]*", rhs)))
}
base_thr    <- num_after("^MDS_PM1_THRESHOLD *<-")
modplus_thr <- num_after("^MDS_PM1_MODPLUS_THRESHOLD *<-")
strong_thr  <- num_after("^MDS_PM1_STRONG_THRESHOLD *<-")
stopifnot(base_thr > modplus_thr, modplus_thr > strong_thr)   # -4 > -8 > -12

# mirrors the server.R Path-4 precedence (strong checked first); the static greps
# above assert that precedence and the emitted tags exist in server.R.
mds_tier <- function(mds) {
  if (is.na(mds) || mds > base_thr) return(NA_character_)
  if (mds <= strong_thr)  return("PM1_strong")
  if (mds <= modplus_thr) return("PM1_moderate_plus")
  "PM1"
}
stopifnot(
  is.na(mds_tier(base_thr + 0.01)),                       # just above -4 -> no fire
  mds_tier(base_thr)           == "PM1",                  # -4  -> Moderate
  mds_tier(modplus_thr + 0.01) == "PM1",                  # -7.99 -> Moderate
  mds_tier(modplus_thr)        == "PM1_moderate_plus",    # -8  -> Moderate-plus
  mds_tier(strong_thr + 0.01)  == "PM1_moderate_plus",    # -11.99 -> Moderate-plus
  mds_tier(strong_thr)         == "PM1_strong",           # -12 -> Strong
  mds_tier(strong_thr - 5)     == "PM1_strong"
)

# ── The lookup on the shipped table reproduces DOLPHIN, and real deltas tier ──
if (file.exists("data/pfam_pssm_human.rds")) {
  source("analyses/lib/pssm_lookup.R")
  t <- pssm_table_load("data/pfam_pssm_human.rds")
  stopifnot(length(t$pssm) > 6000)
  ref <- list(c("SYUA_HUMAN",14,"G","R",-9.161), c("SYUA_HUMAN",51,"G","D",-5.415),
              c("SYUA_HUMAN",53,"A","T",1.165),  c("CASR_HUMAN",143,"G","E",-5.959))
  for (r in ref) {
    d <- pssm_delta(t, r[1], as.integer(r[2]), r[3], r[4])$delta
    stopifnot(!is.na(d), abs(d - as.numeric(r[5])) < 0.3)
  }
  # PM1 fires on the disfavoured set, not on tolerated A53T (delta +1.1 > -4)
  stopifnot(
    pssm_fires_pm1(t, "SYUA_HUMAN", 51, "G", "D"),
    pssm_fires_pm1(t, "CASR_HUMAN", 143, "G", "E"),
    !pssm_fires_pm1(t, "SYUA_HUMAN", 53, "A", "T")
  )
  # real reference deltas land in the expected tier
  d_g14r <- pssm_delta(t, "SYUA_HUMAN", 14, "G", "R")$delta   # ~ -9.16 -> Moderate-plus
  d_g51d <- pssm_delta(t, "SYUA_HUMAN", 51, "G", "D")$delta   # ~ -5.42 -> Moderate
  d_g143e <- pssm_delta(t, "CASR_HUMAN", 143, "G", "E")$delta # ~ -5.96 -> Moderate
  stopifnot(
    mds_tier(d_g14r)  == "PM1_moderate_plus",
    mds_tier(d_g51d)  == "PM1",
    mds_tier(d_g143e) == "PM1"
  )
  cat("mds pm1: wiring + tiers + scores all pass (", length(t$pssm), "families)\n")
} else {
  cat("mds pm1: wiring + tier-boundary checks pass (table not shipped, numeric checks skipped)\n")
}
