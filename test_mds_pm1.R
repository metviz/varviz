# Self-check for the MDS (Missense Disfavour Score) PM1 Path 4 wiring in server.R.
# Run: Rscript test_mds_pm1.R
#
# Two parts: the lookup helper still reproduces DOLPHIN on the shipped table, and
# server.R is wired to use it as Path 4 (source, load, entry-name bridge, threshold
# gate, export column). Static greps guard the wiring; the numeric check guards
# the score. Slices server.R as a trusted in-repo build artifact, never user input.
src <- readLines("server.R", warn = FALSE)

# ── server.R actually wires MDS as Path 4 ──────────────────────────────────
stopifnot(
  any(grepl('source\\("analyses/lib/pssm_lookup.R"', src)),
  any(grepl('MDS_TABLE <- tryCatch\\(pssm_table_load\\("data/pfam_pssm_human.rds"', src)),
  any(grepl('MDS_PM1_THRESHOLD <- -4', src)),
  # entry-name bridge from the UniProt JSON already fetched
  any(grepl('mds_entry_name <-', src)) && any(grepl('pfam_data\\$uniProtkbId', src)),
  # the Path 4 gate, firing branch, and unavailable sentinel
  any(grepl('pm1_pathway_val <- "mds"', src)),
  any(grepl('pm1_pathway_val <- "mds_unavailable"', src)),
  any(grepl('mds_val <= MDS_PM1_THRESHOLD', src)),
  # score is exported
  any(grepl('MDS_Score = ', src)),
  # the retired DOLPHIN call is no longer the runtime Path 4
  !any(grepl('pm1_pathway_val <- "dolphin"', src))
)

# ── The lookup on the shipped table reproduces DOLPHIN ─────────────────────
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
  cat("mds pm1: wiring + scores all pass (", length(t$pssm), "families)\n")
} else {
  cat("mds pm1: wiring checks pass (table not shipped, numeric checks skipped)\n")
}
