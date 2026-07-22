# Self-check for DOLPHIN availability reporting + Pass-Blind PM5 stripping.
# Run: Rscript test_dolphin_availability.R
# ponytail: no network needed — fetch_dolphin is stubbed, so this stays green
# whether or not the DOLPHIN service is up.
source("analyses/lib/dolphin.R")
source("analyses/lib/clinvar_blind.R")

# ── dolphin_pm1_call must separate "unreachable" from "said no" ────────────
# Transport failure (timeout/DNS/TLS) surfaces as fetch_dolphin -> NULL. That
# must return NA, not FALSE: collapsing the two is what let an expired TLS
# certificate masquerade as "no PM1" with no signal to the user.
local({
  fetch_dolphin <<- function(...) NULL
  stopifnot(is.na(dolphin_pm1_call(gene = "CASR", p_notation = "p.Asp433Tyr")))
})
# A valid response carrying no PM1 is a real negative -> FALSE, not NA.
local({
  fetch_dolphin <<- function(...) list(results = "mutation not within domain")
  r <- dolphin_pm1_call(gene = "CASR", p_notation = "p.Asp433Tyr")
  stopifnot(identical(r, FALSE))
})
# A valid response carrying PM1 -> TRUE.
local({
  fetch_dolphin <<- function(...) list(results = list(list(acmg = "PM1;PM2")))
  r <- dolphin_pm1_call(gene = "CASR", p_notation = "p.Gly143Glu")
  stopifnot(isTRUE(r))
})

# ── server.R must record the outage rather than leave the pathway blank ────
src <- readLines("server.R", warn = FALSE)
stopifnot(
  any(grepl('is\\.na\\(dolphin_hit\\)', src)),
  any(grepl('pm1_pathway_val <- "dolphin_unavailable"', src)),
  # The error handler must yield NA too, or a thrown error still reads as "no".
  any(grepl('DOLPHIN call failed', src))
)

# ── Pass-Blind must strip every ClinVar-derived PM5 tier ──────────────────
# PM5_supporting is ClinVar-derived exactly as PM5 is; omitting it would leak
# circular evidence into the arm that exists to be non-circular.
stopifnot(
  all(c("PM5", "PM5_supporting") %in% CLINVAR_DIRECT_TAGS),
  !("PM5_supporting" %in% strip_clinvar_tags(c("PM2", "PP2", "PM5_supporting"))),
  !("PM5" %in% strip_clinvar_tags(c("PM2", "PP2", "PM5"))),
  # every PS1 tier is already stripped — PM5 must match that treatment
  length(strip_clinvar_tags(c("PS1_supporting", "PM5_supporting", "PM2"))) == 1L,
  # non-ClinVar evidence survives untouched
  identical(sort(strip_clinvar_tags(c("PM2", "PP2", "PP3_strong"))),
            sort(c("PM2", "PP2", "PP3_strong")))
)
cat("dolphin availability + PM5 blind-stripping: all checks pass\n")
