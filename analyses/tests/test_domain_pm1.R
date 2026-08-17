# Self-check for the conservation-based domain PM1 scorer.
# Run: Rscript test_domain_pm1.R
# ponytail: pure functions only, no dbNSFP corpus needed, so this stays
# runnable on any checkout without the 110 GB local dataset.
source("analyses/lib/domain_pm1.R")

# ── dbNSFP field parsing ──────────────────────────────────────────────────
# ".;.;.;." is the trap: it is not "." but carries no annotation. Treating it
# as present overcounted domain-resident rows roughly two-fold (84k vs 41k on
# chr21), which would have inflated every coverage figure downstream.
stopifnot(
  is.na(dbnsfp_first_present(".;.;.;.")),
  is.na(dbnsfp_first_present(".")),
  is.na(dbnsfp_first_present("")),
  dbnsfp_first_present(".;Cadherin;.")        == "Cadherin",
  dbnsfp_first_present("GPCR family 3;.;.")   == "GPCR family 3",
  dbnsfp_first_present("  ;  Kinase ;.")      == "Kinase",
  identical(dbnsfp_numeric("5.12;."), 5.12),
  is.na(dbnsfp_numeric(".;.")),
  is.na(dbnsfp_numeric("notanumber"))
)

# ── within-domain percentile ──────────────────────────────────────────────
cons <- c(1, 2, 3, 4, 5)
stopifnot(
  # too few positions to rank honestly -> NA, not a fabricated 100th centile
  all(is.na(domain_position_percentile(cons, min_positions = 20L))),
  # with the gate satisfied, ranks span 0..100 and are monotonic in conservation
  {
    p <- domain_position_percentile(1:25, min_positions = 20L)
    isTRUE(all.equal(min(p), 0)) && isTRUE(all.equal(max(p), 100)) &&
      !is.unsorted(p)
  },
  # Ties share a rank rather than one arbitrarily taking the top slot. The tied
  # group must stay a minority or the discrimination guard (below) voids the
  # domain — which is why this uses 3 tied of 25, not 24 of 25.
  {
    p <- domain_position_percentile(c(rep(30, 3), 1:22), min_positions = 20L)
    length(unique(p[1:3])) == 1L && min(p[1:3]) > max(p[4:25])
  },
  # NAs stay NA and do not shift other positions' ranks
  {
    p <- domain_position_percentile(c(1:24, NA), min_positions = 20L)
    is.na(p[25]) && isTRUE(all.equal(max(p, na.rm = TRUE), 100))
  }
)

# ── PM1 call ──────────────────────────────────────────────────────────────
stopifnot(
  identical(domain_pm1_call(c(95, 89.9, 90), threshold = 90), c(TRUE, FALSE, TRUE)),
  # unscorable stays NA — "not enough data" must not read as "not constrained",
  # the same distinction the DOLPHIN outage showed us matters
  is.na(domain_pm1_call(NA_real_)),
  length(domain_pm1_call(numeric(0))) == 0L
)

# ── end-to-end on a synthetic two-domain table ────────────────────────────
# Domain A: 30 positions, conservation ramps 1..30 -> top 10% must fire.
# Domain B: 5 positions -> below the gate, nothing may fire.
df <- rbind(
  data.frame(gene = "G1", domain = "A", aapos = 1:30, cons = as.numeric(1:30)),
  data.frame(gene = "G1", domain = "B", aapos = 31:35, cons = as.numeric(1:5))
)
s <- score_domain_table(df, min_positions = 20L, threshold = 90)
a <- s[s$domain == "A", ]
b <- s[s$domain == "B", ]
stopifnot(
  all(is.na(b$pctile)), all(is.na(b$pm1)),          # small domain unscorable
  sum(a$pm1, na.rm = TRUE) == 3L,                    # 90th centile of 30 = top 3
  all(a$pm1[a$aapos >= 28], na.rm = TRUE),           # and they are the top ones
  !any(a$pm1[a$aapos <= 27], na.rm = TRUE),
  # A uniformly conserved domain must be reported as unscorable, not blanketed
  # with PM1. Under ties.method="max" every residue ties at the top, so without
  # the discrimination guard this fires at all 30 positions — the mirror of the
  # GERP-saturation bug that made the first build fire at none.
  {
    flat <- data.frame(gene = "G2", domain = "C", aapos = 1:30, cons = rep(7, 30))
    sf <- score_domain_table(flat, min_positions = 20L, threshold = 90)
    all(is.na(sf$pctile)) && all(is.na(sf$pm1)) && !any(sf$pm1, na.rm = TRUE)
  },
  # Guard boundary: a domain where a large minority ties at the ceiling is
  # still usable (real case: GERP saturates on 23% of the CASR ligand-binding
  # domain), but one where the majority ties is not.
  {
    ok  <- domain_position_percentile(c(rep(9, 6), 1:24), min_positions = 20L)
    bad <- domain_position_percentile(c(rep(9, 20), 1:10), min_positions = 20L)
    !all(is.na(ok)) && all(is.na(bad))
  }
)

# ── collapse: one row per residue, independent of alternate-allele count ──
# Three alt alleles at residue 10 must not out-weigh a residue with one.
dup <- data.frame(gene = "G1", domain = "A", aapos = c(10, 10, 10, 11),
                  cons = c(4, 4, 4, 9))
cl <- collapse_domain_positions(dup)
stopifnot(nrow(cl) == 2L, cl$cons[cl$aapos == 10] == 4, cl$cons[cl$aapos == 11] == 9)

cat("domain PM1 scorer: all checks pass\n")
