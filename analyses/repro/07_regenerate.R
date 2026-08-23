#!/usr/bin/env Rscript
# Regenerate every manuscript benchmark number from the authoritative run,
# on the de-duplicated (distinct variant) basis.
#
# Authoritative run : analyses/ps_final/summary.tsv
#   server.R Path 4 ungated -- MDS originates PM1 where no other path fired and
#   corroborates (additive tiers, capped at Strong) where one did; the
#   ClinVar-hotspot upgrade records its pre-upgrade tag so strip_clinvar_tags()
#   demotes it under Pass-Blind.
# Ablation comparator: analyses/ps_nomds_v2/summary.tsv (MDS disabled entirely).
#
# Both runs were scored on the patched server.R with all four data-source
# sentinels active and both pass 10_validate_run.py intrinsic checks. The
# previous pair must not be used: ps_mds_corroborate lost KCNQ1 and TP53 to
# silent UniProt failures, and analyses/ps_nomds (Jul 30) has no run.log at all,
# so a MDS-on/off comparison across them would confound the MDS effect with
# differing contamination.
#
# Basis: distinct (gene, p_notation). The universe file is MaveDB + VariBench and
# MaveDB contributes one row per scoreset, so raw rows are variant x scoreset
# records. Duplicates classify identically (verified by 05_dedup_check.R).
#
# Run from the repository root:  Rscript analyses/repro/07_regenerate.R

RUN <- "analyses/ps_final/summary.tsv"
OFF <- "analyses/ps_nomds_v2/summary.tsv"
OUT <- "analyses/repro/out/regenerated.tsv"

rd <- function(p) read.delim(p, sep = "\t", quote = "", colClasses = "character")
ded <- function(d) d[!duplicated(paste(d$gene, d$p_notation)), ]

d   <- ded(rd(RUN))
off <- ded(rd(OFF))
n   <- nrow(d)

res <- list(); put <- function(k, v) res[[k]] <<- v
VUS <- c("VUS-High", "VUS-Mid", "VUS-Low")
pm1 <- function(t) ifelse(grepl("PM1_strong", t), "PM1_strong",
                   ifelse(grepl("PM1_moderate_plus", t), "PM1_moderate_plus",
                   ifelse(grepl("PM1", t), "PM1", "none")))

put("universe_rows",     nrow(rd(RUN)))
put("universe_variants", n)

f <- d$varviz_classification_full; b <- d$varviz_classification_blind
pct <- function(x, s) 100 * sum(x %in% s) / n
put("full_plp",  pct(f, c("Pathogenic","Likely Pathogenic")))
put("blind_plp", pct(b, c("Pathogenic","Likely Pathogenic")))
put("full_vush", pct(f, "VUS-High"));  put("blind_vush", pct(b, "VUS-High"))
put("full_lbb",  pct(f, c("Likely Benign","Benign")))
put("blind_lbb", pct(b, c("Likely Benign","Benign")))
put("full_vus",  pct(f, VUS));  put("blind_vus", pct(b, VUS))

sh <- f != b
put("shift_n", sum(sh)); put("shift_pct", 100 * sum(sh) / n)
tr <- function(x, y) sum(f == x & b == y)
put("lp_to_vush_n",     tr("Likely Pathogenic","VUS-High"))
put("lp_to_vush_of_lp", 100 * tr("Likely Pathogenic","VUS-High") / sum(f == "Likely Pathogenic"))
put("lp_to_vush_of_sh", 100 * tr("Likely Pathogenic","VUS-High") / sum(sh))
put("p_to_lp_n",        tr("Pathogenic","Likely Pathogenic"))
put("p_to_lp_of_p",     100 * tr("Pathogenic","Likely Pathogenic") / sum(f == "Pathogenic"))
put("p_to_lp_of_sh",    100 * tr("Pathogenic","Likely Pathogenic") / sum(sh))

pf <- pm1(d$tags_full); pb <- pm1(d$tags_blind)
for (k in c("PM1_strong","PM1_moderate_plus","PM1")) {
  put(paste0("pm1_full_",  k), sum(pf == k))
  put(paste0("pm1_blind_", k), sum(pb == k))
}
put("pm1_strong_demoted",     sum(pf == "PM1_strong") - sum(pb == "PM1_strong"))
put("pm1_strong_demoted_pct", 100 * (sum(pf == "PM1_strong") - sum(pb == "PM1_strong")) /
                                     sum(pf == "PM1_strong"))

# "no PM1 pathway" is serialised two ways across runs: older runs write an
# empty string, newer ones write a literal NA (which read.delim turns into a
# real NA). Left unnormalised, `pw == "mds"` propagates NA and every exact-match
# count below silently becomes NA rather than a number.
pw <- d$pm1_pathway
pw[is.na(pw) | pw == "NA"] <- ""
put("mds_corroborate", sum(grepl("\\+mds", pw) & !grepl("hotspot", pw)))
# "Originates" in the manuscript's sense (Supp S3.2: "originates PM1 for N
# where no other pathway had fired") = MDS is the LEADING pathway, which
# includes rows a ClinVar hotspot later upgraded -- the upgrade happens after
# MDS already originated the call. That is fires - corroborate. The strict
# pathway == "mds" count excludes the upgraded rows and is reported separately;
# conflating the two is what made 00_claims.tsv disagree with the document.
put("mds_originate",        sum(grepl("^mds($|\\+)", pw)))
put("mds_originate_strict", sum(pw == "mds"))
put("mds_fires",       sum(grepl("(^|\\+)mds($|\\+)", pw)))
put("hotspot_upgrades", sum(grepl("hotspot_upgrade", pw)))
put("hotspot_only",     sum(pw == "clinvar_hotspot"))

# MDS contribution vs the MDS-off ablation, matched on variant
m <- merge(d[, c("gene","p_notation","varviz_classification_blind")],
           off[, c("gene","p_notation","varviz_classification_blind")],
           by = c("gene","p_notation"), suffixes = c("_on","_off"))
put("ablation_n", nrow(m))
put("mds_resolved_vs_off",
    sum(m$varviz_classification_blind_off %in% VUS &
        !(m$varviz_classification_blind_on %in% VUS)))
put("off_blind_plp", 100 * sum(m$varviz_classification_blind_off %in%
                               c("Pathogenic","Likely Pathogenic")) / nrow(m))

# Panel A metrics, emitted here so 02_reconcile.R checks them. These were
# absent from 00_claims.tsv entirely, so nothing verified Table 1 / Table S2 --
# which is how Pass-Full AUROC 0.985 / MCC 0.473 / spec 0.242 / VUS 0.098
# survived from a superseded run into the submitted draft.
pa <- "analyses/repro/out/panel_a_runs.tsv"
if (file.exists(pa)) {
  m <- read.delim(pa, sep = "\t", quote = "", stringsAsFactors = FALSE)
  m <- m[m$basis == "variants", ]
  pick <- function(run, arm, metric) {
    v <- m$estimate[m$run == run & m$arm == arm & m$metric == metric]
    if (length(v)) round(as.numeric(v[1]), 3) else NA_real_
  }
  put("auroc_full",   pick("MDS corroborate", "full",  "AUROC"))
  put("auroc_blind",  pick("MDS corroborate", "blind", "AUROC"))
  put("auroc_mdsoff", pick("MDS off",         "blind", "AUROC"))
  put("mcc_full",     pick("MDS corroborate", "full",  "MCC"))
  put("mcc_blind",    pick("MDS corroborate", "blind", "MCC"))
  put("spec_full",    pick("MDS corroborate", "full",  "Specificity"))
  put("spec_blind",   pick("MDS corroborate", "blind", "Specificity"))
  put("vus_full",     pick("MDS corroborate", "full",  "VUS_rate"))
  put("vus_blind",    pick("MDS corroborate", "blind", "VUS_rate"))
} else {
  message("[repro] panel_a_runs.tsv absent - run 06_panel_a_runs.R first; ",
          "Table 1 / Table S2 metrics will NOT be reconciled")
}

o <- data.frame(id = names(res), value = round(as.numeric(unlist(res)), 4))
write.table(o, OUT, sep = "\t", row.names = FALSE, quote = FALSE)
for (i in seq_len(nrow(o)))
  cat(sprintf("  %-24s %12s\n", o$id[i], format(o$value[i])))
cat(sprintf("\n[repro] wrote %s\n", OUT))
