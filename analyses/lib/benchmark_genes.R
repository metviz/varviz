# BENCHMARK_GENES — 15-gene set for Pass 2 F'' validation universe.
#
# Per locked Methodological Foundation MF2 (docs/plans/2026-04-26-sk-revisions-pass2.md):
# the benchmark is decoupled from VarViz's user-facing 7-gene panel. The
# expansion (8 genes) is benchmark-only and not exposed in the UI.
#
# Sourced by:
#   analyses/03_build_universe.R           — defensive filter at union step
#   analyses/05_classify_harness.R         — gene loop in classification
#   analyses/06_panel_a_clinical.R         — VariBench panel filter
#   analyses/07_panel_b_dms.R              — MaveDB panel filter
#
# Coverage observed in canonical TSVs (VariBench D2 + MaveDB, 2026-04-27):
#   - With ANY data:   13 / 15 (PPP2R5C and NUDT15 have neither VariBench nor sufficient MaveDB)
#   - With both VariBench AND MaveDB: 7 (TP53, BRCA1, LDLR, PTEN, KRAS, KCNQ1, KCNH2, GCK)
#   - With only VariBench: 5 (CASR, BAP1, DDR2, TSHR, SLC13A5)
#   - With only MaveDB:    1 (NUDT15)

BENCHMARK_GENES <- c(
  # 7 user-facing VarViz panel genes
  "TP53", "PPP2R5C", "BAP1", "DDR2", "TSHR", "SLC13A5", "CASR",
  # O4 expansion — benchmark-only; not exposed in the UI.
  "BRCA1",   # gold-standard ACMG (Findlay 2018 saturating DMS; ENIGMA PS3 precedent)
  "LDLR",    # Matreyek 2021 DMS; deepest VariBench coverage (~307 variants)
  "PTEN",    # Mighell 2018 + Matreyek 2018 DMS
  "KRAS",    # 22 MaveDB scoresets (drug-response screens)
  "KCNQ1",   # LQTS — VariBench + KCNQ1 DMS
  "KCNH2",   # LQTS — saturating DMS in MaveDB
  "GCK",     # MODY/diabetes — DMS available
  "NUDT15"   # pharmacogenomics — Suiter 2020 DMS (no VariBench)
)
