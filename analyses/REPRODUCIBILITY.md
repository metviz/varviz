# VarViz MDS / benchmark reproducibility

Every manuscript number for the Missense Disfavor Score (MDS) and the
triangulated non-circular benchmark, mapped to the script that produces it.
Paths are repo-relative. Committed scripts live under `analyses/`;
`data/pfam_pssm_human.rds` is bundled.

## Raw data sources (external, one-time)

| Source | Local landing | Used for |
|---|---|---|
| Pfam-A 33.1 full alignments | → `data/pfam_pssm_human.rds` | MDS PSSM (6,534 families, 11.4M mapped residues) |
| ClinVar `variant_summary.txt.gz` ([NCBI FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz)) | `analyses/tmp/clinvar/variant_summary.txt.gz` | genome-wide LR calibration |
| VariBench (PON-PS_D2) | `analyses/derived/varibench_canonical.tsv` | clinical AUROC + MDS calibration |
| MaveDB DMS (21 scoresets) | `analyses/derived/mavedb_canonical.tsv` | functional concordance panel |
| gnomAD v4.1 / dbNSFP 4.9a (local) | `analyses/tmp/…` | predictors + CAPS |
| MANE / UniProt maps | `analyses/tmp/mane/uniprot_acc_entry.tsv`, `analyses/derived/gene_mane.tsv` | gene → family residue mapping |

## Pipeline

| # | Script | Inputs | Output | Produces |
|---|---|---|---|---|
| 1 | `lib/pfam_pssm.R` + `20_fix_pssm_dims.R` | Pfam-A alignments | `data/pfam_pssm_human.rds` | MDS score table (int8 quantized) |
| 2 | `03_build_universe.R` | VariBench + gene set | `derived/variant_universe.tsv` | 90,701-variant / 14-gene universe |
| 3 | `01_pull_varibench.R`, `02_pull_mavedb.R` | external | `derived/varibench_canonical.tsv`, `derived/mavedb_canonical.tsv` | clinical + functional truth sets |
| 4 | `25_clinvar_extract.R` | `tmp/clinvar/variant_summary.txt.gz` | `tmp/clinvar/clinvar_missense_2star.tsv` | 2★ ClinVar missense P/LP + B/LB |
| 5a | `ps_final_harness.R` with `VARVIZ_MDS_PM1=FALSE`, `VARVIZ_OUT_DIR=analyses/ps_nomds_v2` | `server.R`, `lib/clinvar_blind.R`, `lib/local_predictors.R`, `lib/harness_guard.R` | `ps_nomds_v2/summary.tsv` | dual-pass, **MDS off** (ablation) |
| 5b | `ps_final_harness.R` (defaults) | same | `ps_final/summary.tsv` | dual-pass, **Moderate-only MDS** — the submission run |
| 5c | `ps_final_harness.R` with `VARVIZ_MDS_TIERED=TRUE`, `VARVIZ_OUT_DIR=analyses/ps_tiered` (or `ps_tiered_harness.R`) | same | `ps_tiered/summary.tsv` | dual-pass, **LR-tiered MDS** — exploratory calibration only |
| 6 | `21_assess_mds.R` | `varibench_canonical.tsv`, `pfam_pssm_human.rds` | stdout | VariBench AUROC 0.785; DOLPHIN ρ=1.000 |
| 7 | `22_mds_benchmark.R` | `variant_universe.tsv`, `varviz_classifications_dolphin.tsv`, pssm | `derived/varviz_classifications_mds.tsv` | per-variant MDS/DOLPHIN deltas |
| 8 | `23_mds_vs_dolphin_mcnemar.py` | `varviz_classifications_mds.tsv`, `ps_nomds/summary.tsv` | stdout | 3,584 vs 375; McNemar χ²=3,121 |
| 9 | `24_clinvar_mds_lr.R` | `clinvar_missense_2star.tsv`, pssm, MANE maps | stdout + `tmp/clinvar/clinvar_mds.rds` | LR+ 5.4/11.2/39.0 (n=18,570/14,646) |

Metric computation (`AUROC`, `Sensitivity`, `Specificity`, `MCC`, `VUS_rate`)
is `lib/metric_suite.R::compute_metric_suite`, applied to each `summary.tsv`
joined to the VariBench-labelled subset of `variant_universe.tsv`.

## Run configurations (single source of truth)

One engine (`server.R` v2.0.0, points-only `classify_acmg`), one harness body,
several option sets. A run is defined by its options, not by which wrapper file
launched it. Regenerating any run requires `--force` (shared harness) or
`VARVIZ_FORCE=1` (`ps_final_harness.R`); the harness refuses to overwrite a
`summary.tsv` otherwise.

| Run directory | `varviz.mds_pm1` | `varviz.mds_tiered` | Role | Status |
|---|---|---|---|---|
| `ps_final` | TRUE | FALSE | **Submission run.** Moderate-only MDS; corroboration takes the stronger of pathway PM1 and MDS | authoritative once regenerated on v2.0.0 |
| `ps_nomds_v2` | FALSE | FALSE | MDS-off ablation (same harness, same sentinels) | authoritative once regenerated on v2.0.0 |
| `ps_tiered` | TRUE | TRUE | Exploratory: +3 / +4 LR tiers at MDS ≤ −8 / −12 | supplementary calibration only; not a manuscript claim |
| `ps_mds_corroborate`, `ps_mds_consfree` | TRUE | TRUE | pre-review builds; `mds_frees_cons` is a no-op since 2.0.0 so the two are identical | **superseded, do not use** (also contaminated by silent fetch failures, see `repro/README.md`) |
| `ps_nomds`, `ps_reval`, `ps_baseline` | – | – | Jul 30 runs, no `run.log`, pre-sentinel engine | **superseded, do not use** |

Every `summary.tsv` on disk that predates v2.0.0 was produced by the rule-ladder
engine and is not comparable to the current one. Regenerate before quoting.

## Evidence-overlap constraint

Each biological signal is scored under exactly one criterion (also stated at
the MDS constants in `server.R`):

- **MDS** (Pfam PSSM, substitution-specific) is counted once, under PM1. It
  never feeds PP3. Where a pathway PM1 (site / CCRS / domain / ClinVar hotspot)
  and MDS both fire, the stronger is taken, never the sum: Moderate + Moderate
  stays Moderate. No joint calibration of "PM1 and MDS together" exists, so no
  points are awarded for the conjunction.
- **Position conservation** (PhyloP / PhastCons / GERP / ConSurf) spent to
  reach PM1_strong is withheld from the PP3 conservation tier
  (`cons_used_for_pm1`).
- **Sequence predictors** (REVEL / CADD / AlphaMissense / meta-predictors) are
  one line of evidence, PP3, at the Pejaver 2022 calibrated strength.
  Agreement among predictors does not upgrade PP3 and does not stand in for
  PS3.
- **Locus evidence** PP1 + PP4 is capped jointly at 5 points (ClinGen SVI).

## Exploratory: MDS tier LR+ on ClinVar 1-star (2026-09-06, not a manuscript claim)

The 2-star Strong tier rests on 8 benign variants at MDS ≤ −12 (LR+ 39.0, 95% CI
18.0–84.9, lower bound below the 18.7 Strong cut-point). To test the tiers out of
sample, `25_clinvar_extract.R --stars=1` writes a separate ≥1-star set and
`24_clinvar_mds_lr.R --in=<that> --review=single` scores the 1-star-ONLY rows,
which share no variant with the 2-star calibration set.

| MDS ≤ | 2-star (calibration) | 1-star-only (held-out, n=25,241 P / 36,379 B) | pooled ≥1-star |
|---|---|---|---|
| −4 | 5.4 [5.1–5.8] | 8.0 [7.6–8.4] | 7.0 [6.7–7.3] |
| −8 | 11.2 [9.5–13.2] | 14.5 [12.9–16.3] | 13.3 [12.1–14.6] |
| −10 | 18.6 [13.5–25.7] | 26.3 [20.7–33.4] | 23.5 [19.4–28.4] |
| −12 | 39.0 [18.0–84.9] | 32.1 [20.6–50.2] | 33.0 [22.5–48.6] |

Strong replicates out of sample (28 benign in the ≤ −12 tail, lower bound 20.6 >
18.7), and ≤ −10 also clears Strong. ≤ −8 stays between Moderate and Strong on
every set. Both sets are ClinVar, so this is not an orthogonal truth source and
does not calibrate PM1 pathway + MDS jointly; the submission configuration
remains Moderate-only. Outputs: `analyses/tmp/clinvar/clinvar_missense_1star.tsv`,
`clinvar_mds_1star_only.rds`, `clinvar_mds_1star_pooled.rds`,
`lr_1star_exploratory.log` (gitignored, regenerable in ~5 min).

## Manuscript number → source

| Number | From |
|---|---|
| LR+ 5.4 / 11.2 / 39.0 at MDS ≤ −4 / −8 / −12 (n=18,570 P / 14,646 B) | `24_clinvar_mds_lr.R` |
| Pass-Blind AUROC 0.938 → 0.959 → 0.960; spec 0.567, sens 1.000, MCC 0.747 | `ps_{nomds_v2,final,tiered}/summary.tsv` × `lib/metric_suite.R` on VariBench (n=1,108). **Pre-2.0.0 values; regenerate.** |
| MDS resolves 3,584 vs DOLPHIN 375; McNemar χ²=3,121, p<10⁻³⁰⁰ | `23_mds_vs_dolphin_mcnemar.py` |
| LR-tiering adds +130 → 3,714 total | `ps_tiered` vs `ps_nomds_v2`. **Exploratory only since 2.0.0; not a manuscript claim.** |
| DOLPHIN concordance ρ=1.000 (CASR G143E −5.95 vs −5.96); VariBench AUROC 0.785 | `21_assess_mds.R` |
| MDS artifact: 6,534 families, 11,355,753 residues, ~50 MB | `data/pfam_pssm_human.rds` (from `lib/pfam_pssm.R`) |

## Notes

- **Engine flags.** Since v2.0.0 `varviz.mds_tiered` defaults to FALSE
  (Moderate-only MDS). `ps_final_harness.R` reads `VARVIZ_MDS_PM1` and
  `VARVIZ_MDS_TIERED` (TRUE/FALSE only; anything else errors). The DOLPHIN
  comparison uses `analyses/raw/dolphin/by_gene/*.tsv`, where `pm1 = NA` means
  the API call failed and is retried on the next `fetch_dolphin_gene()` run.
- **`clinvar_missense_2star.tsv` fidelity.** `25_clinvar_extract.R` regenerates
  the 2★ set to ~99.8% of the committed file (dedup / p.-parse edge cases);
  immaterial to the LR tiers.
- **Non-circularity.** `ps_*_harness.R` produce a Pass-Full and a Pass-Blind
  column per variant; Pass-Blind strips ClinVar-derived tags via
  `lib/clinvar_blind.R::strip_clinvar_tags()`. MDS, VariBench, and CAPS are each
  orthogonal to the stripped ClinVar evidence.
- **Rejected variants** (not in the manuscript): MDS-first ordering; MDS
  strong-upgrade of site/CCRS/domain (DMS-flat, overstates). The DMS enrichment
  at the strict tail is non-monotone, so the Strong tier is anchored on the
  genome-wide ClinVar LR (spec 0.999 at ≤ −12), not on DMS.
- **Deployment / data storage** (engineering only; does not affect any result).
  The CCRS track (~10.7 M rows) and the Pfam PSSM (~11.4 M residues) are the two
  large reference tables. To keep the cloud instance within memory, the deploy
  build (`26_build_varviz_sqlite.R`) serves CCRS from an indexed SQLite store
  (`data/ccrs.sqlite`, queried one gene at a time via `server.R::ccrs_slice()`)
  and ships the PSSM as a pre-baked raw/int8 + factor `.rds` that loads in ~1 s.
  Storage format is invisible to the science: the same tables computed the same
  way (Pipeline above), only lazily loaded. `data/*.sqlite` and the pre-baked rds
  are git-ignored, regeneratable from the canonical `data/pfam_pssm_human.rds`
  and `VarViz.RData` by that script.
