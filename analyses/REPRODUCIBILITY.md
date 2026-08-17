# VarViz MDS / benchmark reproducibility

Every manuscript number for the Missense Disfavour Score (MDS) and the
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
| 5a | `ps_nomds_harness.R` | `server.R`, `lib/clinvar_blind.R`, `lib/local_predictors.R` | `ps_nomds/summary.tsv` | dual-pass, **MDS off** (baseline) |
| 5b | `ps_reval_harness.R` | same | `ps_reval/summary.tsv` | dual-pass, **flat MDS** |
| 5c | `ps_tiered_harness.R` (`options(varviz.mds_tiered=TRUE)`) | same | `ps_tiered/summary.tsv` | dual-pass, **LR-tiered MDS** |
| 6 | `21_assess_mds.R` | `varibench_canonical.tsv`, `pfam_pssm_human.rds` | stdout | VariBench AUROC 0.785; DOLPHIN ρ=1.000 |
| 7 | `22_mds_benchmark.R` | `variant_universe.tsv`, `varviz_classifications_dolphin.tsv`, pssm | `derived/varviz_classifications_mds.tsv` | per-variant MDS/DOLPHIN deltas |
| 8 | `23_mds_vs_dolphin_mcnemar.py` | `varviz_classifications_mds.tsv`, `ps_nomds/summary.tsv` | stdout | 3,584 vs 375; McNemar χ²=3,121 |
| 9 | `24_clinvar_mds_lr.R` | `clinvar_missense_2star.tsv`, pssm, MANE maps | stdout + `tmp/clinvar/clinvar_mds.rds` | LR+ 5.4/11.2/39.0 (n=18,570/14,646) |

Metric computation (`AUROC`, `Sensitivity`, `Specificity`, `MCC`, `VUS_rate`)
is `lib/metric_suite.R::compute_metric_suite`, applied to each `summary.tsv`
joined to the VariBench-labelled subset of `variant_universe.tsv`.

## Manuscript number → source

| Number | From |
|---|---|
| LR+ 5.4 / 11.2 / 39.0 at MDS ≤ −4 / −8 / −12 (n=18,570 P / 14,646 B) | `24_clinvar_mds_lr.R` |
| Pass-Blind AUROC 0.938 → 0.959 → 0.960; spec 0.567, sens 1.000, MCC 0.747 | `ps_{nomds,reval,tiered}/summary.tsv` × `lib/metric_suite.R` on VariBench (n=1,108) |
| MDS resolves 3,584 vs DOLPHIN 375; McNemar χ²=3,121, p<10⁻³⁰⁰ | `23_mds_vs_dolphin_mcnemar.py` |
| LR-tiering adds +130 → 3,714 total | `ps_tiered` vs `ps_nomds` (upward-only PM1 augmentation delta) |
| DOLPHIN concordance ρ=1.000 (CASR G143E −5.95 vs −5.96); VariBench AUROC 0.785 | `21_assess_mds.R` |
| MDS artifact: 6,534 families, 11,355,753 residues, ~50 MB | `data/pfam_pssm_human.rds` (from `lib/pfam_pssm.R`) |

## Notes

- **Engine flag.** LR-tiered MDS is default on (`server.R`, commit `e0227c9`);
  disable with `options(varviz.mds_tiered = FALSE)`. `ps_nomds` and `ps_reval`
  set the flag off / flat for the A/B baselines.
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
