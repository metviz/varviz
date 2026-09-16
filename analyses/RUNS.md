# Canonical classification runs

Thirteen run directories support every number in the manuscript. Everything else
produced during development was moved to `runs_archive_20260916/` rather than
deleted, so a figure can still be traced if a question arises, without those runs
being mistaken for current ones.

Each run is a full dual-pass classification (Pass-Full and Pass-Blind) over the
same variant universe, reading one pinned snapshot of every external source.

## The 1.2.2 release (current)

| directory | configuration | supports |
| --- | --- | --- |
| `ps_final_v231` | shipped defaults | §3.1, §3.5, Tables S5 and S6 |
| `ps_nomds_v231` | MDS PM1 pathway disabled | §3.1, §3.4 |
| `ps_pm2mod_v231` | PM2 at Moderate (2 pt) | §3.4 |
| `ps_pp3proxy_v231` | retired PP3→Strong proxy restored | §3.4 |
| `ps_nometa_v231` | meta-predictor consensus disabled | §3.4 |
| `ps_pp33pt_v231` | three-point PP3 rung enabled | §3.4 |
| `ras_vcep_pm2_1_v231` | RASopathy cohort, PM2 = 1 pt | §3.4, §3.6 |
| `ras_vcep_pm2_2_v231` | RASopathy cohort, PM2 = 2 pt | §3.4 |
| `external163_v231` | external 163 cohort | §3.6 |

## Earlier releases, kept for the comparisons the manuscript makes

| directory | configuration | supports |
| --- | --- | --- |
| `ps_final_v221` | 1.2.0, before AlphaMissense corroboration | §3.1 and §3.5 before/after |
| `ps_final_v230` | 1.2.1, corroboration without the Pathogenic cap | §3.5 dilution evidence |
| `ps_amcalib_v230` | calibrated AlphaMissense intervals, ungated | §3.5 |
| `ps_amoff_v230` | AlphaMissense at developer threshold only | §3.5 |

## Reproducing

```bash
python3 analyses/build_v230_numbers.py          # every manuscript number, with provenance
python3 analyses/validate_runs.py --baseline analyses/ps_final_v231 \
        --ignore PM1,PM2,PP3,BP4 analyses/ps_{nomds,pm2mod,pp3proxy,nometa,pp33pt}_v231
REVAL_DIR=analyses/ps_final_v231/classifications Rscript analyses/panel_c_caps_reval.R
AM_TSV=analyses/humu/build_v231/numbers/am_scores_14genes.tsv \
        python3 analyses/am_local_calibration.py
```

`--ignore` must include PP3 for any ablation that changes points: the
Pathogenic-boundary constraint is decided on a variant's total, so a weight
change can flip it and with it a PP3 tag (§2.5).

## External sources are pinned, not fetched

`analyses/cache/clinvar/`, `analyses/cache/ensembl/` and `analyses/cache/uniprot/`
hold primed copies. Every run above read from them; no run fetched ClinVar live.
This is what makes two runs comparable, and it exists because NCBI eutils
truncates under concurrent load rather than erroring: six concurrent runs once
returned five different ClinVar record counts for *PTEN* and none at all for
*LDLR*, thinning PM5 and PS1 and collapsing the ClinVar density behind the PM1
hotspot pathway, with row counts intact and nothing reported.

Re-prime with `analyses/prime_clinvar.R` and `analyses/prime_exons.R`, both of
which take a gene list through `VARVIZ_PRIME_GENES`.
