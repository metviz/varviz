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
| `ps_am792_v124_refnum` | shipped defaults | §3.1, §3.5, Tables 2, S5, S6 |
| `ps_nomds_v124` | MDS PM1 pathway disabled | §3.1, §3.4 |
| `ps_pm2mod_v124` | PM2 at Moderate (2 pt) | §3.4 |
| `ps_pp3proxy_v124` | retired PP3→Strong proxy restored | §3.4 |
| `ps_nometa_v124` | meta-predictor consensus disabled | §3.4 |
| `ps_pp33pt_v124` | three-point PP3 rung enabled | §3.4 |
| `ps_amcal_nohybrid_v124` | full calibrated ladder, hybrid off | §3.4, §3.5 |
| `ps_am564_v124` | supporting rung at the pre-1.2.3 0.564 | §3.5 fallback sensitivity |
| `ps_rasopathy_v124` | RASopathy cohort, PM2 = 1 pt | §3.6, Table 4 |
| `ps_rasopathy_pm2mod_v124` | RASopathy cohort, PM2 = 2 pt | §3.4 |
| `ps_ras_am564_v124` | RASopathy at the 0.564 rung | §3.5 fallback sensitivity |
| `ps_ext163_moi_v124` | external 65 cohort, per-gene inheritance | §3.6, §S4 |
| `ps_ext163_moi_pm2mod_v124` | external 65, PM2 = 2 pt | §S4.3 |
| `ps_ext163_w10amono_v124` | *WNT10A* held monoallelic, PM2 = 1 pt | §S4.2 headline |
| `ps_ext163_w10amono_pm2mod_v124` | *WNT10A* monoallelic, PM2 = 2 pt | §S4.3 |

All runs above are release 1.2.3. Each carries a per-gene `__manifest.json` recording the
engine flags, the universe path and its sha256, so a directory can be checked against the
configuration it claims rather than against its name.

**The external cohort needs the inheritance-aware universe.** `variant_universe_external163.tsv`
carries no `inh_param` column, so a run against it silently treats all 44 genes as
monoallelic. The reported results use `variant_universe_external163_moi.tsv`, which assigns
26 genes biallelic and 18 monoallelic, with `variant_universe_external163_wnt10a_mono.tsv`
overriding *WNT10A* for the conservative headline.

## Reproducibility check, 2026-10-11

Two canonical runs were repeated against the same universes to test whether the
pinned-snapshot design holds. Kept rather than deleted so the check is auditable;
both directories are gitignored, so this table is the record.

| directory | repeats | result |
| --- | --- | --- |
| `ps_drift_20261011_ras` | `ps_rasopathy_v124` | 0 of 226 variants differ |
| `ps_drift_20261011_ext` | `ps_ext163_moi_v124` | 0 of 65 variants differ |

Identical on classification, points and the full tag string in both passes.
Headline counts unchanged: RASopathy Pass-Full 147 P/LP and Pass-Blind 123;
external cohort 54 of 65 and 47 of 65. Both universes still match the
`universe_sha256` recorded in the originals' manifests, so the inputs are
provably the same ones.

**This doubles as proof that 1.2.4 moves no calls.** The originals were produced
at 1.2.3 and these at 1.2.4, across 291 variants in 60 genes, with no difference
anywhere. The release notes claim a PATCH on the grounds that nothing reaches
the engine; this is the measurement behind that claim.

Why it was run: the *CASR* case-study table drifted between its September session
export and October, by one variant. That table comes from a live application
session. These cohorts come from the harness, which reads one pinned snapshot of
every external source, so they were expected not to move. They did not. Drift
exposure in the manuscript is therefore confined to the two case-study tables,
which now carry their export date.

**The degraded-cache guard blocks both of these runs and should not be overridden
blindly.** It refused to write `summary.tsv`, reporting REVEL missing for four
RASopathy genes and REVEL or AlphaMissense missing for most external ones. What
it detects is an absent *bulk per-chromosome preload*, not absent values: REVEL
and AlphaMissense were fully populated in every flagged gene and identical to the
original. Only eleven REVEL slices exist under `analyses/raw/revel/by_chrom/`
(chr 1, 3, 4, 7, 10, 11, 12, 13, 14, 17, 19), so any cohort reaching other
chromosomes trips it. `FAM20A` carries no REVEL in either run, which is a real
gap and is reproduced faithfully. The diffs above were taken from the per-gene
checkpoints, which are written regardless, so `VARVIZ_ALLOW_DEGRADED=1` was never
needed. Setting that flag cannot distinguish "preload absent, values complete"
from "values genuinely missing", which is the reason to tighten the guard rather
than to habitually override it.

## Superseded runs, kept for provenance

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
