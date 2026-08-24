# External validation cohort — Supplementary §S4.8

Everything needed to verify the external validation reported in
Supplementary §S4.8 and Table S6 of the VarViz manuscript, without downloading
the large Zenodo deposit. 1.4 MB total.

## Source data

The variant list originates with:

> Ghasemnejad T, Liang Y, Jahanian KH, et al. Comprehensive evaluation of
> ACMG/AMP-based variant classification tools. *Bioinformatics*.
> 2026;42(2):btaf623. doi:10.1093/bioinformatics/btaf623

Their manual-curation file, the four benchmarked tool outputs and the LIRICAL
results are released CC0 in a Code Ocean capsule at
**doi:10.24433/CO.6562438.v1**. Those files are *not* mirrored here — fetch
them from the capsule. `variants_parsed.tsv` records the derivation from
`manual_163new.xlsx`, so the chain is auditable without redistributing another
group's data.

## What VarViz did with it

152 records → 143 alleles → 111 coding substitutions → **65 distinct missense
variants across 44 genes and 83 probands**. cDNA notation was resolved against
MANE Select with the Ensembl VEP REST API; the source file carries no
transcript identifier and mixes standard (`c.289C>T`) with legacy reversed
(`c.G379T`) notation.

No gene overlaps the 14-gene development universe, and no variant contributed
to threshold calibration.

## Files

| file | contents |
|---|---|
| `RESULTS.md` | full write-up: method, results, comparisons, limitations |
| `TABLE_external_comparison.tsv` | Supplementary Table S6 source — 65 variants × VarViz full/blind/points/PM1 pathway × InterVar re-run × 4 published tool columns |
| `variant_universe_external163.tsv` | harness input universe |
| `varviz_dualpass_summary.tsv` | VarViz Pass-Full / Pass-Blind output, 65 rows |
| `classifications/*__dual.tsv` | per-gene checkpoints, 44 files |
| `harness_run.log` | full run log — the four data-source sentinels are the audit trail that no fetch degraded silently |
| `variants_parsed.tsv` | all 143 alleles, parsed and categorised from the source file |
| `variants_annotated.tsv` | + VEP consequence, protein position, p. notation |
| `missense_testset.tsv` | the 65 distinct missense variants |
| `missense_by_gene.tsv` | per-gene comma-separated variant lists |
| `unresolved_alt_tx.json` | the six MANE reference-allele mismatches and the transcript scan that resolved them |
| `hgvs_queries.txt`, `vep_raw.json` | VEP queries and raw responses (GRCh38 coordinates) |
| `intervar_raw.json` | wInterVar API responses, per-criterion ACMG flags |
| `inheritance_hints.tsv` | LIRICAL zygosity + OMIM MOI per gene, for a future inheritance-corrected run |

## Reproducing the VarViz run

```bash
VARVIZ_UNIVERSE=analyses/derived/variant_universe_external163.tsv \
VARVIZ_OUT_DIR=analyses/ps_external \
  Rscript analyses/ps_final_harness.R
```

Copy `variant_universe_external163.tsv` to `analyses/derived/` first. The
harness aborts a gene rather than checkpointing it if any of the four
instrumented data sources degrades, so a completed 44/44 run is itself the
integrity check.

## Read this before quoting the numbers

**Sensitivity only.** Every variant is expert-curated causative, so no
specificity, false-positive rate or AUROC can be derived. A classifier
returning Pathogenic for everything scores 100%.

**One parameter set.** All 44 genes ran monoallelic with PM2 < 1×10⁻⁴ although
roughly 19 are recessive, which favours PM2 activation and biases toward the
pathogenic direction. `inheritance_hints.tsv` supports a corrected run.

**The published tool columns are truncated.** The capsule's per-tool outputs
cap at 100 variants per proband (Franklin) and 90 (GeneBe, InterVar, TAPES).
An empty cell means the variant fell outside that export, *not* that the tool
declined to classify it. Only the InterVar column headed `intervar_rerun` is a
like-for-like comparison; it was regenerated here through the public wInterVar
interface so InterVar saw every variant.
