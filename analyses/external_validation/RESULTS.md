# External validation — Ghasemnejad et al. manual-curation set

Source: Code Ocean capsule 4308488, `data/manual_163new.xlsx` (CC0), the
expert-curated causative variants underlying "Comprehensive Evaluation of
ACMG/AMP-based Variant Classification Tools".

Assembled 2026-08-24 for reviewer item P3 (robustness on external variants).

## 1. Test set construction

152 rows → 143 alleles (4 rows carry two variants, split on `,` / `;`).
Two notation styles present, both handled: standard `c.289C>T` and legacy
reversed `c.G379T` (15 rows).

| category | alleles |
|---|---|
| coding substitution | 111 |
| indel | 20 |
| intronic | 11 |
| empty | 1 |

Coding substitutions resolved against MANE Select via Ensembl VEP REST:

| consequence | alleles |
|---|---|
| missense | 76 |
| missense + splice_region | 8 |
| stop_gained | 23 |
| unmappable | 3 |
| ambiguous | 1 |

**Test set: 65 distinct missense variants, 44 genes, 83 probands.**
Zero gene overlap with the 14-gene development universe.

Six MANE reference-allele mismatches resolved by scanning all coding
transcripts of each gene: CRB1 p.C790F and TBX3 p.A549D retained (non-MANE
transcript, flagged); CEP128 c.1657G>T (stop_gained), PCDH15 c.4726C>T
(consequence differs across the 6 transcripts matching the reference allele),
COL1A1 c.11299G>C (beyond CDS on all 19 coding transcripts), FAM20A c.625G>T
(no transcript matches the stated reference) and `ELAINE` c.377C>T (not an
HGNC symbol) dropped.

Sanity check: BRAF c.1799T>A resolves to p.V600E.

## 2. VarViz result

Run: `analyses/ps_external/`, harness `analyses/ps_final_harness.R` with
`VARVIZ_UNIVERSE=analyses/derived/variant_universe_external163.tsv`.
44/44 genes, 65/65 variants, all four data-source sentinels clean.

| | Pathogenic | Likely Path | VUS-High | VUS-Mid | P+LP |
|---|---|---|---|---|---|
| Pass-Full | 36 | 20 | 8 | 1 | **56/65 = 86.2%** |
| Pass-Blind | 7 | 38 | 15 | 5 | **45/65 = 69.2%** |

42/65 (64.6%) change class under ClinVar blinding — a larger penalty than the
14-gene benchmark, because this cohort leans harder on ClinVar lookups.

Excluding COL1A1/COL1A2 (Gly-X-Y motif genes that favour MDS), on the
remaining 50: Pass-Full 84.0%, Pass-Blind 64.0%. Not a collagen artifact.

Per proband: the causative variant reaches P/LP in 70/83 (84.3%) Pass-Full,
58/83 (69.9%) Pass-Blind.

### PM1 pathway distribution (Pass-Full)

| pathway | n |
|---|---|
| mds | 26 |
| clinvar_hotspot | 20 |
| ccrs | 5 |
| mds + hotspot_upgrade | 5 |
| ccrs + mds | 2 |
| ccrs + hotspot_upgrade | 2 |
| uniprot_domain + mds | 1 |
| uniprot_site + mds | 1 |
| mds_unavailable | 1 |
| none | 2 |

**MDS originates or contributes PM1 on 32/65 (49%)** — the largest single
pathway, larger than clinvar_hotspot, on genes with no overlap with the
development set.

### The 9 misses

| gene | variant | VarViz | PM1 pathway | other tools |
|---|---|---|---|---|
| ADAMTS1 | p.E752K | VUS-Mid | none | InterVar VUS |
| COL1A1 | p.S1199G | VUS-High | clinvar_hotspot | — |
| ESRRB | p.E167K | VUS-High | none | — |
| P3H1 | p.Q714R | VUS-High | clinvar_hotspot | Franklin LP |
| TBX3 | p.A549D | VUS-High | mds_unavailable | — |
| WFS1 | p.K193Q | VUS-High | clinvar_hotspot | — |
| WNT10A | p.R70W | VUS-High | mds | — |
| WNT10A | p.R171C | VUS-High | mds | Franklin Benign |
| WNT10A | p.G213S | VUS-High | mds | Franklin Benign |

Three of nine are WNT10A hypomorphic alleles with reduced penetrance and
appreciable population frequency, where VUS-High is arguably better supported
than Franklin's Benign. TBX3 p.A549D uses a non-MANE transcript, so
`mds_unavailable` may reflect the numbering mismatch rather than a real miss.

## 3. Comparison with the benchmarked tools

### 3.1 The published tool outputs are truncated per proband

Franklin caps at 100 variants per proband, Genebe/InterVar/TAPES at 90.
80 of 151 Franklin probands sit exactly at 100; 131 of 151 Genebe probands sit
exactly at 90.

| tool | probands | min rows | median | max | n at max |
|---|---|---|---|---|---|
| Franklin | 151 | 8 | 100 | 100 | 80 |
| Genebe | 151 | 76 | 90 | 90 | 131 |
| InterVar | 151 | 6 | 74 | 90 | 7 |
| TAPES | 151 | 3 | 84 | 90 | 7 |

**Absence of a call means the variant fell outside that export, not that the
tool declined to classify it.** Raw sensitivity across all 65 is therefore not
comparable between VarViz (which saw all 65) and the published tool columns
(which did not). All 83 of our probands appear in all four tool outputs, so
nothing was lost at the proband level.

Of the 28 unmatched Franklin sample-alleles: 23 are genes never in that
proband's top-100; 5 are probands where a different variant of the same gene
was reported instead (mostly probands carrying two of our variants in one
gene).

### 3.2 InterVar re-run on all 63 — a fair comparison

InterVar was re-queried through the wInterVar public API (hg38, no
authentication) so it saw every variant. 63 of 65 have hg38 coordinates; CRB1
p.C790F and TBX3 p.A549D are excluded because they sit on non-MANE
transcripts.

| | P/LP | sensitivity | Benign-side calls |
|---|---|---|---|
| VarViz Pass-Full | 55/63 | **87.3%** | 0 |
| VarViz Pass-Blind | 45/63 | **71.4%** | 0 |
| InterVar | 20/63 | **31.7%** | 1 |

InterVar returns Uncertain significance on 42 of 63 and never reaches
Pathogenic on any variant in this set — Likely Pathogenic is its ceiling here.
VarViz Pass-Blind, with all ClinVar-derived evidence withheld, still more than
doubles InterVar's rate.

**The re-run validates the published data.** All 9 overlapping variants agree
exactly between the paper's truncated InterVar export and the fresh wInterVar
query (8 Likely Pathogenic, 1 VUS, zero discrepancies), confirming both that
the published calls are genuine and that wInterVar runs the same InterVar the
authors used. It also confirms the truncation diagnosis from the other
direction: InterVar's 9 calls in the benchmark were not InterVar declining on
54 variants.

### 3.3 Franklin — head-to-head only

Franklin's API is a paid Genoox service, so it could not be re-run. The
defensible comparison is restricted to variants where both tools produced a
call:

| | n | P/LP |
|---|---|---|
| VarViz Pass-Full | 43 | 40 (93.0%) |
| Franklin | 43 | 39 (90.7%) |

A tie. VarViz made zero Benign/Likely Benign calls on causative variants;
Franklin made two (WNT10A p.R171C, p.G213S).

### 3.4 Genebe and TAPES — not re-run

Genebe requires an account. TAPES has no public API and needs
ANNOVAR-annotated input; ANNOVAR is installed locally but `humandb/` is empty,
so a re-run needs ~50-60 GB of annotation databases. Their published columns
(10 and 2 calls respectively) are reported for completeness but are too sparse
for a head-to-head.

### 3.5 LIRICAL

Different task — phenotype-driven gene ranking, not variant classification.
On the same 83 probands:

| metric | n/83 | % |
|---|---|---|
| top-1 | 39 | 47.0 |
| top-5 | 58 | 69.9 |
| top-10 | 63 | 75.9 |
| top-20 | 65 | 78.3 |
| causative gene absent from list | 15 | 18.1 |

Median rank 1 when found. Our 75.9% top-10 exceeds the paper's headline 68.21%
because this subset is restricted to probands whose causative variant is a
clean missense.

These numbers must not be tabulated alongside classification sensitivity.
LIRICAL searches a whole exome and returns a ranked gene shortlist from HPO
terms; VarViz is handed one variant and asked whether it is pathogenic. The
supportable statement is complementarity: LIRICAL narrows an exome to
candidate genes from phenotype, VarViz supplies the ACMG classification and
evidence trail for the survivors. That is a concrete, data-backed statement
about workflow position without claiming a user study.

## 4. Limitations

1. **No benign arm.** Every variant is causative by construction, so this
   measures sensitivity only — no specificity, no FPR, no AUC. A classifier
   returning Pathogenic for everything scores 100%.
2. **Single parameter set.** All 44 genes ran monoallelic with PM2 < 1e-4,
   while roughly 19 are recessive (ABCA4, USH2A, MYO15A, OTOF, PCDH15,
   TMPRSS3, SLC26A4, GALC, RPE65, TULP1, CNGA3, RDH12, CRB1, PJVK, ESRRB,
   MYO7A, WNT10A, FAM20A, ENAM). Recessive genes tolerate higher allele
   frequencies, so this over-fires PM2 in the pathogenic direction.
   Inheritance hints from LIRICAL zygosity and OMIM are in
   `inheritance_hints.tsv` for a corrected run. EDA is X-linked and was
   separately re-run in the app at prevalence 1 in 25,000, hetA 0.1, hetG 0.8,
   penetrance 1.0 (af_cutoff 1.6e-6); the call did not change because the
   variant is absent from gnomAD.
3. **Tool export truncation** — §3.1.
4. CRB1 p.C790F and TBX3 p.A549D use non-MANE transcripts; their positions do
   not align with VarViz's MANE-based numbering.
5. Franklin, Genebe and TAPES could not be re-run, so only InterVar has a
   like-for-like comparison across the full set.

## 5. Incidental finding — ClinVar esummary 10 MB cap

NCBI's esummary XML→JSON converter rejects payloads over 10 MB with an
HTTP 200 carrying no `$result`:

    {"ERROR": "Input XML size is 12150021 bytes, and cannot be transformed
     to JSON. the max size is 10MB"}

EDA's third batch of 92 IDs is 12.15 MB, so the existing 4-attempt retry
failed identically every time — deterministic, not a rate limit. Fixed in
commit `97d61aa` by replacing the fixed-stride batch loop with a work queue
that halves a refused batch and re-queues both halves. EDA now extracts 250
ClinVar records instead of 216. This affected the live app as well as the
harness: any gene with heavy ClinVar annotation was silently losing its final
batch.

## 6. Files

| file | contents |
|---|---|
| `TABLE_external_comparison.tsv` | **supplementary-ready**: all 65 variants × VarViz full/blind/points/pathway × InterVar re-run × 4 published tool columns |
| `variants_parsed.tsv` | all 143 alleles, parsed and categorised |
| `variants_annotated.tsv` | + VEP consequence, protein position, p. notation |
| `missense_testset.tsv` | the 65 distinct missense variants |
| `missense_by_gene.tsv` | per-gene comma-separated variant lists |
| `inheritance_hints.tsv` | LIRICAL zygosity + OMIM MOI per gene |
| `intervar_raw.json` | wInterVar API responses, per-criterion ACMG flags |
| `../../../derived/variant_universe_external163.tsv` | harness input universe |
| `../../../ps_external/summary.tsv` | VarViz dual-pass output |
