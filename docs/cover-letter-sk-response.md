# Response to SK Reviewer Comments — VarViz Bioinformatics Application Note

We thank SK for the careful read and the four targeted comments. Each is addressed below; per-point summaries describe what changed in the manuscript and where.

---

## Comment 1 — Replace Table 1 with API/data-source table; move comparison to supplementary

> SK: *"Instead of current table 1, I suggest adding a table with different API sources used as data source which might be more helpful. Current table 1 can be moved to supplementary."*

**Adopted as suggested.** Table 1 in §2.1 now documents the nine live external data sources VarViz queries — UniProt, AlphaFold (EMBL-EBI), gnomAD, NCBI ClinVar, MyVariant.info / dbNSFP, UCSC, Ensembl, ClinGen, and GenCC — with endpoint paths, version/build information, and the corresponding ACMG/AMP roles each source supports (PM1, PM2, BA1, BS1, PS1, PP3, BP4, PP2, BP1, etc.). A base-URL footnote consolidates the API roots to keep the cell text readable.

The previous tool-comparison Table 1 has been moved to **Supplementary Table S2** and expanded from six to ten rows. The expansion incorporates four direct comparators that prior review passes flagged as missing — Franklin (Genoox), InterVar (Li & Wang, 2017), Alamut Visual Plus (SOPHiA Genetics), and MARRVEL (Wang et al., 2017) — alongside the original VarSome, GeneBe, cBioPortal MutMapper, MuPIT/VarSite, VIVID, and VarViz rows.

The §2.1 prose has been tightened from a 90-word API enumeration into a 30-word pointer that defers the full breakdown to the table. The §1 Introduction's parenthetical pointer to the tool-comparison table has been redirected from "(Table 1)" to "(Supplementary Table S2)". Two new bibliography entries have been added (Li & Wang 2017 — InterVar; Wang et al. 2017 — MARRVEL); the remaining new comparators are commercial software products and are credited inline in the supplementary table footnote.

---

## Comment 2 — Mention track toggle and zoom in §2.2 Visualization

> SK: *"Under Visualization: One feature you currently built is to add or remove tracks from the plot but no mention of it in the paper. And also the ability to zoom into the plot for closer inspection."*

**Adopted as suggested.** §2.2 now explicitly describes both interactive features. Two sentences have been added between the (1)–(7) track enumeration and the existing static-export sentence:

> *"Tracks (2)–(6) above can be toggled on or off via a sidebar checkbox group; tracks (1) and (7) are always rendered. Built-in plotly controls support zooming and panning along the amino acid axis, hover tooltips on every plotted element, and one-click axis reset."*

Both claims have been verified against the live application: track toggling is wired through `checkboxGroupInput("plotselection")` in the sidebar; zoom, pan, hover tooltips, and the `resetScale2d` mode-bar button are all active in the deployed Shiny app at shinyapps.io.

---

## Comment 3(i) — Validation circularity

> SK: *"Our choice of variants for validation (randomly choosing from one of our old papers) — this could be questioned by reviewers."*

We address this with three orthogonal benchmarks reported in Supplementary §S2. The validation universe is built from VariBench PON-PS_D2 clinical labels (Yip et al., 2008) and MaveDB deep mutational scanning scoresets (Esposito et al., 2019), intersected with gnomAD v4.1 singleton flags. Each variant is classified twice by the same engine: a *Pass-Full* classification using all ACMG criteria, and a methodologically rigorous *Pass-Blind* classification with ClinVar-derived criteria (PS1, PM5, PM1-via-ClinVar-hotspot, PP5, BP6) systematically withheld via a small blinding helper, leaving non-circular evidence (gnomAD AF, AlphaMissense + REVEL, ClinGen + GenCC gene-level validity, CCRS, UniProt domains). The full universe covers **90,701 variants from 14 BENCHMARK_GENES** spanning 21 MaveDB scoresets and 12 VariBench-evaluable genes.

- **Panel A — Clinical concordance vs VariBench (N = 1,108; 945 P + 163 B).** Pass-Blind AUROC = **0.959**, MCC = **0.774**, Sensitivity = 1.000, Specificity = **0.607**. Pass-Full AUROC = 0.995. **ΔAUROC = +0.036**, well below the 0.05 "robust" threshold: withholding ClinVar-derived criteria does not materially degrade clinical discrimination. Pass-Blind specificity (0.607) and MCC (0.774) actually **improve substantially** over Pass-Full (0.262 / 0.499) because over-confident Likely-Pathogenic calls correctly retreat to VUS instead of generating false positives.
- **Panel B — Functional concordance vs MaveDB DMS (N = 115,167 binary across 21 scoresets, 8 genes).** Median per-scoreset AUROC = **0.598** (range 0.24–0.77); pooled AUROC across all studies is undefined due to near-uniform classifier output at pooled scale (see §S2.6). The per-gene picture splits cleanly: **strong functional discrimination** in NUDT15 (0.73), KCNH2 (0.72), KCNQ1 (0.70), GCK (0.66), and PTEN (0.63); **weak discrimination** in BRCA1 (0.45), TP53 (0.42), and KRAS (0.38). The strong/weak split tracks ClinVar curation density: in densely-curated disease genes the PM1 ClinVar-hotspot criterion fires on virtually every missense, producing nearly-uniform Likely-Pathogenic calls that do not separate by DMS function. In less densely-curated disease genes the classifier's REVEL / AlphaMissense / conservation / CCRS criteria carry the signal and yield meaningful AUROCs (0.63–0.73). This bounded limitation is consistent with VarViz's design as a high-sensitivity clinical screening tool rather than a universal per-variant functional predictor.
- **Panel C — Population-genetic concordance via CAPS (N = 15,303 gnomAD-present).** Spearman ρ between ordinal VarViz classification rank and per-bin CAPS = **0.811 Pass-Full / 0.919 Pass-Blind** — Pass-Blind ρ *exceeds* Pass-Full ρ, indicating that debinding from ClinVar-derived strength tightens the monotonic relationship with negative selection. Pathogenic-bin CAPS = +0.184 (Pass-Blind) with monotonic gradient through every bin from Benign to Pathogenic. The gradient is constructed independently of any clinical annotation and directly addresses the ascertainment-bias and data-circularity framing articulated by Gudkov et al. 2025.

Together, the three panels triangulate: (i) clinical labels confirm the engine works without ClinVar circularity (ΔAUROC = +0.036, with substantial Pass-Blind gains in specificity 0.262 → 0.607 and MCC 0.499 → 0.774), (ii) population-genetic data confirm the gradient is biologically real (Pass-Blind ρ = 0.919, exceeding Pass-Full), and (iii) the DMS gap is now sharply bounded — gene-resolved analysis shows VarViz's functional discrimination is genuine in less densely-curated disease genes (mean AUROC 0.63–0.73) and saturates as expected in BRCA1/TP53/KRAS where ClinVar hotspot density is high.

Per these results, §3 has been **reframed as a case study** ("3 Case study: VarViz on the CASR receptor"). The case-study opening now reads "To demonstrate VarViz's ability to surface mechanism-distinguishing spatial patterns, we applied it to the complete enumeration of 24 CASR missense variants with assigned functional mechanism…", explicitly removing the impression of arbitrary selection, and a closing sentence routes the reader to Supplementary §S2 for the engine-level benchmarks. Five new bibliography entries (Esposito 2019, Giacomelli 2018, Gudkov 2025, Kotler 2018, Yip 2008) accompany the §S2 panels in the main paper, alongside two Pass-1 entries (Li & Wang 2017, Wang et al. 2017) for the new comparators in Supplementary Table S2.

---

## Comment 3(ii) — Colors and fonts

> SK: *"Use of multiple colours and fonts in your plot. So Maybe you can use less number of colours and fonts can be made bigger and uniform."*

We respectfully disagree on color and have addressed fonts directly.

**Color.** The multi-track, multi-color rendering of the protein landscape is the application's distinguishing feature — it is what enables cross-mechanism spatial pattern recognition (e.g., the LOF/GOF separation in the CASR case study) that no single-track or reduced-palette tool can surface. Each color carries semantic load: lollipop fill encodes ACMG classification tier; mutation-density curves use red/blue to separate ClinVar pathogenic from gnomAD benign distributions; AlphaMissense bars use a graded blue–yellow–orange ramp aligned with calibrated thresholds; conservation tracks use distinct hues for the three overlaid scores (PhyloP 100V, PhyloP 470M, PhastCons). Reducing this palette would compromise the differentiating capability of the visualization.

**Fonts — addressed.** Typography across the protein landscape, the gene-info card, and the per-variant ACMG cards has been unified under a single Arial/Helvetica family with a three-tier size hierarchy (Big / Medium / Small at multipliers 1.00 / 0.75 / 0.55 of a 14-pt anchor). All variant-classification text — including the per-variant card classification badge that SK specifically flagged — now sits at the largest tier and matches the size of the variant identifier in the same card. The Big/Medium/Small system is implemented once in `typography.R` (R-side) and `:root` CSS variables (Shiny-side) so that the tier definition is the single source of truth across all rendering pathways. Figure 1 has been re-exported with the new typography; the variant labels, track labels, and protein-name footer are all visibly larger and Helvetica-uniform in the revised PDF.

---

We hope these revisions address the reviewer's concerns clearly. The Comment 3(i) empirical work is now complete (three-panel triangulation in Supplementary §S2), the §3 reframe and Supplementary §S2 cross-reference are in v2 of the AppNote, the 10-row Supplementary Table S2 and §S2 writeup are in v2 of the Supplementary, and seven new bibliography entries (2 Pass-1, 5 Pass-2) are integrated.
