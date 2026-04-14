<<<<<<< HEAD
# varviz
Protein-centric Gene Variant Visualization. 
VarViz displays multiple layers of protein-level annotation on a shared amino-acid axis so you can evaluate variant significance at a glance. All data is fetched live from public APIs; no pre-downloaded databases are required beyond a gene-coordinate file and CCRS regions.


Artificial intelligence disclosure: AI coding assistance was used for the code review and improvements of the project
=======
# VarViz

**A protein-centric web application for interactive missense variant visualisation and automated ACMG/AMP pathogenicity classification**

[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![R Shiny](https://img.shields.io/badge/Built%20with-R%20Shiny-blue.svg)](https://shiny.posit.co/)
[![Live App](https://img.shields.io/badge/Live%20App-varviz.shinyapps.io-brightgreen.svg)](https://varviz.shinyapps.io/varviz/)

---

## Overview

Interpreting a panel of missense variants in a disease gene requires simultaneously consulting population databases, structural predictors, conservation resources, and clinical repositories. These tools share no common protein-level view and must be queried and collated by hand for each variant.

VarViz addresses this by fetching **nine evidence streams** through live API calls and rendering them as vertically aligned, interactive tracks on a shared amino acid axis, with automated ACMG/AMP classification combining Richards et al. (2015) combinatorial rules, Tavtigian et al. (2020) Bayesian point scoring, and Pejaver et al. (2022) calibrated PP3/BP4 thresholds.

![VarViz workflow and CASR validation](docs/VarViz_Figure1_AB.png)

**(A)** VarViz workflow architecture. **(B)** 24-variant CASR protein landscape showing spatial separation between loss-of-function (mean AlphaMissense 0.83) and gain-of-function variants (mean 0.39).

---

## Live application

**[https://varviz.shinyapps.io/varviz/](https://varviz.shinyapps.io/varviz/)**

No installation required. A modern web browser is sufficient. No user data are stored between sessions.

---

## Features

### Nine interactive evidence tracks
| Track | Data source |
|---|---|
| AlphaFold pLDDT structural confidence | AlphaFold EBI API |
| gnomAD v4 allele frequency rainfall | gnomAD GraphQL API v4 (GRCh38) |
| ClinVar pathogenic / gnomAD benign mutation density | NCBI E-utilities |
| AlphaMissense mean pathogenicity per residue | AlphaFold EBI API |
| ClinVar variant ticks / PTM sites / CCRS constrained regions | NCBI + local CCRS data |
| ConSurf evolutionary conservation (uploaded grades file) | User upload |
| Multi-species conservation (PhyloP 100V, 470M; PhastCons) | UCSC + Ensembl REST APIs |
| UniProt domain architecture | UniProt REST API |
| Lollipop variant map | Assembled from all above |

### Automated hybrid ACMG/AMP engine
- **Richards et al. (2015)** combinatorial rule-based classification (first pass)
- **Tavtigian et al. (2020)** Bayesian point scoring when no rule matches (≥10 Pathogenic, 6-9 LP, -3 to +5 VUS, -4 to -6 LB, <=-7 Benign)
- **Pejaver et al. (2022)** calibrated PP3/BP4 thresholds across seven predictors (REVEL, CADD, MetaSVM, MetaLR, MetaRNN, DANN, AlphaMissense) with Supporting/Moderate/Strong tiers
- **PM1** via three pathways: CCRS ≥90th percentile, UniProt functional domain, or 15-residue ClinVar hotspot neighbourhood
- **PM2** against a prevalence-derived maximum credible allele frequency (Whiffin et al., 2017)
- **PS1** strength scaled by ClinVar review star count (GeneBe approach)
- **PP2/BP1** from ClinGen gene-disease validity or GeVIR percentile
- **Three evidence-independence guards** prevent double-dipping (PP5 suppressed by PS1; conservation withheld from PP3 when used for PM1_strong; PM1 neighbourhood independent of PM2)

### Per-variant clinical inputs
Each variant card provides real-time dropdowns for:
- **PP1** cosegregation (Supporting / Moderate / Strong; Jarvik & Browning, 2016)
- **PS2 / PM6** de novo status (confirmed / assumed)
- **PM3** compound heterozygosity (Supporting / Moderate / Strong; Li et al., 2025)

Classification and explanatory comment update immediately on any change.

### Reproducible exports
- **TSV download**: all ACMG tags, evidence values, in silico scores, population frequencies, and analysis parameters; any classification can be exactly reproduced from the file alone
- **gnomAD raw data TSV**: full AC/AN/AF and ancestry-specific frequencies for all variants in the queried gene
- **Static figure**: publication-quality PDF/PNG/JPEG via cowplot with all tracks aligned

---

## How to use

### 1. Select a gene
Type or select a HGNC gene symbol in the **Choose a gene** dropdown. VarViz resolves it to a UniProt accession and GRCh38 coordinates via a pre-computed GENCODE v49 lookup table.

### 2. Enter variants
Enter missense variants as a comma-separated list in **p-notation** (e.g. `p.Arg220Trp, p.Gly509Arg, p.Thr676Arg`) in the Variants field, or upload a plain-text file (one variant per line). The `p.` prefix is added automatically if omitted.

### 3. Set analysis parameters (optional)
Click **Edit parameters** to configure:

| Parameter | Default | Notes |
|---|---|---|
| Inheritance | Monoallelic | Switches BS1/PM2 thresholds for biallelic genes |
| Prevalence | 1 in 2,000 | Used to compute PM2 maximum credible AF |
| Allelic heterogeneity | 0.2 | Fraction of cases attributable to this gene |
| Genetic heterogeneity | 1.0 | Fraction attributable to this locus |
| Penetrance | 0.5 | Used in AF cutoff calculation |
| Confidence | 0.95 | Wilson CI for AF threshold |

### 4. Upload ConSurf grades (optional)
Download a ConSurf grades file from [consurfdb.tau.ac.il](https://consurfdb.tau.ac.il) for your protein of interest and upload it using the **ConSurf grades file** browser. The ConSurf conservation track appears automatically.

### 5. Click Go
All eight external APIs are queried in parallel. Results appear within 10-30 seconds depending on gene size and network conditions. API responses are cached in-session; re-querying the same gene is instantaneous.

### 6. Explore the protein landscape
- **Hover** any track element for detailed values
- **Click and drag** to zoom; double-click to reset
- **Toggle tracks** using the Output Tracks checkboxes in the sidebar

### 7. Enter clinical evidence per variant
Scroll to the **Variant Summary** tab. Each variant card has dropdowns for cosegregation, de novo status, and compound heterozygosity. The ACMG score, classification tier, and explanatory comment update in real time.

### 8. Export results
- **Download the plot**: static PDF/PNG/JPEG of all active tracks
- **Download TSV**: full annotation table with all ACMG tags and analysis parameters
- **Download gnomAD data**: raw gnomAD v4 variant table for the queried gene

---

## Validation

Applied to **24 CASR missense variants** with functionally validated loss-of-function (FHH1) and gain-of-function (ADH1) disease mechanisms from a population genomic cohort of 51,289 individuals (Gorvin et al., 2020):

- **Sensitivity 93%** for loss-of-function variants (14/15 classified as P/LP; Wilson 95% CI 70-99%)
- The protein landscape revealed spatial separation between LOF and GOF variants in the AlphaMissense mean pathogenicity track; this pattern was not recoverable from any individual predictor score

---

## Local installation

VarViz is designed as a web application and runs without installation at the link above. To run locally:

```r
# Prerequisites: R >= 4.2, the following packages
install.packages(c(
  "shiny", "shinyjs", "plotly", "ggplot2", "cowplot",
  "httr", "jsonlite", "data.table", "dplyr", "stringr",
  "BiocManager"
))

# Clone and run
git clone https://github.com/metviz/varviz.git
cd varviz
Rscript -e "shiny::runApp('.')"
```

The app will open at `http://127.0.0.1:XXXX` in your browser. Internet access is required for all eight external API calls.

---

## Data sources

| Source | Data retrieved | Reference |
|---|---|---|
| UniProt REST API | Protein domains, PTMs, disulfide bonds, topological regions | UniProt Consortium, 2025 |
| AlphaFold EBI API | Per-residue pLDDT, AlphaMissense mean pathogenicity | Jumper et al., 2021; Cheng et al., 2023 |
| gnomAD GraphQL API v4 | Allele frequencies, allele counts, homozygotes (GRCh38 / MANE Select) | Karczewski et al., 2020 |
| NCBI E-utilities | ClinVar variant assertions, review stars, VCV IDs | Landrum et al., 2018 |
| MyVariant.info | dbNSFP v4.4 (15 in silico predictor scores) | Xin et al., 2016; Liu et al., 2020 |
| UCSC + Ensembl REST | PhyloP 100V, PhyloP 470M, PhastCons per-residue scores | |
| ClinGen / GenCC APIs | Gene-disease validity classification | Rehm et al., 2015 |
| GeVIR | Gene-level missense constraint percentile | Abramovs et al., 2020 |
| CCRS (local) | Constrained coding region scores | Havrilla et al., 2019 |

---

## Citation

If you use VarViz in your research, please cite:

> [Author names] (2025). VarViz: A Protein-Centric Web Application for Interactive Missense Variant Visualization and Automated ACMG/AMP Classification. *Bioinformatics* (submitted).

---

## Licence

MIT. See [LICENSE](LICENSE) for details.

---

## Disclaimer

> **For research use only. VarViz is not a clinical diagnostic tool and must not be used as the sole basis for clinical decision-making.**
>
> The authors and contributors accept no liability for any clinical, diagnostic, or therapeutic decisions made on the basis of outputs produced by this software. Use is entirely at the user's own risk.

VarViz is an open-source research tool intended to support variant triaging, hypothesis generation, and evidence review in research and pre-clinical settings. The following limitations apply to all use of this software:

**Clinical use**
All ACMG/AMP classifications generated by VarViz are automated and must be reviewed and confirmed by a qualified clinical geneticist or certified variant interpretation laboratory before any clinical action is taken. VarViz has not been validated as a clinical diagnostic tool, has not received regulatory clearance in any jurisdiction, and does not constitute medical advice. Variant classifications should not be communicated to patients or used to guide treatment decisions without independent expert review.

**Classification accuracy**
The ACMG/AMP engine implements published guidelines (Richards et al., 2015; Tavtigian et al., 2020) but automated application of any classification framework involves assumptions and approximations. Criteria that require clinical judgement (PS2, PS3, PM3, PP1) are only partially automatable; VarViz provides per-variant dropdowns for these inputs but cannot substitute for direct clinical assessment. Classifications may differ from those produced by accredited diagnostic laboratories.

**Data provenance**
All variant data are fetched live from third-party public databases (gnomAD, ClinVar, UniProt, AlphaFold, MyVariant.info, UCSC, Ensembl, ClinGen, GeVIR) at the time of analysis. VarViz does not control the content, accuracy, or availability of these resources. Data may change between sessions. gnomAD v4 annotates variants using the MANE Select transcript; variants submitted in p-notation derived from an alternative transcript may differ in residue numbering and may be reported as absent from gnomAD when they are in fact present under a different residue number. Users should verify the transcript underlying their p-notation when a variant is reported as absent.

**No data storage**
No user-submitted gene or variant data are stored by VarViz between sessions. All analysis is performed client-server within the active session only.

**Third-party API availability**
VarViz depends on eight external APIs. Service interruptions at any upstream provider will affect the corresponding track or classification criteria. The authors make no guarantee of uninterrupted availability.

**Liability**
The authors, contributors, and institutions associated with VarViz accept no responsibility or liability for any loss, harm, or damage arising from the use or misuse of this software or its outputs. The software is provided "as is" without warranty of any kind, express or implied, including but not limited to warranties of merchantability, fitness for a particular purpose, or non-infringement. In no event shall the authors be liable for any direct, indirect, incidental, special, or consequential damages arising from use of this software.

**Not a substitute for database subscriptions**
VarViz retrieves publicly available data only. It does not access proprietary databases, internal laboratory knowledge bases, or subscription-only resources. Classifications are based solely on publicly deposited evidence at the time of query.

---

## Artificial intelligence disclosure

Consistent with the ISCB acceptable use policy, we disclose that AI-assisted tools (Grammarly and Claude) were used in this project to aid in editing the manuscript, code review, development, and improvements. All scientific content, analyses, interpretations, and conclusions are the work of the authors. The authors take full responsibility for the accuracy and integrity of the submitted work.

---

## Notes

- gnomAD v4 uses the MANE Select transcript; p-notation from older transcripts may differ in residue numbering
- Network dependence on eight external APIs means availability may vary with API outages
- The server will be maintained for a minimum of three years following publication
- Bug reports and feature requests are welcome via the GitHub Issues tab
>>>>>>> 3660da8 (docs: update comprehensive README with workflow figure and usage guide)
