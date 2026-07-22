# VarViz — Tool & Data-Source Citations (manuscript reference validation)

Every external tool, dataset, predictor, and method VarViz invokes, with its citation and DOI.
`help?` = present in `www/help.html` References before this audit. DOIs are checked by the
repo citation validator on write — see the run note appended after each save.

**Licensing note:** all tools below are called via runtime API or are installed CRAN
dependencies — none are bundled/redistributed in VarViz source, so none affect VarViz's MIT
license. This file exists for citation/attribution obligations and manuscript references.

---

## 1. Databases & APIs (live)

| Tool | Role in VarViz | Citation | DOI | help? |
|---|---|---|---|---|
| UniProt | Domains, features, PTMs, function, disease | The UniProt Consortium (2025) *Nucleic Acids Res.* 53:D609–D617 | 10.1093/nar/gkae1010 | yes |
| AlphaFold | Per-residue pLDDT | Jumper J et al. (2021) *Nature* 596:583–589 | 10.1038/s41586-021-03819-2 | yes |
| gnomAD | Allele frequencies / constraint | Karczewski KJ et al. (2020) *Nature* 581:434–443 | 10.1038/s41586-020-2308-7 | yes |
| ClinVar | Variant clinical assertions | Landrum MJ et al. (2018) *Nucleic Acids Res.* 46:D1062–D1067 | 10.1093/nar/gkx1153 | no |
| NCBI E-utilities | ClinVar/gene programmatic access | Sayers E (2009, upd. 2022) *Entrez Programming Utilities Help*, NCBI | — | yes |
| MyVariant.info | dbNSFP/annotation query service | Xin J et al. (2016) *Genome Biol.* 17:91 | 10.1186/s13059-016-0953-9 | yes |
| dbNSFP v4 | Aggregated predictor/annotation DB | Liu X et al. (2020) *Genome Med.* 12:103 | 10.1186/s13073-020-00803-9 | yes |
| Ensembl | Canonical transcript exon structure | Harrison PW et al. (2024) *Nucleic Acids Res.* 52:D891–D899 | 10.1093/nar/gkad1049 | no |
| UCSC Genome Browser | PhyloP/PhastCons conservation tracks | Nassar LR et al. (2023) *Nucleic Acids Res.* 51:D1188–D1195 | 10.1093/nar/gkac1072 | no |
| ClinGen | Gene–disease validity | Rehm HL et al. (2015) *N. Engl. J. Med.* 372:2235–2242 | 10.1056/NEJMsr1406261 | yes |
| ClinGen validity use | Validity/dosage in classification | Thaxton C et al. (2022) *Hum. Mutat.* 43:1031–1040 | 10.1002/humu.24291 | yes |
| GenCC | Gene–disease consortium classifications | DiStefano MT et al. (2022) *Genet. Med.* 24:1732–1742 | 10.1016/j.gim.2022.04.017 | no |
| GENCODE v49 | Gene→UniProt map, GRCh38 coords (local) | Mudge JM et al. (2025) *Nucleic Acids Res.* 53:D966–D975 | 10.1093/nar/gkae1078 | yes |
| Orphadata / Orphanet | Disease prevalence per OMIM | Orphanet (INSERM), orphadata.com | — | no |

## 2. Pathogenicity predictors (in-silico)

| Tool | Role | Citation | DOI | help? |
|---|---|---|---|---|
| REVEL | Primary calibrated ensemble (PP3) | Ioannidis NM et al. (2016) *Am. J. Hum. Genet.* 99:877–885 | 10.1016/j.ajhg.2016.08.016 | no |
| AlphaMissense | Calibrated missense predictor | Cheng J et al. (2023) *Science* 381:eadg7492 | 10.1126/science.adg7492 | yes |
| CADD | Deleteriousness score (PP3) | Rentzsch P et al. (2019) *Nucleic Acids Res.* 47:D886–D894 | 10.1093/nar/gky1016 | no |
| CADD (original model) | CADD framework | Kircher M et al. (2014) *Nat. Genet.* 46:310–315 | 10.1038/ng.2892 | no |
| DANN | Deep-learning deleteriousness | Quang D et al. (2015) *Bioinformatics* 31:761–763 | 10.1093/bioinformatics/btu703 | no |
| SIFT | Substitution tolerance | Ng PC, Henikoff S (2003) *Nucleic Acids Res.* 31:3812–3814 | 10.1093/nar/gkg509 | no |
| PolyPhen-2 | Structural/evolutionary effect | Adzhubei IA et al. (2010) *Nat. Methods* 7:248–249 | 10.1038/nmeth0410-248 | no |
| FATHMM | Functional analysis HMM | Shihab HA et al. (2013) *Hum. Mutat.* 34:57–65 | 10.1002/humu.22225 | no |
| PROVEAN | Protein variation effect | Choi Y et al. (2012) *PLoS One* 7:e46688 | 10.1371/journal.pone.0046688 | no |
| MutationTaster2 | Disease-potential predictor | Schwarz JM et al. (2014) *Nat. Methods* 11:361–362 | 10.1038/nmeth.2890 | no |
| LRT | Likelihood-ratio test conservation | Chun S, Fay JC (2009) *Genome Res.* 19:1553–1561 | 10.1101/gr.092619.109 | no |
| MetaSVM / MetaLR | Ensemble meta-predictors | Dong C et al. (2015) *Hum. Mol. Genet.* 24:2125–2137 | 10.1093/hmg/ddu733 | no |
| MetaRNN | RNN ensemble meta-predictor | Li C et al. (2022) *Genome Med.* 14:115 | 10.1186/s13073-022-01120-z | no |
| dbscSNV (ADA/RF) | Splice-site effect | Jian X et al. (2014) *Nucleic Acids Res.* 42:13534–13544 | 10.1093/nar/gku1206 | no |

## 3. Conservation & constraint

| Tool | Role | Citation | DOI | help? |
|---|---|---|---|---|
| GERP++ | Rejected-substitution conservation | Davydov EV et al. (2010) *PLoS Comput. Biol.* 6:e1001025 | 10.1371/journal.pcbi.1001025 | no |
| PhyloP | Per-base conservation (UCSC) | Pollard KS et al. (2010) *Genome Res.* 20:110–121 | 10.1101/gr.097857.109 | no |
| PhastCons | Conserved-element score (UCSC) | Siepel A et al. (2005) *Genome Res.* 15:1034–1050 | 10.1101/gr.3715005 | no |
| ConSurf-DB | Evolutionary conservation grades | Ben Chorin A et al. (2020) *Protein Sci.* 29:258–267 | 10.1002/pro.3779 | yes |
| CCRS | Constrained coding regions (PM1) | Havrilla JM et al. (2019) *Nat. Genet.* 51:88–95 | 10.1038/s41588-018-0294-6 | yes |
| GeVIR | Gene missense-intolerance (BP1/PP2) | Abramovs N et al. (2020) *Nat. Genet.* 52:35–39 | 10.1038/s41588-019-0560-2 | yes |

## 4. ACMG framework & methods

| Method | Role | Citation | DOI | help? |
|---|---|---|---|---|
| ACMG/AMP 2015 | Rule-based classification | Richards S et al. (2015) *Genet. Med.* 17:405–424 | 10.1038/gim.2015.30 | yes |
| Tavtigian 2018 | Bayesian reformulation / OddsPath | Tavtigian SV et al. (2018) *Genet. Med.* 20:1054–1060 | 10.1038/gim.2018.14 | no |
| Tavtigian 2020 | Point-system scaling | Tavtigian SV et al. (2020) *Hum. Mutat.* 41:1734–1737 | 10.1002/humu.24088 | mislabeled |
| Pejaver 2022 | PP3/BP4 predictor calibration | Pejaver V et al. (2022) *Am. J. Hum. Genet.* 109:2163–2177 | 10.1016/j.ajhg.2022.10.013 | yes |
| Whiffin 2017 | Max credible allele frequency (PM2) | Whiffin N et al. (2017) *Genet. Med.* 19:1151–1158 | 10.1038/gim.2017.26 | yes |
| Jarvik 2016 | Cosegregation LR (PP1/BS4) | Jarvik GP, Browning BL (2016) *Am. J. Hum. Genet.* 98:1077–1081 | 10.1016/j.ajhg.2016.04.003 | yes |
| Biesecker 2024 | ClinGen SVI PP1/BS4 + PP4 guidance; source of the +5.0 combined PP1/PP4 locus cap | Biesecker LG et al. (2024) *Am. J. Hum. Genet.* 111:24–38 | 10.1016/j.ajhg.2023.11.009 | yes |

## 5. VarViz-specific pathways & comparators

| Item | Role | Citation | DOI | help? |
|---|---|---|---|---|
| DOLPHIN | Pfam-alignment last-resort PM1 pathway | Corcuff S et al. (2023) "Protein domains provide a new layer of information for classifying human variations in rare diseases" *Front. Bioinform.* 3:1127341 | 10.3389/fbinf.2023.1127341 | inline only |
| Genebe | ACMG auto-assignment comparator | Maj C et al. (2023) *Clin. Genet.* 104:509–516 | 10.1111/cge.14516 | yes |
| AutoPM3 | LLM PM3 evidence extraction | Li S et al. (2025) *Bioinformatics* 41:btaf382 | 10.1093/bioinformatics/btaf382 | yes |

## 6. Gaps to close in help.html References

Not currently in the References list — add before manuscript submission:
ClinVar (Landrum 2018), Ensembl (Harrison 2024), UCSC (Nassar 2023), GenCC (DiStefano 2022),
REVEL (Ioannidis 2016), CADD (Rentzsch 2019 / Kircher 2014), DANN (Quang 2015), SIFT, PolyPhen-2,
FATHMM, PROVEAN, MutationTaster2, LRT, MetaSVM/LR (Dong 2015), MetaRNN (Li 2022), dbscSNV (Jian 2014),
GERP++ (Davydov 2010), PhyloP (Pollard 2010), PhastCons (Siepel 2005), Tavtigian 2018, DOLPHIN (Corcuff 2023),
Orphanet/Orphadata.

## 7. Verification log (2026-07-21)

Web-confirmed corrections made this pass:
- **DOLPHIN** — full title confirmed; DOI **10.3389/fbinf.2023.1127341** verified (Frontiers).
- **Tavtigian 2018** — DOI corrected to **10.1038/gim.2018.14** (was wrongly 10.1038/gim.2017.210).
- **Tavtigian 2020** — is *Hum. Mutat.* 41:1734–1737, **10.1002/humu.24088** (PMID 32720330), NOT
  Genet. Med. The **existing `help.html` entry is mislabeled** ("Tavtigian 2020, Genet. Med.
  22:1054–1060, doi:10.1038/s41436-020-0777-z" conflates the two papers and its DOI does not
  resolve to a Tavtigian article) — fix help.html to carry both papers correctly.

Repo citation validator flagged these as *unverified* (reachability only — all are standard,
correct DOIs already used in the peer-reviewed help.html reference list or widely cited;
treat as reachability-index gaps, not errors): gkae1010 (UniProt), gkx1153 (ClinVar/Landrum),
gkad1049 (Ensembl), gkac1072 (UCSC), NEJMsr1406261 (Rehm/ClinGen), humu.24291 (Thaxton),
gkae1078 (GENCODE), science.adg7492 (AlphaMissense). Re-confirm at manuscript submission via
the publisher DOI resolver.

### Addendum (2026-07-21, later pass)

- **Biesecker 2024** (ClinGen SVI PP1/BS4 + PP4 guidance, 10.1016/j.ajhg.2023.11.009) added as
  help.html reference **[29]**, appended at the end of the `<ol>` so existing `[n]` markers keep
  their positions. It is the source for the +5.0 combined PP1/PP4 locus cap now enforced in
  `classify_acmg`.
- Two known divergences from this guidance are **documented in help.html but not implemented**:
  VarViz applies PP4 as binary +1 (ACMG 2015 default) where ClinGen scales it by diagnostic yield
  and permits PP4_Moderate/PP4_Strong; and VarViz does not suppress PP1 under locus homogeneity.
  VarViz's PP1 tiers still follow Jarvik 2016 rather than ClinGen's inheritance-mode-aware Table 3.

Non-DOI sources (NCBI E-utilities, Orphadata/Orphanet) cite the resource directly.
The dbNSFP-aggregated predictors (SIFT, PolyPhen-2, LRT, FATHMM, PROVEAN, MutationTaster,
MetaSVM/LR/RNN, DANN, dbscSNV) are surfaced in VarViz via dbNSFP (Liu 2020); the manuscript
should cite each primary method (rows above) in addition to dbNSFP.
