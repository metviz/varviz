# Changelog

VarViz follows [semantic versioning](https://semver.org): MAJOR for changes that
alter classifications, MINOR for new evidence sources or user-facing features,
PATCH for fixes that leave every call unchanged.

Because this tool assigns ACMG classifications, each entry states explicitly
whether it can move a variant's call.

## [2.2.1] - 2026-09-06

### Fixed

- **`PM1_Derivation` was scoped inside the BS1 branch.** `pm1_deriv_val` was
  initialised inside `if (!bs1_fires)`, but the export reads it for every
  variant. When BS1 fires the PM1 ladder is skipped, so the variable is either
  undefined (`object 'pm1_deriv_val' not found`, aborting the gene when the
  first variant of that gene is common) or, worse, still holds the PREVIOUS
  variant's derivation string and reports it as this variant's provenance. Now
  initialised per variant alongside `pm1_pathway_val`. Found when CASR
  p.Ala986Ser (gnomAD AF 0.134) aborted a full-universe run; the same latent
  fault predates 2.2.0.

## [2.2.0] - 2026-09-06

Predictor-lookup correctness and run-integrity guards. **Changes
classifications** wherever the dbNSFP protein-key collision below applied: every
`ps_*` run and every case-study export produced before this release must be
regenerated.

### Fixed

- **dbNSFP protein lookup returned a different substitution.**
  `load_dbnsfp_for_region()` indexed rows by `(aapos, aaalt)` only. A dbNSFP row
  carries a `;`-separated multi-transcript `aapos` list, so one row is indexed at
  several positions and first-write-wins across substitutions that collide.
  SNCA's canonical `G51D` and an alternate-transcript `N51D` both key to `51_D`;
  a `G51D` query was served the `N51D` row (REVEL absent, MetaSVM and MetaLR
  tolerated) instead of its own (REVEL 0.72, both damaging). PP3 stayed
  Supporting instead of reaching Moderate, costing a point and moving the call
  from Pathogenic to Likely Pathogenic. The index now also carries an
  `(aapos, aaref, aaalt)` key; `lookup_dbnsfp_by_aa()` gains `aa_ref` and, when
  given one, requires an exact match, returning `NULL` on mismatch so the caller
  falls back rather than scoring another variant's predictors. All three
  harnesses now pass the reference residue.

### Added

- **Degraded predictor caches are reported and refused.** REVEL, AlphaMissense
  and UCSC conservation load independently of dbNSFP, and a missing one lowers
  PP3 strength without changing the row count. `record_cache_degradation()`
  records any gene missing a cache or below 95% UCSC coverage; the run lists
  them and exits non-zero unless `--allow-degraded` is passed.
- `analyses/tests/test_dbnsfp_aa_lookup.R`, built on the real G51D/N51D
  collision.

### Changed

- **Analysis parameters match the manuscript.** The harnesses hardcoded allelic
  heterogeneity 0.5, penetrance 1.0, 125,748 alleles and `ac_cutoff` 13, giving a
  maximum credible allele frequency of 1.25e-4. Supplementary S1 states 1.0e-4,
  which comes from the app defaults (0.2 / 0.5 / 251,496 / 34). All three
  harnesses now use those, and read per-gene overrides from the universe so SNCA
  runs at its own prevalence 1 in 10,000, allelic heterogeneity 0.5 and
  penetrance 0.6 (max AF 4.17e-5, max AC 16).
- Headless plots save on white. `densityplot()` and `clinvar_ccrsplot()` add only
  `theme()` overrides, so their background is transparent; the app paints white
  through plotly but a PNG device writes RGB without alpha and flattened it to
  black, making half of every headless figure black.
- The CCRS panel's y-axis title read `ClinVar/PTMs/CCRs/PS`. The row label was
  renamed to Hot in 1.1.0 but the axis title was missed.

## [2.1.1] - 2026-09-06

User-facing text only. **No classification changes** — no engine code was
touched, only the strings the app displays.

### Changed

- Landing page, feature list and Variant Summary legend described the engine as
  "hybrid Richards 2015 rule-based + Tavtigian Bayesian scoring". Classification
  has been points-only since 2.0.0, so they now read: Richards et al. (2015)
  criteria, scored on the Tavtigian et al. (2020) point scale with Pejaver et al.
  (2022) PP3/BP4 calibration and ClinGen SVI strength recommendations.
- PP5 and BP6 tooltips say the criterion is shown for context and not scored
  (retired by ClinGen SVI).
- PS3 tooltip says the published mutagenesis is for *this* substitution, and the
  site branch no longer implies a PTM annotation alone establishes function.

## [2.1.0] - 2026-09-06

MDS strength tiers reduced to two. **No change to the default configuration**
(`varviz.mds_tiered = FALSE`, MDS at Moderate); the opt-in tiered configuration
now yields Strong instead of Moderate-plus for MDS ≤ −8 to −12.

### Changed

- **`PM1_moderate_plus` (+3) retired.** The Tavtigian point scale has no
  half-step between Moderate (2) and Strong (4), and the LR+ at MDS ≤ −8 lands
  between the two strength cut-points on every benign arm measured (11.2
  genome-wide 2-star, 14.5 on 1-star-only, 12.1–17.4 against the Kwon et al.
  2026 control sets), so that range is reported as a confident Moderate rather
  than a separate tier. With the exploratory tier enabled, PM1 from MDS is
  Moderate (+2) at ≤ −4 and Strong (+4) at ≤ −12; nothing in between.
- The retired tag is still recognised by `lib/clinvar_blind.R` and
  `repro/07_regenerate.R` so summaries written by earlier versions still
  demote and tabulate correctly.

## [2.0.1] - 2026-09-06

Analysis scripts only. **No classification changes.**

### Fixed

- `24_clinvar_mds_lr.R` had been unrunnable since 13 Aug: the rebuilt
  `data/pfam_pssm_human.rds` stores `family` as a factor and `nzchar(factor)`
  errors. Now `as.character()`. The manuscript LR table reproduces exactly.

### Added

- `25_clinvar_extract.R --stars=1` writes a separate ≥1-star ClinVar missense
  set; `24_clinvar_mds_lr.R --in= --review=single|all --out=` scores it. Used
  for an exploratory out-of-sample check of the MDS tiers (REPRODUCIBILITY.md
  "Exploratory"): Strong replicates on 1-star-only variants (LR+ 32.1,
  95% CI 20.6–50.2 at MDS ≤ −12). Not a manuscript number; defaults unchanged.

## [2.0.0] - 2026-09-06

Evidence-engine revision. **Changes classifications.** Every
`ps_*` run, the manuscript numbers, Figure panels and Supplementary Table S2
must be regenerated against this engine before being quoted.

### Changed

- **Classification is points-only.** `classify_acmg()` sums Tavtigian 2020
  points over one module-level table, `ACMG_TAG_PTS`, and bands them
  (P >= 10, LP 6-9, VUS 0-5 with High/Mid/Low sub-tiers, LB -1..-6, B <= -7),
  with BA1 stand-alone. The Richards 2015 rule ladder that ran first is gone:
  it counted tags by prefix (`PS1_supporting + PM2` reached Likely Pathogenic
  at 3 points; `PS1_moderate + PS2` reached Pathogenic at 6) and matched
  pathogenic combinations before looking at benign evidence (`PS1 + PS2 + BA1`
  returned Pathogenic at 0 points). Reduced-strength tags now score the
  strength they were assigned, and conflicting evidence cancels. The `rule`
  field reads `"score N pts"` or `"BA1 (stand-alone)"`.
- **PM2 is Supporting (+1)**, per the ClinGen SVI recommendation of 4 Sep 2020.
  PVS1 + PM2 = 9 still reaches Likely Pathogenic, as that recommendation
  specifies.
- **PP5 and BP6 no longer score** (retired by ClinGen SVI, 2018). They are still
  emitted and displayed at 0 points.
- **No PP3 upgrade from predictor agreement.** AlphaMissense >= 0.90 with
  REVEL >= 0.773 used to promote PP3 to Strong as a "PS3 proxy". Two in-silico
  predictors are one evidence line; PP3 strength now comes only from the
  Pejaver 2022 calibrated thresholds. The convergence is still noted in the
  narrative.
- **PTM location is not experimental evidence.** A phospho / ubiquitin /
  acetyl / disulfide annotation at the residue used to emit `PS3_supporting`;
  it now carries the non-scoring `PP_PTM` marker and stays visible in
  `PTM_Info` / `PTM_Strength`.
- **Mutagenesis PS3 requires the tested substitution to match.** The UniProt
  parser keeps `alternativeSequence` (new `alt_aa` column); `PS3_supporting`
  fires only when the queried alt residue was the one assayed. S50A "Loss of
  activity" no longer scores for S50F. Unmatched or unknown: note shown,
  nothing scored.
- **MDS corroborates, it does not stack.** Where a pathway PM1 and MDS both
  fire, the stronger of the two is taken (`max`), never the sum; Moderate PM1 +
  Moderate MDS stays Moderate. The "release conservation to PP3" step that
  depended on the sum is removed; `varviz.mds_frees_cons` is accepted as a
  no-op.
- **MDS LR tiers are off by default** (`varviz.mds_tiered = FALSE`). The
  submission configuration is Moderate-only MDS; the +3 / +4 tiers remain
  available as an exploratory option (`ps_tiered_harness.R`, or
  `VARVIZ_MDS_TIERED=TRUE` for `ps_final_harness.R`, whose env flags are now
  validated instead of treating any non-"FALSE" string as TRUE).

### Added

- `analyses/tests/test_acmg_points_only.R`: strength-suffix and conflicting-
  evidence cases, SVI PM2 rule, retired criteria, band edges, PP1/PP4 cap.
  `test_uniprot_sites.R` gains the alt-residue gating cases.

## [1.1.3] - 2026-09-06

Pipeline audit fixes (analysis scripts and harnesses only). **No classification
changes** — `classify_acmg()` and `build_variant_table()` are untouched; every
call the app or a harness produces is identical. What changes is whether a
failed run can masquerade as a finished one.

### Fixed

- **Harnesses no longer checkpoint a gene whose data sources failed.** New
  `analyses/lib/harness_guard.R`: `harness_fetch()` wraps the eight gene-level
  fetches (Pfam, UniProt, gnomAD, ClinVar, CCRS, AlphaFold, mean pathogenicity,
  gene info) plus the Ensembl exon and dbNSFP pre-loads, records each failure in
  `options(varviz.harness_fetch_failed)`, and `run_one()` refuses the checkpoint
  when anything new failed. Applied to the shared `05_classify_harness.R` (all
  six `ps_*` wrappers), `ps_final_harness.R` (which previously watched only
  three sentinel families) and `repro/08_casestudy_harness.R`. A run with any
  aborted gene now exits 1 instead of writing a partial `summary.tsv`.
- **`summary.tsv` is never overwritten silently.** The shared harness needs
  `--force`, `ps_final_harness.R` needs `VARVIZ_FORCE=1`. The manuscript numbers
  were read from those files.
- **`repro/run_all.sh` propagates a claims mismatch.** `02_reconcile.R`'s exit
  status was swallowed by `|| true`, so the driver always exited 0.
- **DOLPHIN failures are no longer cached as answers.** `lib/dolphin.R` stops
  caching NA (a timeout / 429 / 500 lived for the rest of the session);
  `lib/dolphin_bulk.R` writes `pm1 = NA` on a failed fetch and retries NA rows
  on the next run instead of recording `FALSE`. Re-probing the existing
  per-gene caches showed GCK, KRAS, NUDT15 and PTEN held thousands of HTTP-500s
  recorded as "no PM1"; those rows are being refetched.
- **Frozen inputs need `--force` to regenerate.** `01_pull_varibench.R`,
  `02_pull_mavedb.R`, `25_clinvar_extract.R` refuse to overwrite
  `varibench_canonical.tsv`, `mavedb_canonical.tsv`,
  `clinvar_missense_2star.tsv` when they exist, and stop on zero rows.
- **Empty joins fail loud.** `22_mds_benchmark.R` stops if the VariBench join
  is empty, under half-matched, or single-class; `24_clinvar_mds_lr.R` prints
  the rows dropped at each join and stops if either P or B arm is empty.
- **`17_merge_mane_into_rdata.R` checks that its backup of `VarViz.RData`
  succeeded** before overwriting the file.
- **AUROC bootstrap CI is reproducible.** `compute_metric_suite()` gains
  `seed = 1L` (scoped; caller's RNG restored). `seed = NULL` keeps the old
  behaviour.

### Added

- Plain-R offline tests: `test_harness_guard.R`, `test_dolphin_bulk_failure.R`,
  `test_dolphin_failure_cache.R`, `test_metric_suite_seed.R`.

## [1.1.2] - 2026-08-30

Dead code and dependency removal. **No classification changes** — nothing in this
release touches the evidence engine, the point ladder, or any export column.

### Removed

- **Five `output$` blocks that no UI ever rendered.** `multicons_status` (57
  lines), `consurf_status` (88), `mplot_container` (9), `debug_cols` (4) and
  `highlight` (2) had no `uiOutput`/`plotOutput` counterpart anywhere in `ui.R`
  or in any `renderUI`, so none of them was reachable. The dynamic-height
  `mplot_container` in particular was superseded by the fixed-height
  `plotlyOutput("mplot")` that `ui.R` actually renders. The `highlight_data`
  reactive fed only the dead `output$highlight` table and went with it.
- **Three package dependencies.** `DT` (its only call site was the dead
  `output$highlight`), `markdown` (loaded in `ui.R`, never called), and
  `shinycssloaders` from `server.R` (its only `withSpinner` was in the dead
  `mplot_container`; `ui.R` still loads it for the one live spinner).
- **`tm`.** The word-cloud text pipeline ran `tolower` → `removeNumbers` →
  stopwords → `removePunctuation` → `stripWhitespace` on text that the
  preceding `gsub("[^[:alpha:][:space:]]", " ", ...)` had already reduced to
  letters and spaces, so two of those steps were no-ops and the contraction
  stopwords (`don't`, `it's`, …) were unreachable. Replaced with base-R
  tokenisation plus the reachable stopwords, keeping `TermDocumentMatrix`'s
  three-character minimum. `analyses/tests/test_wordcloud_tokens.R` asserts the
  new tokeniser reproduces the `tm` term frequencies exactly, and skips when
  `tm` is not installed.

Together this drops four packages from the deployed bundle, which matters
against the shinyapps.io free-tier memory ceiling.

## [1.1.1] - 2026-08-30

Display and export text only. **No classification changes.**

### Fixed

- **A transcript accession was shown as the ClinVar condition.** When ClinVar
  returns no `trait_set`, the fallback took the title text before the gene
  symbol; on `NM_003042.4(SLC6A1):c.914C>T (p.Ala305Val)` its pattern matched
  only the `4(`, so the capture kept `NM_003042.` and that reached the variant
  card and the `ClinVar_Trait` column as a condition name. The fallback now
  rejects accession-shaped results and trims a trailing accession stem, so
  `Noonan syndrome, NM_002834.5(PTPN11):...` yields `Noonan syndrome`.
- **`clean_trait()` was applied to the narrative comment only**, so the raw
  value still reached the card, the evidence card and the export. It now runs
  on the column, and was hoisted to the top level -- it had been local to
  `generate_acmg_comment()`, where the column could not have called it.
- **A same-codon ClinVar match displayed the other variant's call as this
  variant's.** `p.Ala305Thr` showed a red "Pathogenic" from the `p.Ala305Val`
  record while being absent from ClinVar itself; only a small "(same codon)"
  marked the difference. The row now names the source and drops the red:
  `p.Ala305Val: Pathogenic - same codon, not this variant`. Scoring was already
  correct -- PP5 requires an exact match, so the variant received
  `PM5_supporting`, not PP5.

## [1.1.0] - 2026-08-29

### Fixed — wrong answers

- **Coding and genomic HGVS were coerced into protein positions.** `p.` was
  prepended to whatever was typed, so `c.289C>T` scored residue 289 and
  `NM_000388.4:c.2968G>A` scored 388, the transcript version. Now refused
  outright, all-or-nothing: those notations cannot be placed on the sequence.
- **The reference residue was never checked.** `A2ML1 p.I554N` returned a full
  classification of residue 554, which carries V. Inputs are now compared
  against the UniProt canonical sequence; mismatches are skipped individually
  with a persistent warning, and the run is blocked only when none survive.
- **`highlight()` failed when every variant was filtered out** with "arguments
  imply differing number of rows: 0, 1", because the gate read the reactive
  that had already removed them instead of the raw input.
- **Gene-level APIs fired on rejected submissions** — UniProt, UCSC, ClinVar and
  gnomAD were all called before the notation was validated.

### Added

- **ClinGen gene-disease validity restored.** All three per-gene sources had
  failed silently: `ldh.genome.network` 404s for every gene, GenCC's
  `validity-prop` answers `{"message": ""}`, and the eRepo branch queried
  variant interpretations with a filter that endpoint ignores. Now read from the
  published curation table (~3668 curations, ~3029 genes), fetched once per
  session. **Does not change any classification**: of the 58 genes across both
  validation cohorts, the only three that could newly gain PP2 already have it
  from ClinVar counts, and the sole BP1-suppression candidate is already
  Pathogenic without BP1.
- **PM1 derivation on hover**, recorded in the new `PM1_Derivation` column where
  the points are decided — the ladder caps at 4, so a final `PM1_strong` cannot
  be decomposed after the fact.
- **Run-parameter header on exported TSVs**, stamped when Go is pressed rather
  than at download, including `variants_requested` alongside `variants` so a
  skipped entry is visible in the file.

### Changed

- The plot's **PS** row is now **Hot**. It draws permutation-significant ClinVar
  hotspots, and the old label read as the ACMG Strong tier — an empty row looked
  like a contradiction when a variant scored `PM1_strong` by another pathway.
  Renamed in the plot, the track checkbox, the landing page and help.
- Help documents all five rows of the ClinVar/PTM/CCRS/Hot panel and states that
  they are independent PM1 pathways.
- `ClinGen_Disease` and `ClinGen_MOI` no longer export empty (the call resolved
  to no function).
- Runtime scratch files write to `tempdir()` instead of accumulating in the app
  directory.
- The completion toast waits for the gene panel to render.

## [1.0.0] - 2026-08-26

State of the engine at Human Mutation submission, tagged
`pre-human_mutation-submission`. Every number in the manuscript and
supplementary was produced by this version.
