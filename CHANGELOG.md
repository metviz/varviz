# Changelog

VarViz follows [semantic versioning](https://semver.org): MAJOR for changes that
alter classifications, MINOR for new evidence sources or user-facing features,
PATCH for fixes that leave every call unchanged.

Because this tool assigns ACMG classifications, each entry states explicitly
whether it can move a variant's call.

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
