# Changelog

VarViz follows [semantic versioning](https://semver.org): MAJOR for changes that
alter classifications, MINOR for new evidence sources or user-facing features,
PATCH for fixes that leave every call unchanged.

Because this tool assigns ACMG classifications, each entry states explicitly
whether it can move a variant's call.

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
