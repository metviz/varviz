# Option C — handoff note for follow-up session

## What's done

**Downloads complete and on disk** (gitignored under `analyses/raw/`):

| Path | Size | Source | Coverage |
|---|---|---|---|
| `analyses/raw/revel/revel-v1.3_all_chromosomes.zip` | 637 MB (6.5 GB extracted) | `https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip` (REVEL v1.3, May 2021) | All chromosomes; one row per (chrom, pos, ref, alt) with REVEL score and transcript IDs |
| `analyses/raw/alphamissense/{UID}-F1-aa-substitutions.csv` × 11 | 2.5 MB total | `https://alphafold.ebi.ac.uk/files/AF-{UID}-F1-aa-substitutions.csv` | Per-UniProt: every (residue, alt) → AlphaMissense pathogenicity + class |

UniProt accessions covered: P01116 (KRAS), P01130 (LDLR), P35557 (GCK), P41180 (CASR), P51787 (KCNQ1), P60484 (PTEN), Q12809 (KCNH2), Q13362 (PPP2R5C), Q16832 (DDR2), Q86YT5 (SLC13A5), Q9NV35 (NUDT15).

**Missing AM CSV: TSHR (P16473) returned 404** at AlphaFold-EBI (no AM bulk file exists for this protein). TSHR has 4 universe variants — fall back to MyVariant for those.

## What's left (7 tasks; ~3-4h work)

### 1. Extract + index bulk REVEL  ✅ **DONE 2026-04-30**

REVEL bulk has been unzipped, awk-split into 10 per-chromosome CSVs at
`analyses/raw/revel/by_chrom/chr_{N}.csv` (where N ∈ {1, 3, 7, 10, 11, 12,
13, 14, 17, 19} — all 14 BENCHMARK_GENES covered). The 6.5 GB intermediate
was deleted; the 637 MB source zip is retained as backup. 3.3 GB on disk
across 10 files; 43.5M total rows.

Schema (shared across files): `chr, hg19_pos, grch38_pos, ref, alt, aaref, aaalt, REVEL, Ensembl_transcriptid`

Same (chr, pos, ref, alt) appears multiple times — once per Ensembl
transcript context. Use `dplyr::distinct(chr, grch38_pos, ref, alt, REVEL)`
when consolidating, or pick the MANE Select transcript via Ensembl lookup.

Per-chrom row counts:
  chr1: 8.5M | chr3: 4.8M | chr7: 3.8M | chr10: 3.3M | chr11: 4.9M
  chr12: 4.4M | chr13: 1.4M | chr14: 2.7M | chr17: 4.5M | chr19: 5.2M

Recommended load pattern in R:
```r
library(data.table)
revel_chr <- fread("analyses/raw/revel/by_chrom/chr_17.csv")
setkey(revel_chr, grch38_pos, ref, alt)
# instant lookups by binary search; ~150 MB per chrom in memory
```

### 2. Build `analyses/lib/local_predictors.R`  *(~1h)*

Two functions:
```r
load_revel_for_gene(gene_symbol)           # → data.table indexed by (grch38_pos, ref, alt)
                                            # uses chr from gene_data, filters REVEL bulk
lookup_revel(revel_dt, chrom, pos, ref, alt)  # → numeric REVEL score or NA
load_alphamissense_csv(uniprot_id)         # → data.table indexed by p_notation (one-letter)
                                            # already exists as afs_data in harness
lookup_alphamissense(am_dt, p_notation)    # → list(score, class) or NA
```

The `(gene, p_notation)` → `(chrom, pos, ref, alt)` translation:
- Use `aa_to_genomic()` already in `server.R:2036`
- Needs Ensembl exon table from `fetch_ensembl_exons(gene)` (already available)
- One-time per gene: 3 codons × 3 alt nucleotides per AA → 9 candidate genomic substitutions; match by `aaalt` from REVEL row

### 3. Patch `server.R` for UCSC conservation → PM1_strong  *(~10 lines)*

Around `server.R:4519-4528`, when `dbnsfp$conservation$PhyloP_100V` is NA, fall back to UCSC values from the existing `consurf_data` / per-base track fetch. Specifically:

```r
phylop_v <- if (dbnsfp$has_data) dbnsfp$conservation$PhyloP_100V else NA_real_
if (is.na(phylop_v) && !is.null(consurf_data$ucsc_phylop100)) {
  # consurf_data here is misnamed — it actually carries UCSC tracks
  phylop_v <- consurf_data$ucsc_phylop100[consurf_data$prot_pos == pos]
}
# Same pattern for phastcons_v (PhastCons_100V) and gerp_v (GERP_RS)
```

Verify the exact attribute path of UCSC PhyloP/PhastCons in `consurf_data` first via `fetch_ucsc_window` return shape.

### 4. Patch `analyses/05_classify_harness.R` to monkey-patch `fetch_dbnsfp`

After `source("server.R")`:
```r
fetch_dbnsfp_remote <- fetch_dbnsfp  # save original

fetch_dbnsfp <<- function(gene_name, hgvsp) {
  revel_dt <- gene_revel_caches[[gene_name]]   # pre-loaded per gene
  am_dt    <- gene_am_caches[[gene_name]]
  # Try local lookup
  local <- list(
    has_data = FALSE,
    ensemble = list(REVEL = list(score = lookup_revel(revel_dt, ...))),
    pathogenicity = list(AlphaMissense = lookup_alphamissense(am_dt, hgvsp)),
    conservation = list(PhyloP_100V = NA, PhastCons_100V = NA, GERP_RS = NA, PhyloP_470M = NA)
  )
  if (!is.na(local$ensemble$REVEL$score) || !is.na(local$pathogenicity$AlphaMissense$score)) {
    local$has_data <- TRUE
    return(local)
  }
  # Fall back
  fetch_dbnsfp_remote(gene_name, hgvsp)
}
```

Pre-load `gene_revel_caches` and `gene_am_caches` once per gene at start of `classify_gene()`. Avoids per-variant file I/O.

### 5. Smoke test on LDLR (299 vars; UID P01130)
- Delete `analyses/classifications/LDLR__dual.tsv` (the 0-row stub)
- Run harness; expect ~3-5 min for LDLR with local lookups vs 30+ min via MyVariant
- Spot-check 5 LDLR variants: REVEL score should match MyVariant reference; classification should be similar

### 6. Production run on remaining 10 genes
Order to delete stubs (re-enable harness processing):
```
CASR DDR2 GCK KCNH2 KCNQ1 KRAS NUDT15 PTEN SLC13A5
```
Skip TSHR (no AM CSV); it'll re-run from scratch via MyVariant for its 4 vars. Total: ~3-4h serial.

### 7. Regenerate panels + update writeup numbers
- `Rscript analyses/09_run_all_panels.R`
- `Rscript analyses/10_compose_supp_s2.R` → updated `analyses/reports/supp_s2_skeleton.md`
- Manually edit `analyses/reports/supp_s2.md` with new full-universe numbers
- Re-inject §S2 into Supplementary.docx (use `analyses/manuscript_inject_supp_s2.py`)
- Update cover-letter Comment 3(i) numbers
- Update §3 closing sentence in v2.docx with new ΔAUROC

### 8. Commit + final verification

## Why TSHR is fine to leave on MyVariant
- 4 universe variants total (per Phase 1.4 gnomAD annotation summary)
- 4 × 1.5s API time = 6s additional runtime
- Won't affect aggregate Pass-Blind AUROC or CAPS ρ
- Documented limitation in cover letter (P16473 lacks AlphaMissense bulk file)

## Quick start when resuming

```bash
cd /home/raghu/tools/varviz
ls analyses/raw/revel/   # confirm zip still there
ls analyses/raw/alphamissense/   # 11 CSVs

# Then start at task #53 above (build lookup helper)
```
