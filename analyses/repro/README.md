# Reproducibility harness for the VarViz manuscript benchmark

Regenerates every benchmark statistic quoted in the manuscript and supplement
directly from the derived classification runs, and reports which numbers
reconcile and which do not.

Base R only. No package dependencies. Run everything from the repository root.

```sh
Rscript analyses/repro/01_recompute.R    # recompute from the default run
Rscript analyses/repro/02_reconcile.R    # compare against the manuscript's claims
Rscript analyses/repro/03_which_run.R    # score every run on disk against every claim
```

## Files

| File | Purpose |
|------|---------|
| `00_claims.tsv` | Every benchmark number the manuscript reports, with its location. Machine-readable so claims are data, not prose. |
| `01_recompute.R` | Recomputes the statistics from a classification run. Takes the run path as `argv[1]` (default: `analyses/derived/varviz_classifications_mds.tsv`) and an output path as `argv[2]`. |
| `02_reconcile.R` | Joins claims to recomputed values, flags disagreements. Exits non-zero if any claim fails. |
| `03_which_run.R` | Scores all candidate runs against all claims. Use this to identify the authoritative run. |
| `out/` | Generated. Not intended to be edited by hand. |

## Status as of 20 Aug 2026

`analyses/derived/varviz_classifications_mds.tsv` is the **authoritative run**.
Every benchmark number in the manuscript, supplement and cover letter has been
regenerated from it; `02_reconcile.R` reports 22 / 22 claims reconciling.

Two consequences of that choice are recorded in the paper:

1. **The run is untiered.** Every MDS-augmented row receives PM1 at Moderate
   (+2); the point delta is 0 or 2 in every MDS band, including the <= -12 tail.
   The likelihood-ratio tiers are implemented in `server.R` but are not
   exercised by this benchmark. `analyses/ps_tiered_harness.R` exists and has
   completed 1 of 14 genes; running it to completion is what would let the
   tiering claim be benchmarked.
2. **Three terms are now defined separately** because they differ threefold:
   *fires* (MDS <= -4 at a scorable residue, 18,887 rows / 20.8%), *augments*
   (fires and no prior PM1 tag, 6,492 / 7.2%), *resolves* (augments and leaves
   the VUS band, 4,165 / 4.6%). `01_recompute.R` emits all three for both MDS
   and DOLPHIN.

## Coverage figures

MDS is scorable at 59,136 / 90,701 residues (65.2%) within the 14 benchmark
genes. The "~45% of residues" figure in Supp S4.6 is the proteome-wide Pfam-A
rate; both now appear in the Discussion.

## Panels B and C

`04_panel_a.R` regenerates Panel A from the authoritative run using
`analyses/lib/metric_suite.R` unchanged. Panels B (MaveDB) and C (CAPS) have
not yet been regenerated the same way -- `analyses/07_panel_b_dms.R` and
`analyses/08_panel_c_caps.R` still read `varviz_classifications.tsv`, which
shares the Pass-Full and Pass-Blind columns with the authoritative run but
carries no MDS arm. The MDS-augmented CAPS figure (rho 0.847 -> 0.991) is the
one number in the paper that this harness has not yet verified.

## RESOLVED 2026-08-23 — authoritative runs

**Authoritative:** `analyses/ps_final/summary.tsv` (MDS on) and
`analyses/ps_nomds_v2/summary.tsv` (MDS-off ablation). Both are 90,701 rows /
54,454 distinct variants / 14 genes and both pass
`10_validate_run.py` intrinsic checks. All 24 claims in `00_claims.tsv`
reconcile against them.

**Superseded, do not use:** `ps_mds_corroborate`, `ps_mds_consfree`, and
`ps_nomds` (Jul 30, no `run.log`, unauditable). The first sharded `ps_final`
attempt is also superseded; its artifacts are quarantined under
`analyses/ps_final/_bogus/`.

### What was actually wrong

Four fetches returned degraded data instead of raising. Every one left the row
count correct, so row-count checks passed while the numbers were wrong. The
earlier diagnosis in this file blamed a single missing timeout; that was one
symptom of the first vector, not the cause.

| # | fetch | failure | cost |
|---|-------|---------|------|
| 1 | UniProt `gene_info` | `req_retry(max_tries=3)` does not retry connection errors without `retry_on_failure=TRUE`, so the retry never once fired. KCNH2 (Q12809, 245 KB) also needs ~130 s against a 30 s timeout, so it could never succeed | KCNQ1, TP53, KCNH2, KRAS, CASR, DDR2 across runs |
| 2 | UniProt Pfam | failed with a different message carrying no accession, so one grep could not catch both vectors | 3 genes in `ps_mds_consfree` (HTTP 500) |
| 3 | Ensembl exons | a failed lookup leaves `chrom = NA`, which silently disables the dbNSFP *and* REVEL preloads; conservation then reads `GERP=NA` and flips `cons_strong` | BRCA1: 1,326 wrong classifications, no error line anywhere |
| 4 | ClinVar esummary | a throttled batch has no `$result` and was skipped by a bare `next` | LDLR extracted 203 of 1,731 records from 4,734 IDs |

Vectors 3 and 4 were provoked by running 10 gene shards concurrently: NCBI
eutils allows 3 req/s unauthenticated and `rest.ensembl.org` began timing out
under the same load. **Run one gene at a time.** tabix on the local dbNSFP is
concurrency-safe (verified: 14,286 rows returned identically by 10 simultaneous
readers); the external APIs are not.

### What now prevents it

- All four record the affected gene in a sentinel option
  (`varviz.uniprot_failed`, `varviz.localpred_failed`,
  `varviz.clinvar_batch_failed`); `run_one()` in the harness refuses to write a
  checkpoint for any gene that trips one, so a corrupt gene is re-run rather
  than cached.
- Timeouts raised to 180 s; `retry_on_failure = TRUE` on every UniProt fetch;
  ClinVar batches retry 4x with backoff before being counted as lost.
- Out-of-band caches (`analyses/cache/uniprot`, `analyses/cache/ensembl`) take
  the flaky remotes off the critical path. `Q12809.json` is primed there because
  that payload cannot be fetched inside any reasonable timeout.
- `10_validate_run.py` rewritten: intrinsic checks decide pass/fail (gene count
  vs universe, checkpoint count, all four failure strings across *all* logs,
  `dbNSFP=NA` local caches, ClinVar extracted/found ratio below 8%); reference
  comparison is advisory only, because a reference run is not ground truth.
- `analyses/tests/test_uniprot_sentinel.R` asserts both UniProt fetches record
  into the sentinel, by pointing them at a dead host.

### Reproducing

```
Rscript analyses/ps_final_harness.R                     # one gene at a time via VARVIZ_GENES
VARVIZ_OUT_DIR=analyses/ps_nomds_v2 VARVIZ_MDS_PM1=FALSE Rscript analyses/ps_final_harness.R
python3 analyses/repro/10_validate_run.py analyses/ps_final
Rscript analyses/repro/06_panel_a_runs.R                # must run before 07
Rscript analyses/repro/07_regenerate.R
Rscript analyses/repro/02_reconcile.R                   # expect 24/24
```

`06_panel_a_runs.R` must precede `07_regenerate.R`: the Panel A metrics
(AUROC / MCC / specificity / VUS rate) are read from its output. Those metrics
were absent from `00_claims.tsv` until 2026-08-23, which is how Table 1 and
Table S2 carried Pass-Full AUROC 0.985 / MCC 0.473 from a superseded run.
