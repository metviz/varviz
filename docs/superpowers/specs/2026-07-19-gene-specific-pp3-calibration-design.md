# Gene-specific PP3 calibration (Evidence Strength §4.2) — Design

**Status:** Approved design, pre-implementation.
**Source:** `evidence_strength_review.md` §4.2 ("Gene-specific PP3 recalibration — the differentiator").
**Depends on:** §4.1 posterior readout (shipped, `acmg_posterior`, commit c49bcf7).

---

## 1. Purpose

VarViz applies genome-wide Pejaver 2022 REVEL/CADD/AM cut-points to every gene. Those
thresholds are mis-calibrated for some genes (the premise behind gene-specific VCEP
thresholds, e.g. MEN1). VarViz already pulls per-gene ClinVar and dbNSFP data live, so it
can compute the gene's *own* Evidence Strength for each predictor threshold at query time
and show where the global default disagrees — a capability static-threshold tools lack.

**Evidence Strength = LR+ = OddsPath, prior-free.** The prior enters exactly once, at the
posterior (§4.1). This spec never multiplies LR+ by a prior.

## 2. Goals / non-goals

**Goals**
- For the queried gene, compute per-predictor (REVEL, AM, CADD) LR+ from its ClinVar
  P/LP vs LB/B variants.
- Two outputs per predictor: (a) validate each fixed Pejaver cut-point, (b) report a
  gene-optimal threshold.
- Report as an advisory badge; allow gated opt-in override of the PP3 tier.
- Never double-count the queried variant's own ClinVar assertion against its calibration.

**Non-goals (v1)**
- No auto-override of the ACMG call.
- No persistence of overrides beyond the session.
- No calibration of predictors outside REVEL/AM/CADD.
- No local dbNSFP file dependency (hosting can't carry the 73 GB set).

## 3. Data flow (per gene, cached)

New cache namespace `api_cache$gene_calib`, keyed by gene name.

1. **P/LP arm** — existing `extract_clinvar(gene)` (server.R:824), already P/LP-only.
2. **Benign arm** — new sibling esearch: clone the query at server.R:835 with
   `("benign"[clinsig] OR "likely benign"[clinsig])`, same esummary parse. One extra
   esearch + esummary round trip.
3. **Scores for both arms** — extract each variant's protein change by parsing
   `p.Arg130Gly` out of the ClinVar `name` string, normalize 3-letter → 1-letter (reuse
   the existing `variants_to_try` hgvsp normalization used by `fetch_dbnsfp`), then issue
   **one MyVariant.info batch POST**:
   `POST https://myvariant.info/v1/query`, body `q=[hgvsp...]`,
   `scopes=dbnsfp.hgvsp`, `fields=dbnsfp.revel,dbnsfp.alphamissense,cadd.phred`.
   A gene's combined arms are < 1000 variants → a single POST. Replaces the per-variant
   GET loop (server.R:1944, `Sys.sleep(0.2)` each) for this bulk path.
4. Cache the two labelled, scored arms. LR+ and the gene-optimal sweep recompute cheaply
   from the cached arms when the min-sensitivity param changes — no re-fetch.

## 4. The statistic

For a predictor and a threshold `t`, positives = P/LP arm, negatives = LB/B arm:

- `sens = P(score >= t | P/LP)`, `spec = P(score < t | LB/B)`
- `LR+ = sens / (1 - spec)`
- **Zero-cell guard:** Haldane-Anscombe +0.5 continuity correction on the 2x2 cells so a
  clean-separation threshold yields a finite LR+, not infinity.
- **95% CI (Katz log method):**
  `Var(ln LR+) = (1 - sens)/(sens * n_pos) + spec/((1 - spec) * n_neg)`,
  `CI = exp(ln LR+ +/- 1.96 * sqrt(Var))`.
- **LR+ -> Evidence-Strength tier:** map via the Tavtigian OddsPath boundaries reused from
  §4.1 — Supporting >= 2.0813, Moderate >= 4.33, Strong >= 18.7, Very strong >= 350
  (C = 2.0813). LR+ *is* the OddsPath; no prior involved.

## 5. Two outputs, per predictor

Predictor Pejaver cut-points as encoded in VarViz today:
- REVEL: 0.644 (Supporting) / 0.773 (Moderate) / 0.932 (Strong)
- CADD: 28.1 (Supporting) / 35 (Moderate)
- AM: 0.564 (Supporting)   *(VarViz encodes only the supporting cut-point for AM)*

1. **Validation table** — at each fixed Pejaver cut-point: gene LR+ (95% CI), gene tier,
   and agree / disagree vs the genome-wide Pejaver label.
2. **Gene-optimal** — threshold that maximizes LR+ subject to `sens >= floor`, where
   `floor` is a user parameter (default 0.90) exposed in the parameters panel (same UX as
   prevalence/penetrance). Report threshold + LR+ (CI) + tier. Sweep over the observed
   score grid for the gene (candidate thresholds = the distinct scores present), not a
   fixed step — exact and cheap at these N.

## 6. Circularity guard + the variant's own ClinVar evidence (kept separate)

Two independent roles of the queried variant's ClinVar record, never merged:

- **Calibration (leave-one-out):** if the queried variant is present in ClinVar, drop it
  from its arm before computing LR+. Prevents the variant from helping define the
  threshold that then scores it. Badge notes when LOO fires. On small-N genes, dropping
  one variant visibly widens the CI — correct.
- **Independent assertion line:** show the queried variant's own ClinVar record —
  significance + review stars (`clinvar_goldstar`, parsed at server.R:863) + submitted
  criteria when available — as a distinct "This variant in ClinVar" row. It informs
  confidence as ClinVar-assertion evidence, tracked separately from the gene-specific PP3
  LR+, never folded in.
- **No double-count:** the two lines are labelled as independent evidence streams. This is
  the PP4/PS4 double-dip Bindra et al. committed; the UI states the separation explicitly.

## 7. Gating / confidence

- Compute only if both arms have >= 1 variant (Haldane covers per-cell zeros at a
  threshold). Below that, show "insufficient gene data", fall back to Pejaver.
- **N >= 20 per arm** -> high confidence. Override toggle enabled normally.
- **N < 20 per arm** -> still show LR+ + CI + N, flagged "low confidence (n = X/Y)".
  Override available but requires explicit user confirmation.
- CI is always displayed regardless of confidence tier.

## 8. Authority & UI

- **Default: advisory.** New collapsible "Gene-specific calibration" section on the
  variant card. Never mutates `acmg_res` or the PP3 tag on its own.
- **Opt-in override toggle.** When flipped (and gating per §7 satisfied), applies the
  gene-specific tier to the PP3 tag for the current session only. High confidence =
  one-click; low confidence = explicit confirm.
- **Min-sensitivity input** added to the parameters panel, default 0.90.
- Fall-back and LOO states are surfaced in-line, not hidden.

## 9. Files touched

- `server.R` — benign-arm esearch (clone ~835); `fetch_dbnsfp_batch()` (new, MyVariant.info
  POST); hgvsp parse/normalize + join; `gene_lr_plus()` / CI / tier mapping;
  gene-optimal sweep; `api_cache$gene_calib` namespace; card section render; override wiring.
- `ui.R` — min-sensitivity parameter input.
- This is an app feature; server.R/ui.R edits are in scope (not a Pass-2 analysis script).
- Back up server.R/ui.R to timestamped `.bak` before edits (project rule).

## 10. Risks

- **Small N** — the dominant failure mode (Bindra had 15 LB/B). Mitigated by the N>=20
  confidence gate, always-on CI, and LOO.
- **hgvsp format mismatch** ClinVar `name` vs dbNSFP hgvsp — mitigated by reusing the
  existing normalizer; unmatched variants are dropped from the arm and counted (report the
  match rate so silent loss is visible).
- **MyVariant.info coverage** — variants with no dbNSFP hit are excluded from that
  predictor's arm; report per-predictor N so a sparse predictor is obvious.
- **Circularity** — handled by LOO (§6).
- **Double-count** — handled by separating the assertion line from calibration (§6).

## 11. Testing

- Unit (standalone `test_gene_calib.R`, pattern of `test_acmg_posterior.R`):
  - LR+ on a hand-built 2x2 matches `sens/(1-spec)`.
  - Haldane correction: clean-separation arm -> finite LR+, not Inf.
  - Katz CI: known arm -> CI matches closed-form.
  - LOO: removing the query variant changes LR+ in the expected direction; large-N ->
    negligible, small-N -> visible.
  - Tier mapping: LR+ 11.83 -> Moderate; 2.5 -> Supporting; 20 -> Strong (reuses §4.1 C).
  - Gene-optimal: on a synthetic gene, returns the max-LR+ threshold at/above the floor.
- Live smoke test against a well-populated gene (both arms >= 20) and a sparse gene
  (confirm low-confidence fallback + explicit-confirm override path).

## 12. Open items

None blocking. AM currently contributes only a supporting-tier validation (single
cut-point); its gene-optimal output is still meaningful. Revisit adding AM moderate/strong
cut-points if Pejaver-style AM calibration is later encoded in the PP3 ladder.
