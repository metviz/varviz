# Evidence Strength (calibrated pathogenicity odds) — relevance to VarViz

**Scope.** What Evidence Strength is, a worked cautionary error from a current JCEM MEN1
paper (Bindra et al., 2026, doi:10.1210/clinem/dgag239 — DOI as printed in the source PDF),
how the concept already lives inside the VarViz engine, and where it could add value.
Written for the VarViz codebase (`server.R`), not as a general primer.

> **Terminology.** The canonical published term is **OddsPath** (odds of pathogenicity,
> Tavtigian 2018). Throughout VarViz we use the plain-language label **"Evidence Strength"**
> for the tier and **"odds of pathogenicity"** for the numeric value, citing `(OddsPath;
> Tavtigian 2018)` once so clinical geneticists and reviewers map it to the standard.
> Rationale: "OddsPath" reads like a file path and hides the meaning; "Evidence Strength"
> matches the ACMG tier vocabulary (Supporting / Moderate / Strong) users already see.

---

## 1. What Evidence Strength (OddsPath) actually is

Evidence Strength = the **odds of pathogenicity contributed by one piece of evidence** —
i.e. the likelihood ratio of that evidence, independent of the prior.

```
posterior_odds = prior_odds × evidence_strength
posterior_prob = posterior_odds / (1 + posterior_odds)
```

Evidence Strength is the middle term. It is a property of the *evidence* (how much a
REVEL≥0.7 call, a functional assay, a segregation result shifts the odds), **not** of the
patient's prior. That decoupling is the whole point of the ACMG/AMP Bayesian reformulation
(Tavtigian 2018): separate evidence strength from prior so one tool calibration is reusable
across genes and patients.

Tavtigian's strength tiers are just fixed Evidence-Strength cut-points at prior_P = 0.10,
with points = log_C(strength), C ≈ 2.08:

| Strength tier | Points | Odds of pathogenicity (≈) |
|---|---|---|
| Supporting | 1 | 2.08 |
| Moderate | 2 | 4.33 |
| Strong | 4 | 18.7 |
| Very strong | 8 | 350 |

A total point score is therefore a **log of the combined odds of pathogenicity**, and the
final posterior probability is one formula away. VarViz already computes that point sum
(`server.R:3476–3560`).

---

## 2. The JCEM paper's error — cautionary, and directly instructive

Bindra et al. dichotomized REVEL at ≥0.70, measured **LR+ = 11.83** against ClinVar, then
wrote:

> "Using a fixed Prior_P of 0.1 (prior odds = 0.11), the resulting OddsPath ... was
> approximately 1.30 and was therefore applied as supporting (PP3_Supporting)."

That number is `prior_odds × LR+ = 0.11 × 11.83 ≈ 1.30`. But that product is the
**posterior odds**, not the evidence strength. The Evidence Strength *is* the LR ≈ **11.83**,
prior-independent.

Consequences:
- **11.83 sits between their own quoted tiers** (moderate ≥4.33, strong ≥18.7) → REVEL≥0.7
  should be **Moderate**, not Supporting. They understated their own evidence.
- Their next paragraph proves the bug: assuming prior 0.5 they "get OddsPath 11.83 =
  moderate." Evidence Strength cannot depend on the prior — the fact that it moved when
  they changed the prior *is* the tell that they computed posterior odds and mislabeled it.

**The one-line guard for VarViz:** when reporting *evidence strength*, never multiply by the
prior. Prior enters exactly once, at the final posterior step. Two separate outputs:
- **Evidence Strength** = odds of pathogenicity (drives the tier) — prior-free
- **variant posterior probability** = prior_odds × ∏(evidence strengths) → prob — prior used once

VarViz's Tavtigian pathway is already structured this way (points accumulate prior-free,
threshold applied at the end), so it does **not** have this bug. Worth keeping that
separation explicit in any new calibration code.

---

## 3. Where Evidence Strength already lives in VarViz

VarViz does not currently expose an Evidence Strength number, but it is implicitly baked in:

- **Pejaver 2022 PP3/BP4 ladder** (`server.R:4882–4909`): the REVEL cut-points
  0.644 / 0.773 / 0.932 for supporting / moderate / strong are Pejaver's *local-posterior*
  Evidence-Strength calibrations. Every threshold in that block is a strength boundary that
  someone else pre-computed on a genome-wide set.
- **Tavtigian point sum** (`server.R:3476–3560`): a running log-strength accumulator. The
  final tier (≥10 P, 6–9 LP, …) is a strength→posterior mapping.

So the machinery is present; what is missing is (a) surfacing the calibrated posterior as a
number, and (b) letting the calibration be *gene-specific* rather than genome-global.

---

## 4. Opportunities for VarViz (ranked by value / effort)

### 4.1 Report a calibrated posterior probability, not just a tier — *high value, low effort*
The Tavtigian point sum already computed is `log_2.08(combined odds of pathogenicity)`. One
function turns it into a probability the user can read:

```r
# points -> calibrated posterior probability of pathogenicity
acmg_posterior <- function(points, prior_p = 0.10, C = 2.0813) {
  prior_odds <- prior_p / (1 - prior_p)
  post_odds  <- prior_odds * C^points     # C^points == combined evidence strength
  post_odds / (1 + post_odds)
}
```

Show it on each variant card next to the tier: *"LP — posterior P(path) ≈ 0.94."* More
honest than a bare tier and makes the prior assumption visible. Ties directly to the
prevalence/penetrance prior already collected in the parameters panel.

### 4.2 Gene-specific PP3 recalibration — *high value, medium effort, the differentiator*
This is exactly what the JCEM paper did by hand for MEN1, and VarViz is already positioned to
do it **live**: it pulls per-gene ClinVar LP/P and gnomAD/benign sets for every query. So for
the queried gene it can:
1. Extract REVEL (and AM) for ClinVar LP/P vs LB/B in that gene.
2. Compute LR+ / Evidence Strength at candidate thresholds.
3. Report the gene-specific Evidence-Strength tier for the threshold used, alongside the
   genome-wide Pejaver default.

Payoff: flags when the global Pejaver threshold is mis-calibrated for a specific gene (the
whole premise that MEN1 needs its own VCEP thresholds). Even a "gene-level calibration
confidence" badge (n benign available, LR+ CI) would be novel in a browser tool.

**Guardrail — copy from the JCEM mistakes:**
- Small benign N makes LR+ unstable (they had 15 LB/B). Require a minimum N (e.g. ≥10 per
  arm) before showing a gene-specific tier; otherwise fall back to Pejaver global and say so.
- Do **not** let the ClinVar-LP set that trained the threshold also be the set you score —
  circularity. If VarViz computes gene-specific thresholds *and* the user's variants are
  ClinVar LP/P, exclude-self or cross-validate.
- Report the LR+ as Evidence Strength, prior-free (§2).

### 4.3 Continuous per-predictor Evidence Strength instead of hard tiers — *medium value, medium effort*
Pejaver publishes continuous local-posterior curves, not just 3 cut-points. VarViz could map
a raw REVEL/AM score to its continuous Evidence Strength and show "REVEL 0.81 → strength 6.2×
→ Moderate" on hover. More information than snapping to a tier, and it exposes borderline
scores honestly (0.775 vs 0.771 currently flip moderate/supporting at a cliff).

### 4.4 Make the prior explicit and gene-aware — *low effort, conceptual tidy*
VarViz already derives a prevalence-based max-credible-AF (Whiffin) for PM2. The same
prevalence/penetrance/heterogeneity inputs define a defensible Prior_P for §4.1 rather than
the default 0.10. Surfacing "prior used = X" closes the loop and pre-empts the JCEM confusion
by showing the prior only enters at the posterior.

### 4.5 Evidence Strength for the clinical dropdowns (PP1/PS3/PM3) — *lower priority*
The per-variant segregation/functional/comphet dropdowns map to fixed strength tiers. Those
tiers are themselves Evidence-Strength values; if VarViz shows the running posterior (§4.1),
each dropdown change visibly moves the probability — a strong teaching/QC feature and a
sanity check against over-counting (the PP4/PS4 double-dip problem the JCEM paper had).

---

## 5. Minimum concrete next step

Implement §4.1 (posterior-probability readout) first — a pure function over a number VarViz
already has, no new API call, and it immediately makes the existing Bayesian engine legible.
§4.2 (live gene-specific Evidence Strength from the ClinVar/gnomAD data already in memory) is
the genuinely novel, publishable extension and the natural VarViz differentiator over
static-threshold tools.

**Invariant to enforce in any of this:** Evidence Strength (OddsPath) is prior-free; prior_P
enters exactly once, at the posterior. That single rule is what the JCEM paper got wrong.
