#!/usr/bin/env python3
"""Calibrate AlphaMissense on this benchmark's own ClinVar labels.

Adopting the Bergquist et al. 2025 intervals wholesale moves 19.11% of Pass-Full
calls, almost all upward, and nearly all of that lands on variants ClinVar has
never reviewed (+10.6 points of actionable share) rather than on ones it has
(+0.1). That pattern is consistent with thresholds calibrated at one prior being
applied at a much lower one, but the strata split alone cannot demonstrate it.

This measures it directly. For each score threshold, the positive likelihood
ratio

    LR+ = P(score >= t | pathogenic) / P(score >= t | benign)

is computed on the ClinVar-labeled variants of the 14 benchmark genes, and
compared against the evidence-strength odds of Tavtigian et al. 2020: 2.08 for
supporting, 4.33 for moderate, 18.7 for strong. If the published thresholds
transfer to this universe, the LR+ at each should land near its intended
strength. If they were calibrated on a set with a higher prior, the LR+ here
will fall short of it.

Scores come from the AlphaMissense bulk release, not from the classification
runs, so this is independent of anything the engine did.

  python3 analyses/am_local_calibration.py
"""
import csv, math, os, sys

AM   = os.environ.get("AM_TSV", "")
CV   = "analyses/tmp/clinvar_14genes_all.tsv"
OUT  = "analyses/humu/build_v230"
BERG = [("supporting", 0.792), ("moderate", 0.906), ("3-point", 0.972), ("strong", 0.990)]
DEV  = [("developer default", 0.564)]
TAVT = [("supporting", 2.08), ("moderate", 4.33), ("strong", 18.7), ("very strong", 350.0)]


def wilson(k, n, z=1.96):
    if n == 0: return (float("nan"),)*2
    p = k/n; d = 1 + z*z/n
    c = (p + z*z/(2*n))/d
    h = z*math.sqrt(p*(1-p)/n + z*z/(4*n*n))/d
    return (max(0.0, c-h), min(1.0, c+h))


def main():
    if not AM or not os.path.isfile(AM):
        sys.exit("set AM_TSV to the extracted gene/variant/score table")
    am = {}
    for line in open(AM, encoding="utf8"):
        p = line.rstrip("\n").split("\t")
        if len(p) < 3: continue
        try: am[(p[0], p[1])] = float(p[2])
        except ValueError: pass

    # Best-reviewed ClinVar assertion per variant, restricted to unambiguous calls.
    lab = {}
    for r in csv.DictReader(open(CV, newline="", encoding="utf8"), delimiter="\t"):
        s = r["sig"].lower()
        if "conflicting" in s or "uncertain" in s: continue
        l = 1 if ("pathogenic" in s and "benign" not in s) else (0 if "benign" in s and "pathogenic" not in s else None)
        if l is None: continue
        st = 0 if r["review"].startswith("no ") else (3 if "expert" in r["review"] else (2 if "multiple" in r["review"] else 1))
        k = (r["gene"], r["p_notation"].replace("p.", ""))
        if k not in lab or st > lab[k][1]: lab[k] = (l, st)

    for min_star in (1, 2):
        pos = [am[k] for k, v in lab.items() if k in am and v[0] == 1 and v[1] >= min_star]
        neg = [am[k] for k, v in lab.items() if k in am and v[0] == 0 and v[1] >= min_star]
        print(f"\n{'='*82}")
        print(f"AlphaMissense LR+ on {len(pos):,} pathogenic and {len(neg):,} benign "
              f"ClinVar >={min_star}-star variants, 14 benchmark genes")
        print(f"{'='*82}")
        if len(neg) < 20:
            print("  too few benign labels to estimate LR+"); continue
        print(f"{'threshold':<26}{'sens':>8}{'1-spec':>9}{'LR+':>9}{'95% CI on LR+':>22}{'implies':>14}")
        for name, t in DEV + BERG:
            tp = sum(1 for s in pos if s >= t); fp = sum(1 for s in neg if s >= t)
            sens = tp/len(pos); fpr = fp/len(neg)
            lr = sens/fpr if fpr > 0 else float("inf")
            # CI via the Wilson bounds of each rate, taken in the worst/best pairing.
            sl, sh = wilson(tp, len(pos)); fl, fh = wilson(fp, len(neg))
            lo = sl/fh if fh > 0 else float("inf"); hi = sh/fl if fl > 0 else float("inf")
            implies = "below supporting"
            for nm, thr in TAVT:
                if lr >= thr: implies = nm
            ci = f"[{lo:.1f}, {hi:.1f}]" if hi != float("inf") else f"[{lo:.1f}, inf)"
            lrs = f"{lr:.1f}" if lr != float("inf") else "inf"
            print(f"{name+' (>='+str(t)+')':<26}{sens:>8.3f}{fpr:>9.3f}{lrs:>9}{ci:>22}{implies:>14}")

        # Where would each strength threshold actually fall on this cohort?
        print(f"\n  empirical thresholds on this cohort (lowest score reaching each LR+):")
        xs = sorted({round(s, 3) for s in pos + neg})
        for nm, need in TAVT:
            hit = None
            for t in xs:
                fp = sum(1 for s in neg if s >= t)
                if fp == 0: continue
                lr = (sum(1 for s in pos if s >= t)/len(pos))/(fp/len(neg))
                if lr >= need: hit = t; break
            print(f"    {nm:<12} LR+ >= {need:<6} at AlphaMissense >= "
                  + (f"{hit:.3f}" if hit is not None else "not reached below the point where no benign variant scores higher"))

    os.makedirs(f"{OUT}/numbers", exist_ok=True)
    print(f"\n(scores from the AlphaMissense bulk release; labels from {CV})")


if __name__ == "__main__":
    main()
