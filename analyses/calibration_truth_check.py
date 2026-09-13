#!/usr/bin/env python3
"""Does the Bergquist AlphaMissense calibration improve discrimination, or just shift calls up?

The calibration moves 19.11% of Pass-Full calls, almost all upward. Two readings
fit that: the thresholds are right and the engine was under-calling, or they are
calibrated at a prior this universe does not have and over-call here.

The two readings differ in a measurable way. If the calibration finds true
pathogenic variants it was missing, sensitivity rises faster than the false
positive rate and MCC improves. If it is simply shifting mass upward, both rise
together and MCC falls.

Pass-Blind is the pass to measure on: it withholds every ClinVar-derived
criterion, so ClinVar's own labels are an independent truth set rather than an
input. Pass-Full would be circular.

The prior hypothesis is tested separately, by splitting the universe into the
variants ClinVar has reviewed (a prior resembling the calibration set) and those
it has never seen (a far lower prior), and comparing how far each moves.

  python3 analyses/calibration_truth_check.py
"""
import collections, csv, math, sys

RUNS = [("release", "analyses/ps_final_v221"),
        ("calibrated (no override)", "analyses/ps_amdefer_v221"),
        ("calibrated + override", "analyses/ps_amcalib_v221"),
        ("hybrid", "analyses/ps_amhybrid_v221")]
CLINVAR = "analyses/tmp/clinvar_14genes_all.tsv"
ACTIONABLE = {"Pathogenic", "Likely Pathogenic"}

STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, single submitter": 1,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classification provided": 0,
}


def stars(review):
    r = (review or "").strip().lower()
    return STARS.get(r, 1 if r.startswith("criteria provided") else 0)


def truth(sig):
    """Binary label, or None when ClinVar does not take a position."""
    s = (sig or "").lower()
    if "conflicting" in s or "uncertain" in s:
        return None
    path = "pathogenic" in s
    ben = "benign" in s
    if path and not ben:
        return 1
    if ben and not path:
        return 0
    return None


def load(path):
    d = {}
    with open(path, newline="", encoding="utf8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            d.setdefault((r["gene"], r["p_notation"]), r)
    return d


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return ((tp * tn - fp * fn) / den) if den else float("nan")


def main():
    runs = {n: load(f"{d}/summary.tsv") for n, d in RUNS}
    cv = {}
    with open(CLINVAR, newline="", encoding="utf8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r["gene"], r["p_notation"])
            st = stars(r["review"])
            # Keep the best-reviewed assertion per variant.
            if k not in cv or st > cv[k][1]:
                cv[k] = (truth(r["sig"]), st, r["sig"])

    universe = set(runs["release"])
    for min_star in (1, 2):
        labeled = {k: v for k, v in cv.items()
                   if k in universe and v[0] is not None and v[1] >= min_star}
        npos = sum(1 for v in labeled.values() if v[0] == 1)
        nneg = len(labeled) - npos
        print(f"\n{'='*78}")
        print(f"Pass-Blind discrimination against ClinVar >={min_star}-star labels "
              f"({npos:,} pathogenic, {nneg:,} benign)")
        print(f"{'='*78}")
        if nneg < 20:
            print("  too few benign labels at this star level to score specificity")
            continue
        print(f"{'run':<26}{'sens':>8}{'spec':>8}{'FPR':>8}{'MCC':>8}"
              f"{'TP':>7}{'FP':>7}{'TN':>7}{'FN':>7}")
        for name, _ in RUNS:
            r = runs[name]
            tp = fp = tn = fn = 0
            for k, (lab, _st, _s) in labeled.items():
                pred = r[k]["varviz_classification_blind"] in ACTIONABLE
                if lab == 1:
                    tp += pred; fn += not pred
                else:
                    fp += pred; tn += not pred
            sens = tp / (tp + fn) if tp + fn else float("nan")
            spec = tn / (tn + fp) if tn + fp else float("nan")
            print(f"{name:<26}{sens:>8.3f}{spec:>8.3f}{1-spec:>8.3f}"
                  f"{mcc(tp, fp, tn, fn):>8.3f}{tp:>7,}{fp:>7,}{tn:>7,}{fn:>7,}")

    # Prior test: how far does the calibration move each stratum?
    reviewed = {k for k in universe if k in cv and cv[k][1] >= 1}
    unseen = universe - set(cv)
    print(f"\n{'='*78}")
    print("Prior test: actionable share by ClinVar exposure (Pass-Blind)")
    print(f"{'='*78}")
    print(f"{'stratum / mode':<50}{'n':>8}{'release':>9}{'mode':>10}{'shift':>8}")
    for label, keys in (("ClinVar-reviewed (>=1 star)", reviewed),
                        ("never in ClinVar", unseen)):
        if not keys:
            continue
        def share(run):
            r = runs[run]
            return sum(1 for k in keys
                       if r[k]["varviz_classification_blind"] in ACTIONABLE) / len(keys)
        a = share("release")
        for run in ("calibrated (no override)", "hybrid"):
            b = share(run)
            print(f"{label+' / '+run:<50}{len(keys):>8,}{100*a:>9.1f}%{100*b:>10.1f}%{100*(b-a):>+8.1f}")


if __name__ == "__main__":
    main()
