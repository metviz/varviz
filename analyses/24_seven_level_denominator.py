#!/usr/bin/env python3
"""Seven-level distribution computed per record vs per distinct variant.

analyses/derived/variant_universe_gnomad.tsv holds 90,701 variant x scoreset
records that collapse to 54,454 distinct (gene, p_notation). Supplementary
Table S3 currently tabulates all 90,701 records, so a variant assayed in five
MaveDB scoresets contributes five times to the distribution, while the S4.2
body text describes the same table as covering "each of the 54,454 distinct
missense variants".

This prints both denominators side by side so the choice can be made on the
size of the difference rather than on principle alone.

The harness classifies distinct variants, so ps_final/summary.tsv already has
one row per distinct variant. Per-record counts are obtained by weighting each
distinct variant by how many universe rows carry it.

Usage:
  python3 analyses/24_seven_level_denominator.py
"""

import collections
import csv
import sys

UNIVERSE = "analyses/derived/variant_universe_gnomad.tsv"
SUMMARY = "analyses/ps_final/summary.tsv"

ORDER = [
    "Pathogenic",
    "Likely Pathogenic",
    "VUS-High",
    "VUS-Mid",
    "VUS-Low",
    "Likely Benign",
    "Benign",
]

# Table S3 as published, for comparison.
PUBLISHED = {
    "Pathogenic": (17656, 7868),
    "Likely Pathogenic": (57942, 38142),
    "VUS-High": (7824, 32820),
    "VUS-Mid": (6682, 11253),
    "VUS-Low": (568, 589),
    "Likely Benign": (12, 12),
    "Benign": (17, 17),
}


def main():
    weight = collections.Counter()
    for r in csv.DictReader(open(UNIVERSE), delimiter="\t"):
        weight[(r["gene"], r["p_notation"])] += 1
    n_rec, n_var = sum(weight.values()), len(weight)
    print(
        f"universe: {n_rec:,} records -> {n_var:,} distinct variants "
        f"({n_rec / n_var:.2f} records per variant)\n"
    )

    rows = list(csv.DictReader(open(SUMMARY), delimiter="\t"))
    seen = set()
    per_var = {"full": collections.Counter(), "blind": collections.Counter()}
    per_rec = {"full": collections.Counter(), "blind": collections.Counter()}
    missing = 0
    for r in rows:
        k = (r["gene"], r["p_notation"])
        if k in seen:
            continue
        seen.add(k)
        w = weight.get(k)
        if w is None:
            missing += 1
            continue
        for arm in ("full", "blind"):
            c = r[f"varviz_classification_{arm}"]
            per_var[arm][c] += 1
            per_rec[arm][c] += w
    if missing:
        print(f"WARNING {missing} classified variants absent from the universe\n")

    tv = sum(per_var["full"].values())
    tr = sum(per_rec["full"].values())
    print(f"classified: {tv:,} distinct variants  /  {tr:,} records\n")

    hdr = (
        f"{'subclass':18s} | {'PER RECORD (as published)':>28s} | "
        f"{'PER DISTINCT VARIANT':>26s} | {'shift':>12s}"
    )
    print(hdr)
    print("-" * len(hdr))
    for arm, label in (("full", "Pass-Full"), ("blind", "Pass-Blind")):
        print(f"\n{label}")
        for c in ORDER:
            rec, var = per_rec[arm][c], per_var[arm][c]
            pr = 100 * rec / tr if tr else 0
            pv = 100 * var / tv if tv else 0
            pub = PUBLISHED[c][0 if arm == "full" else 1]
            flag = "" if rec == pub else f"  (published {pub:,})"
            print(
                f"  {c:16s} | {rec:>10,} {pr:>7.2f}% | {var:>10,} {pv:>7.2f}% | "
                f"{pv - pr:+7.2f} pp{flag}"
            )
        print(
            f"  {'Total':16s} | {tr:>10,} {100.0:>7.2f}% | {tv:>10,} {100.0:>7.2f}% |"
        )

    print("\nlargest absolute shift, percentage points:")
    for arm, label in (("full", "Pass-Full"), ("blind", "Pass-Blind")):
        d = max(
            ORDER,
            key=lambda c: abs(100 * per_var[arm][c] / tv - 100 * per_rec[arm][c] / tr),
        )
        pv = 100 * per_var[arm][d] / tv
        pr = 100 * per_rec[arm][d] / tr
        print(f"  {label:11s} {d:18s} {pr:6.2f}% -> {pv:6.2f}%  ({pv - pr:+.2f} pp)")

    for arm, label in (("full", "Pass-Full"), ("blind", "Pass-Blind")):
        plp_r = sum(per_rec[arm][c] for c in ("Pathogenic", "Likely Pathogenic"))
        plp_v = sum(per_var[arm][c] for c in ("Pathogenic", "Likely Pathogenic"))
        print(
            f"  {label:11s} P+LP combined      "
            f"{100 * plp_r / tr:6.2f}% -> {100 * plp_v / tv:6.2f}%  "
            f"({100 * plp_v / tv - 100 * plp_r / tr:+.2f} pp)"
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
