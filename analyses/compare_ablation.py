#!/usr/bin/env python3
"""Compare two classification runs that differ by exactly one criterion weight.

Every ablation in this project is run as a full re-classification of the same
variant universe with one option changed, so the difference between two
summary.tsv files is attributable to that option and nothing else. This script
reports that difference in the shape §3.4 of the manuscript quotes it: the
share of calls that move, which transitions carry the movement, and how many
variants cross the actionable boundary in each direction.

Rows are keyed on (gene, p_notation). A summary holds more rows than distinct
variants, because a variant enumerated on more than one transcript position
appears once per checkpoint row; the key collapses those. A variant missing
from either side is reported rather than silently dropped, since a missing gene
is the signature of an incomplete run, not of a criterion change.

  python3 analyses/compare_ablation.py BASE_DIR ABLATION_DIR [--label NAME]
"""
import argparse, collections, csv, os, sys

ACTIONABLE = {"Pathogenic", "Likely Pathogenic"}
BENIGN     = {"Benign", "Likely Benign"}


def load(path):
    """Read summary.tsv into {(gene, p_notation): row}, refusing duplicates."""
    if not os.path.isfile(path):
        sys.exit(f"missing: {path}")
    rows, dupes = {}, 0
    with open(path, newline="", encoding="utf8") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            k = (r["gene"], r["p_notation"])
            if k in rows:
                dupes += 1
                # Identical re-enumeration is expected; a genuine disagreement
                # between two rows for the same variant is not, and would make
                # every downstream count depend on read order.
                if (r["varviz_classification_full"] != rows[k]["varviz_classification_full"]
                        or r["varviz_classification_blind"] != rows[k]["varviz_classification_blind"]):
                    sys.exit(f"conflicting duplicate rows for {k} in {path}")
                continue
            rows[k] = r
    return rows, dupes


def compare(base, abl, pas):
    """Transitions for one pass ('full' or 'blind'), base -> ablation."""
    col = f"varviz_classification_{pas}"
    moves = collections.Counter()
    changed = 0
    gained = lost = 0          # actionable status crossing, ablation vs base
    for k, b in base.items():
        a = abl.get(k)
        if a is None:
            continue
        bc, ac = b[col], a[col]
        if bc == ac:
            continue
        changed += 1
        moves[(bc, ac)] += 1
        if bc not in ACTIONABLE and ac in ACTIONABLE:
            gained += 1
        elif bc in ACTIONABLE and ac not in ACTIONABLE:
            lost += 1
    return changed, moves, gained, lost


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base"); ap.add_argument("ablation")
    ap.add_argument("--label", default="")
    ap.add_argument("--top", type=int, default=8)
    a = ap.parse_args()

    bp = os.path.join(a.base, "summary.tsv") if os.path.isdir(a.base) else a.base
    xp = os.path.join(a.ablation, "summary.tsv") if os.path.isdir(a.ablation) else a.ablation
    base, bd = load(bp)
    abl,  xd = load(xp)

    shared = base.keys() & abl.keys()
    only_b = base.keys() - abl.keys()
    only_x = abl.keys() - base.keys()

    print(f"\n{'='*74}")
    print(f"{a.label or 'ablation'}\n  base     {bp}\n  ablation {xp}")
    print(f"{'='*74}")
    print(f"variants: base {len(base):,}  ablation {len(abl):,}  shared {len(shared):,}"
          f"  (collapsed duplicate rows: {bd:,} / {xd:,})")
    if only_b or only_x:
        print(f"  !! present in only one run: base-only {len(only_b):,}, "
              f"ablation-only {len(only_x):,}")
        for k in list(only_b)[:3] + list(only_x)[:3]:
            print(f"     {k[0]} {k[1]}")
        print("     an incomplete run cannot be compared; check checkpoint counts")

    n = len(shared)
    for pas in ("full", "blind"):
        changed, moves, gained, lost = compare(
            {k: base[k] for k in shared}, abl, pas)
        pct = 100.0 * changed / n if n else 0.0
        print(f"\nPass-{pas.capitalize()}: {changed:,} of {n:,} calls change ({pct:.2f}%)")
        print(f"  actionable status: {gained:,} gained, {lost:,} lost")
        for (frm, to), c in moves.most_common(a.top):
            print(f"    {frm:>18}  ->  {to:<18} {c:7,}")
        tail = sum(moves.values()) - sum(c for _, c in moves.most_common(a.top))
        if tail:
            print(f"    {'other transitions':>18}      {'':<18} {tail:7,}")


if __name__ == "__main__":
    main()
