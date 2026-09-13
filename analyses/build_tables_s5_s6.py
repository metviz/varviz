#!/usr/bin/env python3
"""Regenerate Supplementary Tables S5 and S6 from a classification run.

S5 is the seven-level distribution under both passes; S6 is the Pass-Full to
Pass-Blind transition matrix. Both were computed ad hoc when first written, so
neither could be checked after the underlying run was re-made. This rebuilds
them from summary.tsv and the variant universe, and compares the result against
the values currently in the supplement so a drift is visible rather than
silent.

Two denominators are in play. The harness classifies distinct (gene,
p_notation) variants, so summary.tsv has one row per distinct variant; the
universe holds one row per variant x scoreset record, so a variant assayed in
five MaveDB scoresets weighs five times in a per-record count. S5 reports
distinct counts with per-record counts in parentheses, which is what this
computes.

  python3 analyses/build_tables_s5_s6.py [--run analyses/ps_final_v221]
"""
import argparse, collections, csv, os, sys

ORDER = ["Pathogenic", "Likely Pathogenic", "VUS-High", "VUS-Mid",
         "VUS-Low", "Likely Benign", "Benign"]
SHORT = ["P", "LP", "VUS-H", "VUS-M", "VUS-L", "LB", "B"]
RANK = {c: i for i, c in enumerate(ORDER)}      # 0 = most pathogenic

# The values currently in analyses/humu/supp_humu.md, for drift detection.
PUBLISHED_S5 = {
    "Pathogenic":        (2431, 4240, 1254, 2152),
    "Likely Pathogenic": (22402, 41714, 15812, 32069),
    "VUS-High":          (25231, 34216, 9043, 15697),
    "VUS-Mid":           (3777, 8856, 27061, 38272),
    "VUS-Low":           (592, 1651, 1263, 2487),
    "Likely Benign":     (5, 7, 5, 7),
    "Benign":            (16, 17, 16, 17),
}
PUBLISHED_MOVED = 31786


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="analyses/ps_final_v221")
    ap.add_argument("--universe", default="analyses/derived/variant_universe_gnomad.tsv")
    a = ap.parse_args()

    summary = os.path.join(a.run, "summary.tsv")
    for p in (summary, a.universe):
        if not os.path.isfile(p):
            sys.exit(f"missing: {p}")

    # One row per distinct variant.
    rows = {}
    for r in csv.DictReader(open(summary, newline="", encoding="utf8"), delimiter="\t"):
        rows.setdefault((r["gene"], r["p_notation"]), r)

    # Per-record weight: how many universe rows carry each distinct variant.
    weight = collections.Counter()
    for r in csv.DictReader(open(a.universe, newline="", encoding="utf8"), delimiter="\t"):
        weight[(r["gene"], r["p_notation"])] += 1

    unweighted = [k for k in rows if k not in weight]
    if unweighted:
        print(f"!! {len(unweighted):,} classified variants absent from the universe "
              f"(per-record counts would undercount); e.g. {unweighted[:2]}")

    dn = {p: collections.Counter() for p in ("full", "blind")}
    rn = {p: collections.Counter() for p in ("full", "blind")}
    for k, r in rows.items():
        w = weight.get(k, 1)
        for p in ("full", "blind"):
            c = r[f"varviz_classification_{p}"]
            dn[p][c] += 1
            rn[p][c] += w

    N = len(rows)
    NR = sum(rn["full"].values())
    print(f"\nTable S5 — seven-level distribution ({N:,} distinct / {NR:,} records)")
    print(f"{'Subclass':<20}{'Pass-Full n':>20}{'%':>8}{'Pass-Blind n':>20}{'%':>8}  check")
    drift = 0
    for c in ORDER:
        got = (dn["full"][c], rn["full"][c], dn["blind"][c], rn["blind"][c])
        exp = PUBLISHED_S5[c]
        ok = got == exp
        drift += (not ok)
        print(f"{c:<20}{f'{got[0]:,} ({got[1]:,})':>20}{100*got[0]/N:>8.2f}"
              f"{f'{got[2]:,} ({got[3]:,})':>20}{100*got[2]/N:>8.2f}  "
              f"{'ok' if ok else 'DRIFT was ' + str(exp)}")
    act = lambda p: 100 * (dn[p]["Pathogenic"] + dn[p]["Likely Pathogenic"]) / N
    print(f"{'actionable %':<20}{act('full'):>20.1f}{'':>8}{act('blind'):>20.1f}")

    # S6: transitions.
    mat = collections.Counter()
    for r in rows.values():
        mat[(r["varviz_classification_full"], r["varviz_classification_blind"])] += 1
    moved = sum(v for (f, b), v in mat.items() if f != b)
    toward_path = sum(v for (f, b), v in mat.items() if RANK[b] < RANK[f])

    print(f"\nTable S6 — Pass-Full to Pass-Blind transitions")
    print(f"{'Full \\\\ Blind':<20}" + "".join(f"{s:>9}" for s in SHORT) + f"{'% moved':>10}")
    for f in ORDER:
        tot = dn["full"][f]
        row = "".join(f"{mat[(f, b)]:>9,}" if mat[(f, b)] else f"{'':>9}" for b in ORDER)
        pm = 100 * (tot - mat[(f, f)]) / tot if tot else 0.0
        print(f"{f:<20}{row}{pm:>10.1f}")

    print(f"\nmoved: {moved:,} (published {PUBLISHED_MOVED:,}) "
          f"{'ok' if moved == PUBLISHED_MOVED else 'DRIFT'}")
    print(f"toward pathogenic: {toward_path:,} "
          f"{'ok — blinding only removes evidence' if toward_path == 0 else 'UNEXPECTED'}")
    drift += (moved != PUBLISHED_MOVED) + (toward_path != 0)
    print(f"\n{'S5/S6 reproduce the supplement exactly' if not drift else f'{drift} DRIFTS — supplement needs updating'}")
    return 1 if drift else 0


if __name__ == "__main__":
    sys.exit(main())
