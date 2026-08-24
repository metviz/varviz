#!/usr/bin/env python3
"""Compare the 44-gene external cohort before and after the inheritance correction.

Original run (analyses/ps_external) classified all 44 genes as monoallelic,
PM2 < 1e-4. Corrected run (analyses/ps_external_moi) sets inh_param per gene
from observed proband zygosity plus curated gene-disease MOI, so the 26
recessive genes get the ACMG recessive thresholds instead (PM2 < 0.01,
BS1 > 0.05; server.R:5136-5144).

Sensitivity is expected to move for the biallelic genes only. Monoallelic genes
are a control: their parameters did not change, so any difference there would
indicate run-to-run instability rather than the correction, and is reported
separately for exactly that reason.

Usage:
  python3 analyses/compare_external163_moi.py
"""

import collections
import csv
import os
import sys

OLD = "analyses/ps_external/summary.tsv"
NEW = "analyses/ps_external_moi/summary.tsv"
ASSIGN = "analyses/tmp/compare/out/inheritance_assignment.tsv"
OUT = "analyses/tmp/compare/out/moi_comparison.tsv"

PLP = ("Pathogenic", "Likely Pathogenic")
ORDER = [
    "Pathogenic",
    "Likely Pathogenic",
    "VUS-High",
    "VUS-Mid",
    "VUS-Low",
    "Likely Benign",
    "Benign",
]


def load(p):
    return {
        (r["gene"], r["p_notation"]): r for r in csv.DictReader(open(p), delimiter="\t")
    }


def dist(rows, field):
    c = collections.Counter(r[field] for r in rows)
    return " ".join(f"{k}={c[k]}" for k in ORDER if c[k])


def main():
    for p in (OLD, NEW, ASSIGN):
        if not os.path.exists(p):
            sys.exit(f"missing {p} - has the corrected run finished?")

    inh = {r["gene"]: r for r in csv.DictReader(open(ASSIGN), delimiter="\t")}
    old, new = load(OLD), load(NEW)

    keys = sorted(set(old) & set(new))
    only_old, only_new = set(old) - set(new), set(new) - set(old)
    if only_old or only_new:
        print(
            f"WARNING asymmetric: {len(only_old)} only in old, {len(only_new)} only in new"
        )
        for k in sorted(only_old)[:5]:
            print(f"   old only: {k}")
        for k in sorted(only_new)[:5]:
            print(f"   new only: {k}")

    print(f"variants compared: {len(keys)}\n")

    rows = []
    for k in keys:
        o, n = old[k], new[k]
        g = k[0]
        rows.append(
            {
                "gene": g,
                "p_notation": k[1],
                "inh_param": inh.get(g, {}).get("inh_param", "?"),
                "confidence": inh.get(g, {}).get("confidence", ""),
                "full_old": o["varviz_classification_full"],
                "full_new": n["varviz_classification_full"],
                "pts_old": o["varviz_pts_full"],
                "pts_new": n["varviz_pts_full"],
                "blind_old": o["varviz_classification_blind"],
                "blind_new": n["varviz_classification_blind"],
                "tags_old": o["tags_full"],
                "tags_new": n["tags_full"],
            }
        )

    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    for label, sel in (
        ("ALL", lambda r: True),
        (
            "biallelic genes (parameters changed)",
            lambda r: r["inh_param"] == "biallelic",
        ),
        (
            "monoallelic genes (control - unchanged)",
            lambda r: r["inh_param"] == "monoallelic",
        ),
    ):
        sub = [r for r in rows if sel(r)]
        if not sub:
            continue
        n = len(sub)
        print(f"=== {label}: n={n} ===")
        for arm in ("full", "blind"):
            a = sum(1 for r in sub if r[f"{arm}_old"] in PLP)
            b = sum(1 for r in sub if r[f"{arm}_new"] in PLP)
            print(
                f"  Pass-{arm.capitalize():5s} P/LP  {a}/{n} ({100 * a / n:.1f}%)"
                f"  ->  {b}/{n} ({100 * b / n:.1f}%)   delta {b - a:+d}"
            )
        print(f"  old dist (full): {dist(sub, 'full_old')}")
        print(f"  new dist (full): {dist(sub, 'full_new')}")
        ch = [r for r in sub if r["full_old"] != r["full_new"]]
        print(f"  classification changed: {len(ch)}/{n}\n")

    changed = [
        r
        for r in rows
        if r["full_old"] != r["full_new"] or r["blind_old"] != r["blind_new"]
    ]
    print(f"=== per-variant changes ({len(changed)}) ===")
    if not changed:
        print("  none")
    for r in sorted(changed, key=lambda x: (x["inh_param"], x["gene"])):
        print(
            f"  {r['gene']:9s} {r['p_notation']:10s} [{r['inh_param'][:4]}/{r['confidence'][:4]}]"
        )
        if r["full_old"] != r["full_new"]:
            print(
                f"      Full : {r['full_old']:18s} ({r['pts_old']:>3s} pts)  ->  "
                f"{r['full_new']:18s} ({r['pts_new']:>3s} pts)"
            )
        if r["blind_old"] != r["blind_new"]:
            print(
                f"      Blind: {r['blind_old']:18s}          ->  {r['blind_new']:18s}"
            )
        to = set(r["tags_old"].replace(" ", "").split(","))
        tn = set(r["tags_new"].replace(" ", "").split(","))
        if to != tn:
            gone, came = sorted(to - tn), sorted(tn - to)
            print(f"      tags : -{gone}  +{came}")

    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
