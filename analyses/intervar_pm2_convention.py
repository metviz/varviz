#!/usr/bin/env python3
"""Infer which PM2 weight InterVar applies, from its own per-criterion flags.

VarViz scores PM2 at Supporting, following the ClinGen SVI 2020 recommendation.
InterVar implements Richards 2015, which predates it and scores PM2 at
Moderate. That difference is worth more than a footnote in any head-to-head:
on the external cohort every one of InterVar's P/LP calls fires PM2, so the
comparison silently measures a criterion-weight convention alongside whatever
else it is meant to measure.

The wInterVar API returns a per-criterion flag for every ACMG code, so the
question is decidable rather than assumed. This re-derives InterVar's own
verdict from its flags under both weights and reports which one reproduces it.

Usage:
    python3 analyses/intervar_pm2_convention.py \
        analyses/external_validation/intervar_raw.json
"""
import argparse
import json
import sys
from collections import Counter

PVS = ["PVS1"]
PS = [f"PS{i}" for i in (1, 2, 3, 4)]
PM = [f"PM{i}" for i in range(1, 7)]
PP = [f"PP{i}" for i in range(1, 6)]
BA = ["BA1"]
BS = [f"BS{i}" for i in (1, 2, 3, 4)]
BP = [f"BP{i}" for i in range(1, 8)]


def richards(flags, pm2_supporting):
    """The Richards 2015 combining rules, with PM2 counted either way."""
    pvs = sum(flags[k] for k in PVS)
    ps = sum(flags[k] for k in PS)
    pm = sum(flags[k] for k in PM)
    pp = sum(flags[k] for k in PP)
    if pm2_supporting and flags["PM2"]:
        pm -= 1
        pp += 1
    ba = sum(flags[k] for k in BA)
    bs = sum(flags[k] for k in BS)
    bp = sum(flags[k] for k in BP)

    pathogenic = (
        (pvs >= 1 and (ps >= 1 or pm >= 2 or (pm >= 1 and pp >= 1) or pp >= 2))
        or ps >= 2
        or (ps == 1 and (pm >= 3 or (pm >= 2 and pp >= 2) or (pm >= 1 and pp >= 4)))
    )
    likely_path = (
        (pvs == 1 and pm == 1)
        or (ps == 1 and 1 <= pm <= 2)
        or (ps == 1 and pp >= 2)
        or pm >= 3
        or (pm >= 2 and pp >= 2)
        or (pm >= 1 and pp >= 4)
    )
    benign = ba >= 1 or bs >= 2
    likely_ben = (bs >= 1 and bp >= 1) or bp >= 2

    if benign:
        return "Benign"
    if pathogenic:
        return "Pathogenic"
    if likely_path:
        return "Likely pathogenic"
    if likely_ben:
        return "Likely benign"
    return "Uncertain significance"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("raw", help="wInterVar responses, keyed 'GENE|p.notation'")
    args = ap.parse_args()
    data = json.load(open(args.raw))

    agree = {"moderate": 0, "supporting": 0}
    conflicts = []
    for key, flags in data.items():
        got = flags["Intervar"].strip().lower()
        mod = richards(flags, pm2_supporting=False)
        sup = richards(flags, pm2_supporting=True)
        agree["moderate"] += mod.lower() == got
        agree["supporting"] += sup.lower() == got
        if mod != sup:
            conflicts.append((key, flags["Intervar"].strip(), mod, sup))

    n = len(data)
    print(f"variants: {n}")
    print(f"PM2 as Moderate   reproduces InterVar's verdict on {agree['moderate']}/{n}")
    print(f"PM2 as Supporting reproduces InterVar's verdict on {agree['supporting']}/{n}")

    pl = {k: f for k, f in data.items()
          if f["Intervar"].strip().lower() in ("pathogenic", "likely pathogenic")}
    firing = sum(f["PM2"] == 1 for f in pl.values())
    survives = [k for k, f in pl.items()
                if richards(f, pm2_supporting=True).lower() in
                ("pathogenic", "likely pathogenic")]
    print(f"\nInterVar P/LP calls: {len(pl)}; firing PM2: {firing}; "
          f"surviving PM2 at Supporting: {len(survives)}")
    for k in sorted(survives):
        print(f"  survives: {k}")

    crit = PVS + PS + PM + PP
    fired = Counter()
    for f in data.values():
        for k in crit + BA + BS + BP:
            fired[k] += f[k]
    per = [sum(f[k] for k in crit) for f in data.values()]
    print(f"\npathogenic criteria fired per variant: mean {sum(per)/len(per):.1f}, "
          f"max {max(per)}")
    print("most-fired criteria:",
          ", ".join(f"{k}={v}" for k, v in fired.most_common(6)))

    if conflicts:
        print(f"\nverdict differs between the two weights on {len(conflicts)} variants:")
        for key, got, mod, sup in conflicts:
            print(f"  {key:24s} InterVar={got:20s} moderate={mod:20s} supporting={sup}")

    # The whole point is that the answer is decidable; say so if it is not.
    if agree["moderate"] == agree["supporting"]:
        print("\nINCONCLUSIVE: both weights reproduce the same number of verdicts.")
        sys.exit(1)


if __name__ == "__main__":
    main()
