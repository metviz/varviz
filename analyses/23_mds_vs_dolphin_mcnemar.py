#!/usr/bin/env python3
"""Exact current-engine McNemar: MDS vs DOLPHIN Pass-Blind VUS resolution.

Both pathways are post-hoc +2 augmentations on the SAME current split-roles
MDS-off baseline (analyses/ps_nomds), using per-variant MDS/DOLPHIN deltas
(sequence properties, from varviz_classifications_mds.tsv). Reproduces the
manuscript's ~3,596 MDS figure and yields the paired McNemar statistic.

Augmentation rule (from 11_merge_dolphin.R): add +2 to Pass-Blind iff the
pathway fires AND pm1_pathway is '' or 'clinvar_hotspot'. VUS band = pts in
[-3, 5] (pts_to_bin, 22_mds_benchmark.R). Resolution is upward-only (+2).
"""

import csv
from math import erfc, sqrt

BASELINE = "analyses/ps_nomds/summary.tsv"  # current engine, MDS OFF
DELTAS = "analyses/derived/varviz_classifications_mds.tsv"  # per-variant deltas
GATE = {"", "clinvar_hotspot"}


def isvus(p):  # VUS-Low..VUS-High
    return -3 <= p <= 5


# per-variant MDS / DOLPHIN firing (constant across duplicate label rows)
delta = {}
for r in csv.DictReader(open(DELTAS), delimiter="\t"):
    delta[(r["gene"], r["p_notation"])] = (
        r["mds_pm1"] == "TRUE",
        r["dolphin_pm1"] == "TRUE" if r["dolphin_pm1"] else False,
    )

res_mds = res_dol = both = mds_only = dol_only = neither = n = 0
for r in csv.DictReader(open(BASELINE), delimiter="\t"):
    k = (r["gene"], r["p_notation"])
    if k not in delta:
        continue
    try:
        pts = float(r["varviz_pts_blind"])
    except ValueError:
        continue
    n += 1
    gate = r["pm1_pathway"] in GATE
    rm = isvus(pts) and not isvus(pts + 2 * (delta[k][0] and gate))
    rd = isvus(pts) and not isvus(pts + 2 * (delta[k][1] and gate))
    res_mds += rm
    res_dol += rd
    if rm and rd:
        both += 1
    elif rm:
        mds_only += 1
    elif rd:
        dol_only += 1
    else:
        neither += 1

chi = (abs(mds_only - dol_only) - 1) ** 2 / (mds_only + dol_only)
p = erfc(sqrt(chi / 2))  # underflows to 0 for large chi2
print(f"rows: {n}")
print(f"resolved  MDS={res_mds}  DOLPHIN={res_dol}  ({res_mds / max(res_dol, 1):.1f}x)")
print(
    f"McNemar 2x2: both={both} MDS-only={mds_only} DOLPHIN-only={dol_only} neither={neither}"
)
print(f"McNemar chi2(1)={chi:.1f}  p={p:.3e}  (p<1e-300 when 0)")
