#!/usr/bin/env python3
"""Build a VarViz universe for the SVC v4.0 pilot missense variants.

Unlike every previous VarViz cohort, this one carries the pilot's own curated
disease-model parameters per variant, so the maximum credible allele frequency
is computed rather than assumed.

af_cutoff reproduces server.R:6313-6322 exactly:
    monoallelic  0.5 * (1/prev) * hetA * hetG * (1/pen)
    biallelic    sqrt(1/prev) * hetA * sqrt(hetG) * (1/sqrt(pen))

MOI maps to VarViz's two inheritance settings: AD -> monoallelic,
AR -> biallelic, semidominant -> monoallelic (the stricter of the two).
X-linked variants are dropped: VarViz offers no X-linked mode, and the pilot
sheets for those variants give a binning approach instead of numeric
prevalence, so there is nothing to compute from.

Usage:
  python3 analyses/build_svc_pilot_universe.py [in.tsv] [out.tsv]
"""

import csv
import math
import re
import sys

AA3 = {
    "Ala": "A",
    "Arg": "R",
    "Asn": "N",
    "Asp": "D",
    "Cys": "C",
    "Gln": "Q",
    "Glu": "E",
    "Gly": "G",
    "His": "H",
    "Ile": "I",
    "Leu": "L",
    "Lys": "K",
    "Met": "M",
    "Phe": "F",
    "Pro": "P",
    "Ser": "S",
    "Thr": "T",
    "Trp": "W",
    "Tyr": "Y",
    "Val": "V",
}

P3 = re.compile(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$")


def to_one_letter(p):
    m = P3.match(p)
    if not m:
        return ""
    a, pos, b = m.groups()
    if a not in AA3 or b not in AA3:
        return ""
    return f"p.{AA3[a]}{pos}{AA3[b]}"


def prevalence(s):
    """'1 in 5000', '1 in 1,000,000', '1/1000' -> float fraction."""
    s = (s or "").strip().replace(",", "")
    m = re.search(r"1\s*(?:in|/)\s*([\d.]+)", s, re.I)
    return 1.0 / float(m.group(1)) if m else None


def moi_to_inh(moi):
    m = (moi or "").lower()
    if m.startswith("x-link") or "x-link" in m:
        return "x-linked"
    if "recessive" in m or m.strip() == "ar":
        return "biallelic"
    if "semidom" in m:
        return "monoallelic"
    if "dominant" in m or m.strip() == "ad":
        return "monoallelic"
    return ""


def af_cutoff(inh, prev, hetA, hetG, pen):
    if inh == "monoallelic":
        return 0.5 * prev * hetA * hetG * (1.0 / pen)
    return math.sqrt(prev) * hetA * math.sqrt(hetG) * (1.0 / math.sqrt(pen))


def num(s):
    try:
        return float(str(s).strip())
    except (TypeError, ValueError):
        return None


UNIVERSE_COLS = [
    "gene",
    "p_notation",
    "p_notation_raw",
    "source",
    "label",
    "score_raw",
    "study",
    "set",
    "hgvs_g",
    "clinvar_alleleid",
    "hgvs_pro_raw",
    "gnomad_ac",
    "gnomad_an",
    "gnomad_af",
    "gnomad_singleton",
    "gnomad_present",
    "inh_param",
    "af_cutoff",
    "prevalence_1_in_n",
    "allelic_het",
    "genetic_het",
    "penetrance",
]


def main():
    src = sys.argv[1] if len(sys.argv) > 1 else "analyses/tmp/svc_practice_set.tsv"
    out = (
        sys.argv[2]
        if len(sys.argv) > 2
        else "analyses/derived/variant_universe_svc_pilot.tsv"
    )

    rows = list(csv.DictReader(open(src), delimiter="\t"))
    kept, skipped = [], []

    for r in rows:
        # numbered pilot variants only -- drop the two FBN1 worked examples and
        # the unnumbered PKP2 sheet, which are teaching material.
        if not re.match(r"^v\d+ ", r["sheet"]):
            skipped.append((r["sheet"], "not a numbered pilot variant"))
            continue
        if r["vtype"] != "missense":
            skipped.append((r["sheet"], r["vtype"]))
            continue

        p1 = to_one_letter(r["p_notation"])
        if not p1:
            skipped.append((r["sheet"], f"cannot convert {r['p_notation']}"))
            continue

        inh = moi_to_inh(r["moi"])
        if inh == "x-linked":
            skipped.append(
                (r["sheet"], "X-linked: no VarViz mode, no numeric prevalence")
            )
            continue
        if not inh:
            skipped.append((r["sheet"], f"unmapped MOI {r['moi']!r}"))
            continue

        prev = prevalence(r["prevalence"])
        hetA, hetG, pen = (
            num(r["allelic_het"]),
            num(r["genetic_het"]),
            num(r["penetrance"]),
        )
        if None in (prev, hetA, hetG, pen) or 0 in (hetA, hetG, pen):
            skipped.append((r["sheet"], "incomplete parameters"))
            continue

        cut = af_cutoff(inh, prev, hetA, hetG, pen)
        kept.append(
            {
                "gene": r["gene"],
                "p_notation": p1,
                "p_notation_raw": p1,
                "source": "SVC_v4.0_pilot",
                "label": "pilot",
                "score_raw": "NA",
                "study": r["sheet"],
                "set": "svc_practice",
                "hgvs_g": "NA",
                "clinvar_alleleid": "NA",
                "hgvs_pro_raw": r["p_notation"],
                "gnomad_ac": "NA",
                "gnomad_an": "NA",
                "gnomad_af": "NA",
                "gnomad_singleton": "NA",
                "gnomad_present": "NA",
                "inh_param": inh,
                "af_cutoff": f"{cut:.6g}",
                "prevalence_1_in_n": f"{1 / prev:.0f}",
                "allelic_het": hetA,
                "genetic_het": hetG,
                "penetrance": pen,
            }
        )

    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=UNIVERSE_COLS, delimiter="\t")
        w.writeheader()
        w.writerows(kept)

    print(f"{len(kept)} variants -> {out}\n")
    hdr = f"{'study':12s} {'gene':8s} {'variant':10s} {'inh':12s} {'prev':>9s} {'hetA':>5s} {'hetG':>5s} {'pen':>5s} {'af_cutoff':>11s}"
    print(hdr)
    print("-" * len(hdr))
    for k in kept:
        print(
            f"{k['study']:12s} {k['gene']:8s} {k['p_notation']:10s} {k['inh_param']:12s} "
            f"{'1/' + k['prevalence_1_in_n']:>9s} {k['allelic_het']:>5} {k['genetic_het']:>5} "
            f"{k['penetrance']:>5} {k['af_cutoff']:>11s}"
        )
    print(f"\ndefault harness cutoff for comparison: 1e-04")
    print("\nskipped:")
    for s, why in skipped:
        print(f"   {s:16s} {why}")


if __name__ == "__main__":
    main()
