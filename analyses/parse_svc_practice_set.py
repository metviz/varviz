#!/usr/bin/env python3
"""Parse the SVC v4.0 pilot practice variant set into a tidy table.

Each worksheet is one variant, laid out as labelled blocks rather than columns:
Variant Data (gene, NM/NP HGVS, CAID, gnomAD link), Disease Data (disease name,
MONDO id, mechanism statement, MOI), Pop Frequency Data (prevalence estimate,
allelic heterogeneity, genetic heterogeneity, penetrance), Variant Impact Data
(REVEL, splicing) and Clinical Data (per-proband narrative).

The population-frequency block is the reason this file matters to VarViz: it
supplies, per variant, the four disease-model parameters that VarViz otherwise
requires the user to supply and that feed the maximum credible allele frequency
and hence PM2 -- prevalence, allelic heterogeneity, genetic heterogeneity and
penetrance. Every other VarViz cohort so far has run on assumed values.

Usage:
  python3 analyses/parse_svc_practice_set.py <Practice Variants Set.xlsx> [out.tsv]
"""

import csv
import re
import sys

import openpyxl

FIELDS = [
    "sheet",
    "gene",
    "transcript",
    "c_notation",
    "protein_acc",
    "p_notation",
    "caid",
    "disease",
    "mondo",
    "moi",
    "mechanism",
    "prevalence",
    "allelic_het",
    "genetic_het",
    "penetrance",
    "revel",
    "splicing",
    "n_probands",
    "gnomad_link",
    "vtype",
]

LABEL = {
    "prevalence estimate": "prevalence",
    "prevalence": "prevalence",
    "allelic heterogeneity": "allelic_het",
    "genetic (locus) heterogeneity": "genetic_het",
    "genetic heterogeneity": "genetic_het",
    "penetrance": "penetrance",
    "caid": "caid",
    "gene": "gene",
}


def cells(ws):
    for row in ws.iter_rows(values_only=True):
        vals = [v for v in row if v not in (None, "")]
        if vals:
            yield [str(v).strip() for v in vals]


def parse_sheet(ws):
    r = {k: "" for k in FIELDS}
    r["sheet"] = ws.title
    probands = 0
    for vals in cells(ws):
        head = vals[0]
        key = head.rstrip(":").strip().lower()

        if key in LABEL and len(vals) > 1:
            r[LABEL[key]] = vals[1]
            continue

        m = re.match(r"^(NM_[\w.]+):(c\..+)$", head)
        if m:
            r["transcript"], r["c_notation"] = m.group(1), m.group(2)
            continue
        m = re.match(r"^(NP_[\w.]+):(p\..+)$", head)
        if m:
            r["protein_acc"], r["p_notation"] = m.group(1), m.group(2)
            continue
        # "Variant | NM_...:c...." puts the HGVS in the second cell
        if key == "variant" and len(vals) > 1:
            m = re.match(r"^(NM_[\w.]+):(c\..+)$", vals[1])
            if m:
                r["transcript"], r["c_notation"] = m.group(1), m.group(2)
            continue

        if head.startswith("MONDO:"):
            r["mondo"] = head
            continue
        if head.startswith("MOI:"):
            r["moi"] = head.split(":", 1)[1].strip()
            continue
        if "mechanism of disease" in head.lower():
            r["mechanism"] = head
            continue
        if head.lower().startswith(("revel score", "revel:")):
            r["revel"] = head.split(":", 1)[1].strip() if ":" in head else head
            continue
        if "splicing" in head.lower() and len(head) < 90:
            r["splicing"] = head
            continue
        if head.startswith("http") and "gnomad" in head.lower():
            r["gnomad_link"] = head
            continue
        if re.match(r"^Proband\s*\d+\s*:", head):
            probands += 1
            continue
        # disease name: the line directly under the Disease Data header, before MONDO
        if (
            not r["disease"]
            and r["mondo"] == ""
            and r["gene"]
            and head
            not in (
                "Disease Data",
                "Variant Data",
                "Pop Frequency Data",
                "Variant Impact Data",
                "Clinical Data",
            )
            and not head.startswith("http")
            and len(head) < 90
            and ":" not in head
        ):
            r["disease"] = head
    r["n_probands"] = probands
    r["vtype"] = vtype(r["c_notation"], r["p_notation"])
    return r


def vtype(c, p):
    """Classify from HGVS. VarViz classifies missense only, so the pilot set has
    to be split before any comparison."""
    if re.search(r"\d[+-]\d", c) or re.search(r"[+-]\d+[ACGT]>", c):
        return "splice/intronic"
    if p.endswith("="):
        return "synonymous"
    if "fsTer" in p or "fs" in p:
        return "frameshift"
    if p.endswith("Ter") or "Ter" in p and "fs" not in p and "del" not in p:
        return "nonsense"
    if "del" in p or "dup" in p or "ins" in p or "del" in c or "dup" in c:
        return "indel"
    if re.fullmatch(r"p\.[A-Z][a-z]{2}\d+[A-Z][a-z]{2}", p):
        return "missense"
    if p in ("", "p.?"):
        return "non-coding/unknown"
    return "other"


def main():
    src = (
        sys.argv[1]
        if len(sys.argv) > 1
        else "/home/raghu/Downloads/Practice Variants Set.xlsx"
    )
    out = sys.argv[2] if len(sys.argv) > 2 else "analyses/tmp/svc_practice_set.tsv"

    wb = openpyxl.load_workbook(src, data_only=True)
    rows = [parse_sheet(ws) for ws in wb]

    with open(out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    print(f"{len(rows)} sheets -> {out}\n")
    hdr = f"{'sheet':16s} {'gene':9s} {'p_notation':22s} {'MOI':6s} {'prev':12s} {'hetA':5s} {'hetG':5s} {'pen':5s}"
    print(hdr)
    print("-" * len(hdr))
    for r in rows:
        print(
            f"{r['sheet']:16s} {r['gene']:9s} {r['p_notation'][:22]:22s} {r['moi'][:6]:6s} "
            f"{r['prevalence'][:12]:12s} {r['allelic_het'][:5]:5s} {r['genetic_het'][:5]:5s} "
            f"{r['penetrance'][:5]:5s}"
        )

    import collections

    print("\nvariant types:", dict(collections.Counter(r["vtype"] for r in rows)))
    mis = [r for r in rows if r["vtype"] == "missense"]
    print(f"missense (VarViz-classifiable): {len(mis)}")
    for r in mis:
        print(
            f"   {r['sheet']:16s} {r['gene']:8s} {r['p_notation']:20s} "
            f"prev={r['prevalence'] or '-':12s} hetA={r['allelic_het'] or '-':5s} "
            f"hetG={r['genetic_het'] or '-':5s} pen={r['penetrance'] or '-':5s} MOI={r['moi'][:22]}"
        )

    missing = [r["sheet"] for r in rows if not r["prevalence"]]
    print(f"\nsheets with no prevalence parsed: {missing if missing else 'none'}")
    prot = sum(1 for r in rows if r["p_notation"].startswith("p."))
    print(f"sheets with a p. notation: {prot}/{len(rows)}")


if __name__ == "__main__":
    main()
