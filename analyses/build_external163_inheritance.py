#!/usr/bin/env python3
"""Assign per-gene inheritance for the 44-gene external cohort and rebuild its universe.

The published §S4.8 run classified all 44 genes as monoallelic (PM2 < 1e-4)
although many are recessive. server.R:5136-5144 makes the thresholds
inheritance-aware -- biallelic uses PM2 < 0.01 and BS1 > 0.05 -- so running a
recessive gene as monoallelic makes PM2 fire ~100x too readily, biasing the
cohort toward the pathogenic direction. That is the limitation recorded in
Supplementary S4.8.

Two evidence sources, combined explicitly:

  zygosity  Observed genotype of the pathogenic-scored variants in each gene's
            own proband, taken from the LIRICAL output of the source capsule.
            Direct observation in the very cases being classified, so it is the
            primary evidence. hom-only -> biallelic; het-only -> monoallelic.

  moi       Curated gene-disease mode of inheritance for the disease the cohort
            is about. Used where zygosity is absent or ambiguous, and recorded
            for every gene so conflicts are visible rather than silent.

Only inh_param is set. No prevalence estimate is available per gene here, so
af_cutoff stays at the ACMG default 1e-4 for monoallelic genes; for biallelic
genes server.R ignores af_cutoff for PM2 and applies the 0.01 recessive
standard, so no prevalence is needed.

Usage:
  python3 analyses/build_external163_inheritance.py
"""

import csv
import sys

HINTS = "analyses/tmp/compare/out/inheritance_hints.tsv"
TESTSET = "analyses/tmp/compare/out/missense_testset.tsv"
OUT_MAP = "analyses/tmp/compare/out/inheritance_assignment.tsv"
OUT_UNI = "analyses/derived/variant_universe_external163_moi.tsv"

# Curated mode of inheritance for the disease each gene is implicated in within
# this cohort. AD/AR/XL; "AD/AR" where both are established and the zygosity
# evidence decides.
MOI = {
    "ABCA4": "AR",
    "ADAMTS1": "AD",
    "AXIN2": "AD",
    "BMP1": "AR",
    "BMP4": "AD",
    "BRAF": "AD",
    "CNGA3": "AR",
    "CNNM4": "AR",
    "COL1A1": "AD",
    "COL1A2": "AD",
    "CRB1": "AR",
    "DSPP": "AD",
    "EDA": "XL",
    "EDARADD": "AD/AR",
    "ENAM": "AD/AR",
    "ESRRB": "AR",
    "FAM20A": "AR",
    "FGFR2": "AD",
    "GALC": "AR",
    "KCNE1": "AD/AR",
    "KREMEN1": "AR",
    "LRP6": "AD",
    "MYO15A": "AR",
    "MYO7A": "AR",
    "OTOF": "AR",
    "P3H1": "AR",
    "PAX9": "AD",
    "PCDH15": "AR",
    "PITX2": "AD",
    "PJVK": "AR",
    "RPE65": "AR", "PTCH1": "AD",
    "PTH1R": "AD",
    "RDH12": "AR",
    "RGS9": "AR",
    "SLC26A4": "AR",
    "TBX3": "AD",
    "TMPRSS3": "AR",
    "TP63": "AD",
    "TULP1": "AR",
    "TYR": "AR",
    "USH2A": "AR",
    "WFS1": "AR",
    "WNT10A": "AR",
}

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
]


def decide(gene, hom, het):
    """Return (inh_param, basis, confidence)."""
    moi = MOI.get(gene, "")
    zyg = ""
    if hom and not het:
        zyg = "biallelic"
    elif het and not hom:
        zyg = "monoallelic"
    elif hom and het:
        zyg = "mixed"

    moi_inh = {"AR": "biallelic", "AD": "monoallelic", "XL": "monoallelic"}.get(moi, "")

    # EDA is X-linked; VarViz has no X-linked mode and hemizygous males read as
    # homozygous, so the zygosity signal is misleading. Monoallelic is the
    # closer approximation -- an X allele is exposed in males the way a dominant
    # allele is, and the recessive PM2 threshold of 0.01 would be far too loose.
    if moi == "XL":
        return (
            "monoallelic",
            f"X-linked (zyg {hom}hom/{het}het, not informative)",
            "medium",
        )

    if zyg in ("biallelic", "monoallelic"):
        if moi_inh and moi_inh != zyg:
            return moi_inh, f"MOI {moi} over zygosity {zyg} ({hom}hom/{het}het)", "low"
        if moi == "AD/AR":
            return zyg, f"zygosity {hom}hom/{het}het decides AD/AR gene", "medium"
        return zyg, f"zygosity {hom}hom/{het}het agrees with MOI {moi}", "high"

    if moi_inh:
        return moi_inh, f"MOI {moi}; zygosity {'mixed' if zyg else 'absent'}", "medium"
    return "monoallelic", "no evidence; ACMG default", "low"


def main():
    hints = {r["gene"]: r for r in csv.DictReader(open(HINTS), delimiter="\t")}
    variants = list(csv.DictReader(open(TESTSET), delimiter="\t"))
    genes = sorted({v["gene"] for v in variants})

    unknown = [g for g in genes if g not in MOI]
    if unknown:
        sys.exit(f"no curated MOI for: {unknown}")

    assign = {}
    rows = []
    for g in genes:
        h = hints.get(g, {})
        hom = int(h.get("hom_path") or 0)
        het = int(h.get("het_path") or 0)
        inh, basis, conf = decide(g, hom, het)
        assign[g] = inh
        rows.append(
            {
                "gene": g,
                "inh_param": inh,
                "curated_moi": MOI[g],
                "hom_path": hom,
                "het_path": het,
                "basis": basis,
                "confidence": conf,
            }
        )

    with open(OUT_MAP, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    with open(OUT_UNI, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=UNIVERSE_COLS, delimiter="\t")
        w.writeheader()
        for v in variants:
            w.writerow(
                {
                    "gene": v["gene"],
                    "p_notation": v["p_notation"],
                    "p_notation_raw": v["p_notation"],
                    "source": "Ghasemnejad2026_manual",
                    "label": "Causative",
                    "score_raw": "NA",
                    "study": "NA",
                    "set": "manual_163new",
                    "hgvs_g": "NA",
                    "clinvar_alleleid": "NA",
                    "hgvs_pro_raw": v["cnot"],
                    "gnomad_ac": "NA",
                    "gnomad_an": "NA",
                    "gnomad_af": "NA",
                    "gnomad_singleton": "NA",
                    "gnomad_present": "NA",
                    "inh_param": assign[v["gene"]],
                }
            )

    bi = [r for r in rows if r["inh_param"] == "biallelic"]
    print(
        f"{len(genes)} genes: {len(bi)} biallelic, {len(genes) - len(bi)} monoallelic"
    )
    print(f"  -> {OUT_MAP}\n  -> {OUT_UNI}\n")
    hdr = (
        f"{'gene':9s} {'inh':12s} {'MOI':6s} {'hom':>4s} {'het':>4s} {'conf':7s} basis"
    )
    print(hdr)
    print("-" * 96)
    for r in sorted(rows, key=lambda x: (x["inh_param"], x["gene"])):
        print(
            f"{r['gene']:9s} {r['inh_param']:12s} {r['curated_moi']:6s} "
            f"{r['hom_path']:>4d} {r['het_path']:>4d} {r['confidence']:7s} {r['basis']}"
        )

    low = [r["gene"] for r in rows if r["confidence"] == "low"]
    print(
        f"\nlow-confidence assignments (zygosity contradicts curated MOI, or no evidence): {low}"
    )
    nvar = sum(1 for v in variants if assign[v["gene"]] == "biallelic")
    print(
        f"variants affected by the change: {nvar}/{len(variants)} now biallelic "
        f"(PM2 < 0.01 instead of < 1e-4)"
    )


if __name__ == "__main__":
    main()
