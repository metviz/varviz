#!/usr/bin/env python3
"""Recompute Supplementary Section S4 from a current external-163 harness run.

The cohort is the Ghasemnejad et al. 2026 manual-curation set (65 missense
variants, 44 genes, 83 probands). Its first pass was run on 2026-08-24, before
release 1.0.0, and its tag strings still carry PP3_strong from the retired
AlphaMissense+REVEL proxy. Those figures therefore do not describe the release
the manuscript reports, so every number in S4 is regenerated here from a run of
the current engine.

Two inputs are release-independent and are carried forward from the original
assembly, keyed on gene + p. notation: the per-variant proband counts and
transcript/cDNA notation, and the `intervar_rerun` column obtained from the
public wInterVar API. Neither depends on the VarViz release.

WNT10A is one of four genes whose recessive inheritance rests on curated
gene-disease MOI that the het-only observed zygosity does not corroborate, and
it is the only one of the four that moves a call. The conservative reading
substitutes its three variants from a monoallelic pass; both readings are
reported.

Usage:
    python3 analyses/build_external163_s4.py \
        --run analyses/ps_external163_v122 \
        --wnt10a analyses/ps_external163_v122_wnt10a_mono \
        --out analyses/external_validation
"""
import argparse
import csv
import json
import sys
from collections import Counter, OrderedDict

PL = ("Pathogenic", "Likely Pathogenic")
# Excluded from the InterVar comparison: both sit on non-MANE transcripts, so
# they have no hg38 coordinate under VarViz's MANE-based numbering.
NO_HG38 = {("CRB1", "p.C790F"), ("TBX3", "p.A549D")}
COLLAGEN = {"COL1A1", "COL1A2"}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def key(row):
    return (row["gene"], row["p_notation"])


def load_run(run_dir, wnt10a_dir=None):
    """Summary rows keyed by (gene, p.notation); WNT10A optionally overridden."""
    rows = OrderedDict((key(r), r) for r in read_tsv(f"{run_dir}/summary.tsv"))
    if len(rows) != 65:
        sys.exit(f"{run_dir}/summary.tsv holds {len(rows)} variants, expected 65")
    if wnt10a_dir:
        sub = read_tsv(f"{wnt10a_dir}/summary.tsv")
        if len(sub) != 3:
            sys.exit(f"{wnt10a_dir}/summary.tsv holds {len(sub)} variants, expected 3")
        for r in sub:
            if key(r) not in rows:
                sys.exit(f"WNT10A pass carries {key(r)}, absent from the main pass")
            rows[key(r)] = r
    return rows


def sens(rows, field):
    hit = [r for r in rows if r[field] in PL]
    return len(hit), len(rows), 100.0 * len(hit) / len(rows)


def report(rows_by_key, label):
    rows = list(rows_by_key.values())
    out = {"label": label, "n": len(rows)}
    for pas, field in (("full", "varviz_classification_full"),
                       ("blind", "varviz_classification_blind")):
        n, tot, pct = sens(rows, field)
        out[f"{pas}_pl"] = n
        out[f"{pas}_pct"] = round(pct, 1)
        out[f"{pas}_dist"] = dict(Counter(r[field] for r in rows))
    out["disagree"] = sum(
        r["varviz_classification_full"] != r["varviz_classification_blind"]
        for r in rows)
    out["disagree_pct"] = round(100.0 * out["disagree"] / len(rows), 1)
    out["pm1_pathways"] = dict(Counter(r["pm1_pathway"] for r in rows).most_common())
    out["pm1_mds_any"] = sum("mds" in r["pm1_pathway"] and
                             r["pm1_pathway"] != "mds_unavailable" for r in rows)
    # Collagen genes carry Gly-X-Y motifs that the Pfam PSSM scores strongly, so
    # the non-collagen subset shows the result is not a motif artefact.
    nc = [r for r in rows if r["gene"] not in COLLAGEN]
    out["non_collagen_n"] = len(nc)
    out["non_collagen_full_pct"] = round(sens(nc, "varviz_classification_full")[2], 1)
    out["non_collagen_blind_pct"] = round(sens(nc, "varviz_classification_blind")[2], 1)
    out["misses"] = [
        {"gene": r["gene"], "variant": r["p_notation"],
         "full": r["varviz_classification_full"], "pm1_pathway": r["pm1_pathway"]}
        for r in rows if r["varviz_classification_full"] not in PL]
    return out


def intervar_block(rows_by_key, prior):
    """Like-for-like VarViz vs InterVar on the 63 variants with hg38 coordinates."""
    sub = [r for k, r in rows_by_key.items() if k not in NO_HG38]
    iv = Counter()
    for k, r in rows_by_key.items():
        if k in NO_HG38:
            continue
        call = (prior.get(k) or {}).get("intervar_rerun", "").strip()
        if not call:
            sys.exit(f"no intervar_rerun call carried forward for {k}")
        iv[call] += 1
    iv_pl = sum(n for c, n in iv.items()
                if c.lower() in ("pathogenic", "likely pathogenic"))
    n = len(sub)
    return {
        "n": n,
        "varviz_full_pl": sens(sub, "varviz_classification_full")[0],
        "varviz_full_pct": round(sens(sub, "varviz_classification_full")[2], 1),
        "varviz_blind_pl": sens(sub, "varviz_classification_blind")[0],
        "varviz_blind_pct": round(sens(sub, "varviz_classification_blind")[2], 1),
        "intervar_pl": iv_pl,
        "intervar_pct": round(100.0 * iv_pl / n, 1),
        "intervar_dist": dict(iv.most_common()),
        "varviz_benign_side": sum(
            r["varviz_classification_full"].startswith(("Benign", "Likely Benign"))
            for r in sub),
    }


def proband_block(rows_by_key, parsed_path, testset_path):
    """Per-proband reach: a proband counts as reached when at least one of its
    test-set variants is classified P/LP.

    Summing the per-variant proband counts overcounts, because one proband
    carries two of the 65 variants: the 65 variants span 84 proband-variant
    pairs but only 83 distinct probands. The sample identifiers live in the
    parsing intermediate, so the link runs variant -> (gene, cDNA) -> sample.
    """
    testset = read_tsv(testset_path)
    cnot = {(t["gene"], t["cnot"]): (t["gene"], t["p_notation"]) for t in testset}
    by_sample = {}
    for r in read_tsv(parsed_path):
        k = cnot.get((r["gene"], r["cnot"]))
        if k:
            by_sample.setdefault(r["sample_id"], set()).add(k)
    reached = {"full": 0, "blind": 0}
    for variants in by_sample.values():
        for pas, field in (("full", "varviz_classification_full"),
                           ("blind", "varviz_classification_blind")):
            if any(rows_by_key[k][field] in PL for k in variants):
                reached[pas] += 1
    n = len(by_sample)
    return {"probands": n,
            "full": reached["full"], "full_pct": round(100.0 * reached["full"] / n, 1),
            "blind": reached["blind"],
            "blind_pct": round(100.0 * reached["blind"] / n, 1)}


def write_table(rows_by_key, prior, path):
    """Table S7 source: current VarViz columns beside the carried-forward ones."""
    cols = ["gene", "p_notation", "c_notation", "transcript", "n_probands",
            "inheritance", "inh_confidence", "varviz_full", "varviz_full_pts",
            "varviz_blind", "varviz_blind_pts", "pm1_pathway", "intervar_rerun",
            "franklin_published", "genebe_published", "intervar_published",
            "tapes_published"]
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for k, r in rows_by_key.items():
            p = prior.get(k, {})
            w.writerow({
                "gene": r["gene"], "p_notation": r["p_notation"],
                "c_notation": p.get("c_notation", ""),
                "transcript": p.get("transcript", ""),
                "n_probands": p.get("n_probands", ""),
                "inheritance": p.get("inheritance", ""),
                "inh_confidence": p.get("inh_confidence", ""),
                "varviz_full": r["varviz_classification_full"],
                "varviz_full_pts": r["varviz_pts_full"],
                "varviz_blind": r["varviz_classification_blind"],
                "varviz_blind_pts": r["varviz_pts_blind"],
                "pm1_pathway": r["pm1_pathway"],
                "intervar_rerun": p.get("intervar_rerun", ""),
                "franklin_published": p.get("franklin_published", ""),
                "genebe_published": p.get("genebe_published", ""),
                "intervar_published": p.get("intervar_published", ""),
                "tapes_published": p.get("tapes_published", ""),
            })


def markdown_tables(rows_by_key, prior):
    """The two data tables of §S4 as markdown, ready to paste into the supplement."""
    miss = [(k, r) for k, r in rows_by_key.items()
            if r["varviz_classification_full"] not in PL]
    lines = ["| Gene | Variant | Pass-Full | Pts | PM1 pathway | InterVar (re-run) |",
             "|---|---|---|---|---|---|"]
    for k, r in miss:
        iv = (prior.get(k) or {}).get("intervar_rerun", "") or "n/a"
        path = r["pm1_pathway"] or "none"
        lines.append(f"| *{r['gene']}* | {r['p_notation']} | "
                     f"{r['varviz_classification_full']} | {r['varviz_pts_full']} | "
                     f"{path} | {iv} |")
    misses_md = "\n".join(lines)

    lines = ["| Gene | Variant | cDNA | Inh. | Pass-Full | Pts | Pass-Blind | Pts | "
             "PM1 pathway | InterVar |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for k, r in rows_by_key.items():
        p = prior.get(k, {})
        iv = p.get("intervar_rerun", "") or "n/a"
        path = r["pm1_pathway"] or "none"
        lines.append(
            f"| *{r['gene']}* | {r['p_notation']} | {p.get('c_notation','')} | "
            f"{ {'biallelic': 'bi', 'monoallelic': 'mono'}.get(p.get('inheritance',''), '?') } | {r['varviz_classification_full']} | "
            f"{r['varviz_pts_full']} | {r['varviz_classification_blind']} | "
            f"{r['varviz_pts_blind']} | {path} | {iv} |")
    full_md = "\n".join(lines)
    return misses_md, full_md


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", required=True, help="current harness run directory")
    ap.add_argument("--wnt10a", help="WNT10A monoallelic run directory")
    ap.add_argument("--prior", default="analyses/external_validation/"
                                       "TABLE_external_comparison.tsv",
                    help="source of the release-independent columns")
    ap.add_argument("--out", default="analyses/external_validation")
    ap.add_argument("--parsed",
                    default="analyses/external_validation/variants_parsed.tsv")
    ap.add_argument("--testset",
                    default="analyses/external_validation/missense_testset.tsv")
    args = ap.parse_args()

    prior = {key(r): r for r in read_tsv(args.prior)}

    assigned = load_run(args.run)
    numbers = {"assigned_moi": report(assigned, "inheritance as assigned")}
    numbers["intervar_assigned"] = intervar_block(assigned, prior)

    primary = assigned
    if args.wnt10a:
        conservative = load_run(args.run, args.wnt10a)
        numbers["wnt10a_monoallelic"] = report(
            conservative, "WNT10A held monoallelic (conservative)")
        numbers["intervar_conservative"] = intervar_block(conservative, prior)
        numbers["probands"] = proband_block(conservative, args.parsed, args.testset)
        primary = conservative
    else:
        numbers["probands"] = proband_block(assigned, args.parsed, args.testset)

    write_table(primary, prior, f"{args.out}/TABLE_external_comparison.tsv")
    misses_md, full_md = markdown_tables(primary, prior)
    with open(f"{args.out}/S4_tables.md", "w") as fh:
        fh.write("<!-- Supplementary Table S8 -->\n" + misses_md +
                 "\n\n<!-- Supplementary Table S10 -->\n" + full_md + "\n")
    with open(f"{args.out}/S4_numbers.json", "w") as fh:
        json.dump(numbers, fh, indent=2)
    print(json.dumps(numbers, indent=2))


if __name__ == "__main__":
    main()
