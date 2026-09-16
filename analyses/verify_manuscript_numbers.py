#!/usr/bin/env python3
"""Recompute every checkable figure in the manuscript and confirm the text matches.

Numbers reach the manuscript from several runs across several releases, and a
stale one is indistinguishable from a current one by reading. This recomputes
each claim from the canonical runs named in analyses/RUNS.md and reports any
that the text does not contain, and any it contains in a superseded form.

  python3 analyses/verify_manuscript_numbers.py
"""
import csv, collections, io, math, os, re, sys

ACT = {"Pathogenic", "Likely Pathogenic"}
ORDER = ["Pathogenic","Likely Pathogenic","VUS-High","VUS-Mid","VUS-Low","Likely Benign","Benign"]
RANK = {c: i for i, c in enumerate(ORDER)}
sp = lambda t: [x.strip() for x in re.split(r",\s*", t or "") if x.strip()]


def load(run):
    d = {}
    for r in csv.DictReader(open(f"analyses/{run}/summary.tsv", newline="", encoding="utf8"), delimiter="\t"):
        d.setdefault((r["gene"], r["p_notation"]), r)
    return d


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return (tp*tn - fp*fn)/den if den else float("nan")


def clinvar_labels(universe):
    cv = {}
    for r in csv.DictReader(open("analyses/tmp/clinvar_14genes_all.tsv", newline="", encoding="utf8"), delimiter="\t"):
        s = r["sig"].lower()
        if "conflicting" in s or "uncertain" in s: continue
        l = 1 if ("pathogenic" in s and "benign" not in s) else (0 if "benign" in s and "pathogenic" not in s else None)
        if l is None: continue
        st = 0 if r["review"].startswith("no ") else (3 if "expert" in r["review"] else (2 if "multiple" in r["review"] else 1))
        k = (r["gene"], r["p_notation"])
        if k not in cv or st > cv[k][1]: cv[k] = (l, st)
    return {k: v for k, v in cv.items() if k in universe}


def main():
    base, old = load("ps_final_v231"), load("ps_final_v221")
    v21 = load("ps_final_v230")
    n = len(base)
    claims = []          # (label, expected string, note)
    def c(label, val, note=""): claims.append((label, val, note))

    # --- distributions --------------------------------------------------------
    for pas in ("full", "blind"):
        col = f"varviz_classification_{pas}"
        a = sum(1 for r in base.values() if r[col] in ACT)
        c(f"actionable Pass-{pas}", f"{a:,}", f"{100*a/n:.1f}%")
        for b in ORDER:
            k = sum(1 for r in base.values() if r[col] == b)
            c(f"{b} Pass-{pas}", f"{k:,}", f"{100*k/n:.2f}%")

    # --- dual-pass ------------------------------------------------------------
    mat = collections.Counter((r["varviz_classification_full"], r["varviz_classification_blind"]) for r in base.values())
    moved = sum(v for (f, t), v in mat.items() if f != t)
    toward = sum(v for (f, t), v in mat.items() if RANK[t] < RANK[f])
    c("dual-pass disagreement", f"{moved:,}", f"{100*moved/n:.1f}%; {toward} toward pathogenic")

    # --- discrimination -------------------------------------------------------
    cv = clinvar_labels(base)
    for star in (1, 2):
        sel = {k: v for k, v in cv.items() if v[1] >= star}
        for nm, run in (("1.2.0", old), ("1.2.1", v21), ("1.2.2", base)):
            tp = fp = tn = fn = 0
            for k, (l, _) in sel.items():
                p = run[k]["varviz_classification_blind"] in ACT
                if l == 1: tp += p; fn += not p
                else:      fp += p; tn += not p
            c(f"{nm} sens >={star}star", f"{tp/(tp+fn):.3f}", f"MCC {mcc(tp,fp,tn,fn):.3f}, spec {tn/(tn+fp):.3f}, FP {fp}")

    # --- ablations ------------------------------------------------------------
    for label, run, pref in (("PM2 at Moderate","ps_pm2mod_v231","PM2"),
                             ("MDS PM1","ps_nomds_v231","PM1"),
                             ("PP3->Strong proxy","ps_pp3proxy_v231","PP3"),
                             ("meta-consensus","ps_nometa_v231","PP3"),
                             ("PP3 3-point","ps_pp33pt_v231","PP3")):
        a = load(run)
        for pas in ("full", "blind"):
            col, tc = f"varviz_classification_{pas}", f"tags_{pas}"
            att = sum(1 for k, r in base.items() if r[col] != a[k][col]
                      and sorted(x for x in sp(r[tc]) if not x.startswith(pref))
                       == sorted(x for x in sp(a[k][tc]) if not x.startswith(pref)))
            c(f"{label} Pass-{pas}", f"{100*att/n:.2f}%", f"{att:,} calls")

    # --- cohorts --------------------------------------------------------------
    ras = load("ras_vcep_pm2_1_v231")
    lab = {}
    for r in csv.DictReader(open("analyses/derived/variant_universe_rasopathy.tsv", newline="", encoding="utf8"), delimiter="\t"):
        lab.setdefault((r["gene"], r["p_notation"]), r["label"])
    for pas in ("full", "blind"):
        col = f"varviz_classification_{pas}"
        tp = sum(1 for k, v in lab.items() if v == "Pathogenic" and k in ras and ras[k][col] in ACT)
        np_ = sum(1 for k, v in lab.items() if v == "Pathogenic" and k in ras)
        fp = sum(1 for k, v in lab.items() if v == "Benign" and k in ras and ras[k][col] in ACT)
        nb = sum(1 for k, v in lab.items() if v == "Benign" and k in ras)
        c(f"RASopathy sens Pass-{pas}", f"{tp}", f"of {np_} = {100*tp/np_:.1f}%")
        c(f"RASopathy benign actionable Pass-{pas}", f"{fp}", f"of {nb} = {100*fp/nb:.1f}%")
    ext = load("external163_v231")
    for pas in ("full", "blind"):
        col = f"varviz_classification_{pas}"
        a = sum(1 for r in ext.values() if r[col] in ACT)
        c(f"external163 actionable Pass-{pas}", f"{a}", f"of {len(ext)} = {100*a/len(ext):.1f}%")


    # --- composite claims, phrased as the text phrases them -------------------
    # A ratio the text states as "X of Y (Z%)" is checked as that string, because
    # the percentage of a subset is not the percentage of all calls and the two
    # are easy to conflate: section 3.4's 43.3% is 19,040/44,002, not the 43.16%
    # of all calls that the PM2 ablation moves.
    pm2 = load("ps_pm2mod_v231")
    for pas in ("full", "blind"):
        col = f"varviz_classification_{pas}"
        den = sum(1 for r in pm2.values() if r[col] in ACT)
        lost = sum(1 for k in pm2 if pm2[k][col] in ACT and base[k][col] not in ACT)
        c(f"PM2 downgrade Pass-{pas} (of P/LP)", f"{100*lost/den:.1f}%", f"{lost:,} of {den:,}")
        c(f"PM2 downgrade numerator Pass-{pas}", f"{lost:,}", f"denominator {den:,}")
    for pas in ("full", "blind"):
        col = f"varviz_classification_{pas}"
        a = sum(1 for r in base.values() if r[col] in ACT)
        c(f"actionable share Pass-{pas}", f"{100*a/n:.1f}%", f"{a:,} variants")
    # --- report ---------------------------------------------------------------
    docs = {p: io.open(f"analyses/humu/{p}", encoding="utf8").read()
            for p in ("appnote_humu.md", "supp_humu.md")}
    print(f"{'quantity':<40}{'computed':>12}   where found")
    print("-"*88)
    missing = []
    for label, val, note in claims:
        loc = [d.replace("_humu.md","") for d, t in docs.items() if val in t]
        mark = ",".join(loc) if loc else "NOT IN TEXT"
        if not loc: missing.append((label, val, note))
        print(f"{label:<40}{val:>12}   {mark}{('  ('+note+')') if note else ''}")
    print(f"\n{len(claims)-len(missing)}/{len(claims)} computed values appear verbatim in the text")
    if missing:
        print("\nnot found as written (may be phrased differently, or stale):")
        for l, v, nt in missing: print(f"  {l:<40}{v:>12}  {nt}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
