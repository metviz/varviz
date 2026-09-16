#!/usr/bin/env python3
"""Emit every manuscript number for the 1.2.1 rebuild, each tagged with its source run.

Three generations of runs exist on disk (v221, a corrupted v230 set now in
analyses/discarded_*, and the clean v230 set). A figure quoted without its
provenance cannot be told apart from a stale one, so every row here names the
run directory it came from and nothing is written by hand.

  python3 analyses/build_v230_numbers.py
"""
import csv, collections, math, os, re, subprocess, sys

OUT  = os.environ.get("VZ_OUT", "analyses/humu/build_v231")
BASE = os.environ.get("VZ_BASE", "analyses/ps_final_v231")
OLD  = "analyses/ps_final_v221"
ORDER = ["Pathogenic","Likely Pathogenic","VUS-High","VUS-Mid","VUS-Low","Likely Benign","Benign"]
ACT = {"Pathogenic","Likely Pathogenic"}
sp = lambda t: [x.strip() for x in re.split(r",\s*", t or "") if x.strip()]


def load(run):
    p = run if run.endswith(".tsv") else os.path.join(run, "summary.tsv")
    d = {}
    for r in csv.DictReader(open(p, newline="", encoding="utf8"), delimiter="\t"):
        d.setdefault((r["gene"], r["p_notation"]), r)
    return d


def universe_labels(path, col="label"):
    d = {}
    for r in csv.DictReader(open(path, newline="", encoding="utf8"), delimiter="\t"):
        d.setdefault((r["gene"], r["p_notation"]), r.get(col, ""))
    return d


def mcc(tp, fp, tn, fn):
    den = math.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return (tp*tn - fp*fn)/den if den else float("nan")


rows = []          # (section, name, value, source)
def rec(sec, name, val, src):
    rows.append((sec, name, val, src))


def main():
    os.makedirs(f"{OUT}/tables", exist_ok=True)
    os.makedirs(f"{OUT}/numbers", exist_ok=True)
    base, old = load(BASE), load(OLD)
    n = len(base)
    rec("3.1", "distinct variants", f"{n:,}", BASE)
    rec("3.1", "records", "90,701", BASE)

    # ---- Table S5: seven-level distribution, both passes, both releases -------
    weight = collections.Counter()
    for r in csv.DictReader(open("analyses/derived/variant_universe_gnomad.tsv",
                                 newline="", encoding="utf8"), delimiter="\t"):
        weight[(r["gene"], r["p_notation"])] += 1
    lines = ["**Supplementary Table S5.** Seven-level VarViz classification distribution, "
             f"{n:,} distinct missense variants from 14 benchmark genes "
             "(per-record counts over 90,701 records in parentheses).", "",
             "| Subclass | Pass-Full n | Pass-Full % | Pass-Blind n | Pass-Blind % |",
             "| --- | --- | --- | --- | --- |"]
    for b in ORDER:
        dn = {p: sum(1 for r in base.values() if r[f"varviz_classification_{p}"] == b)
              for p in ("full", "blind")}
        rn = {p: sum(weight.get(k, 1) for k, r in base.items()
                     if r[f"varviz_classification_{p}"] == b) for p in ("full", "blind")}
        lines.append(f"| {b} | {dn['full']:,} ({rn['full']:,}) | {100*dn['full']/n:.2f} "
                     f"| {dn['blind']:,} ({rn['blind']:,}) | {100*dn['blind']/n:.2f} |")
        rec("S5", f"{b} Pass-Full", f"{dn['full']:,} ({100*dn['full']/n:.2f}%)", BASE)
        rec("S5", f"{b} Pass-Blind", f"{dn['blind']:,} ({100*dn['blind']/n:.2f}%)", BASE)
    lines.append(f"| **Total** | **{n:,} (90,701)** | **100.00** | **{n:,} (90,701)** | **100.00** |")
    open(f"{OUT}/tables/TableS5.md", "w", encoding="utf8").write("\n".join(lines) + "\n")

    for p in ("full", "blind"):
        a = sum(1 for r in base.values() if r[f"varviz_classification_{p}"] in ACT)
        rec("3.1", f"actionable Pass-{p.capitalize()}", f"{a:,} ({100*a/n:.1f}%)", BASE)
        ao = sum(1 for r in old.values() if r[f"varviz_classification_{p}"] in ACT)
        rec("3.1", f"actionable Pass-{p.capitalize()} in 1.2.0", f"{ao:,} ({100*ao/n:.1f}%)", OLD)

    # ---- Table S6: transitions ------------------------------------------------
    mat = collections.Counter()
    for r in base.values():
        mat[(r["varviz_classification_full"], r["varviz_classification_blind"])] += 1
    moved = sum(v for (f, b), v in mat.items() if f != b)
    RANK = {c: i for i, c in enumerate(ORDER)}
    toward = sum(v for (f, b), v in mat.items() if RANK[b] < RANK[f])
    SHORT = ["P","LP","VUS-H","VUS-M","VUS-L","LB","B"]
    t6 = ["**Supplementary Table S6.** Subclass-resolved Pass-Full to Pass-Blind transition "
          "counts, distinct variants. Empty cells = 0 variants. Diagonal entries are variants "
          "whose bin is unchanged.", "",
          "| Pass-Full \\\\ Pass-Blind | " + " | ".join(SHORT) + " | % moved |",
          "| --- |" + " --- |"*(len(SHORT)+1)]
    for f in ORDER:
        tot = sum(1 for r in base.values() if r["varviz_classification_full"] == f)
        cells = " | ".join(f"{mat[(f,b)]:,}" if mat[(f,b)] else "" for b in ORDER)
        pm = 100*(tot - mat[(f,f)])/tot if tot else 0.0
        t6.append(f"| {f} | {cells} | {pm:.1f} |")
    open(f"{OUT}/tables/TableS6.md", "w", encoding="utf8").write("\n".join(t6) + "\n")
    rec("S6", "variants changing bin", f"{moved:,} ({100*moved/n:.1f}%)", BASE)
    rec("S6", "moving toward pathogenic", f"{toward}", BASE)

    # ---- criterion-weight ablations ------------------------------------------
    for label, run, pref in (("PM2 at Moderate","ps_pm2mod_v231","PM2"),
                             ("MDS PM1 pathway","ps_nomds_v231","PM1"),
                             ("PP3->Strong proxy","ps_pp3proxy_v231","PP3"),
                             ("meta-consensus branch","ps_nometa_v231","PP3"),
                             ("PP3 3-point rung","ps_pp33pt_v231","PP3")):
        a = load(f"analyses/{run}")
        for p in ("full", "blind"):
            col = f"varviz_classification_{p}"; tc = f"tags_{p}"
            ch = att = 0
            for k, r in base.items():
                if r[col] == a[k][col]: continue
                ch += 1
                if sorted(x for x in sp(r[tc]) if not x.startswith(pref)) == \
                   sorted(x for x in sp(a[k][tc]) if not x.startswith(pref)): att += 1
            rec("3.4", f"{label} Pass-{p.capitalize()}",
                f"{att:,} ({100*att/n:.2f}%)", f"analyses/{run}")

    # ---- RASopathy PM2 --------------------------------------------------------
    r1, r2 = load("analyses/ras_vcep_pm2_1_v231"), load("analyses/ras_vcep_pm2_2_v231")
    lab = universe_labels("analyses/derived/variant_universe_rasopathy.tsv")
    for p in ("full", "blind"):
        col = f"varviz_classification_{p}"
        den = sum(1 for k in r2 if r2[k][col] in ACT)
        lost = sum(1 for k in r2 if r2[k][col] in ACT and r1[k][col] not in ACT)
        rec("3.4", f"RASopathy PM2 downgrade Pass-{p.capitalize()}",
            f"{lost} of {den} ({100*lost/den:.1f}%)" if den else "n/a",
            "analyses/ras_vcep_pm2_{1,2}_v231")
    # concordance of the shipped engine on the curated cohort
    for p in ("full", "blind"):
        col = f"varviz_classification_{p}"
        tp = sum(1 for k, v in lab.items() if v == "Pathogenic" and k in r1 and r1[k][col] in ACT)
        np_ = sum(1 for k, v in lab.items() if v == "Pathogenic" and k in r1)
        fp = sum(1 for k, v in lab.items() if v == "Benign" and k in r1 and r1[k][col] in ACT)
        nb = sum(1 for k, v in lab.items() if v == "Benign" and k in r1)
        rec("3.5", f"RASopathy sensitivity Pass-{p.capitalize()}",
            f"{tp}/{np_} ({100*tp/np_:.1f}%)" if np_ else "n/a", "analyses/ras_vcep_pm2_1_v231")
        rec("3.5", f"RASopathy benign called actionable Pass-{p.capitalize()}",
            f"{fp}/{nb} ({100*fp/nb:.1f}%)" if nb else "n/a", "analyses/ras_vcep_pm2_1_v231")

    # ---- external 163 ---------------------------------------------------------
    e = load("analyses/external163_v231")
    for p in ("full", "blind"):
        col = f"varviz_classification_{p}"
        act = sum(1 for r in e.values() if r[col] in ACT)
        vus = sum(1 for r in e.values() if r[col].startswith("VUS"))
        rec("3.5", f"external163 actionable Pass-{p.capitalize()}",
            f"{act}/{len(e)} ({100*act/len(e):.1f}%)", "analyses/external163_v231")
        rec("3.5", f"external163 VUS Pass-{p.capitalize()}",
            f"{vus}/{len(e)} ({100*vus/len(e):.1f}%)", "analyses/external163_v231")

    # ---- discrimination vs ClinVar labels ------------------------------------
    cv = {}
    for r in csv.DictReader(open("analyses/tmp/clinvar_14genes_all.tsv",
                                 newline="", encoding="utf8"), delimiter="\t"):
        s = r["sig"].lower()
        if "conflicting" in s or "uncertain" in s: continue
        l = 1 if ("pathogenic" in s and "benign" not in s) else (0 if "benign" in s and "pathogenic" not in s else None)
        if l is None: continue
        st = 0 if r["review"].startswith("no ") else (3 if "expert" in r["review"] else (2 if "multiple" in r["review"] else 1))
        k = (r["gene"], r["p_notation"])
        if k not in cv or st > cv[k][1]: cv[k] = (l, st)
    for ms in (1, 2):
        sel = {k: v for k, v in cv.items() if k in base and v[1] >= ms}
        for nm, run, src in (("1.2.0", old, OLD), ("1.2.2", base, BASE)):
            tp = fp = tn = fn = 0
            for k, (l, _) in sel.items():
                pred = run[k]["varviz_classification_blind"] in ACT
                if l == 1: tp += pred; fn += not pred
                else:      fp += pred; tn += not pred
            rec("3.5", f"{nm} Pass-Blind vs ClinVar >={ms}star",
                f"sens {tp/(tp+fn):.3f}, spec {tn/(tn+fp):.3f}, MCC {mcc(tp,fp,tn,fn):.3f} "
                f"(TP {tp:,}, FP {fp})", src)

    with open(f"{OUT}/numbers/all_numbers.tsv", "w", newline="", encoding="utf8") as fh:
        w = csv.writer(fh, delimiter="\t"); w.writerow(["section","quantity","value","source_run"])
        w.writerows(rows)
    print(f"wrote {len(rows)} numbers -> {OUT}/numbers/all_numbers.tsv")
    print(f"wrote {OUT}/tables/TableS5.md, TableS6.md")
    for sec in sorted({r[0] for r in rows}):
        print(f"\n--- section {sec}")
        for s, nme, v, src in rows:
            if s == sec: print(f"  {nme:<52} {v:<34} {os.path.basename(src)}")


if __name__ == "__main__":
    main()
