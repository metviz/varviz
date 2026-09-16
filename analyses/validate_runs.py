#!/usr/bin/env python3
"""Flag genes whose evidence collapsed in one run but not the others.

A harness run fetches evidence per gene. When a fetch fails partway and
VARVIZ_ALLOW_DEGRADED is set, the gene is still checkpointed: the row count is
correct, the file looks complete, and the criteria that depended on the failed
source are simply absent. Nothing in the run reports it. One such gene was found
only because its bin distribution looked wrong three steps downstream.

Comparing runs catches it cheaply. An ablation changes one criterion, so for
every OTHER criterion a gene's tagged count should match the baseline almost
exactly. A gene where an untargeted criterion moves by more than a few percent
did not change because of the option; it lost evidence.

  python3 analyses/validate_runs.py --baseline analyses/ps_final_v230 RUN [RUN ...]
"""
import argparse, collections, csv, os, re, sys

PREFIXES = ["PM1", "PM2", "PM5", "PP2", "PP3", "PP5", "PS1", "BP4", "BS1", "BA1", "PM4", "PP1", "PP4"]
split = lambda t: [x.strip() for x in re.split(r",\s*", t or "") if x.strip()]


def profile(run):
    """{gene: {prefix: count}} over distinct variants, plus per-gene variant count."""
    f = os.path.join(run, "summary.tsv")
    if not os.path.isfile(f):
        return None, None
    seen, counts, n = set(), collections.defaultdict(collections.Counter), collections.Counter()
    for r in csv.DictReader(open(f, newline="", encoding="utf8"), delimiter="\t"):
        k = (r["gene"], r["p_notation"])
        if k in seen:
            continue
        seen.add(k)
        n[r["gene"]] += 1
        tags = split(r["tags_full"])
        for p in PREFIXES:
            if any(t.startswith(p) for t in tags):
                counts[r["gene"]][p] += 1
    return counts, n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--baseline", required=True)
    ap.add_argument("runs", nargs="+")
    # NOTE on coupling: since 1.2.2 the Pathogenic-boundary cap is decided on a
    # variant's total points, so any ablation that changes points can flip a cap
    # decision and with it a PP3 tag. That is a real downstream effect, not
    # contamination, and PP3 must be added to --ignore for such ablations. On the
    # RASopathy cohort, moving PM2 from 1 to 2 points changed 11 PP3 tags this
    # way: at 11 points reverting the raise drops the variant below Pathogenic so
    # the cap fires, at 12 it does not and the cap correctly declines.
    ap.add_argument("--tolerance", type=float, default=2.0,
                    help="percent of a gene's variants an untargeted criterion may move (default 2)")
    ap.add_argument("--ignore", default="",
                    help="comma-separated prefixes the run legitimately targets, e.g. PP3,BP4")
    a = ap.parse_args()

    base, bn = profile(a.baseline)
    if base is None:
        sys.exit(f"baseline has no summary.tsv: {a.baseline}")
    ignore = {p.strip() for p in a.ignore.split(",") if p.strip()}
    bad = 0
    for run in a.runs:
        prof, n = profile(run)
        name = os.path.basename(run.rstrip("/"))
        if prof is None:
            print(f"{name:<22} no summary.tsv (still running?)")
            continue
        hits = []
        for gene in sorted(base):
            if n.get(gene, 0) != bn[gene]:
                hits.append(f"    {gene}: {n.get(gene,0):,} variants vs {bn[gene]:,} in baseline")
                continue
            for p in PREFIXES:
                if p in ignore:
                    continue
                d = prof[gene][p] - base[gene][p]
                if bn[gene] and abs(d) > a.tolerance / 100.0 * bn[gene]:
                    hits.append(f"    {gene} {p}: {prof[gene][p]:,} vs {base[gene][p]:,} "
                                f"({d:+,}, {100.0*d/bn[gene]:+.1f}% of the gene)")
        if hits:
            bad += 1
            print(f"{name:<22} SUSPECT")
            for h in hits:
                print(h)
        else:
            print(f"{name:<22} ok")
    print(f"\n{bad} run(s) with untargeted criteria moving more than {a.tolerance}% of a gene")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
