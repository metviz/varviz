#!/usr/bin/env python3
"""How many clinically-pathogenic variants does each DMS scoreset call WT-like?

Why this exists. Panel B treats a scoreset's quartile split as pathogenicity truth:
the damaging tail is LOF and the other tail is WT-like, which the metrics read as
Benign. A protein-level assay cannot see every route to disease, so some clinically
Pathogenic variants land in the WT-like tail and are scored as benign truth. That caps
the AUROC any classifier can reach, independently of the classifier.

This measures the cap using VariBench clinical labels, which are independent of the DMS
score and of every predictor.

Where each one sits decides whether any rescoring could recover it, so the script
reports the distance from the damaging end as well as the count.

  far extreme, >=90th percentile   the assay scores it beyond wild type in the opposite
                                   direction, which unsigned scoring would recover
  near wild type, <90th            the assay cannot separate it from wild type at all,
                                   and no rescoring recovers it

An earlier version of this file asserted that PTEN R130 and BRCA1 R71 were near-wild-type
blind spots, on the mechanistic reasoning that a catalytically dead but stably folded
variant has normal abundance, and that a splicing defect is invisible to an assay on the
substituted protein. The measurement contradicts that: PTEN R130G/L/Q and R159S sit at the
90th to 99th percentile of abundance and BRCA1 R71G/K/T at the 91st to 100th, so the assay
scores them ABOVE wild type rather than merely normal. The reasoning was plausible and
wrong, and only the conclusion survives.

That does not make unsigned scoring the answer. The script also estimates its cost, and on
the labels that exist it is unfavourable: recovering 14 of 264 Pathogenic variants (5.3%)
while dragging 7 of 38 Benign variants (18.4%) into the damaging group. Where the apparent
gain is largest, PTEN with 6 recovered and 0 moved, there are ZERO benign labels, so the
cost is unmeasured rather than absent.

And scoring above wild type in an assay is not established gain-of-function disease. Only
GCK W99R has that verified: PMID 37101203 names it as a hyperinsulinemic hypoglycemia
variant that "grew faster than WT". For PTEN R130 and BRCA1 R71 the readings that suggest
themselves, stabilised-but-dead and splice-disrupting, are not gain of function, and
neither has been verified.

So the panel keeps a single signed quartile split for every gene, and this ceiling is
declared rather than transformed away. None of these is a defect in a scoreset, and none is
grounds for excluding one.

This is NOT the criterion used to exclude a scoreset. A low count here means the assay
agrees with the clinical labels; a high count can mean either a blind spot or a
bidirectional mechanism. Eligibility is decided on what the assay measures, recorded in
analyses/derived/mavedb_scoreset_polarity.tsv.

Run from project root:
  python3 analyses/check_dms_truth_ceiling.py
"""
import io, sys, collections

UNIV = "analyses/derived/variant_universe_gnomad.tsv"
POL  = "analyses/derived/mavedb_scoreset_polarity.tsv"
PATHO  = {"pathogenic", "likely pathogenic"}
BENIGN = {"benign", "likely benign"}


def main():
    pol = {}
    try:
        with io.open(POL, encoding="utf-8", newline="") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            si, di = hdr.index("study"), hdr.index("damaging_direction")
            ei = hdr.index("ineligible_reason")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(si, di, ei):
                    pol[f[si]] = (f[di], f[ei])
    except (FileNotFoundError, ValueError):
        print("polarity table missing; run analyses/determine_mavedb_polarity.py first")
        return 1

    dms = collections.defaultdict(dict)
    lab = {}
    with io.open(UNIV, encoding="utf-8", newline="") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in
              ("gene", "p_notation", "source", "study", "score_raw", "label")}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ix.values()):
                continue
            g, pn, src, stu, sr, lb = (f[ix[k]] for k in
                ("gene", "p_notation", "source", "study", "score_raw", "label"))
            if src == "MaveDB":
                try:
                    dms[stu][(g, pn)] = float(sr)
                except (TypeError, ValueError):
                    pass
            else:
                l = lb.strip().lower()
                if l in PATHO:
                    lab[(g, pn)] = "P"
                elif l in BENIGN:
                    lab[(g, pn)] = "B"

    print(f"{'gene':<8}{'study':<15}{'dir':<7}{'nP':>4}{'called WT-like':>16}{'pct':>7}"
          f"{'far':>6}{'nearWT':>7}")
    print("-" * 73)
    rows = []
    for stu in sorted(dms):
        if stu not in pol:
            continue
        direction, inel = pol[stu]
        if inel:
            continue                        # excluded on what it measures
        d = dms[stu]
        v = sorted(d.values())
        n = len(v)
        if n < 8:
            continue
        q25, q75 = v[n // 4], v[3 * n // 4]
        nP, bad, far, nearwt = 0, [], [], []
        nB, movedB = 0, 0
        for k, val in d.items():
            l = lab.get(k)
            if l is None:
                continue
            pc = 100 * sum(1 for y in v if y < val) / n
            dist = pc if direction == "lower" else 100 - pc
            if l == "B":
                nB += 1
                if dist >= 90:
                    movedB += 1        # unsigned scoring would call this damaging
                continue
            nP += 1
            wtlike = (val >= q75) if direction == "lower" else (val <= q25)
            if wtlike:
                bad.append(k[1])
                (far if dist >= 90 else nearwt).append((k[1], dist))
        if not nP:
            continue
        gene = next(iter(d))[0]
        rows.append((gene, stu, direction, nP, sorted(bad), far, nearwt, nB, movedB))
        print(f"{gene:<8}{stu:<15}{direction:<7}{nP:>4}{len(bad):>16}{100*len(bad)/nP:>6.0f}%"
              f"{len(far):>6}{len(nearwt):>7}")
        for v, dist in sorted(far):
            print(f"          {v:<10} {dist:>3.0f}%  far extreme, unsigned scoring would recover")
        for v, dist in sorted(nearwt):
            print(f"          {v:<10} {dist:>3.0f}%  near wild type, no rescoring recovers")

    print("-" * 73)
    tot_p = sum(r[3] for r in rows)
    tot_b = sum(len(r[4]) for r in rows)
    n_far = sum(len(r[5]) for r in rows)
    n_near = sum(len(r[6]) for r in rows)
    lab_b = sum(r[7] for r in rows)
    mov_b = sum(r[8] for r in rows)
    if tot_p:
        print(f"across {len(rows)} eligible scoresets: {tot_b} of {tot_p} "
              f"clinically-Pathogenic variants ({100*tot_b/tot_p:.1f}%) are scored as "
              f"WT-like truth")
        print("Every classifier is measured against that, so it is a ceiling on Panel B,")
        print("not a property of the classifier.\n")
        print(f"of those {tot_b}: {n_far} sit at the far opposite extreme, {n_near} near wild type")
        print(f"unsigned scoring would recover {n_far} of {tot_p} Pathogenic "
              f"({100*n_far/tot_p:.1f}%)")
        if lab_b:
            print(f"and would move {mov_b} of {lab_b} Benign ({100*mov_b/lab_b:.1f}%) "
                  f"into the damaging group")
            print("The cost rate exceeds the recovery rate, so the panel keeps a single")
            print("signed split and declares this ceiling instead.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
