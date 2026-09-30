#!/usr/bin/env python3
"""How many clinically-pathogenic variants does each DMS scoreset call WT-like?

Why this exists. Panel B treats a scoreset's quartile split as pathogenicity truth:
the damaging tail is LOF and the other tail is WT-like, which the metrics read as
Benign. A protein-level assay cannot see every route to disease, so some clinically
Pathogenic variants land in the WT-like tail and are scored as benign truth. That caps
the AUROC any classifier can reach, independently of the classifier.

This measures the cap using VariBench clinical labels, which are independent of the DMS
score and of every predictor.

The reason matters more than the count, and it differs per case:

  GCK activity        W99R is a gain-of-function variant causing hyperinsulinemic
                      hypoglycemia. PMID 37101203 names it as an HH variant that "grew
                      faster than WT", and the activity distribution there spans "from
                      total loss of function to above-wildtype activity". The assay is
                      right; a single signed quartile split cannot represent a gene with
                      both loss- and gain-of-function disease.
  PTEN VAMP-seq       R130 is the catalytic residue. VAMP-seq measures abundance, and a
                      catalytically dead but stably folded variant has normal abundance.
  BRCA1 HDR and M2H   R71 substitutions disrupt splicing. An assay on the substituted
                      protein cannot see a splicing defect.

None of these is a defect in the scoreset, and none is grounds for excluding one. They
are the ceiling of DMS as pathogenicity truth, and worth stating rather than absorbing
into an unexplained AUROC.

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

    print(f"{'gene':<8}{'study':<15}{'dir':<7}{'nP':>4}{'called WT-like':>16}{'pct':>7}")
    print("-" * 60)
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
        nP, bad = 0, []
        for k, val in d.items():
            if lab.get(k) != "P":
                continue
            nP += 1
            wtlike = (val >= q75) if direction == "lower" else (val <= q25)
            if wtlike:
                bad.append(k[1])
        if not nP:
            continue
        gene = next(iter(d))[0]
        rows.append((gene, stu, direction, nP, sorted(bad)))
        print(f"{gene:<8}{stu:<15}{direction:<7}{nP:>4}{len(bad):>16}{100*len(bad)/nP:>6.0f}%")
        if bad:
            print(f"          {', '.join(sorted(bad))}")

    print("-" * 60)
    tot_p = sum(r[3] for r in rows)
    tot_b = sum(len(r[4]) for r in rows)
    if tot_p:
        print(f"across {len(rows)} eligible scoresets: {tot_b} of {tot_p} "
              f"clinically-Pathogenic variants ({100*tot_b/tot_p:.1f}%) are scored as "
              f"WT-like truth")
        print("Every classifier is measured against that, so it is a ceiling on Panel B,")
        print("not a property of the classifier.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
