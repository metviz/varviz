#!/usr/bin/env python3
"""Per-scoreset numbering offsets, and what correcting them does to the
DMS-versus-predictor agreement.

analyses/test_universe_numbering.py established that 8,048 of 87,974 universe
variants (9.15%) are not numbered on the UniProt canonical sequence, and that two
MaveDB sets are each reconciled by one constant shift: TP53 +101 (fixes
2,754/2,754) and BRCA1 +1 (fixes 3,420/3,429).

Two things follow that this script measures.

First, the offset belongs to a scoreset, not a gene, so it is resolved per MaveDB
study URN here. A gene with several scoresets can mix numbering.

Second, if a scoreset's variants are numbered on a different sequence, then
joining them to a predictor by p_notation pairs a DMS measurement for one residue
with a prediction for another, which drives any rank agreement toward zero. The
test is therefore the absolute Spearman correlation between the DMS score and
AlphaMissense, before and against the shifted position. Absolute, because the
score polarity of these scoresets is not recorded anywhere in this repository, and
the magnitude is what the offset affects; the sign is a separate question.

AlphaMissense is read from the per-accession AlphaFold CSV, which is keyed by
UniProt accession and has no transcript dimension, so this measurement is
independent of the MANE transcript work.

Run from project root:
  python3 analyses/test_mavedb_offset_effect.py
"""
import io, re, sys, collections, math

UID = {"BAP1":"Q92560","BRCA1":"P38398","CASR":"P41180","DDR2":"Q16832",
       "GCK":"P35557","KCNH2":"Q12809","KCNQ1":"P51787","KRAS":"P01116",
       "LDLR":"P01130","NUDT15":"Q9NV35","PTEN":"P60484","SLC13A5":"Q86YT5",
       "TP53":"P04637","TSHR":"P16473"}
AM_DIR = "analyses/raw/alphamissense"
UNIV   = "analyses/derived/variant_universe_gnomad.tsv"
VAR_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

def load_am(uid):
    """(ref, pos, alt) -> score, plus pos -> ref."""
    try:
        raw = io.open(f"{AM_DIR}/{uid}-F1-aa-substitutions.csv",
                      encoding="utf-8", newline="").read()
    except FileNotFoundError:
        return None, None
    sc, ref = {}, {}
    for line in re.split(r"[\r\n]+", raw):
        parts = line.split(",")
        if len(parts) < 2:
            continue
        m = VAR_RE.match(parts[0])
        if not m:
            continue
        try:
            v = float(parts[1])
        except ValueError:
            continue
        r, p, a = m.group(1), int(m.group(2)), m.group(3)
        sc[(r, p, a)] = v
        ref[p] = r
    return (sc, ref) if sc else (None, None)

def spearman(xs, ys):
    n = len(xs)
    if n < 10:
        return None
    def ranks(v):
        order = sorted(range(n), key=lambda i: v[i])
        r = [0.0] * n
        i = 0
        while i < n:
            j = i
            while j + 1 < n and v[order[j + 1]] == v[order[i]]:
                j += 1
            avg = (i + j) / 2.0 + 1.0
            for k in range(i, j + 1):
                r[order[k]] = avg
            i = j + 1
        return r
    rx, ry = ranks(xs), ranks(ys)
    mx, my = sum(rx) / n, sum(ry) / n
    num = sum((a - mx) * (b - my) for a, b in zip(rx, ry))
    dx = math.sqrt(sum((a - mx) ** 2 for a in rx))
    dy = math.sqrt(sum((b - my) ** 2 for b in ry))
    return num / (dx * dy) if dx and dy else None

def main():
    with io.open(UNIV, encoding="utf-8", newline="") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        ix = {k: hdr.index(k) for k in ("gene", "p_notation", "source", "study", "score_raw")}
        sets = collections.defaultdict(list)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(ix.values()):
                continue
            if f[ix["source"]] != "MaveDB":
                continue
            sets[(f[ix["gene"]], f[ix["study"]])].append((f[ix["p_notation"]], f[ix["score_raw"]]))

    print(f"{'gene':<8}{'study':<18}{'n':>6}{'ref_ok%':>9}{'shift':>7}{'fixed%':>8}"
          f"{'|rho| as-is':>13}{'|rho| shifted':>15}")
    print("-" * 92)
    flagged = []
    stats = []
    for (gene, study) in sorted(sets):
        uid = UID.get(gene)
        if not uid:
            continue
        am, ref = load_am(uid)
        if am is None:
            continue
        recs = []
        for p_not, sraw in sets[(gene, study)]:
            m = VAR_RE.match(p_not[2:] if p_not.startswith("p.") else p_not)
            if not m:
                continue
            try:
                sc = float(sraw)
            except (TypeError, ValueError):
                sc = None
            recs.append((m.group(1), int(m.group(2)), m.group(3), sc))
        if not recs:
            continue
        ok = sum(1 for r, p, a, _ in recs if ref.get(p) == r)
        bad = [(r, p, a) for r, p, a, _ in recs if ref.get(p) != r]
        shift, fixed = 0, 0
        if bad:
            best, best_n = 0, -1
            for k in range(-1200, 1201):
                n = sum(1 for r, p, a in bad if ref.get(p + k) == r)
                if n > best_n:
                    best, best_n = k, n
            shift, fixed = best, best_n

        def rho_at(k):
            xs, ys = [], []
            for r, p, a, sc in recs:
                if sc is None:
                    continue
                v = am.get((r, p + k, a))
                if v is None:
                    continue
                xs.append(sc); ys.append(v)
            rr = spearman(xs, ys)
            return (rr, len(xs))

        r0, n0 = rho_at(0)
        r1, n1 = rho_at(shift) if shift else (r0, n0)
        f0 = f"{abs(r0):.3f} (n={n0})" if r0 is not None else f"-- (n={n0})"
        f1 = f"{abs(r1):.3f} (n={n1})" if r1 is not None else f"-- (n={n1})"
        print(f"{gene:<8}{study:<18}{len(recs):>6}{100*ok/len(recs):8.1f}%"
              f"{shift:>7d}{(100*fixed/len(bad) if bad else 0):7.1f}%{f0:>13}{f1:>15}")
        keep = shift if (bad and fixed / len(bad) > 0.9 and shift != 0) else 0
        stats.append((gene, study, keep, len(recs), ok, shift, len(bad), fixed))
        if keep:
            flagged.append((gene, study, shift, len(recs), ok, r0, r1, n0, n1))

    import os
    os.makedirs("analyses/derived", exist_ok=True)
    out = "analyses/derived/mavedb_scoreset_offsets.tsv"
    # offset is the shift to APPLY (0 where none is warranted). best_shift and
    # best_shift_fixes record what the search found even where it was rejected,
    # so a set like KRAS -- a tenth off UniProt numbering with no single shift
    # that reconciles it -- is visible rather than silently recorded as clean.
    with io.open(out, "w", encoding="utf-8") as fh:
        fh.write("gene\tstudy\toffset\tn_variants\ton_uniprot_asis\tpct_asis"
                 "\tn_mismatched\tbest_shift\tbest_shift_fixes\n")
        for gene, study, keep, n, ok, shift, nbad, fixed in sorted(stats):
            fh.write(f"{gene}\t{study}\t{keep}\t{n}\t{ok}\t{100*ok/n:.1f}"
                     f"\t{nbad}\t{shift}\t{fixed}\n")
    print(f"\nwrote {out}")

    if not flagged:
        print("No scoreset needs a shift.")
        return 0
    print("\nScoresets numbered on a different sequence:")
    for gene, study, shift, n, ok, r0, r1, n0, n1 in flagged:
        print(f"\n  {gene} {study}")
        print(f"    variants                     : {n}")
        print(f"    on UniProt numbering as-is   : {ok} ({100*ok/n:.1f}%)")
        print(f"    constant shift that fixes it : {shift:+d}")
        if r0 is not None and r1 is not None:
            print(f"    |Spearman| vs AlphaMissense  : {abs(r0):.3f} (n={n0}) "
                  f"-> {abs(r1):.3f} (n={n1}) after the shift")
            if abs(r0) < 0.2 <= abs(r1):
                print(f"    the as-is agreement is near zero and the shift restores it,")
                print(f"    so the join was pairing different residues")
    print(f"\nFAIL: {len(flagged)} scoreset(s) are not on UniProt numbering")
    return 1

if __name__ == "__main__":
    sys.exit(main())
