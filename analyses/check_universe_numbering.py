#!/usr/bin/env python3
"""Does the variant universe's p_notation use UniProt canonical numbering?

Why this exists. 10,570 of the 90,701 rows in
analyses/ps_am792_v123_final/master_all_variants.tsv carry no AlphaMissense score
from either source. Spot-checking three of them showed why: BRCA1 p.A101H,
p.A101L and p.A107C are absent from the P38398 AlphaMissense CSV because UniProt
P38398 position 101 is Tyr and position 107 is Phe, not Ala. So those variants are
not numbered on the UniProt canonical sequence the CSV, the AlphaFold model and
the plotted landscape all use.

All three come from MaveDB, study urn:mavedb:00000003-a-2, and MaveDB scoresets
are numbered on the scoreset's own target sequence, which may be a fragment or a
different isoform.

The test compares, per gene and per source, the reference amino acid in each
variant's p_notation against the reference amino acid at the same position in the
gene's UniProt canonical sequence. The AlphaMissense per-accession CSV supplies
the reference sequence: it enumerates every substitution of every residue, so the
reference residue at each position is read straight off the variant identifiers.

A mismatching set is then searched for a constant offset, because a fragment
construct numbered from its own residue 1 shows up as one shift applied to every
variant, whereas a genuinely wrong gene or accession does not.

Run from project root:
  python3 analyses/test_universe_numbering.py
"""
import io, re, sys, collections

UID = {"BAP1":"Q92560","BRCA1":"P38398","CASR":"P41180","DDR2":"Q16832",
       "GCK":"P35557","KCNH2":"Q12809","KCNQ1":"P51787","KRAS":"P01116",
       "LDLR":"P01130","NUDT15":"Q9NV35","PTEN":"P60484","SLC13A5":"Q86YT5",
       "TP53":"P04637","TSHR":"P16473"}
AM_DIR = "analyses/raw/alphamissense"
UNIV   = "analyses/derived/variant_universe_gnomad.tsv"
VAR_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

def uniprot_ref(uid):
    """position -> reference residue, from the AlphaMissense substitution CSV."""
    path = f"{AM_DIR}/{uid}-F1-aa-substitutions.csv"
    try:
        # The files use bare CR line endings in places, so split on either.
        raw = io.open(path, encoding="utf-8", newline="").read()
    except FileNotFoundError:
        return None
    ref = {}
    for line in re.split(r"[\r\n]+", raw):
        m = VAR_RE.match(line.split(",", 1)[0])
        if m:
            ref[int(m.group(2))] = m.group(1)
    return ref or None

def best_offset(pairs, ref, lo=-1200, hi=1200):
    """The single shift that reconciles the most variants, and its hit rate."""
    best, best_n = None, -1
    for k in range(lo, hi + 1):
        n = sum(1 for p, r in pairs if ref.get(p + k) == r)
        if n > best_n:
            best, best_n = k, n
    return best, best_n

def main():
    with io.open(UNIV, encoding="utf-8", newline="") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        gi, pi, si = hdr.index("gene"), hdr.index("p_notation"), hdr.index("source")
        rows = collections.defaultdict(list)
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) <= max(gi, pi, si):
                continue
            rows[f[gi]].append((f[pi], f[si]))

    print(f"{'gene':<9}{'source':<12}{'n':>7}{'ref_ok':>8}{'ref_bad':>8}"
          f"{'past_end':>9}{'%ok':>7}   offset fix")
    print("-" * 78)
    grand = collections.Counter()
    suspect = []
    for gene in sorted(rows):
        ref = uniprot_ref(UID.get(gene, ""))
        if ref is None:
            print(f"{gene:<9}{'(no AlphaMissense CSV -- cannot test)':<40}")
            continue
        L = max(ref)
        by_src = collections.defaultdict(list)
        for p_not, src in rows[gene]:
            by_src[src].append(p_not)
        for src in sorted(by_src):
            ok = bad = past = 0
            badpairs = []
            for p_not in by_src[src]:
                m = VAR_RE.match(p_not[2:] if p_not.startswith("p.") else p_not)
                if not m:
                    continue
                r, pos = m.group(1), int(m.group(2))
                if pos > L:
                    past += 1
                    continue
                if ref.get(pos) == r:
                    ok += 1
                else:
                    bad += 1
                    badpairs.append((pos, r))
            n = ok + bad + past
            if not n:
                continue
            grand["ok"] += ok; grand["bad"] += bad; grand["past"] += past
            note = ""
            if badpairs and bad > 5:
                k, hit = best_offset(badpairs, ref)
                note = (f"  shift {k:+d} fixes {hit}/{bad} "
                        f"({100*hit/bad:.1f}%)")
                if hit / bad > 0.9:
                    suspect.append((gene, src, k, bad, hit, L))
            print(f"{gene:<9}{src:<12}{n:>7}{ok:>8}{bad:>8}{past:>9}"
                  f"{100*ok/n:6.1f}%{note}")

    tot = grand["ok"] + grand["bad"] + grand["past"]
    print("-" * 78)
    print(f"universe total: {tot}  ref_ok {grand['ok']} ({100*grand['ok']/tot:.1f}%)  "
          f"ref_bad {grand['bad']} ({100*grand['bad']/tot:.1f}%)  "
          f"past_end {grand['past']} ({100*grand['past']/tot:.1f}%)")

    if suspect:
        print("\nSets a single constant shift reconciles (>90% of mismatches):")
        print(f"  {'gene':<9}{'source':<12}{'shift':>7}{'mismatched':>12}"
              f"{'fixed':>8}{'uniprot_len':>13}")
        for gene, src, k, bad, hit, L in suspect:
            print(f"  {gene:<9}{src:<12}{k:>+7d}{bad:>12}{hit:>8}{L:>13}")
        print("\nA constant shift means the set is numbered on a different sequence,")
        print("not that individual identifiers are wrong.")
    else:
        print("\nNo mismatching set is explained by a single constant shift.")

    # Fails when a material share of the universe is off-sequence, so this can
    # gate a re-run rather than being read by eye.
    bad_share = grand["bad"] / tot if tot else 0
    print(f"\n{'FAIL' if bad_share > 0.01 else 'PASS'}: "
          f"{100*bad_share:.2f}% of the universe does not match UniProt numbering")
    return 1 if bad_share > 0.01 else 0

if __name__ == "__main__":
    sys.exit(main())
