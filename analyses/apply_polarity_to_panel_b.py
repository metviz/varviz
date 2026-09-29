#!/usr/bin/env python3
"""What the Panel B AUROC table looks like once scoreset polarity is corrected.

analyses/determine_mavedb_polarity.py establishes that 9 of the 21 scoresets are
scored in the opposite direction to the blanket lower-is-damaging convention that
analyses/07_panel_b_dms.R applies. Reversing the damaging direction swaps which
quartile is labelled loss-of-function and which WT-like, over the same variants,
so the AUROC becomes 1 - AUROC.

That correction is valid only where the variant set is unchanged. TP53
00000059-a-1 is also numbered on a domain construct offset by +101, so correcting
it changes which variants enter the join at all; its value has to be recomputed by
a re-run, not flipped, and it is reported here as pending rather than given a
number.

Run from project root:
  python3 analyses/apply_polarity_to_panel_b.py
"""
import io, re, sys, statistics as st

PANEL = "analyses/reports/panel_b_auroc_v123final.txt"
POL   = "analyses/derived/mavedb_scoreset_polarity.tsv"
OFF   = "analyses/derived/mavedb_scoreset_offsets.tsv"

def main():
    pol, off = {}, {}
    with io.open(POL, encoding="utf-8", newline="") as fh:
        hdr = fh.readline().rstrip("\n").split("\t")
        si, di = hdr.index("study"), hdr.index("damaging_direction")
        oi = hdr.index("offset_applied")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > max(si, di, oi):
                pol[f[si]] = f[di]
                off[f[si]] = int(f[oi] or 0)

    rows = []
    for line in io.open(PANEL, encoding="utf-8"):
        m = re.match(r"\s*(\S+)\s+(\d{8}-\S+)\s+(\d+)\s+([\d.]+)\s+([\d.]+)\s*$", line)
        if m:
            rows.append([m.group(1), m.group(2), int(m.group(3)),
                         float(m.group(4)), float(m.group(5))]) 
    if not rows:
        print("could not parse", PANEL); return 1

    print(f"{'gene':<8}{'study':<15}{'n':>6}{'full':>7}{'blind':>7}"
          f"{'direction':>11}{'full_corr':>11}{'blind_corr':>12}  note")
    print("-" * 94)
    full_c, blind_c, pending = [], [], []
    for gene, stu, n, fu, bl in rows:
        d = pol.get(stu, "lower")
        k = off.get(stu, 0)
        if d == "higher" and k:
            note = "PENDING re-run: polarity AND numbering offset"
            fc = bc = None
            pending.append((gene, stu))
        elif d == "higher":
            fc, bc = 1 - fu, 1 - bl
            note = "flipped: higher score is damaging"
        else:
            fc, bc = fu, bl
            note = ""
        if fc is not None:
            full_c.append(fc); blind_c.append(bc)
        fct = "--" if fc is None else f"{fc:.3f}"
        bct = "--" if bc is None else f"{bc:.3f}"
        print(f"{gene:<8}{stu:<15}{n:>6}{fu:>7.3f}{bl:>7.3f}{d:>11}{fct:>11}{bct:>12}  {note}")

    print("-" * 94)
    print(f"\nas published : Pass-Full median {st.median([r[3] for r in rows]):.3f}  "
          f"range {min(r[3] for r in rows):.3f}-{max(r[3] for r in rows):.3f}  (n={len(rows)})")
    print(f"               Pass-Blind median {st.median([r[4] for r in rows]):.3f}  "
          f"range {min(r[4] for r in rows):.3f}-{max(r[4] for r in rows):.3f}")
    print(f"\npolarity-corrected, excluding {len(pending)} pending:")
    print(f"               Pass-Full median {st.median(full_c):.3f}  "
          f"range {min(full_c):.3f}-{max(full_c):.3f}  (n={len(full_c)})")
    print(f"               Pass-Blind median {st.median(blind_c):.3f}  "
          f"range {min(blind_c):.3f}-{max(blind_c):.3f}")
    below = [r for r in rows if r[3] < 0.5]
    below_c = [f for f in full_c if f < 0.5]
    print(f"\nscoresets below 0.5 : {len(below)} as published -> {len(below_c)} after correction")
    if pending:
        print("\npending re-run (cannot be corrected by flipping):")
        for gene, stu in pending:
            print(f"  {gene} {stu}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
