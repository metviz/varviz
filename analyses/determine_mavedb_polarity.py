#!/usr/bin/env python3
"""Score polarity for every MaveDB scoreset in the benchmark universe.

Why this exists. analyses/07_panel_b_dms.R:14-19 applies one blanket convention
to all 21 scoresets -- bottom quartile = loss of function, top quartile = WT-like
-- and says so explicitly: "Convention: lower score = more damaging, the dominant
convention across MaveDB function/fitness assays. Inverted assays (e.g.
abundance-as-stability where higher = more stable) would require per-URN sign
correction." It also notes "there is no URN -> paper mapping in the ingest yet".
Nothing in the repository records a polarity per scoreset, so an inverted assay
reports its AUROC as 1 - the true value.

Nine of the 21 scoresets sit below 0.5 in analyses/reports/panel_b_auroc_v123final.txt,
which is what an inverted scale looks like.

Polarity is decided here from CLINICAL labels, not from agreement with a
computational predictor. Orienting a truth scale by whether it agrees with the
predictor being evaluated, and then reporting that predictor's AUROC against it,
guarantees a good number; it would manufacture the result. The universe carries
608 Pathogenic and 136 Benign VariBench labels that are independent of both the
DMS measurements and every predictor, and many of them also have a DMS score.

For each scoreset the DMS score is asked to separate those labels under the
declared convention that LOWER is more damaging. An AUROC below 0.5 means the
scale runs the other way.

The established numbering offsets are applied before the join. Without them the
TP53 scoreset matches only 5 clinical variants, by coincidence of a construct
position equalling a UniProt position, rather than the 100+ it really shares.

Publication identifiers were each verified against NCBI esummary, and the BRCA1
preprint DOI against Crossref, in the session that wrote this file.

The six KRAS scoresets are the one group whose direction could not be settled from
the MaveDB record alone: they list no publication there, and a free-energy change
has no self-evident sign. The paper was identified by an exact match between the
MaveDB experiment title and PMID 38109937, corroborated by its abstract and by the
six interaction partners matching the six scoresets one to one, and its full text
then settles the sign directly: mutations that increase the binding free energy
are reported as enriched at the binding interface, and distal mutations that
increase it are the inhibitory ones the paper set out to map. Increase therefore
means weaker binding. Where the
assay semantics alone settle the direction, that is recorded too, so a scoreset
with too few clinical anchors still gets a documented answer rather than a guess.

Run from project root:
  python3 analyses/determine_mavedb_polarity.py
"""
import io, sys, collections

UNIV    = "analyses/derived/variant_universe_gnomad.tsv"
OFFSETS = "analyses/derived/mavedb_scoreset_offsets.tsv"
OUT     = "analyses/derived/mavedb_scoreset_polarity.tsv"

# Curated from the MaveDB API and the verified publications. assay_high_means is
# what a HIGH score means physically, which is what fixes the direction.
#   pub      : identifier verified this session (PubMed via esummary, DOI via Crossref)
#   role     : assay_paper, method_paper or scale_definition -- a method or scale
#              citation does not describe the assay and cannot settle polarity
CURATED = {
 "00000003-a-2":  ("BRCA1", "PubMed:28784151", "method_paper",
                   "BRCA1 RING autoubiquitination (E3 ligase) DMS, Enrich2 scores",
                   "more ligase activity", "lower",
                   "assay readout: more E3 ligase activity is healthier"),
 "00000081-a-1":  ("BRCA1", "PubMed:30219179", "assay_paper",
                   "homology-directed DNA repair, DEPLETION score",
                   "more depletion, i.e. less repair", "higher",
                   "title states DEPLETION; clinical labels agree (AUROC 0.104, 6P/4B)"),
 "00000081-a-2":  ("BRCA1", "PubMed:30219179", "assay_paper",
                   "homology-directed DNA repair, mean HDR FLUORESCENCE",
                   "more repair", "lower",
                   "title states FLUORESCENCE; clinical labels agree (AUROC 0.776, 13P/11B)"),
 "00000093-a-1":  ("BRCA1", "DOI:10.1101/092619", "assay_paper",
                   "BRCA1-BARD1 RING mammalian 2-hybrid, reported as BASE_ACMG POINTS",
                   "stronger ACMG evidence FOR pathogenicity", "higher",
                   "score column documented as BASE_ACMG points; positive points are evidence FOR pathogenicity"),
 "00000096-a-1":  ("GCK",   "PubMed:37101203", "assay_paper",
                   "GCK activity by yeast functional complementation",
                   "more enzyme activity", "lower",
                   "assay readout: more enzyme activity is healthier"),
 "00000096-b-1":  ("GCK",   "PubMed:38627865", "assay_paper",
                   "GCK abundance by DHFR-PCA",
                   "more protein abundance", "lower",
                   "assay readout: more protein abundance is healthier"),
 "00001216-a-1":  ("KCNH2", "PubMed:39315434", "assay_paper",
                   "KCNH2 cell-surface trafficking MAVE",
                   "more trafficking to the surface", "lower",
                   "assay readout: more surface trafficking is healthier; clinical labels agree (AUROC 0.959, 41P/6B)"),
 "00000094-a-14": ("KCNQ1", "PubMed:35760561", "assay_paper",
                   "KCNQ1 normalized currents, homozygous state",
                   "more potassium current", "lower",
                   "assay readout: more potassium current is healthier"),
 "00000094-a-15": ("KCNQ1", "PubMed:35760561", "assay_paper",
                   "KCNQ1 normalized currents, homozygous state",
                   "more potassium current", "lower",
                   "assay readout: more potassium current is healthier"),
 "00000055-0-1":  ("NUDT15","PubMed:32094176", "assay_paper",
                   "NUDT15 combined activity and stability",
                   "more activity and stability", "lower",
                   "assay readout: more activity and stability is healthier"),
 "00000055-a-1":  ("NUDT15","PubMed:29785012", "method_paper",
                   "NUDT15 protein stability by VAMP-seq",
                   "more protein abundance", "lower",
                   "assay readout: more protein abundance is healthier"),
 "00000055-b-1":  ("NUDT15","PubMed:32094176", "assay_paper",
                   "NUDT15 activity by thiopurine-cytotoxicity screen",
                   "more enzyme activity", "lower",
                   "assay readout: more enzyme activity is healthier"),
 "00000101-a-1":  ("PTEN",  "PubMed:34649609", "assay_paper",
                   "PTEN variant abundance by VAMP-seq (fill-in)",
                   "more protein abundance", "lower",
                   "assay readout: more protein abundance is healthier"),
 "00000102-0-1":  ("PTEN",  "PubMed:34649609", "assay_paper",
                   "PTEN variant abundance by VAMP-seq (combined)",
                   "more protein abundance", "lower",
                   "assay readout: more protein abundance is healthier"),
 "00000059-a-1":  ("TP53",  "PubMed:29979965", "assay_paper",
                   "p53 relative fitness by cell growth selection, DNA-binding domain",
                   "faster growth, which for a tumour suppressor means LOST function",
                   "higher",
                   "tumour suppressor under growth selection, so LOF enriches; PMID 29979965 states the readout is relative growth rate"),
 # Six binding free-energy scoresets from one experiment, one partner each.
 **{u: ("KRAS", "PubMed:38109937", "assay_paper",
        f"KRAS binding free energy change (ddG) against {p}",
        "larger ddG, i.e. weaker binding", "higher",
        "PMID 38109937 full text: mutations that increase the binding free energy are "
        "enriched at the binding interface, so increase means inhibition")
    for u, p in (("00000115-a-1", "DARPin K55"), ("00000115-a-2", "DARPin K27"),
                 ("00000115-a-3", "SOS1"),       ("00000115-a-4", "PIK3CG RBD"),
                 ("00000115-a-5", "RALGDS RBD"), ("00000115-a-6", "RAF1 RBD"))},
}

PATHO = {"pathogenic", "likely pathogenic", "p", "lp"}
BENIGN = {"benign", "likely benign", "b", "lb"}


def auroc(pos, neg):
    """P(a random positive scores below a random negative) under lower-is-damaging."""
    if not pos or not neg:
        return None
    vals = sorted([(v, 1) for v in pos] + [(v, 0) for v in neg], key=lambda t: t[0])
    n = len(vals)
    i = 0
    ranksum = 0.0
    while i < n:
        j = i
        while j + 1 < n and vals[j + 1][0] == vals[i][0]:
            j += 1
        avg = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            if vals[k][1] == 1:
                ranksum += avg
        i = j + 1
    n1, n2 = len(pos), len(neg)
    u = ranksum - n1 * (n1 + 1) / 2.0
    # Pathogenic is the positive class and should rank LOW, so invert the statistic
    # to read as "the declared lower-is-damaging convention holds".
    return 1.0 - u / (n1 * n2)


def load_offsets(path=OFFSETS):
    off = {}
    try:
        with io.open(path, encoding="utf-8", newline="") as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            si, oi = hdr.index("study"), hdr.index("offset")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) > max(si, oi):
                    off[f[si]] = int(f[oi])
    except (FileNotFoundError, ValueError):
        pass
    return off


def main():
    off = load_offsets()
    if not off:
        print("NOTE: no offset table; run python3 analyses/check_mavedb_offsets.py first.\n")
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
            elif lb.strip().lower() in PATHO | BENIGN:
                lab[(g, pn)] = lb.strip().lower()

    def shift(pn, k):
        """p.X123Y -> p.X<123+k>Y, so the join runs in the reference's numbering."""
        if not k:
            return pn
        body = pn[2:] if pn.startswith("p.") else pn
        i = 1
        while i < len(body) and body[i].isdigit():
            i += 1
        j = i
        while j < len(body) and body[j].isdigit():
            j += 1
        ref, num, alt = body[0], body[1:len(body) - 1], body[-1]
        if not num.isdigit():
            return pn
        return f"p.{ref}{int(num) + k}{alt}"

    rows = []
    print(f"{'gene':<8}{'study':<15}{'off':>5}{'nP':>4}{'nB':>4}"
          f"{'AUROC_lower_dmg':>17}{'clinical call':>16}  {'assay call':<12} agree")
    print("-" * 104)
    for stu in sorted(dms):
        cur = CURATED.get(stu)
        gene = cur[0] if cur else (next(iter(dms[stu]))[0] if dms[stu] else "?")
        k = off.get(stu, 0)
        pos, neg = [], []
        for (g, pn), v in dms[stu].items():
            key = (g, shift(pn, k))
            lb = lab.get(key)
            if lb in PATHO:
                pos.append(v)
            elif lb in BENIGN:
                neg.append(v)
        a = auroc(pos, neg)
        # Two or fewer of either class cannot orient anything.
        if a is None or min(len(pos), len(neg)) < 3:
            clin, a_txt = "insufficient", "--"
        else:
            clin = "lower" if a >= 0.5 else "higher"
            a_txt = f"{a:.3f}"
        assay = cur[5] if cur else "unknown"
        agree = "-" if clin == "insufficient" else ("yes" if clin == assay else "NO <<<")
        print(f"{gene:<8}{stu:<15}{k:>5}{len(pos):>4}{len(neg):>4}"
              f"{a_txt:>17}{clin:>16}  {assay:<12} {agree}")
        rows.append((gene, stu, cur, k, len(pos), len(neg), a, clin, assay))

    import os
    os.makedirs("analyses/derived", exist_ok=True)
    with io.open(OUT, "w", encoding="utf-8") as fh:
        fh.write("gene\tstudy\tpublication\tpublication_role\tassay\thigh_score_means"
                 "\tdamaging_direction\tbasis\tdirection_evidence\tn_pathogenic\tn_benign"
                 "\tauroc_lower_is_damaging\toffset_applied\n")
        for gene, stu, cur, k, npos, nneg, a, clin, assay in rows:
            if cur:
                pub, role, desc, high, ev = cur[1], cur[2], cur[3], cur[4], cur[6]
            else:
                pub = role = desc = high = ev = "UNKNOWN"
            basis = ("clinical_labels_and_assay" if clin == assay and clin != "insufficient"
                     else "assay_semantics" if clin == "insufficient"
                     else "CONFLICT_review_required")
            fh.write(f"{gene}\t{stu}\t{pub}\t{role}\t{desc}\t{high}\t{assay}\t{basis}"
                     f"\t{ev}\t{npos}\t{nneg}\t{'' if a is None else f'{a:.4f}'}\t{k}\n")
    print("-" * 104)
    print(f"wrote {OUT}")

    inverted = [r for r in rows if r[8] == "higher"]
    conflicts = [r for r in rows if r[7] not in ("insufficient",) and r[7] != r[8]]
    print(f"\nScoresets where HIGHER score is the damaging direction, i.e. the blanket")
    print(f"lower-is-damaging convention inverts them ({len(inverted)} of {len(rows)}):")
    for gene, stu, cur, k, npos, nneg, a, clin, assay in inverted:
        print(f"  {gene:<8}{stu:<15}{cur[3] if cur else ''}")
    if conflicts:
        print(f"\nCONFLICT between clinical labels and assay semantics ({len(conflicts)}):")
        for gene, stu, cur, k, npos, nneg, a, clin, assay in conflicts:
            print(f"  {gene:<8}{stu:<15}clinical={clin} assay={assay} "
                  f"AUROC={a:.3f} nP={npos} nB={nneg}")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
