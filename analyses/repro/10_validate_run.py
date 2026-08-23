#!/usr/bin/env python3
"""Detect silent data-source corruption in a harness run.

Four fetches degrade instead of raising, and every one of them leaves the row
count correct, so row-count checks see nothing:

  1. UniProt gene_info   "[GeneInfo] FAILED for <ACC>"
  2. UniProt Pfam        "[Pfam] FAILED for <ACC>"  (older builds:
                         "UniProt API request failed" -- no accession, which is
                         why a single grep never caught both)
  3. Ensembl exons       a NA chrom silently disables the dbNSFP and REVEL
                         preloads; visible only as "local caches: dbNSFP=NA"
  4. ClinVar esummary    a throttled batch has no $result and was skipped by a
                         bare `next`; LDLR extracted 203 of 1731 records this way

The earlier version of this script read only run.log, matched only vector 1,
could not tell a 13-gene summary from a 14-gene one (it once passed a summary
containing a single gene), and called any deviation from the reference runs a
failure -- when the reference runs were themselves the contaminated ones.

So: INTRINSIC checks decide pass/fail. Reference comparison is advisory only,
because a reference is not a ground truth -- ps_mds_corroborate lost KCNQ1 and
ps_mds_consfree lost KCNH2/KRAS/KCNQ1/PTEN/TSHR.

Usage:
  python3 analyses/repro/10_validate_run.py RUN_DIR [REF_SUMMARY ...]
"""

import csv, collections, glob, os, re, sys

UNIVERSE = "analyses/derived/variant_universe_gnomad.tsv"
DBNSFP = "analyses/tmp/dbNSFPv4.9a_hg38_custombuild.bgz"

FAIL_PATTERNS = [
    ("uniprot_geneinfo", re.compile(r"\[GeneInfo\] FAILED for")),
    ("uniprot_pfam", re.compile(r"\[Pfam\] FAILED for|UniProt API request failed")),
    ("localpred", re.compile(r"LOCAL PREDICTORS FAILED TO LOAD")),
    ("clinvar_batch", re.compile(r"\[ClinVar\] BATCH FAILED")),
    (
        "abort",
        re.compile(r"ABORT: required data source failed|ABORT: UniProt fetch failed"),
    ),
]
CACHE_NA = re.compile(r"\[(\w+)\] local caches: dbNSFP=NA")
EXTRACT = re.compile(r"\[ClinVar\] Extracted (\d+) variants \((\d+) missense")
FOUND = re.compile(r"\[ClinVar\] Found (\d+) variant IDs")
# Healthy extracted/found ratios across a verified 14-gene run span 25.6%-63.4%
# (SLC13A5 lowest, GCK highest). The LDLR run whose esummary batches were
# throttled sat at 4.3% (203 of 1731 records from 4734 IDs). 15% separates them
# with margin. This is the only ClinVar-truncation signal available in logs
# written before the BATCH FAILED marker existed; it catches severe truncation
# only -- PTEN's milder loss (965 vs 1123) ran at 29.7% and is invisible here,
# which is what the BATCH FAILED instrumentation is for going forward.
CLINVAR_MIN_RATIO = 15.0
TOL = 5.0


def load(p):
    d = {}
    for r in csv.DictReader(open(p, errors="replace"), delimiter="\t"):
        d.setdefault((r["gene"], r["p_notation"]), r)
    return d


def shares(d):
    tot, hot, una = collections.Counter(), collections.Counter(), collections.Counter()
    for (g, _), r in d.items():
        tot[g] += 1
        b = (r.get("pm1_pathway") or "").split("+")[0]
        if b == "clinvar_hotspot":
            hot[g] += 1
        if b == "mds_unavailable":
            una[g] += 1
    return {g: (100 * hot[g] / tot[g], 100 * una[g] / tot[g], tot[g]) for g in tot}


def scan_logs(run):
    """Every log in the run dir, not just run.log. Sharded runs put each gene's
    failures in its own shard_<GENE>.log, and the concat pass overwrites run.log."""
    logs = sorted(glob.glob(os.path.join(run, "*.log")))
    hits = collections.defaultdict(list)
    cache_na, extracted, found = [], {}, {}
    for p in logs:
        if os.sep + "_bogus" + os.sep in p:
            continue
        for line in open(p, errors="replace"):
            for name, rx in FAIL_PATTERNS:
                if rx.search(line):
                    hits[name].append((os.path.basename(p), line.strip()[:100]))
            m = CACHE_NA.search(line)
            if m:
                cache_na.append((os.path.basename(p), m.group(1)))
            g = (
                os.path.basename(p)[6:-4]
                if os.path.basename(p).startswith("shard_")
                else None
            )
            m = FOUND.search(line)
            if m and g:
                found[g] = int(m.group(1))
            m = EXTRACT.search(line)
            if m and g:
                extracted[g] = (int(m.group(1)), int(m.group(2)))
    return logs, hits, cache_na, extracted, found


def expected_universe():
    if not os.path.exists(UNIVERSE):
        return None, None
    genes, rows = set(), 0
    for r in csv.DictReader(open(UNIVERSE, errors="replace"), delimiter="\t"):
        genes.add(r["gene"])
        rows += 1
    return genes, rows


def main(run, refs):
    fail = []
    truncated = []
    print("=== INTRINSIC CHECKS (decide pass/fail) ===")

    s = os.path.join(run, "summary.tsv")
    if not os.path.exists(s):
        print("  summary.tsv MISSING - run incomplete")
        return 1
    cur_rows = load(s)
    cur = shares(cur_rows)

    exp_genes, exp_rows = expected_universe()
    if exp_genes:
        missing = exp_genes - set(cur)
        extra = set(cur) - exp_genes
        print("  genes: %d/%d present" % (len(cur), len(exp_genes)))
        if missing:
            fail.append("missing genes: %s" % ", ".join(sorted(missing)))
        if extra:
            fail.append("unexpected genes: %s" % ", ".join(sorted(extra)))
        n_raw = sum(1 for _ in open(s, errors="replace")) - 1
        print(
            "  rows: %d raw / %d distinct (universe %d)"
            % (n_raw, len(cur_rows), exp_rows)
        )
        if n_raw != exp_rows:
            fail.append("row count %d != universe %d" % (n_raw, exp_rows))
    else:
        print("  universe file absent - cannot verify completeness")

    ck = os.path.join(run, "classifications")
    if os.path.isdir(ck):
        n_ck = len(glob.glob(os.path.join(ck, "*__dual.tsv")))
        print("  checkpoints: %d" % n_ck)
        if exp_genes and n_ck != len(exp_genes):
            fail.append("checkpoints %d != genes %d" % (n_ck, len(exp_genes)))

    logs, hits, cache_na, extracted, found = scan_logs(run)
    print("  logs scanned: %d" % len(logs))
    if not logs:
        fail.append("no logs found - cannot audit data sources")
    for name, _ in FAIL_PATTERNS:
        if hits[name]:
            fail.append("%s failures: %d" % (name, len(hits[name])))
            print("  %-18s %d hit(s)" % (name, len(hits[name])))
            for f, l in hits[name][:3]:
                print("       %s: %s" % (f, l))
        else:
            print("  %-18s clean" % name)

    if os.path.exists(DBNSFP):
        if cache_na:
            fail.append(
                "local predictors not loaded for: %s"
                % ", ".join(sorted({g for _, g in cache_na}))
            )
            print(
                "  local caches      dbNSFP=NA for: %s"
                % ", ".join(sorted({g for _, g in cache_na}))
            )
        else:
            print("  local caches      all loaded")

    print(
        "\n=== ADVISORY: vs reference runs (references may themselves be corrupt) ==="
    )
    ref = {
        os.path.basename(os.path.dirname(p)): shares(load(p))
        for p in refs
        if os.path.exists(p)
    }
    print("%-9s %8s %9s %9s   %s" % ("gene", "n", "hotspot%", "mds_una%", "vs refs"))
    for g in sorted(cur, key=lambda x: -cur[x][2]):
        h, u, n = cur[g]
        cmp = []
        for nm, rs in ref.items():
            if g in rs:
                cmp.append("%s %+.1f" % (nm[:16], h - rs[g][0]))
        note = ""
        if ref and cmp and all(abs(float(c.split()[-1])) > TOL for c in cmp):
            note = "   <- differs from all refs; confirm which run is correct"
        print("%-9s %8d %8.1f%% %8.1f%%   %s%s" % (g, n, h, u, "  ".join(cmp), note))

    if extracted:
        print("\nClinVar records extracted per gene (truncation check):")
        for g in sorted(extracted):
            n_ex = extracted[g][0]
            n_fo = found.get(g)
            ratio = (100.0 * n_ex / n_fo) if n_fo else None
            mark = ""
            if ratio is not None and ratio < CLINVAR_MIN_RATIO:
                mark = "  <== TRUNCATED"
                truncated.append("%s (%d of %d IDs, %.1f%%)" % (g, n_ex, n_fo, ratio))
            print("  %-9s %5d variants (%d missense)%s%s"
                  % (g, n_ex, extracted[g][1],
                     "" if ratio is None else "  from %d IDs = %.1f%%" % (n_fo, ratio),
                     mark))

    if truncated:
        fail.append("ClinVar truncation: %s" % "; ".join(truncated))

    print()
    if fail:
        print("VALIDATION FAILED - do not use this run")
        for f in fail:
            print("  - %s" % f)
        return 1
    print("VALIDATION PASSED (intrinsic checks clean; review advisory section)")
    return 0


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    sys.exit(main(sys.argv[1], sys.argv[2:]))
