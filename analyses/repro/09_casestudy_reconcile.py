#!/usr/bin/env python3
"""Reconcile VarViz app exports against the case-study values in the manuscript.

The app writes p-notation in one-letter form (p.E46K) even when three-letter is
entered, so both forms are normalised before matching.

Usage:
    python3 analyses/repro/09_casestudy_reconcile.py EXPORT.tsv [EXPORT2.tsv ...]
"""
import csv, sys, os

AA = {'Ala':'A','Arg':'R','Asn':'N','Asp':'D','Cys':'C','Gln':'Q','Glu':'E','Gly':'G',
      'His':'H','Ile':'I','Leu':'L','Lys':'K','Met':'M','Phe':'F','Pro':'P','Ser':'S',
      'Thr':'T','Trp':'W','Tyr':'Y','Val':'V','Ter':'*'}

def norm(v):
    """p.Ala53Thr and p.A53T both -> A53T"""
    s = v.strip()
    if s.startswith('p.'): s = s[2:]
    for three, one in AA.items():
        s = s.replace(three, one)
    return s.upper()

CLASS_ALIAS = {'Likely P':'Likely Pathogenic', 'Likely B':'Likely Benign',
               'LP':'Likely Pathogenic', 'P':'Pathogenic'}

def ncls(c):
    """Table S1 abbreviates 'Likely P'; the engine writes it in full."""
    c = (c or '').strip()
    return CLASS_ALIAS.get(c, c)

def load_claims(p):
    return list(csv.DictReader(open(p), delimiter='\t'))

def load_exports(paths):
    out = {}
    for p in paths:
        for r in csv.DictReader(open(p), delimiter='\t'):
            v = r.get('Variant') or r.get('variant') or ''
            if not v: continue
            out[norm(v)] = r
    return out

def main(paths):
    claims = load_claims('analyses/repro/00_casestudy_claims.tsv')
    exp = load_exports(paths)
    print('export rows: %d distinct variants across %d file(s)\n' % (len(exp), len(paths)))
    hdr = '%-8s %-14s %-19s %-19s %6s %6s  %-28s %s'
    print(hdr % ('gene','variant','manuscript','engine','pts','pts','pm1_pathway','status'))
    print('-'*140)
    miss = mism = ok = absent = 0
    for c in claims:
        key = norm(c['variant'])
        r = exp.get(key)
        if r is None:
            print(hdr % (c['gene'], c['variant'], c['classification'] or '-', 'NOT EXPORTED',
                         c['points'] or '-', '-', '-', 'missing')); absent += 1; continue
        eng_c = r.get('Final_Classification','')
        eng_p = r.get('Final_Points','')
        pw    = r.get('ACMG_PM1_Pathway','')
        want_c, want_p = c['classification'], c['points']
        if not want_c and not want_p:
            st = 'reference'; ok += 1
        elif ncls(eng_c) == ncls(want_c) and str(eng_p) == str(want_p):
            st = 'match'; ok += 1
        else:
            st = 'CHANGED' if str(eng_p) != str(want_p) else 'class-label only'
            mism += 1 if st == 'CHANGED' else 0
            ok += 0 if st == 'CHANGED' else 1
        print(hdr % (c['gene'], c['variant'], want_c or '-', eng_c, want_p or '-', eng_p, pw[:28], st))
    print('-'*140)
    print('%d match/reference, %d CHANGED, %d not in export' % (ok, mism, absent))
    # pathway summary -- did the ungated MDS actually engage?
    pws = [r.get('ACMG_PM1_Pathway','') for r in exp.values()]
    corr = sum(1 for p in pws if '+mds' in p and 'hotspot' not in p)
    orig = sum(1 for p in pws if p == 'mds')
    upg  = sum(1 for p in pws if 'hotspot_upgrade' in p)
    print('\npathway check: %d corroborating (+mds), %d originating (mds), %d hotspot upgrades' % (corr, orig, upg))
    if corr == 0 and upg == 0:
        print('  WARNING: no composite pathways -- the export may come from the pre-change engine.')

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(2)
    for p in sys.argv[1:]:
        if not os.path.exists(p): print('no such file:', p); sys.exit(2)
    main(sys.argv[1:])
