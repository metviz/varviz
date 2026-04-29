"""Insert 7 bib entries (2 Pass 1 + 5 Pass 2) into SK-style v2 docx, sort + renumber."""
import re, sys

DOC = "/tmp/v2_unpacked/word/document.xml"

NEW_ENTRIES = [
    "Li M.M. and Wang K. (2017) InterVar: clinical interpretation of genetic variants by the 2015 ACMG-AMP guidelines. Am. J. Hum. Genet., 100, 267-280.",
    "Wang J. et al. (2017) MARRVEL: integration of human and model organism genetic resources to facilitate functional annotation of the human genome. Am. J. Hum. Genet., 100, 843-853.",
    "Esposito,D. et al. (2019) MaveDB: an open-source platform to distribute and interpret data from multiplexed assays of variant effect. Genome Biol., 20, 223.",
    "Giacomelli,A.O. et al. (2018) Mutational processes shape the landscape of TP53 mutations in human cancer. Nat. Genet., 50, 1381-1387.",
    "Gudkov,M. et al. (2025) Confounder-aware estimation of constraint with CAPS. Bioinform. Adv., online ahead of print.",
    "Kotler,E. et al. (2018) A systematic p53 mutation library links differential functional impact to cancer mutation pattern and evolutionary conservation. Mol. Cell, 71, 178-190.e8.",
    "Yip,Y.L. et al. (2008) The Swiss-Prot variant page and the ModSNP database in the era of computational human genetics. Hum. Mutat., 29, 361-366.",
]

xml = open(DOC, encoding="utf-8").read()

# Match each reference paragraph with this SK-doc style: <w:p ...><w:pPr>...spacing before=40 after=40 line=276...</w:pPr><w:r><w:rPr><w:color w:val="000000"/></w:rPr><w:t>NN. Author...</w:t></w:r></w:p>
pattern = re.compile(
    r'<w:p\b[^>]*>\s*<w:pPr>\s*<w:spacing w:before="40" w:after="40" w:line="276" w:lineRule="exact"/>\s*<w:jc w:val="both"/>\s*</w:pPr>\s*<w:r>\s*<w:rPr>\s*<w:color w:val="000000"/>\s*</w:rPr>\s*<w:t>(\d+)\.\s+([^<]+)</w:t>\s*</w:r>\s*</w:p>',
    re.DOTALL,
)
matches = list(pattern.finditer(xml))
print(f"Found {len(matches)} reference paragraphs", file=sys.stderr)
if not matches:
    sys.exit("ERROR: pattern did not match")

existing = [m.group(2).strip() for m in matches]
all_sorted = sorted(existing + NEW_ENTRIES, key=lambda s: s.lower())

def make_p(num, body):
    return (
        '<w:p w14:paraId="00000000" w14:textId="77777777" w:rsidR="00A10375" w:rsidRDefault="00000000">'
        '<w:pPr><w:spacing w:before="40" w:after="40" w:line="276" w:lineRule="exact"/><w:jc w:val="both"/></w:pPr>'
        f'<w:r><w:rPr><w:color w:val="000000"/></w:rPr><w:t>{num}. {body}</w:t></w:r>'
        '</w:p>'
    )

first, last = matches[0], matches[-1]
new_block = "\n    ".join(make_p(i+1, b) for i, b in enumerate(all_sorted))
new_xml = xml[:first.start()] + new_block + xml[last.end():]
open(DOC, "w", encoding="utf-8").write(new_xml)

import xml.etree.ElementTree as ET
ET.parse(DOC)
print(f"OK. {len(matches)} -> {len(all_sorted)} entries", file=sys.stderr)
print("Inserted positions of new entries:", file=sys.stderr)
for new_e in NEW_ENTRIES:
    pos = all_sorted.index(new_e) + 1
    print(f"  #{pos:2d}: {new_e[:60]}...", file=sys.stderr)
