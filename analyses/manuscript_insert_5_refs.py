"""
Insert 5 new bibliography entries into v2 document.xml in alphabetical order
and renumber surrounding entries. Operates on the unpacked /tmp/v2_unpacked/.
"""

import re
import sys

DOC = "/tmp/v2_unpacked/word/document.xml"

NEW_ENTRIES = [
    "Esposito,D. et al. (2019) MaveDB: an open-source platform to distribute and interpret data from multiplexed assays of variant effect. Genome Biol., 20, 223.",
    "Giacomelli,A.O. et al. (2018) Mutational processes shape the landscape of TP53 mutations in human cancer. Nat. Genet., 50, 1381-1387.",
    "Gudkov,M. et al. (2025) Confounder-aware estimation of constraint with CAPS. Bioinform. Adv., online ahead of print.",
    "Kotler,E. et al. (2018) A systematic p53 mutation library links differential functional impact to cancer mutation pattern and evolutionary conservation. Mol. Cell, 71, 178-190.e8.",
    "Yip,Y.L. et al. (2008) The Swiss-Prot variant page and the ModSNP database in the era of computational human genetics. Hum. Mutat., 29, 361-366.",
]

with open(DOC, "r", encoding="utf-8") as f:
    xml = f.read()

# Match each reference paragraph: <w:p>...<w:t>N. Author,X. et al. ... </w:t></w:r></w:p>
# Use a regex that captures the whole <w:p>..</w:p> block with a numbered entry.
# Reference pattern: <w:t>NN. <Author>... — we extract those blocks in order.
pattern = re.compile(
    r'<w:p>\s*<w:pPr>\s*<w:pStyle w:val="Normal"/>\s*<w:spacing w:lineRule="exact" w:line="276" w:before="40" w:after="40"/>\s*<w:jc w:val="both"/>\s*<w:rPr>\s*<w:sz w:val="24"/>\s*<w:szCs w:val="24"/>\s*</w:rPr>\s*</w:pPr>\s*<w:r>\s*<w:rPr>\s*<w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/>\s*<w:color w:val="000000"/>\s*<w:sz w:val="24"/>\s*<w:szCs w:val="24"/>\s*</w:rPr>\s*<w:t>(\d+)\.\s+([^<]+)</w:t>\s*</w:r>\s*</w:p>',
    re.DOTALL,
)

matches = list(pattern.finditer(xml))
print(f"Found {len(matches)} reference paragraphs", file=sys.stderr)
if len(matches) == 0:
    sys.exit("ERROR: pattern did not match any reference paragraph")

# Existing entries: list of (number, body_text)
entries = [(int(m.group(1)), m.group(2).strip()) for m in matches]
sort_key = lambda body: body.lower()

# Build full sorted list with new entries inserted
all_entries = [body for (_n, body) in entries] + NEW_ENTRIES
all_entries_sorted = sorted(all_entries, key=sort_key)


# Build replacement block — same paragraph styling as the originals
def make_p(num, body):
    return (
        "<w:p>"
        "<w:pPr>"
        '<w:pStyle w:val="Normal"/>'
        '<w:spacing w:lineRule="exact" w:line="276" w:before="40" w:after="40"/>'
        '<w:jc w:val="both"/>'
        "<w:rPr>"
        '<w:sz w:val="24"/>'
        '<w:szCs w:val="24"/>'
        "</w:rPr>"
        "</w:pPr>"
        "<w:r>"
        "<w:rPr>"
        '<w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/>'
        '<w:color w:val="000000"/>'
        '<w:sz w:val="24"/>'
        '<w:szCs w:val="24"/>'
        "</w:rPr>"
        f"<w:t>{num}. {body}</w:t>"
        "</w:r>"
        "</w:p>"
    )


# Replace the entire span from first match to last match with the new sorted block.
first, last = matches[0], matches[-1]
prefix = xml[: first.start()]
suffix = xml[last.end() :]
new_block = "\n    ".join(
    make_p(i + 1, body) for i, body in enumerate(all_entries_sorted)
)

new_xml = prefix + new_block + suffix
with open(DOC, "w", encoding="utf-8") as f:
    f.write(new_xml)

print(
    f"Wrote {DOC}: {len(entries)} -> {len(all_entries_sorted)} entries", file=sys.stderr
)
print("First 5 entries:", file=sys.stderr)
for i, body in enumerate(all_entries_sorted[:5]):
    print(f"  {i + 1}. {body[:80]}...", file=sys.stderr)
print("Last 3 entries:", file=sys.stderr)
for i, body in enumerate(all_entries_sorted[-3:]):
    n = len(all_entries_sorted) - 3 + i + 1
    print(f"  {n}. {body[:80]}...", file=sys.stderr)
