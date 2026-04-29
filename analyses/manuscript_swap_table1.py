"""
Replace the old Table 1 (comparator table) in v2 document.xml with the new
API/data-source Table 1 per implementation doc Tasks 5+6 / design doc 1.1.
"""

import re
import sys

DOC = "/tmp/v2_unpacked/word/document.xml"

with open(DOC, "r", encoding="utf-8") as f:
    xml = f.read()

# ---- Locate the OLD Table 1 block: caption paragraph + table ----
# Caption paragraph contains the literal "Table 1. Comparison of VarViz with related variant"
# Find that <w:p>...</w:p> block.
caption_marker = "Table 1. Comparison of VarViz with related variant"
caption_idx = xml.find(caption_marker)
if caption_idx < 0:
    sys.exit("ERROR: caption marker not found")

# Walk back to the most recent <w:p> open
p_open_re = re.compile(r"<w:p\b[^>]*>")
p_open = None
for m in p_open_re.finditer(xml[:caption_idx]):
    p_open = m
caption_start = p_open.start()

# Find the </w:p> after caption + the <w:tbl>...</w:tbl> after that
caption_end_re = re.compile(r"</w:p>")
m_caption_close = caption_end_re.search(xml, caption_idx)
caption_end = m_caption_close.end()

m_tbl_open = re.search(r"<w:tbl\b", xml, flags=0)
# Find next <w:tbl after caption_end
m_tbl_open = re.search(r"<w:tbl\b", xml[caption_end:])
tbl_start_abs = caption_end + m_tbl_open.start()

m_tbl_close = re.search(r"</w:tbl>", xml[tbl_start_abs:])
tbl_end_abs = tbl_start_abs + m_tbl_close.end()

print(f"Old Table 1 block: bytes {caption_start} .. {tbl_end_abs}", file=sys.stderr)
print(f"  Block length: {tbl_end_abs - caption_start} chars", file=sys.stderr)

# ---- Build NEW Table 1 ----
ROWS = [
    (
        "Protein sequence + domains",
        "UniProt",
        "REST /uniprotkb/{acc}.json",
        "live release",
        "PM1 (domain pathway); visualization",
    ),
    (
        "Structural confidence (pLDDT)",
        "AlphaFold (EMBL-EBI)",
        "/files/AF-{acc}-F1-confidence_v4.json",
        "AFDB v4",
        "Visualization track",
    ),
    (
        "Population allele frequency",
        "gnomAD",
        "GraphQL",
        "v4.1, GRCh38, MANE Select",
        "PM2, BA1, BS1, BS2",
    ),
    (
        "Clinical assertions",
        "NCBI ClinVar",
        "E-utilities (esearch + esummary)",
        "live",
        "PS1, PM5, PP5, BP6; PM1 (hotspot pathway)",
    ),
    (
        "In silico predictors + AlphaMissense",
        "MyVariant.info / dbNSFP",
        "REST /v1/query",
        "dbNSFP v4.4 (REVEL, AlphaMissense, MetaSVM, …)",
        "PP3, BP4",
    ),
    (
        "Multi-species conservation",
        "UCSC Genome Browser",
        "REST /getData/track",
        "hg38; PhyloP 100/470-way; PhastCons",
        "PP3 (supporting)",
    ),
    (
        "Transcript / coordinate lookup",
        "Ensembl",
        "REST /lookup/symbol, /lookup/id",
        "live (GRCh38)",
        "MANE Select resolution",
    ),
    (
        "Gene–disease validity",
        "ClinGen",
        "KB lookup by HGNC ID",
        "live",
        "PP2, BP1 (gene-level prior)",
    ),
    (
        "Gene–disease validity (broad)",
        "GenCC",
        "/api/v1/validity-prop",
        "live (v1 API)",
        "PP2, BP1; ClinGen fallback",
    ),
]
HEADERS = (
    "Evidence layer",
    "Source",
    "Endpoint / API",
    "Version / Build",
    "ACMG/AMP role",
)

# Per implementation doc: total 9360 DXA across 1800/1300/2360/1700/2200
COL_DXA = (1800, 1300, 2360, 1700, 2200)


def cell(text, dxa, bold=False):
    bprops = '<w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:sz w:val="20"/><w:szCs w:val="20"/>{}</w:rPr>'
    rpr = bprops.format("<w:b/><w:bCs/>" if bold else "")
    return (
        f"<w:tc>"
        f'<w:tcPr><w:tcW w:w="{dxa}" w:type="dxa"/>'
        f'<w:tcMar><w:top w:w="60" w:type="dxa"/><w:left w:w="100" w:type="dxa"/><w:bottom w:w="60" w:type="dxa"/><w:right w:w="100" w:type="dxa"/></w:tcMar>'
        f"</w:tcPr>"
        f'<w:p><w:pPr><w:pStyle w:val="Normal"/><w:spacing w:lineRule="exact" w:line="240" w:before="0" w:after="0"/><w:rPr>{rpr[8:-8] if bold else ""}<w:sz w:val="20"/></w:rPr></w:pPr>'
        f'<w:r>{rpr}<w:t xml:space="preserve">{text}</w:t></w:r>'
        f"</w:p>"
        f"</w:tc>"
    )


def row(cells_text, bold=False):
    return (
        "<w:tr>"
        + "".join(cell(t, dxa, bold) for t, dxa in zip(cells_text, COL_DXA))
        + "</w:tr>"
    )


new_table_xml = (
    "<w:tbl>"
    "<w:tblPr>"
    '<w:tblW w:w="9360" w:type="dxa"/>'
    '<w:tblLayout w:type="fixed"/>'
    "<w:tblBorders>"
    '<w:top w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
    '<w:left w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
    '<w:bottom w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
    '<w:right w:val="single" w:sz="6" w:space="0" w:color="000000"/>'
    '<w:insideH w:val="single" w:sz="4" w:space="0" w:color="808080"/>'
    '<w:insideV w:val="single" w:sz="4" w:space="0" w:color="808080"/>'
    "</w:tblBorders>"
    '<w:tblCellMar><w:top w:w="60" w:type="dxa"/><w:left w:w="100" w:type="dxa"/><w:bottom w:w="60" w:type="dxa"/><w:right w:w="100" w:type="dxa"/></w:tblCellMar>'
    '<w:tblLook w:val="04A0" w:firstRow="1" w:lastRow="0" w:firstColumn="1" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/>'
    "</w:tblPr>"
    "<w:tblGrid>"
    + "".join(f'<w:gridCol w:w="{d}"/>' for d in COL_DXA)
    + "</w:tblGrid>"
    + row(HEADERS, bold=True)
    + "".join(row(r) for r in ROWS)
    + "</w:tbl>"
)

# Caption paragraph (above table)
caption_xml = (
    "<w:p>"
    '<w:pPr><w:pStyle w:val="Normal"/><w:spacing w:lineRule="exact" w:line="276" w:before="240" w:after="60"/>'
    '<w:rPr><w:b/><w:sz w:val="20"/></w:rPr></w:pPr>'
    '<w:r><w:rPr><w:b/><w:bCs/><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr>'
    '<w:t xml:space="preserve">Table 1. </w:t></w:r>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr>'
    '<w:t xml:space="preserve">Live API data sources queried by VarViz, with corresponding ACMG/AMP roles.</w:t></w:r>'
    "</w:p>"
)


# Footnote paragraph (below table) — base URLs with superscript letters
def super_run(letter, url):
    return (
        f'<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="18"/><w:szCs w:val="18"/><w:vertAlign w:val="superscript"/></w:rPr><w:t xml:space="preserve">{letter}</w:t></w:r>'
        f'<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr><w:t xml:space="preserve">{url}  </w:t></w:r>'
    )


footnote_xml = (
    "<w:p>"
    '<w:pPr><w:pStyle w:val="Normal"/><w:spacing w:lineRule="exact" w:line="240" w:before="60" w:after="120"/>'
    '<w:rPr><w:sz w:val="18"/></w:rPr></w:pPr>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="18"/><w:szCs w:val="18"/><w:b/></w:rPr><w:t xml:space="preserve">Base URLs: </w:t></w:r>'
    + super_run("a", "https://rest.uniprot.org")
    + super_run("b", "https://alphafold.ebi.ac.uk")
    + super_run("c", "https://gnomad.broadinstitute.org/api")
    + super_run("d", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/")
    + super_run("e", "https://myvariant.info")
    + super_run("f", "https://api.genome.ucsc.edu")
    + super_run("g", "https://rest.ensembl.org")
    + super_run("h", "https://search.clinicalgenome.org/kb/genes/{HGNC:N}")
    + super_run("i", "https://thegencc.org (query: ?hgnc_id=HGNC:{N})")
    + "</w:p>"
)

new_block = caption_xml + new_table_xml + footnote_xml

# Splice
new_xml = xml[:caption_start] + new_block + xml[tbl_end_abs:]
with open(DOC, "w", encoding="utf-8") as f:
    f.write(new_xml)

# Validate
import xml.etree.ElementTree as ET

ET.parse(DOC)
print(f"OK. Wrote {DOC}", file=sys.stderr)
print(
    f"Document size: {len(xml)} -> {len(new_xml)} chars (diff {len(new_xml) - len(xml):+d})",
    file=sys.stderr,
)
