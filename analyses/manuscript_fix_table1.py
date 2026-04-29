"""Replace the broken new Table 1 with a corrected version (no tblBorders, clean cell rPr)."""
import re, sys

DOC = "/tmp/v2_unpacked/word/document.xml"
xml = open(DOC, encoding="utf-8").read()

# Find the broken table block: from caption "Live API data sources queried by VarViz"
# back to nearest <w:p> open, forward through </w:tbl> + footnote </w:p>.
caption_idx = xml.find("Live API data sources queried by VarViz")
if caption_idx < 0:
    sys.exit("ERROR: caption marker not found")

# Walk back to caption <w:p> open
p_open_re = re.compile(r"<w:p[> ]")
p_open = None
for m in p_open_re.finditer(xml[:caption_idx]):
    p_open = m
caption_start = p_open.start()

# Walk forward to find footnote close (after table). Footnote contains "Base URLs:"
footnote_idx = xml.find("Base URLs:", caption_idx)
m_footnote_close = re.search(r"</w:p>", xml[footnote_idx:])
footnote_end = footnote_idx + m_footnote_close.end()

print(f"Old (broken) Table 1 block: bytes {caption_start} .. {footnote_end}", file=sys.stderr)

ROWS = [
    ("Protein sequence + domains",        "UniProt",                  "REST /uniprotkb/{acc}.json",                       "live release",                                  "PM1 (domain pathway); visualization"),
    ("Structural confidence (pLDDT)",     "AlphaFold (EMBL-EBI)",     "/files/AF-{acc}-F1-confidence_v4.json",            "AFDB v4",                                       "Visualization track"),
    ("Population allele frequency",       "gnomAD",                   "GraphQL",                                          "v4.1, GRCh38, MANE Select",                     "PM2, BA1, BS1, BS2"),
    ("Clinical assertions",               "NCBI ClinVar",             "E-utilities (esearch + esummary)",                 "live",                                          "PS1, PM5, PP5, BP6; PM1 (hotspot pathway)"),
    ("In silico predictors + AlphaMissense", "MyVariant.info / dbNSFP", "REST /v1/query",                                "dbNSFP v4.4 (REVEL, AlphaMissense, MetaSVM)",   "PP3, BP4"),
    ("Multi-species conservation",        "UCSC Genome Browser",      "REST /getData/track",                              "hg38; PhyloP 100/470-way; PhastCons",           "PP3 (supporting)"),
    ("Transcript / coordinate lookup",    "Ensembl",                  "REST /lookup/symbol, /lookup/id",                  "live (GRCh38)",                                 "MANE Select resolution"),
    ("Gene-disease validity",             "ClinGen",                  "KB lookup by HGNC ID",                             "live",                                          "PP2, BP1 (gene-level prior)"),
    ("Gene-disease validity (broad)",     "GenCC",                    "/api/v1/validity-prop",                            "live (v1 API)",                                 "PP2, BP1; ClinGen fallback"),
]
HEADERS = ("Evidence layer", "Source", "Endpoint / API", "Version / Build", "ACMG/AMP role")
COL_DXA = (1800, 1300, 2360, 1700, 2200)

def cell(text, dxa, bold):
    rpr = '<w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/>'
    if bold:
        rpr += '<w:b/><w:bCs/>'
    return (
        f'<w:tc>'
        f'<w:tcPr><w:tcW w:w="{dxa}" w:type="dxa"/></w:tcPr>'
        f'<w:p>'
        f'<w:pPr><w:spacing w:lineRule="exact" w:line="240" w:before="0" w:after="0"/><w:jc w:val="left"/></w:pPr>'
        f'<w:r><w:rPr>{rpr}</w:rPr><w:t xml:space="preserve">{text}</w:t></w:r>'
        f'</w:p>'
        f'</w:tc>'
    )

def row(cells_text, bold):
    return '<w:tr>' + ''.join(cell(t, dxa, bold) for t, dxa in zip(cells_text, COL_DXA)) + '</w:tr>'

# Match the original SK comparator table's tblPr structure (no tblBorders block)
new_table_xml = (
    '<w:tbl>'
    '<w:tblPr>'
    '<w:tblW w:w="9360" w:type="dxa"/>'
    '<w:tblLayout w:type="fixed"/>'
    '<w:tblCellMar><w:top w:w="60" w:type="dxa"/><w:left w:w="100" w:type="dxa"/><w:bottom w:w="60" w:type="dxa"/><w:right w:w="100" w:type="dxa"/></w:tblCellMar>'
    '<w:tblLook w:val="04A0" w:firstRow="1" w:lastRow="0" w:firstColumn="1" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/>'
    '</w:tblPr>'
    '<w:tblGrid>' + ''.join(f'<w:gridCol w:w="{d}"/>' for d in COL_DXA) + '</w:tblGrid>'
    + row(HEADERS, bold=True)
    + ''.join(row(r, bold=False) for r in ROWS)
    + '</w:tbl>'
)

caption_xml = (
    '<w:p>'
    '<w:pPr><w:spacing w:lineRule="exact" w:line="276" w:before="240" w:after="60"/><w:jc w:val="left"/></w:pPr>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/><w:b/><w:bCs/></w:rPr><w:t xml:space="preserve">Table 1. </w:t></w:r>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr><w:t xml:space="preserve">Live API data sources queried by VarViz, with corresponding ACMG/AMP roles.</w:t></w:r>'
    '</w:p>'
)

def super_run(letter, url):
    return (
        f'<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="18"/><w:szCs w:val="18"/><w:vertAlign w:val="superscript"/></w:rPr><w:t xml:space="preserve">{letter}</w:t></w:r>'
        f'<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr><w:t xml:space="preserve">{url}  </w:t></w:r>'
    )

footnote_xml = (
    '<w:p>'
    '<w:pPr><w:spacing w:lineRule="exact" w:line="240" w:before="60" w:after="120"/><w:jc w:val="left"/></w:pPr>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="18"/><w:szCs w:val="18"/><w:b/><w:bCs/></w:rPr><w:t xml:space="preserve">Base URLs: </w:t></w:r>'
    + super_run("a", "https://rest.uniprot.org")
    + super_run("b", "https://alphafold.ebi.ac.uk")
    + super_run("c", "https://gnomad.broadinstitute.org/api")
    + super_run("d", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/")
    + super_run("e", "https://myvariant.info")
    + super_run("f", "https://api.genome.ucsc.edu")
    + super_run("g", "https://rest.ensembl.org")
    + super_run("h", "https://search.clinicalgenome.org/kb/genes/{HGNC:N}")
    + super_run("i", "https://thegencc.org (query: ?hgnc_id=HGNC:{N})")
    + '</w:p>'
)

new_block = caption_xml + new_table_xml + footnote_xml
new_xml = xml[:caption_start] + new_block + xml[footnote_end:]
open(DOC, "w", encoding="utf-8").write(new_xml)

import xml.etree.ElementTree as ET
ET.parse(DOC)
print(f"OK. Wrote {DOC}", file=sys.stderr)
print(f"Document: {len(xml)} -> {len(new_xml)} chars", file=sys.stderr)
