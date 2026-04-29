"""Insert Supp Table S2 (10-row comparator) into supplementary docx after Supp Table S1."""
import re, sys

DOC = "/tmp/supp_unpacked/word/document.xml"
xml = open(DOC, encoding="utf-8").read()

# Find S1 table close: search for </w:tbl> after the "Supplementary Table S1:" caption
caption_idx = xml.find("Supplementary Table S1:")
m_close = re.search(r'</w:tbl>', xml[caption_idx:])
insert_pos = caption_idx + m_close.end()

print(f"Insert position (after S1 close): byte {insert_pos}", file=sys.stderr)

ROWS = [
    # (Tool, ACMG, Multi-variant panel, Population freq, Key distinction)
    ("VarSome (Kopanos et al., 2019)",          "Full criteria",            "No",       "Yes",  "Single-variant lookup"),
    ("GeneBe (Maj et al., 2023)",                "Full criteria",            "No",       "Yes",  "API-based; no landscape"),
    ("Franklin (Genoox)",                        "Yes",                      "No",       "Yes",  "Commercial; clinical workflow integration"),
    ("InterVar (Li and Wang, 2017)",             "Full criteria",            "No",       "No",   "Open-source ACMG; command-line"),
    ("Alamut Visual Plus (SOPHiA Genetics)",     "Yes",                      "No",       "Yes",  "Commercial NGS analysis workstation"),
    ("cBioPortal MutMapper (Cerami et al., 2012)", "None",                   "Yes",      "No",   "Somatic / cancer-genomics only"),
    ("MARRVEL (Wang et al., 2017)",              "None",                     "No",       "Yes",  "Cross-organism resource integration"),
    ("MuPIT / VarSite (Niknafs 2013; Laskowski 2020)", "None",               "No",       "No",   "3D structure context only"),
    ("VIVID (Tichkule et al., 2022)",            "None",                     "No",       "No",   "Structural context only"),
    ("VarViz (this work)",                       "Hybrid ACMG + Bayesian",   "Yes",      "Yes",  "Protein landscape + 9 live tracks + per-variant inputs"),
]
HEADERS = ("Tool", "ACMG classification", "Multi-variant panel", "Population freq", "Key distinction")
COL_DXA = (2400, 2000, 1600, 1300, 2660)  # total 9960 to match existing supp tables

def cell(text, dxa, bold):
    rpr = '<w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/>'
    if bold:
        rpr += '<w:b/><w:bCs/>'
    return (
        f'<w:tc>'
        f'<w:tcPr><w:tcW w:w="{dxa}" w:type="dxa"/></w:tcPr>'
        f'<w:p><w:pPr><w:spacing w:lineRule="exact" w:line="240" w:before="0" w:after="0"/><w:jc w:val="left"/></w:pPr>'
        f'<w:r><w:rPr>{rpr}</w:rPr><w:t xml:space="preserve">{text}</w:t></w:r>'
        f'</w:p></w:tc>'
    )

def row(cells_text, bold):
    return '<w:tr>' + ''.join(cell(t, dxa, bold) for t, dxa in zip(cells_text, COL_DXA)) + '</w:tr>'

# Caption
caption_xml = (
    '<w:p>'
    '<w:pPr><w:spacing w:lineRule="exact" w:line="276" w:before="240" w:after="60"/><w:jc w:val="left"/></w:pPr>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/><w:b/><w:bCs/></w:rPr><w:t xml:space="preserve">Supplementary Table S2. </w:t></w:r>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr><w:t xml:space="preserve">Comparison of VarViz with related variant visualization and classification tools.</w:t></w:r>'
    '</w:p>'
)

table_xml = (
    '<w:tbl>'
    '<w:tblPr>'
    '<w:tblW w:w="9960" w:type="dxa"/>'
    '<w:tblLayout w:type="fixed"/>'
    '<w:tblCellMar><w:top w:w="60" w:type="dxa"/><w:left w:w="100" w:type="dxa"/><w:bottom w:w="60" w:type="dxa"/><w:right w:w="100" w:type="dxa"/></w:tblCellMar>'
    '<w:tblLook w:val="04A0" w:firstRow="1" w:lastRow="0" w:firstColumn="1" w:lastColumn="0" w:noHBand="0" w:noVBand="1"/>'
    '</w:tblPr>'
    '<w:tblGrid>' + ''.join(f'<w:gridCol w:w="{d}"/>' for d in COL_DXA) + '</w:tblGrid>'
    + row(HEADERS, bold=True)
    + ''.join(row(r, bold=(i == 9)) for i, r in enumerate(ROWS))  # bold the last (VarViz) row
    + '</w:tbl>'
)

# Footnote
footnote_xml = (
    '<w:p>'
    '<w:pPr><w:spacing w:lineRule="exact" w:line="240" w:before="60" w:after="120"/><w:jc w:val="left"/></w:pPr>'
    '<w:r><w:rPr><w:rFonts w:eastAsia="Times New Roman" w:cs="Times New Roman"/><w:color w:val="000000"/><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr><w:t xml:space="preserve">Franklin is a commercial product of Genoox; Alamut Visual Plus is a commercial product of SOPHiA Genetics.</w:t></w:r>'
    '</w:p>'
)

new_block = caption_xml + table_xml + footnote_xml
new_xml = xml[:insert_pos] + new_block + xml[insert_pos:]
open(DOC, "w", encoding="utf-8").write(new_xml)

import xml.etree.ElementTree as ET
ET.parse(DOC)
print(f"OK. Inserted {len(new_block)} chars after S1 table close", file=sys.stderr)
