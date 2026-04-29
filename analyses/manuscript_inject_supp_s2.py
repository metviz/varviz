"""Strip pandoc-style attrs + jc inside tblPr; inject."""
import re, sys
SUPP = "/tmp/supp_unpacked/word/document.xml"
FRAG = "/tmp/s2_frag/word/document.xml"
frag_xml = open(FRAG, encoding="utf-8").read()
supp_xml = open(SUPP, encoding="utf-8").read()
m_body_open = re.search(r"<w:body[^>]*>", frag_xml)
m_sectpr = re.search(r"<w:sectPr\b", frag_xml[m_body_open.end():])
fragment = frag_xml[m_body_open.end():m_body_open.end() + m_sectpr.start()]

def strip_self_closing(name, body):
    return re.sub(rf'<w:{name}\b[^/>]*/>', '', body)
def strip_open_close(name, body):
    return re.sub(rf'<w:{name}\b[^>]*>.*?</w:{name}>', '', body, flags=re.DOTALL)

fragment = strip_self_closing('pStyle', fragment)
fragment = strip_open_close('numPr', fragment)
fragment = strip_self_closing('bookmarkStart', fragment)
fragment = strip_self_closing('bookmarkEnd', fragment)
fragment = re.sub(r'<w:hyperlink\b[^>]*>', '', fragment)
fragment = re.sub(r'</w:hyperlink>', '', fragment)
fragment = strip_self_closing('tblStyle', fragment)
fragment = strip_open_close('tblPrEx', fragment)

# Strip <w:jc .../> only when inside <w:tblPr>...</w:tblPr> contexts
def strip_jc_in_tblpr(s):
    def fix(m):
        body = m.group(0)
        return re.sub(r'<w:jc\b[^/>]*/>', '', body)
    return re.sub(r'<w:tblPr>.*?</w:tblPr>', fix, s, flags=re.DOTALL)
fragment = strip_jc_in_tblpr(fragment)

m_supp_sectpr = list(re.finditer(r"<w:sectPr\b", supp_xml))
insert_pos = m_supp_sectpr[-1].start()
new_xml = supp_xml[:insert_pos] + fragment + supp_xml[insert_pos:]
open(SUPP, "w", encoding="utf-8").write(new_xml)
import xml.etree.ElementTree as ET
ET.parse(SUPP)
print(f"OK. {len(supp_xml)} -> {len(new_xml)} chars", file=sys.stderr)
