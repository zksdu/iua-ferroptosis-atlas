# -*- coding: utf-8 -*-
"""41_build_docx_chain.py — rebuild P2 manuscript docx + colleague-review preview docx.
Step 1: HTML (from 37_build_html.py) -> submission docx via html4docx.
Step 2: preview docx = submission docx + journal banner + PREVIEW banner +
        embedded figures under legends + embedded supplementary tables S1-S6.
"""
import os, shutil
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.oxml.ns import qn
from docx.oxml import OxmlElement
from htmldocx import HtmlToDocx

BASE = r"C:/Users/Administrator/WorkBuddy/2026-09-27-09-39-22/iua-exosome-bioinfo"
HTML = BASE + r"/output/p2en/stage2/formatted-paper2-en.html"
SUB = BASE + r"/results/paper2_manuscript_EN_v1.0.docx"
PREV = BASE + r"/results/paper2_manuscript_EN_v1.0_WITH_FIGURES_preview.docx"
FIGD = BASE + r"/P2_FINAL_SUBMISSION/03_Figures"
TABD = BASE + r"/P2_FINAL_SUBMISSION/04_SupplementaryTables"
results = BASE + "/results"

# ---------- Step 1 ----------
doc = Document()
with open(HTML, encoding="utf-8") as fh:
    html_text = fh.read()
# body only
import re
m = re.search(r"<body[^>]*>(.*)</body>", html_text, re.S)
body_html = m.group(1) if m else html_text
HtmlToDocx().add_html_to_document(body_html, doc)
doc.save(SUB)
print("submission docx saved:", SUB)

# ---------- Step 2 ----------
d = Document(SUB)

def shade(p, hexcolor):
    pPr = p._p.get_or_add_pPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear"); shd.set(qn("w:fill"), hexcolor)
    pPr.append(shd)

def add_run(p, text, bold=False, size=10.5, color=None):
    r = p.add_run(text)
    r.bold = bold
    r.font.size = Pt(size)
    if color:
        r.font.color.rgb = RGBColor.from_string(color)
    return r

# journal banner at very top
pj = d.paragraphs[0]
p1 = d.add_paragraph(); p2 = d.add_paragraph()
add_run(p1, "TARGET JOURNAL: ", bold=True, size=12, color="1F4E79")
add_run(p1, "Frontiers in Immunology", bold=True, size=12, color="1F4E79")
add_run(p2, "IF 5.9 (2024) | JCR Q1 | CAS Zone 2 | ISSN 1664-2242 | PREVIEW build 2026-09-29 (v1.1 robustness edition)", size=10, color="444444")
p2.paragraph_format.space_after = Pt(6)
pj._p.addprevious(p1._p)
pj._p.addprevious(p2._p)   # order now: [p1, p2, pj, ...]

# PREVIEW banner after journal banner
pb = d.add_paragraph()
add_run(pb, "PREVIEW COPY with embedded figures and supplementary tables — for colleague review only, not for submission.", bold=True, size=10, color="7F1D1D")
shade(pb, "FDECEC")
p2._p.addnext(pb._p)
for pp in [p1, p2]:
    shade(pp, "EAF1F8")

# figures mapping
FIGMAP = {
 "Figure 1.": ["Figure1.png"],
 "Figure 2.": ["Figure2.png"],
 "Figure 3.": ["Figure3.png"],
 "Figure 4.": ["Figure4.png"],
 "Figure 5.": ["Figure5.png"],
 "Figure 6.": ["Figure6.png"],
 "Figure 7.": ["Figure7A.png", "Figure7B.png"],
 "Figure S1.": ["FigureS1.png"],
 "Figure S2.": ["FigureS2.png"],
 "Figure S3.": ["figS3_loo_sensitivity.png"],
 "Figure S4.": ["figS4_deconv_gsea.png"],
 "Figure S5.": ["figS5_wnt_ferro_celltype.png"],
 "Figure S6.": ["figS6_myeloid.png"],
}
inserted = 0
for p in list(d.paragraphs):
    t = p.text.strip()
    for key, files in FIGMAP.items():
        if t.startswith(key):
            anchor = p._p
            for fn in files:
                fp = os.path.join(FIGD, fn)
                if not os.path.exists(fp):
                    fp = os.path.join(BASE, "figures", fn)
                if os.path.exists(fp):
                    np_para = d.add_paragraph()
                    np_para.add_run().add_picture(fp, width=Inches(6.0))
                    anchor.addnext(np_para._p)
                    anchor = np_para._p
                    inserted += 1
print("figures embedded:", inserted)

# supplementary tables S1-S6 appended
TABS = [
 ("TableS1_gene_panel.csv", "Table S1. Gene panels (bulk + single-cell layers)", 9),
 ("TableS2_gse224093_DEG.csv", "Table S2. GSE224093 differentially expressed genes (padj < 0.05)", 9),
 ("TableS3_gse215968_celltype.csv", "Table S3. GSE215968 cell-type module statistics", 8),
 ("TableS4_gse311899_pseudobulk.csv", "Table S4. GSE311899 pseudobulk expression matrix", 6),
 ("TableS5_hmgb1_panel.csv", "Table S5. 20-gene HMGB1 axis panel", 9),
 ("TableS6_robustness_bundle.csv", "Table S6. Robustness and confounding analyses (LOO, external set, coupling, myeloid, composition)", 7),
]
def shade_cell(cell, hexcolor):
    tcPr = cell._tc.get_or_add_tcPr()
    shd = OxmlElement("w:shd")
    shd.set(qn("w:val"), "clear"); shd.set(qn("w:fill"), hexcolor)
    tcPr.append(shd)

import csv
for fname, caption, fs in TABS:
    src = os.path.join(TABD, fname)
    if not os.path.exists(src):
        src = os.path.join(results, "TableS6_robustness_bundle.csv") if "S6" in fname else src
    with open(src, encoding="utf-8-sig") as fh:
        rows = list(csv.reader(fh))
    p = d.add_paragraph()
    r = p.add_run(caption)
    r.bold = True; r.font.size = Pt(10.5)
    r.font.color.rgb = RGBColor.from_string("1F4E79")
    t = d.add_table(rows=len(rows), cols=len(rows[0]))
    t.style = "Table Grid"
    for ri, row in enumerate(rows):
        for ci, val in enumerate(row):
            cell = t.cell(ri, ci)
            cell.text = ""
            cp = cell.paragraphs[0]
            cr = cp.add_run(val)
            cr.font.size = Pt(fs)
            if ri == 0:
                cr.bold = True
                shade_cell(cell, "D9E2F3")
    d.add_paragraph().paragraph_format.space_after = Pt(6)
print("supplementary tables embedded:", len(TABS))

d.save(PREV)
print("preview docx saved:", PREV)

# verify
d2 = Document(PREV)
print("tables:", len(d2.tables), "| paragraphs:", len(d2.paragraphs))
print("top:", [p.text[:60] for p in d2.paragraphs[:3]])
