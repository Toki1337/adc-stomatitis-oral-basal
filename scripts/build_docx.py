import os, re
from datetime import datetime, timezone
from pathlib import Path
import docx
from docx.shared import Inches, Pt, RGBColor
from docx.enum.section import WD_ORIENT, WD_SECTION
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml import OxmlElement, parse_xml
from docx.oxml.ns import nsdecls, qn

ROOT = str(Path(__file__).resolve().parents[1])
MD_PATH = f'{ROOT}/manuscript/manuscript.md'
DOCX_PATH = f'{ROOT}/manuscript/manuscript.docx'
FIG_DIR = f'{ROOT}/figures'

FIG_FILES = {
    1: f'{FIG_DIR}/Fig1_cohort_and_compartments.png',
    2: f'{FIG_DIR}/Fig2_absolute_vs_ratio.png',
    3: f'{FIG_DIR}/Fig3_resolution_and_specificity.png',
    4: f'{FIG_DIR}/Fig4_calibration_and_faers.png',
    5: f'{FIG_DIR}/Fig5_oscc_surfaceome_reranking.png',
}

doc = docx.Document()
now_utc = datetime.now(timezone.utc).replace(tzinfo=None)
doc.core_properties.author = ''
doc.core_properties.comments = ''
doc.core_properties.last_modified_by = ''
doc.core_properties.title = 'Absolute Target Expression in Normal Oral Basal Keratinocytes Predicts Antibody-Drug Conjugate Stomatitis and Re-Ranks Targets for Oral Cavity Carcinoma'
doc.core_properties.created = now_utc
doc.core_properties.modified = now_utc

# Page margins
for s in doc.sections:
    s.top_margin = Inches(0.85)
    s.bottom_margin = Inches(0.85)
    s.left_margin = Inches(0.9)
    s.right_margin = Inches(0.9)

# Base style & Heading 1-3 styles (so Word Navigation Pane & submission systems parse structure)
style = doc.styles['Normal']
style.font.name = 'Calibri'
style.font.size = Pt(10.5)
style.paragraph_format.line_spacing = 1.18
style.paragraph_format.space_after = Pt(5)

for h_name, h_sz, h_rgb, sp_bef, sp_aft in [
    ('Heading 1', 15.0, RGBColor(0x11, 0x22, 0x44), 0, 6),
    ('Heading 2', 12.5, RGBColor(0x1B, 0x4F, 0x72), 12, 4),
    ('Heading 3', 11.0, RGBColor(0x28, 0x74, 0xA6), 9, 3),
]:
    hs = doc.styles[h_name]
    hs.font.name = 'Calibri'
    hs.font.size = Pt(h_sz)
    hs.font.bold = True
    hs.font.italic = False
    hs.font.color.rgb = h_rgb
    hs.paragraph_format.space_before = Pt(sp_bef)
    hs.paragraph_format.space_after = Pt(sp_aft)
    hs.paragraph_format.keep_with_next = True

def add_formatted_runs(p, text):
    # Parse ^superscript^, ***bold-italic***, **bold**, *italic*, `code` sequentially
    tokens = re.split(r'(\^[0-9,–\-†‡§*]+\^|\*\*\*[^*]+\*\*\*|\*\*[^*]+\*\*|\*[^*]+\*|`[^`]+`)', text)
    for tok in tokens:
        if not tok:
            continue
        if tok.startswith('^') and tok.endswith('^'):
            r = p.add_run(tok[1:-1])
            r.font.superscript = True
        elif tok.startswith('***') and tok.endswith('***'):
            r = p.add_run(tok[3:-3])
            r.bold = True
            r.italic = True
        elif tok.startswith('**') and tok.endswith('**'):
            r = p.add_run(tok[2:-2])
            r.bold = True
        elif tok.startswith('*') and tok.endswith('*'):
            r = p.add_run(tok[1:-1])
            r.italic = True
        elif tok.startswith('`') and tok.endswith('`'):
            r = p.add_run(tok[1:-1])
            r.font.name = 'Consolas'
            r.font.size = Pt(9.5)
        else:
            p.add_run(tok)

def set_cell_shading(cell, fill_hex):
    shading = parse_xml(f'<w:shd {nsdecls("w")} w:fill="{fill_hex}"/>')
    cell._tc.get_or_add_tcPr().append(shading)

def set_table_borders(table):
    tblPr = table._tbl.tblPr
    borders = parse_xml(f'''
        <w:tblBorders {nsdecls("w")}>
            <w:top w:val="single" w:sz="8" w:space="0" w:color="333333"/>
            <w:bottom w:val="single" w:sz="8" w:space="0" w:color="333333"/>
            <w:insideH w:val="single" w:sz="4" w:space="0" w:color="DDDDDD"/>
            <w:left w:val="none"/>
            <w:right w:val="none"/>
            <w:insideV w:val="none"/>
        </w:tblBorders>
    ''')
    tblPr.append(borders)

lines = open(MD_PATH, encoding='utf-8').read().splitlines()

i = 0
in_figures_section = False
in_main_tables = False
while i < len(lines):
    ln = lines[i].rstrip()
    if not ln or ln == '---':
        i += 1
        continue
    if ln.startswith('# '):
        p = doc.add_paragraph(style='Heading 1')
        p.add_run(ln[2:])
        i += 1
        continue
    if ln.startswith('## '):
        h_text = ln[3:]
        in_figures_section = (h_text.strip() == 'Figure Legends')
        if h_text.strip() == 'Main Tables':
            in_main_tables = True
            sec = doc.add_section(WD_SECTION.NEW_PAGE)
            sec.orientation = WD_ORIENT.LANDSCAPE
            sec.page_width = Inches(11.0)
            sec.page_height = Inches(8.5)
            sec.left_margin = Inches(0.60)
            sec.right_margin = Inches(0.60)
            sec.top_margin = Inches(0.55)
            sec.bottom_margin = Inches(0.55)
        p = doc.add_paragraph(style='Heading 2')
        p.add_run(h_text)
        i += 1
        continue
    if ln.startswith('### '):
        if in_main_tables and (ln.startswith('### Table 2.') or ln.startswith('### Table 3.')):
            doc.add_page_break()
        p = doc.add_paragraph(style='Heading 3')
        add_formatted_runs(p, ln[4:])
        i += 1
        continue
    # Equation block
    if ln.startswith('$$') and ln.endswith('$$'):
        p = doc.add_paragraph()
        p.alignment = WD_ALIGN_PARAGRAPH.CENTER
        p.paragraph_format.space_before = Pt(4)
        p.paragraph_format.space_after = Pt(4)
        eq_parts = [
            ("logit(y", True, False), ("i", True, True), (") = β", True, False), ("0", True, True),
            (" + β", True, False), ("tub", True, True), (" · I[tubulin", True, False), ("i", True, True),
            ("] + β", True, False), ("DNA", True, True), (" · I[DNA", True, False), ("i", True, True),
            ("] + β", True, False), ("x", True, True), (" · x", True, False), ("i", True, True),
            (" + ε", True, False), ("i", True, True), (",   ε", True, False), ("i", True, True),
            (" ~ N(0, σ²)", True, False)
        ]
        for txt, ital, sub in eq_parts:
            r = p.add_run(txt)
            r.italic = ital
            r.font.subscript = sub
            r.font.size = Pt(10)
        i += 1
        continue
    # Markdown table
    if ln.startswith('|'):
        tbl_lines = []
        while i < len(lines) and lines[i].strip().startswith('|'):
            tbl_lines.append(lines[i].strip())
            i += 1
        # Parse header and data rows (skip separator row)
        rows = []
        for idx_t, tl in enumerate(tbl_lines):
            if idx_t == 1 and re.match(r'^\|[\s:\-\|]+\|$', tl):
                continue
            cells = [c.strip() for c in tl.strip('|').split('|')]
            rows.append(cells)
        if rows:
            ncols = len(rows[0])
            table = doc.add_table(rows=len(rows), cols=ncols)
            table.alignment = WD_TABLE_ALIGNMENT.CENTER
            table.autofit = False
            set_table_borders(table)
            avail_in = 9.80 if in_main_tables else 6.70
            preset_widths = {
                9: [2.05, 1.52, 0.70, 0.85, 0.63, 0.62, 0.85, 1.44, 1.14],
                7: [2.05, 2.25, 0.90, 1.20, 1.20, 1.05, 1.15],
                11: [0.84, 2.08, 0.48, 0.98, 0.74, 0.74, 0.76, 0.86, 1.26, 0.64, 0.42],
            }
            if in_main_tables and ncols in preset_widths:
                col_widths_in = preset_widths[ncols]
            else:
                clean_rows = [[r[c].replace('**', '').replace('*', '').replace('`', '') for r in rows] for c in range(ncols)]
                max_word_lens = [max([max([len(w) for w in cell.split()] + [0]) for cell in clean_rows[c]]) for c in range(ncols)]
                max_cell_lens = [max(len(cell) for cell in clean_rows[c]) for c in range(ncols)]
                weights = [max(max_word_lens[c] + (6.2 if (ncols >= 10 and c == 0) else 2.8), min(max_cell_lens[c], 25.0 if c == 1 else 15.0)) for c in range(ncols)]
                total_w = sum(weights)
                col_widths_in = [avail_in * w / total_w for w in weights]
            col_widths = [Inches(w) for w in col_widths_in]
            col_dxa = [int(round(w * 1440)) for w in col_widths_in]
            total_dxa = sum(col_dxa)
            tblPr = table._tbl.tblPr
            for old_tag in ('w:tblW', 'w:tblLayout', 'w:tblCellMar'):
                old_el = tblPr.find(qn(old_tag))
                if old_el is not None:
                    tblPr.remove(old_el)
            tblPr.append(parse_xml(f'<w:tblW {nsdecls("w")} w:w="{total_dxa}" w:type="dxa"/>'))
            tblPr.append(parse_xml(f'<w:tblLayout {nsdecls("w")} w:type="fixed"/>'))
            tblPr.append(parse_xml(
                f'<w:tblCellMar {nsdecls("w")}>'
                f'<w:top w:w="24" w:type="dxa"/><w:bottom w:w="24" w:type="dxa"/>'
                f'<w:left w:w="50" w:type="dxa"/><w:right w:w="50" w:type="dxa"/>'
                f'</w:tblCellMar>'
            ))
            old_grid = table._tbl.find(qn('w:tblGrid'))
            if old_grid is not None:
                table._tbl.remove(old_grid)
            new_grid = parse_xml(
                f'<w:tblGrid {nsdecls("w")}>'
                + ''.join(f'<w:gridCol w:w="{d}"/>' for d in col_dxa)
                + '</w:tblGrid>'
            )
            tblPr.addnext(new_grid)
            for c_idx, w in enumerate(col_widths):
                table.columns[c_idx].width = w
            if in_main_tables:
                fsize = Pt(7.8 if ncols >= 11 else (8.0 if ncols >= 9 else 8.5))
                vpad = Pt(0.6 if ncols >= 11 else 1.2)
            else:
                fsize = Pt(7.0 if ncols >= 11 else (7.3 if ncols >= 10 else (7.6 if ncols >= 9 else 8.4)))
                vpad = Pt(1.0 if ncols >= 10 else 1.5)
            for r_idx, row_cells in enumerate(rows):
                trPr = table.rows[r_idx]._tr.get_or_add_trPr()
                trPr.append(parse_xml(f'<w:cantSplit {nsdecls("w")}/>'))
                if r_idx == 0:
                    trPr.append(parse_xml(f'<w:tblHeader {nsdecls("w")}/>'))
                for c_idx, val in enumerate(row_cells):
                    cell = table.cell(r_idx, c_idx)
                    cell.width = col_widths[c_idx]
                    tcPr = cell._tc.get_or_add_tcPr()
                    old_tcW = tcPr.find(qn('w:tcW'))
                    if old_tcW is not None:
                        tcPr.remove(old_tcW)
                    tcPr.append(parse_xml(f'<w:tcW {nsdecls("w")} w:w="{col_dxa[c_idx]}" w:type="dxa"/>'))
                    p = cell.paragraphs[0]
                    p.paragraph_format.space_before = vpad
                    p.paragraph_format.space_after = vpad
                    p.paragraph_format.line_spacing = 1.04
                    if c_idx >= 2:
                        p.alignment = WD_ALIGN_PARAGRAPH.RIGHT
                    # Render ***bold-italic***, **bold**, and *italic* cleanly inside table cells
                    tokens = re.split(r'(\*\*\*[^*]+\*\*\*|\*\*[^*]+\*\*|\*[^*]+\*)', val.replace('`', ''))
                    for tok in tokens:
                        if not tok:
                            continue
                        if tok.startswith('***') and tok.endswith('***'):
                            run = p.add_run(tok[3:-3])
                            run.bold = True
                            run.italic = True
                        elif tok.startswith('**') and tok.endswith('**'):
                            run = p.add_run(tok[2:-2])
                            run.bold = True
                        elif tok.startswith('*') and tok.endswith('*'):
                            run = p.add_run(tok[1:-1])
                            run.italic = True
                        else:
                            run = p.add_run(tok)
                        run.font.size = fsize
                        if r_idx == 0:
                            run.bold = True
                    if r_idx == 0:
                        set_cell_shading(cell, "EAF2F8")
            if not (i < len(lines) and lines[i].strip().startswith('*Note:*')):
                doc.add_paragraph().paragraph_format.space_after = Pt(2)
        continue

    # Compact table footnote immediately following a Markdown table (*Note:* ...)
    if ln.startswith('*Note:*'):
        p = doc.add_paragraph()
        p.paragraph_format.space_before = Pt(2)
        p.paragraph_format.space_after = Pt(6)
        p.paragraph_format.line_spacing = 1.05
        tokens = re.split(r'(\*\*\*[^*]+\*\*\*|\*\*[^*]+\*\*|\*[^*]+\*)', ln)
        for tok in tokens:
            if not tok:
                continue
            if tok.startswith('***') and tok.endswith('***'):
                run = p.add_run(tok[3:-3])
                run.bold = True
                run.italic = True
            elif tok.startswith('**') and tok.endswith('**'):
                run = p.add_run(tok[2:-2])
                run.bold = True
            elif tok.startswith('*') and tok.endswith('*'):
                run = p.add_run(tok[1:-1])
                run.italic = True
            else:
                run = p.add_run(tok)
            run.font.size = Pt(8.5)
        i += 1
        continue
    # Check if Figure legend line -> insert corresponding PNG right before legend!
    m_fig = re.match(r'^\*\*Figure\s+(\d+)\.', ln)
    if in_figures_section and m_fig:
        fnum = int(m_fig.group(1))
        if fnum in FIG_FILES and os.path.exists(FIG_FILES[fnum]):
            p_img = doc.add_paragraph()
            p_img.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p_img.paragraph_format.space_before = Pt(8)
            p_img.paragraph_format.space_after = Pt(2)
            p_img.paragraph_format.keep_with_next = True
            p_img.add_run().add_picture(FIG_FILES[fnum], width=Inches(6.55))

    p = doc.add_paragraph()
    if re.match(r'^\d+\.\s', ln):
        p.paragraph_format.left_indent = Inches(0.22)
        p.paragraph_format.first_line_indent = Inches(-0.22)
        p.paragraph_format.space_after = Pt(3)
    add_formatted_runs(p, ln)
    i += 1

doc.save(DOCX_PATH)
print("Saved", DOCX_PATH, "size:", os.path.getsize(DOCX_PATH))
