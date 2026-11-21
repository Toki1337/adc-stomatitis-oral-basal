"""Automated integrity gate for the public code and data archive:
verifies dual-blinded clinical TSVs, processed results/numbers.json,
Supplementary_Tables.xlsx sheets & cross-sheet CI consistency,
Supplementary_References.md, and generated figure files.
"""

import csv
import json
import os
from pathlib import Path
import sys
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
NUM_PATH = ROOT / 'results/numbers.json'
SUPP_XLSX = ROOT / 'manuscript/Supplementary_Tables.xlsx'
SUPP_REFS_MD = ROOT / 'manuscript/Supplementary_References.md'

errors = []

# 1. Verify numbers.json structure and core metrics
if not NUM_PATH.exists():
    errors.append(f'MISSING: {NUM_PATH}')
    nums = {}
else:
    nums = json.loads(NUM_PATH.read_text(encoding='utf-8'))
    for k in ['primary_model', 'target_inference', 'calibration', 'faers', 'joint_models', 'surfaceome', 'hpa_dab', 'clinical_dual_agreement']:
        if k not in nums:
            errors.append(f'NUMBERS FAIL: Missing key "{k}" in results/numbers.json')

# 2. Verify clinical dual-extraction TSVs
for label_ext, rel_p in [
    ('Extractor A', 'clinical/adc_stomatitis_extractor_A.tsv'),
    ('Extractor B', 'clinical/adc_stomatitis_extractor_B.tsv'),
    ('Trials Final', 'clinical/adc_stomatitis_trials.tsv'),
    ('Reconciled', 'clinical/adc_stomatitis_reconciled.tsv'),
]:
    p = ROOT / rel_p
    if not p.exists():
        errors.append(f'CLINICAL TSV FAIL: Missing {rel_p}')
        continue
    with open(p, encoding='utf-8') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    if len(rows) != 66:
        errors.append(f'CLINICAL TSV FAIL: {label_ext} has {len(rows)} rows (expected 66)')
    if label_ext in {'Extractor A', 'Extractor B'}:
        if any(not (r.get('extractor', '').strip() and r.get('verified_by', '').strip()) for r in rows):
            errors.append(f'CLINICAL TSV FAIL: {label_ext} has blank extractor or verified_by cell')

# 3. Verify Supplementary_Tables.xlsx & Supplementary_References.md
if not SUPP_XLSX.exists() or not SUPP_REFS_MD.exists():
    errors.append('SUPPLEMENT FAIL: Missing Supplementary_Tables.xlsx or Supplementary_References.md')
else:
    if '`' in SUPP_REFS_MD.read_text(encoding='utf-8'):
        errors.append('FORMAT FAIL: Found backticks in Supplementary_References.md')
    wb = openpyxl.load_workbook(SUPP_XLSX, read_only=True)
    required_sheets = [
        'README',
        'Table_S1_Clinical_66ADCs',
        'Table_S1a_Extractor_A_66ADCs',
        'Table_S1b_Extractor_B_66ADCs',
        'Table_S1c_Three_Dataset_Results',
        'Table_S1d_Reconciled_66ADCs',
        'Table_S1e_Dual_Agreement',
        'Table_S2a_Features',
        'Table_S2b_Resolution_Ladder',
        'Table_S2c_HPA_DAB_Models',
        'Table_S2d_HPA_DAB_Genes',
        'Table_S2e_HPA_DAB_ParamGrid',
        'Table_S3a_Target_Inference',
        'Table_S3b_Sensitivity_Subgroups',
        'Table_S3c_Protocol_Alignment',
        'Table_S4_CellType_Tissue_Scan',
        'Table_S5a_Axes_Single',
        'Table_S5b_Joint_Nested',
        'Table_S5c_Ratio_Scales',
        'Table_S5d_Axes_Per_ADC',
        'Table_S6a_LOTO_LOO_Calib',
        'Table_S6b_LOTO_Exact_Endpoints',
        'Table_S6c_Calib_Summary',
        'Table_S7_OSCC_420_Surfaceome',
        'Table_S8_FAERS_Models',
        'Supplementary_References',
    ]
    for sh in required_sheets:
        if sh not in wb.sheetnames:
            errors.append(f'SUPPLEMENT FAIL: Missing sheet {sh}')
    if 'Table_S1_Clinical_66ADCs' in wb.sheetnames and 'Supplementary_References' in wb.sheetnames:
        ws_s1 = wb['Table_S1_Clinical_66ADCs']
        hdr_s1 = [c.value for c in next(ws_s1.iter_rows(min_row=1, max_row=1))]
        ref_idx = hdr_s1.index('supp_ref_id')
        s1_refs = [r[ref_idx].value for r in ws_s1.iter_rows(min_row=2)]
        ws_sr = wb['Supplementary_References']
        sr_ids = {r[0].value for r in ws_sr.iter_rows(min_row=2) if r[0].value}
        if len(s1_refs) != 66 or any(not v for v in s1_refs):
            errors.append(f'SUPPLEMENT FAIL: Table S1 supp_ref_id incomplete ({len(s1_refs)} rows)')
        for cell_val in s1_refs:
            if cell_val:
                for part in str(cell_val).split(';'):
                    if part.strip() not in sr_ids:
                        errors.append(f'SUPPLEMENT FAIL: Table S1 supp_ref_id "{part.strip()}" missing from Supplementary_References')
    if 'Table_S2b_Resolution_Ladder' in wb.sheetnames and 'Table_S2c_HPA_DAB_Models' in wb.sheetnames:
        ws_2b = wb['Table_S2b_Resolution_Ladder']
        h2b = [c.value for c in next(ws_2b.iter_rows(min_row=1, max_row=1))]
        m2b = {r[0].value: (float(r[h2b.index('dLR_boot_ci_lo')].value), float(r[h2b.index('dLR_boot_ci_hi')].value)) for r in ws_2b.iter_rows(min_row=2)}
        ws_2c = wb['Table_S2c_HPA_DAB_Models']
        h2c = [c.value for c in next(ws_2c.iter_rows(min_row=1, max_row=1))]
        m2c = {r[0].value: (r[h2c.index('dLR_ci_lo')].value, r[h2c.index('dLR_ci_hi')].value) for r in ws_2c.iter_rows(min_row=2)}
        for k2b, k2c in [('HPA_esoph_DAB', 'esoph_dab_basal'), ('HPA_IHC_oral', 'HPA_IHC_oral_summary'), ('HPA_oral_DAB', 'oral_dab_basal'), ('HPA_skin_DAB', 'skin_dab_basal')]:
            lo2b, hi2b = m2b[k2b]
            lo2c, hi2c = float(m2c[k2c][0]), float(m2c[k2c][1])
            if abs(lo2b - lo2c) > 1e-6 or abs(hi2b - hi2c) > 1e-6:
                errors.append(f'SUPPLEMENT FAIL: Delta-LR CI mismatch between Table S2b ({k2b}) and Table S2c ({k2c})')

# 4. Verify HPA image manifest and generated figures
manifest_p = ROOT / 'dl/hpa_ihc/manifest_509.tsv'
if not manifest_p.exists():
    errors.append('MISSING: dl/hpa_ihc/manifest_509.tsv')
else:
    with open(manifest_p, encoding='utf-8') as f:
        m_rows = list(csv.DictReader(f, delimiter='\t'))
    if len(m_rows) != 509:
        errors.append(f'MANIFEST FAIL: Expected 509 HPA sections, got {len(m_rows)}')

for stem in [
    'Fig1_cohort_and_compartments',
    'Fig2_absolute_vs_ratio',
    'Fig3_resolution_and_specificity',
    'Fig4_calibration_and_faers',
    'Fig5_oscc_surfaceome_reranking',
]:
    for ext in ('pdf', 'png'):
        fp = ROOT / f'figures/{stem}.{ext}'
        if not fp.exists() or fp.stat().st_size == 0:
            errors.append(f'FIGURE FAIL: Missing or empty figures/{stem}.{ext}')

if errors:
    print('\nARCHIVE VERIFICATION FAILED:')
    for e in errors:
        print(' -', e)
    sys.exit(1)

print('ALL ARCHIVE VERIFICATION CHECKS PASSED (clinical TSVs, numbers.json, 26 supplementary sheets, 509-section manifest, and 5 figures verified).')
