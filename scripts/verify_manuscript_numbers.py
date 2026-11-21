"""Automated gate verifying manuscript text, DOCX compliance, citations, supplementary workbook, and consistency with results/numbers.json."""

import json
import os
from pathlib import Path
import re
import sys
import docx
import openpyxl

ROOT = Path(__file__).resolve().parents[1]
MD_PATH = f'{ROOT}/manuscript/manuscript.md'
DOCX_PATH = f'{ROOT}/manuscript/manuscript.docx'
NUM_PATH = f'{ROOT}/results/numbers.json'
SUPP_XLSX = f'{ROOT}/manuscript/Supplementary_Tables.xlsx'
SUPP_REFS_MD = f'{ROOT}/manuscript/Supplementary_References.md'

md = open(MD_PATH, encoding='utf-8').read()
nums = json.load(open(NUM_PATH, encoding='utf-8'))

errors = []

# 1. Formatting compliance: zero raw backticks or code-identifier underscores across md and docx
n_bt_md = md.count('`')
if n_bt_md != 0:
    errors.append(f'FORMAT FAIL: Found {n_bt_md} backticks in manuscript.md')

doc = docx.Document(DOCX_PATH)
docx_texts = [p.text for p in doc.paragraphs]
for tbl in doc.tables:
    for row in tbl.rows:
        for cell in row.cells:
            docx_texts.append(cell.text)
full_docx_text = '\n'.join(docx_texts)
n_bt_docx = full_docx_text.count('`')
if n_bt_docx != 0:
    errors.append(f'FORMAT FAIL: Found {n_bt_docx} backticks in manuscript.docx')

md_no_math = re.sub(r'\$\$.*?\$\$', '', re.sub(r'https?://\S+|\bfile://\S+|\S+@\S+', '', md))
underscore_hits = re.findall(r'\b[A-Za-z][A-Za-z0-9]*_[A-Za-z0-9_]+\b', md_no_math)
if underscore_hits:
    errors.append(f'FORMAT FAIL: Found underscore variable names in manuscript.md: {sorted(set(underscore_hits))[:10]}')

# 2. Structured Abstract word count <= 300
abs_m = re.search(r'## Abstract\s+(.*?)\s+## Introduction', md, re.S)
if not abs_m:
    errors.append('ABSTRACT FAIL: Could not locate Abstract section')
else:
    abs_clean = re.sub(r'\*\*([^*]+)\*\*', r'\1', abs_m.group(1)).strip()
    abs_words = len(abs_clean.split())
    print(f'Abstract word count: {abs_words} words')
    if abs_words > 300:
        errors.append(f'ABSTRACT FAIL: Abstract has {abs_words} words (> 300 limit)')

# 3. Bidirectional bibliography and Supplementary_Tables.xlsx / Supplementary_References.md completeness
ref_block = md.split('## References')[1].split('## Figure Legends')[0]
ref_lines = [ln.strip() for ln in ref_block.splitlines() if re.match(r'^\d+\.\s', ln.strip())]
ref_keys = []
for ln in ref_lines:
    m_yr = re.search(r'\((\d{4})\)', ln)
    yr = m_yr.group(1) if m_yr else ''
    body = ln.split('.', 1)[1].strip()
    if body.startswith('Cancer Genome Atlas Network'):
        surname = 'Cancer Genome Atlas Network'
    elif body.startswith('GTEx Consortium'):
        surname = 'GTEx Consortium'
    elif body.startswith('R Core Team'):
        surname = 'R Core Team'
    elif body.startswith('U.S. Food and Drug Administration'):
        surname = 'U.S. Food and Drug Administration'
    elif body.startswith('UniProt Consortium'):
        surname = 'UniProt Consortium'
    else:
        surname = re.match(r'^([^,]+),', body).group(1).strip()
    ref_keys.append((surname, yr, ln))

body_before_refs = md.split('## References')[0]
body_plus_legends_tables = body_before_refs + '\n' + md.split('## References')[1].split('\n---\n', 1)[1]
ref_pair_set = {(s, yr) for s, yr, _ in ref_keys}
for surname, yr, full_ref in ref_keys:
    # Match surname and year in close proximity in body/legends/tables
    s_pat = re.escape(surname).replace(r'\.', r'\.?')
    if not re.search(rf'{s_pat}[^()\[\]\n]{{0,35}}\b{yr}\b|\b{yr}\b[^()\[\]\n]{{0,15}}{s_pat}', body_plus_legends_tables):
        errors.append(f'CITATION FAIL: Orphan reference in bibliography not cited in text: {surname} ({yr})')

# Reverse check (body/legends/tables -> References): every author-year citation must resolve in ## References
cite_pats = [
    r'\b(Cancer Genome Atlas Network|GTEx Consortium|R Core Team|U\.S\. Food and Drug Administration|UniProt Consortium),?\s*\(?((?:19|20)\d{2})\)?',
    r'\b(Seabold and Perktold|Seabold & Perktold|Ruifrok and Johnston|Ruifrok & Johnston),?\s*\(?((?:19|20)\d{2})\)?',
    r'\b([A-Z][A-Za-zÀ-ÿ-]+)\s+et al\.,?\s*\(?((?:19|20)\d{2})\)?',
    r'\b(Wickham|Slowikowski|Wilke),?\s*\(?((?:19|20)\d{2})\)?',
]
for pat in cite_pats:
    for m_c in re.finditer(pat, body_plus_legends_tables):
        auth_raw, yr_c = m_c.group(1), m_c.group(2)
        auth_norm = 'Seabold' if auth_raw.startswith('Seabold') else ('Ruifrok' if auth_raw.startswith('Ruifrok') else auth_raw)
        if (auth_norm, yr_c) not in ref_pair_set:
            errors.append(f'CITATION FAIL: In-text citation "{auth_raw} ({yr_c})" missing from ## References')

# Check Supplementary_Tables.xlsx supp_ref_id mapping and Supplementary_References.md
if not os.path.exists(SUPP_XLSX) or not os.path.exists(SUPP_REFS_MD):
    errors.append('SUPPLEMENT FAIL: Missing Supplementary_Tables.xlsx or Supplementary_References.md')
else:
    sr_md_txt = open(SUPP_REFS_MD, encoding='utf-8').read()
    if '`' in sr_md_txt:
        errors.append('FORMAT FAIL: Found backticks in Supplementary_References.md')
    wb = openpyxl.load_workbook(SUPP_XLSX, read_only=True)
    required_s1_sheets = [
        'Table_S1a_Extractor_A_66ADCs',
        'Table_S1b_Extractor_B_66ADCs',
        'Table_S1c_Three_Dataset_Results',
        'Table_S1d_Reconciled_66ADCs',
        'Table_S1e_Dual_Agreement',
        'Table_S2e_HPA_DAB_ParamGrid',
        'Table_S3c_Protocol_Alignment',
        'Table_S6a_LOTO_LOO_Calib',
        'Table_S6b_LOTO_Exact_Endpoints',
        'Table_S6c_Calib_Summary',
        'Table_S7_OSCC_420_Surfaceome',
        'Table_S8_FAERS_Models',
    ]
    for sh in required_s1_sheets:
        if sh not in wb.sheetnames:
            errors.append(f'SUPPLEMENT FAIL: Missing sheet {sh} in Supplementary_Tables.xlsx')
    if 'Supplementary_References' not in wb.sheetnames:
        errors.append('SUPPLEMENT FAIL: Missing Supplementary_References sheet in Supplementary_Tables.xlsx')
    else:
        ws_s1 = wb['Table_S1_Clinical_66ADCs']
        hdr_s1 = [c.value for c in next(ws_s1.iter_rows(min_row=1, max_row=1))]
        if 'supp_ref_id' not in hdr_s1:
            errors.append('SUPPLEMENT FAIL: Table_S1_Clinical_66ADCs missing supp_ref_id column')
        else:
            ref_idx = hdr_s1.index('supp_ref_id')
            s1_refs = [r[ref_idx].value for r in ws_s1.iter_rows(min_row=2)]
            if len(s1_refs) != 66 or any(not v for v in s1_refs):
                errors.append(f'SUPPLEMENT FAIL: Table S1 supp_ref_id incomplete ({len(s1_refs)} rows)')
            ws_sr = wb['Supplementary_References']
            sr_ids = {r[0].value for r in ws_sr.iter_rows(min_row=2) if r[0].value}
            if f'S1–S{len(sr_ids)}' not in md:
                errors.append(f'SUPPLEMENT FAIL: Manuscript range does not match actual Supplementary_References count (S1–S{len(sr_ids)})')
            for cell_val in s1_refs:
                if cell_val:
                    for part in str(cell_val).split(';'):
                        if part.strip() not in sr_ids:
                            errors.append(f'SUPPLEMENT FAIL: Table S1 supp_ref_id "{part.strip()}" missing from Supplementary_References')
        if 'Table_S7_OSCC_420_Surfaceome' in wb.sheetnames:
            ws_s7 = wb['Table_S7_OSCC_420_Surfaceome']
            hdr_s7 = [c.value for c in next(ws_s7.iter_rows(min_row=1, max_row=1))]
            for req_dm in ['depmap_common_essential', 'depmap_median_chronos', 'depmap_frac_lt_neg05']:
                if req_dm not in hdr_s7:
                    errors.append(f'SUPPLEMENT FAIL: Table S7 missing {req_dm} column')
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
                    errors.append(f'SUPPLEMENT FAIL: Delta-LR CI mismatch between Table S2b ({k2b}: {lo2b:.4f}, {hi2b:.4f}) and Table S2c ({k2c}: {lo2c:.4f}, {hi2c:.4f})')

# 4. Every italicized uppercase gene symbol in Discussion appears in Results or Main Tables
disc_sec = md.split('## Discussion')[1].split('## Methods')[0]
res_sec = md.split('## Results')[1].split('## Discussion')[0]
tables_sec = md.split('## Main Tables')[1]
res_and_tables = res_sec + '\n' + tables_sec
disc_genes = sorted(set(re.findall(r'\*([A-Z][A-Z0-9-]{1,10})\*', disc_sec)))
for g in disc_genes:
    if f'*{g}*' not in res_and_tables and f'**{g}**' not in res_and_tables:
        errors.append(f'GENE CHECK FAIL: Gene *{g}* in Discussion not introduced in Results or Main Tables')

# 5. Verify standalone clinical extraction TSVs
import csv
for label_ext, ext_p in [
    ('Extractor A', f'{ROOT}/clinical/adc_stomatitis_extractor_A.tsv'),
    ('Extractor B', f'{ROOT}/clinical/adc_stomatitis_extractor_B.tsv'),
]:
    if not os.path.exists(ext_p):
        errors.append(f'CLINICAL TSV FAIL: Missing {ext_p}')
    else:
        with open(ext_p, encoding='utf-8') as f:
            rows_ext = list(csv.DictReader(f, delimiter='\t'))
        if len(rows_ext) != 66:
            errors.append(f'CLINICAL TSV FAIL: {label_ext} has {len(rows_ext)} rows (expected 66)')
        elif any(not (r.get('extractor', '').strip() and r.get('verified_by', '').strip()) for r in rows_ext):
            errors.append(f'CLINICAL TSV FAIL: {label_ext} has blank extractor or verified_by cell')

# 6. All Supplementary Table sheets cited in manuscript.md
for tok in ['Table S1', 'Table S1a', 'Table S1b', 'Table S1c', 'Table S1d', 'Table S1e', 'Table S2a', 'Table S2b', 'Table S2c', 'Table S2d', 'Table S2e', 'Table S3a', 'Table S3b', 'Table S3c', 'Table S4', 'Table S5a', 'Table S5b', 'Table S5c', 'Table S5d', 'Table S6a', 'Table S6b', 'Table S6c', 'Table S7', 'Table S8']:
    if tok in {'Table S1a', 'Table S1b', 'Table S1c', 'Table S1d', 'Table S1e'} and (tok in md or 'Table S1a–e' in md or 'Tables S1a–e' in md):
        continue
    if tok in {'Table S6a', 'Table S6b', 'Table S6c'} and (tok in md or 'Table S6a–c' in md or 'Tables S6a–c' in md):
        continue
    if tok == 'Table S2a' and ('Table S2a' in md or 'Table S2a–b' in md or 'Table S2a–e' in md):
        continue
    if tok == 'Table S2b' and ('Table S2b' in md or 'Table S2a–b' in md or 'Table S2a–e' in md):
        continue
    if tok == 'Table S2c' and ('Table S2c' in md or 'Table S2c–e' in md or 'Table S2c–d' in md or 'Table S2a–e' in md):
        continue
    if tok == 'Table S2d' and ('Table S2d' in md or 'Table S2c–e' in md or 'Table S2c–d' in md or 'Table S2a–e' in md):
        continue
    if tok == 'Table S2e' and ('Table S2e' in md or 'Table S2c–e' in md or 'Table S2a–e' in md):
        continue
    if tok == 'Table S5b' and ('Table S5b' in md or 'Table S5b–c' in md):
        continue
    if tok not in md:
        errors.append(f'SUPPLEMENT FAIL: Missing citation to {tok} in manuscript.md')

# 7. Outcome heterogeneity count consistency ("seven" == numbers.json)
n_oh = nums['n_outcome_heterogeneity']
if n_oh != 7 or 'seven post hoc outcome-heterogeneity checks' not in res_sec or 'seven outcome-heterogeneity sensitivity checks' not in md:
    errors.append(f'COUNT FAIL: Outcome-heterogeneity count mismatch (numbers.json={n_oh})')

# 8. Programmatic numerical spot-checks against numbers.json
ag_map = {r['metric']: r['value'] for r in nums['clinical_dual_agreement']}
dab_map = {r['feature']: r for r in nums['hpa_dab']['models']}
t3_map = {r['gene']: r for r in nums['table3_targets']}
expected_tokens = [
    str(nums['primary_model']['W_basal_LR']),
    f"{nums['primary_model']['W_basal_spearman']:.3f}",
    f"{nums['primary_model']['floor_x0_topoI_pct']:.1f}",
    f"{nums['primary_model']['floor_x0_tubulin_pct']:.1f}",
    f"{nums['primary_model']['floor_x0_dna_pct']:.1f}",
    f"{nums['target_inference']['perm_p']:.4f}",
    f"{nums['target_inference']['boot_beta_ci_lo']:.2f}",
    f"{nums['target_inference']['boot_beta_ci_hi']:.2f}",
    f"{nums['posthoc_4payload']['basal_LR']:.2f}",
    f"{nums['posthoc_4payload']['basal_beta']:.2f}",
    f"{nums['posthoc_4payload']['basal_perm_p']:.4f}",
    f"{nums['posthoc_4payload']['payload_4v3_LR']:.2f}",
    f"{nums['weighted_tobit']['n_weight_beta']:.2f}",
    f"{nums['weighted_tobit']['n_weight_boot_ci'][0]:.2f}",
    f"{nums['weighted_tobit']['n_weight_boot_ci'][1]:.2f}",
    f"{nums['weighted_tobit']['sqrtn_weight_beta']:.2f}",
    f"{nums['weighted_tobit']['sqrtn_weight_boot_ci'][0]:.2f}",
    f"{nums['weighted_tobit']['sqrtn_weight_boot_ci'][1]:.2f}",
    f"{nums['calibration'][0]['cens_below_x']}/{nums['calibration'][0]['n_censored']}",
    f"{nums['calibration'][0]['c_index_all44_x']:.3f}",
    f"{nums['calibration'][0]['p_perm_c_index']:.4f}",
    f"{nums['calibration'][0]['auc_ge15_x']:.3f}",
    f"{nums['calibration'][0]['p_mw_auc_ge15']:.4f}",
    f"{nums['calibration'][0]['spearman_x']:.3f}",
    f"{nums['calibration'][0]['mae_x']:.1f}",
    f"{nums['calibration'][1]['c_index_all44_x']:.3f}",
    f"{nums['calibration'][1]['c_index_all44_payload']:.3f}",
    f"{nums['calibration'][1]['auc_ge15_x']:.3f}",
    f"{nums['calibration'][1]['auc_ge15_payload']:.3f}",
    f"{nums['faers'][0]['beta']:.2f}",
    f"{nums['faers'][0]['LR']:.2f}",
    f"{nums['faers'][0]['p_LR']:.4f}",
    f"{nums['faers'][0]['p_perm']:.4f}",
    f"{nums['faers'][0]['ols_LR']:.2f}",
    f"{nums['faers'][0]['p_ols_perm']:.4f}",
    f"{nums['joint_models'][0]['added_LR_ratio_over_abs']:.2f}",
    f"{nums['joint_models'][0]['added_LR_abs_over_ratio']:.2f}",
    f"{nums['surfaceome']['spearman_choi_basal']:.3f}",
    str(nums['surfaceome']['n_oscc_expressed']),
    str(nums['surfaceome']['n_winners_4of4']),
    str(nums['surfaceome']['n_winners_vital15']),
    str(nums['surfaceome']['n_winners_ext15']),
    str(nums['surfaceome']['n_robust_7of9']),
    str(nums['surfaceome']['n_robust_6of9']),
    str(nums['specificity_scan']['hpa_n_q05']),
    str(nums['specificity_scan']['gtex_n_q05']),
    f"{nums['specificity_scan']['hpa_min_q_noTROP2']:.3f}",
    f"{nums['specificity_scan']['gtex_min_q_noTROP2']:.3f}",
    f"{ag_map['Status (yes / not_listed / no_data)']:.3f}",
    f"{ag_map['Numeric stom_any ICC(2,1)']:.3f}",
    f"{nums['clinical_three_datasets'][0]['LR']:.2f}",
    f"{nums['clinical_three_datasets'][1]['LR']:.2f}",
    str(nums['hpa_dab']['n_sections_total']),
    f"{dab_map['esoph_dab_basal']['LR']:.2f}",
    f"{dab_map['esoph_dab_basal']['perm_p']:.4f}",
    f"{dab_map['esoph_dab_basal']['holm_p_LR']:.3f}",
    f"{dab_map['esoph_dab_p75']['LR']:.2f}",
    f"{dab_map['oral_dab_basal']['LR']:.2f}",
    f"{dab_map['skin_dab_basal']['LR']:.2f}",
    f"{nums['hpa_dab']['ptprz1_esoph_dab']:.3f}",
    f"{nums['hpa_dab']['low_rna_esoph_dab_median']:.3f}",
    f"{nums['hpa_dab']['param_grid_81']['esoph_n_p05']}/{nums['hpa_dab']['param_grid_81']['n_grid']}",
    f"{nums['hpa_dab']['param_grid_81']['esoph_median_LR']:.2f}",
    f"{t3_map['FAT2']['ext_max_tpm']:.2f}",
    f"{t3_map['PTPRZ1']['ext_max_tpm']:.2f}",
    f"{t3_map['CDH3']['ext_max_tpm']:.2f}",
    f"{t3_map['PTK7']['ext_max_tpm']:.2f}",
    f"{t3_map['DSG2']['ext_max_tpm']:.2f}",
    f"{t3_map['ATP2B1']['ext_max_tpm']:.2f}",
]
for tok in expected_tokens:
    if tok not in md:
        errors.append(f'NUMERIC SYNC FAIL: Expected exact statistic token "{tok}" from numbers.json not found in manuscript.md')

# 9. Context-anchored numerical assertions (prevent wrong-sentence numeric matches)
pm = nums['primary_model']
cal0 = nums['calibration'][0]
cal1 = nums['calibration'][1]
loto_num_map = {r['drug']: r for r in nums['loto_numeric_predictions']}
es_dab = dab_map['esoph_dab_basal']
sf = nums['surfaceome']
t3_tac = t3_map['TACSTD2']
t3_cdh3 = t3_map['CDH3']

# Verify Table 1 rows are strictly non-increasing in LR
t1_block = tables_sec.split('### Table 1.')[1].split('### Table 2.')[0]
t1_lrs = []
for ln_t1 in t1_block.splitlines():
    if ln_t1.startswith('|') and not ln_t1.startswith('|---') and 'Modality / Resolution' not in ln_t1:
        cells_t1 = [c.strip().strip('*') for c in ln_t1.strip('|').split('|')]
        if len(cells_t1) >= 6:
            t1_lrs.append(float(cells_t1[5]))
if len(t1_lrs) != 21 or any(t1_lrs[i] < t1_lrs[i + 1] for i in range(len(t1_lrs) - 1)):
    errors.append(f'TABLE 1 FAIL: Table 1 must have 21 rows sorted in descending LR order (got {len(t1_lrs)} rows: {t1_lrs})')

d1t = nums['drop_one_target']
if d1t['n_targets_tested'] != 28 or f"{d1t['min_LR']:.2f}" not in res_sec or f"{d1t['alltissue_noTROP2_LR']:.2f}" not in res_sec:
    errors.append(f'SENSITIVITY SYNC FAIL: Drop-one-target or pseudobulk noTROP2 summary mismatch in numbers.json: {d1t}')

context_checks = [
    (
        'Primary basal model in Abstract & Results',
        rf"Oral basal RNA predicted stomatitis \(LR = {pm['W_basal_LR']:.2f}, β = {pm['W_basal_beta']:.2f}, target-cluster bootstrap 95% CI",
        md,
    ),
    (
        'Primary Tobit slope & LR in Results',
        rf"predicted clinical stomatitis incidence \(LR = {pm['W_basal_LR']:.2f}.*β = {pm['W_basal_beta']:.3f}.*ρ = {pm['W_basal_spearman']:.3f}",
        res_sec,
    ),
    (
        'Choi normal oral epithelium replication in Results',
        rf"Choi et al\. \(2023\).*LR = {pm['Choi_normal_epi_LR']:.2f}, P = {pm['Choi_normal_epi_p']:.4f}, ρ = {pm['Choi_normal_epi_spearman']:.3f}.*LR = {pm['Choi_panHNSCC_LR']:.2f}",
        res_sec,
    ),
    (
        'Zero-antigen-overlap LOTO metrics in Results',
        rf"LOTO out-of-target predictions improved.*Spearman ρ = {cal0['spearman_x']:.3f} \(P = {cal0['p_spearman_x']:.4f}, MAE = {cal0['mae_x']:.1f} vs\. {cal0['mae_payload']:.1f}|Spearman ρ = {cal0['spearman_x']:.3f} \(P = {cal0['p_spearman_x']:.4f}, MAE = {cal0['mae_x']:.1f}",
        res_sec,
    ),
    (
        'LOTO residual predictions (SG, iza-bren, sigvotatug, gemtuzumab) in Results',
        rf"sacituzumab govitecan \(observed 19\.0% vs\. {loto_num_map['sacituzumab govitecan']['pred_loto']:.1f}%.*izalontamab brengitecan \(observed 44\.6% vs\. {loto_num_map['izalontamab brengitecan']['pred_loto']:.1f}%\).*sigvotatug vedotin \(observed 22\.0% vs\. {loto_num_map['sigvotatug vedotin']['pred_loto']:.1f}%\).*gemtuzumab ozogamicin \(observed 21\.0% vs\. {loto_num_map['gemtuzumab ozogamicin']['pred_loto']:.1f}%\)",
        res_sec,
    ),
    (
        'TACSTD2 LOTO prediction & 95% PI in Results',
        rf"LOTO {t3_tac['loto_pred_stom_pct']:.1f}%, 95% PI {t3_tac['loto_pi_lo']:.1f}%–{t3_tac['loto_pi_hi']:.1f}%",
        res_sec,
    ),
    (
        'LOTO censored & log-score comparison in Results',
        rf"below-threshold placement was identical.*{cal0['cens_below_x']}/{cal0['n_censored']}.*overall log-score \([−-]{abs(cal0['mean_logscore_all_x']):.2f} vs\. [−-]{abs(cal0['mean_logscore_all_payload']):.2f}\)",
        res_sec,
    ),
    (
        'LOTO all-44 C-index & exact >=15% ROC-AUC in Results',
        rf"C-index = {cal0['c_index_all44_x']:.3f}.*P = {cal0['p_perm_c_index']:.4f}.*ROC-AUC = {cal0['auc_ge15_x']:.3f}.*P = {cal0['p_mw_auc_ge15']:.4f}",
        res_sec,
    ),
    (
        'LOO-ADC all-44 C-index & exact >=15% ROC-AUC in Results',
        rf"Exploratory LOO-ADC metrics yielded C-index = {cal1['c_index_all44_x']:.3f} versus {cal1['c_index_all44_payload']:.3f} across all 44 ADCs, exact ≥15% ROC-AUC = {cal1['auc_ge15_x']:.3f} versus {cal1['auc_ge15_payload']:.3f}",
        res_sec,
    ),
    (
        'FAERS negative-binomial oral replication in Results',
        rf"negative-binomial count regression replicated.*β = \+{nums['faers'][0]['beta']:.2f}.*LR = {nums['faers'][0]['LR']:.2f}, P = {nums['faers'][0]['p_LR']:.4f}.*P = {nums['faers'][0]['p_perm']:.4f}",
        res_sec,
    ),
    (
        'FAERS equal-weight OLS 12-target permutation in Results',
        rf"equal-weight empirical log-odds OLS.*β = \+{nums['faers'][0]['ols_beta']:.2f}, LR = {nums['faers'][0]['ols_LR']:.2f}, target-permutation P = {nums['faers'][0]['p_ols_perm']:.4f}",
        res_sec,
    ),
    (
        'Esophageal algorithm-defined basal DAB + Holm LR + 81/81 grid in Results',
        rf"mean basal DAB \(LR = {es_dab['LR']:.2f}, P = {es_dab['p_LR']:.4f}, 9-test LR Holm P = {es_dab['holm_p_LR']:.3f}, target-permutation P = {es_dab['perm_p']:.4f}",
        res_sec,
    ),
    (
        'OSCC surfaceome count (420), coupling (0.714), & winner funnel (47/12/4)',
        rf"Across these {sf['n_oscc_expressed']} OSCC surface proteins.*ρ = {sf['spearman_choi_basal']:.3f}.*Forty-seven of {sf['n_oscc_expressed']} targets.*and {sf['n_winners_vital15']} additionally.*{sf['n_winners_ext15']} of the {sf['n_winners_vital15']} candidates",
        res_sec,
    ),
    (
        'Table 3 coherent full-cohort TACSTD2 estimate & PI',
        rf"\|\s*\*{{2,3}}TACSTD2\*{{2,3}}\s*\|[^|\n]+\|[^|\n]+\|\s*{t3_tac['pred_stom_pct']:.1f} \({t3_tac['pred_pi_lo']:.1f}–{t3_tac['pred_pi_hi']:.1f}\)",
        tables_sec,
    ),
    (
        'Table 3 coherent full-cohort CDH3 estimate & PI',
        rf"\|\s*\*{{2,3}}CDH3\*{{2,3}}\s*\|[^|\n]+\|[^|\n]+\|\s*{t3_cdh3['pred_stom_pct']:.1f} \({t3_cdh3['pred_pi_lo']:.1f}–{t3_cdh3['pred_pi_hi']:.1f}\)",
        tables_sec,
    ),
]
for cname, cpat, ctext in context_checks:
    if not re.search(cpat, ctext, re.S):
        errors.append(f'CONTEXT SYNC FAIL ({cname}): pattern not matched in target section')

if errors:
    print('\nVERIFICATION FAILED:')
    for e in errors:
        print(' -', e)
    sys.exit(1)

print('ALL VERIFICATION CHECKS PASSED (0 mismatches, 0 formatting violations).')
