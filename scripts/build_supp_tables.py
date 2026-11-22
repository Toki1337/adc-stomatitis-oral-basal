"""Build multi-sheet Excel supplementary workbook and numbered Supplementary References."""

import glob
import json
import os
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import pandas as pd

ROOT = str(Path(__file__).resolve().parents[1])
R = f'{ROOT}/results/'
M = f'{ROOT}/manuscript/'
os.makedirs(M, exist_ok=True)

clin = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_trials.tsv', sep='\t')
ext_a = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_extractor_A.tsv', sep='\t')
ext_b = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_extractor_B.tsv', sep='\t')
reconciled = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_reconciled.tsv', sep='\t')
ds3_sbs = pd.read_csv(R + 'part1/clinical_three_dataset_side_by_side.tsv', sep='\t')
dual_ag = pd.read_csv(R + 'part1/clinical_dual_agreement.tsv', sep='\t')
adc_feat = pd.read_csv(R + 'part1/adc_table_with_features.tsv', sep='\t')
feat_comp = pd.read_csv(R + 'part1/feature_comparison.tsv', sep='\t')
res_ladder = pd.read_csv(R + 'part1/resolution_ladder.tsv', sep='\t')
dab_models = pd.read_csv(R + 'part1/hpa_ihc_dab_models.tsv', sep='\t')
dab_genes = pd.read_csv(R + 'part1/hpa_ihc_dab_genes.tsv', sep='\t')
dab_grid = pd.read_csv(R + 'part1/hpa_ihc_dab_param_grid.tsv', sep='\t')
tgt_inf = pd.read_csv(R + 'part1/target_inference.tsv', sep='\t')
sens = pd.read_csv(R + 'part1/sensitivity_subgroups.tsv', sep='\t')
scan = pd.read_csv(R + 'part1/celltype_tissue_scan.tsv', sep='\t').sort_values(['source', 'rank_all'])
axes_comp = pd.read_csv(R + 'part3_axes_comparison.tsv', sep='\t')
joint_comp = pd.read_csv(R + 'part3_joint_models.tsv', sep='\t')
ratio_sens = pd.read_csv(R + 'part3_ratio_sensitivity.tsv', sep='\t')
p3_axes = pd.read_csv(R + 'part3_axes.tsv', sep='\t')
p5 = pd.read_csv(R + 'part5_loto_predictions.tsv', sep='\t')
p5_exact = pd.read_csv(R + 'part5_loto_exact_endpoints.tsv', sep='\t')
faers_sum = pd.read_csv(R + 'fig_faers_summary.tsv', sep='\t')
faers_disprop = pd.read_csv(R + 'faers_per_drug_disprop.tsv', sep='\t')
faers_sens = pd.read_csv(R + 'faers_sensitivity.tsv', sep='\t')
faers_readus = pd.read_csv(R + 'faers_readus_baseline.tsv', sep='\t')
faers_pairs = pd.read_csv(R + 'fig_faers_pairs.tsv', sep='\t')
faers_pair_rows = pd.DataFrame([
    dict(
        category='Pairwise Oral ROR Contrast',
        item=f"{r.drug1} vs. {r.drug2} ({r.pair_type})",
        value=f"ROR {r.ror:.2f} (95% CI {r.lo:.2f}–{r.hi:.2f}); {int(r.k1)}/{int(r.n1):,} ({r.pct1:.2f}%) vs. {int(r.k2)}/{int(r.n2):,} ({r.pct2:.2f}%)",
    )
    for _, r in faers_pairs.iterrows()
])
faers_readus_full = pd.concat([faers_readus, faers_pair_rows], ignore_index=True)
p4 = pd.read_csv(R + 'part4_final_ranking.tsv', sep='\t', index_col=0).reset_index().rename(columns={'index': 'gene'})
p4['pass_vital_15tpm'] = p4['vital_max_tpm'] <= 15
p4 = p4.sort_values(['n_pass', 'grid_pass', 'pass_vital_15tpm', 'tum_choi'], ascending=[False, False, False, False])

# 1. Build numbered Supplementary References from local FDA label files, PMC XMLs, and PubMed abstract entries
JOURNAL_ISO_MAP = {
    'Blood': 'Blood',
    'Blood Advances': 'Blood Adv.',
    'Breast Cancer Research : BCR': 'Breast Cancer Res.',
    'Cancer': 'Cancer',
    'Cancer Chemotherapy and Pharmacology': 'Cancer Chemother. Pharmacol.',
    'Cancer Communications': 'Cancer Commun.',
    'Cancer Discovery': 'Cancer Discov.',
    'Cancer Research Communications': 'Cancer Res. Commun.',
    'Cell Reports Medicine': 'Cell Rep. Med.',
    'Clinical Cancer Research': 'Clin. Cancer Res.',
    'Clinical cancer research : an official journal of the American Association for Cancer Research': 'Clin. Cancer Res.',
    'Clinical lymphoma, myeloma & leukemia': 'Clin. Lymphoma Myeloma Leuk.',
    'Haematologica': 'Haematologica',
    'International Journal of Gynaecology and Obstetrics': 'Int. J. Gynaecol. Obstet.',
    'Investigational New Drugs': 'Invest. New Drugs',
    'JAMA Oncology': 'JAMA Oncol.',
    'JTO Clinical and Research Reports': 'JTO Clin. Res. Rep.',
    'Journal of Clinical Oncology': 'J. Clin. Oncol.',
    'NPJ Breast Cancer': 'npj Breast Cancer',
    'NPJ Precision Oncology': 'npj Precis. Oncol.',
    'Nature Communications': 'Nat. Commun.',
    'Nature Medicine': 'Nat. Med.',
    'Neuro-Oncology': 'Neuro Oncol.',
    'Signal Transduction and Targeted Therapy': 'Signal Transduct. Target. Ther.',
    'Targeted Oncology': 'Target. Oncol.',
    'The BMJ': 'BMJ',
    'The Lancet. Haematology': 'Lancet Haematol.',
    'The Lancet. Oncology': 'Lancet Oncol.',
}
PMC_NO_FULLTEXT_PMIDS = {'29196412', '31689155', '28817371'}


def normalize_pmc_entry(pmid, entry):
    e = dict(entry)
    e['journal'] = JOURNAL_ISO_MAP.get(e.get('journal', ''), e.get('journal', ''))
    if e.get('auth'):
        e['auth'] = e['auth'].replace('Trn\u0115n\u00fd', 'Trn\u011bn\u00fd')
    if pmid in PMC_NO_FULLTEXT_PMIDS:
        e['source_type'] = 'PubMed abstract (PMC record without full text)'
    return e


META_CACHE = f'{R}supp_refs_metadata.json'
label_files = sorted(glob.glob(f'{ROOT}/dl/labels/*.json'))
pmc_files = sorted(glob.glob(f'{ROOT}/dl/pmc/*.xml'))
if label_files and pmc_files:
    label_meta = {}
    for jf in label_files:
        data = json.load(open(jf, encoding='utf-8'))['results'][0]
        of = data.get('openfda', {})
        eff = str(data.get('effective_time', ''))
        eff_fmt = f'{eff[:4]}-{eff[4:6]}-{eff[6:8]}' if len(eff) == 8 else eff
        brands = of.get('brand_name') or []
        brand = brands[0].split()[0].upper() if brands else ('ZYNLONTA' if 'lonca' in jf else os.path.basename(jf))
        appls = of.get('application_number') or []
        appl = appls[0] if appls else ('BLA761196' if brand == 'ZYNLONTA' else '')
        set_id = data.get('set_id', '')
        label_meta[brand] = dict(brand=brand, eff_fmt=eff_fmt, appl=appl, set_id=set_id)

    pmc_meta = {}
    for xf in pmc_files:
        root = ET.parse(xf).getroot()
        am = root.find('.//article-meta')
        if am is None:
            am = root
        pmid = next((x.text for x in am.findall('./article-id[@pub-id-type="pmid"]') if x.text), '')
        if not pmid:
            continue
        doi = next((x.text for x in am.findall('./article-id[@pub-id-type="doi"]') if x.text), '')
        pmcid = os.path.basename(xf).replace('.xml', '')
        jtitle = next((x.text for x in root.findall('.//journal-meta//journal-title') if x.text), '')
        at_el = am.find('.//article-title')
        atitle = re.sub(r'\s+', ' ', ''.join(at_el.itertext()).strip()) if at_el is not None else ''
        yr = (
            next((x.text for x in am.findall('./pub-date[@pub-type="ppub"]/year') if x.text), '')
            or next((x.text for x in am.findall('./pub-date[@pub-type="collection"]/year') if x.text), '')
            or next((x.text for x in am.findall('./pub-date/year') if x.text), '')
        )
        vol = next((x.text for x in am.findall('./volume') if x.text), '')
        fpage = next((x.text for x in am.findall('./fpage') if x.text), '')
        lpage = next((x.text for x in am.findall('./lpage') if x.text), '')
        eloc = next((x.text for x in am.findall('./elocation-id') if x.text), '')
        auths = []
        for c in am.findall('.//contrib[@contrib-type="author"]'):
            sn = c.findtext('name/surname') or ''
            gn = c.findtext('name/given-names') or ''
            if sn:
                inits = '. '.join([p[0] for p in gn.replace('-', ' ').split() if p]) + ('.' if gn else '')
                auths.append(f'{sn}, {inits}'.strip())
        auth_str = (auths[0] + ' et al.') if len(auths) > 1 else (auths[0] if auths else '')
        loc = eloc if eloc else (f'{fpage}–{lpage}' if fpage and lpage else fpage)
        body_el = root.find('.//body')
        has_body = body_el is not None and len(''.join(body_el.itertext()).strip()) >= 200
        src_type = 'PMC full text' if has_body else 'PubMed abstract (PMC record without full text)'
        pmc_meta[pmid] = normalize_pmc_entry(pmid, dict(
            pmid=pmid, pmcid=pmcid, doi=doi, auth=auth_str,
            title=atitle, journal=jtitle, year=yr, vol=vol, pages=loc,
            source_type=src_type,
        ))
    with open(META_CACHE, 'w', encoding='utf-8') as f:
        json.dump({'label_meta': label_meta, 'pmc_meta': pmc_meta}, f, indent=2, sort_keys=True)
elif os.path.exists(META_CACHE):
    with open(META_CACHE, encoding='utf-8') as f:
        _meta_cache = json.load(f)
    label_meta = _meta_cache['label_meta']
    pmc_meta = {k: normalize_pmc_entry(k, v) for k, v in _meta_cache['pmc_meta'].items()}
else:
    raise FileNotFoundError(
        f'Missing both raw label/PMC downloads ({ROOT}/dl/) and cached metadata ({META_CACHE})'
    )

# Non-PMC PubMed trial publications verified directly via PubMed
NON_PMC_PUBMED = {
    '25925619': dict(auth='Palanca-Wessels, M. C. et al.', title='Safety and activity of the anti-CD79B antibody-drug conjugate polatuzumab vedotin in relapsed or refractory B-cell non-Hodgkin lymphoma and chronic lymphocytic leukaemia: a phase 1 study', journal='Lancet Oncol.', year='2015', vol='16', pages='704–715', doi='10.1016/S1470-2045(15)70128-2', source_type='PubMed abstract (no PMC full text)'),
    '31540980': dict(auth='Gerber, D. E. et al.', title='Phase Ia study of anti-NaPi2b antibody-drug conjugate lifastuzumab vedotin DNIB0600A in patients with non-small cell lung cancer and platinum-resistant ovarian cancer', journal='Clin. Cancer Res.', year='2020', vol='26', pages='364–372', doi='10.1158/1078-0432.CCR-18-3965', source_type='PubMed abstract (no PMC full text)'),
    '29401246': dict(auth='Banerjee, S. et al.', title='Anti-NaPi2b antibody-drug conjugate lifastuzumab vedotin (DNIB0600A) compared with pegylated liposomal doxorubicin in patients with platinum-resistant ovarian cancer in a randomized, open-label, phase II study', journal='Ann. Oncol.', year='2018', vol='29', pages='917–923', doi='10.1093/annonc/mdy023', source_type='PubMed abstract (no PMC full text)'),
    '39442070': dict(auth='Turner, N. et al.', title='Trastuzumab duocarmazine in pretreated human epidermal growth factor receptor 2-positive advanced or metastatic breast cancer: an open-label, randomized, phase III trial (TULIP)', journal='J. Clin. Oncol.', year='2025', vol='43', pages='513–523', doi='10.1200/JCO.24.00529', source_type='PubMed abstract (no PMC full text)'),
    '31257177': dict(auth='Banerji, U. et al.', title='Trastuzumab duocarmazine in locally advanced and metastatic solid tumours and HER2-expressing breast cancer: a phase 1 dose-escalation and dose-expansion study', journal='Lancet Oncol.', year='2019', vol='20', pages='1124–1135', doi='10.1016/S1470-2045(19)30328-6', source_type='PubMed abstract (no PMC full text)'),
    '30930134': dict(auth='Jagannath, S. et al.', title='Indatuximab ravtansine (BT062) monotherapy in patients with relapsed and/or refractory multiple myeloma', journal='Clin. Lymphoma Myeloma Leuk.', year='2019', vol='19', pages='372–380', doi='10.1016/j.clml.2019.02.006', source_type='PubMed abstract (no PMC full text)'),
    '18203652': dict(auth='Riechelmann, H. et al.', title='Phase I trial with the CD44v6-targeting immunoconjugate bivatuzumab mertansine in head and neck squamous cell carcinoma', journal='Oral Oncol.', year='2008', vol='44', pages='823–829', doi='10.1016/j.oraloncology.2007.10.009', source_type='PubMed abstract (no PMC full text)'),
    '17062682': dict(auth='Tijink, B. M. et al.', title='A phase I dose escalation study with anti-CD44v6 bivatuzumab mertansine in patients with incurable squamous cell carcinoma of the head and neck or esophagus', journal='Clin. Cancer Res.', year='2006', vol='12', pages='6064–6072', doi='10.1158/1078-0432.CCR-06-0910', source_type='PubMed abstract (no PMC full text)'),
    '17351401': dict(auth='Rupp, U. et al.', title='Safety and pharmacokinetics of bivatuzumab mertansine in patients with CD44v6-positive metastatic breast cancer: final results of a phase I study', journal='Anticancer Drugs', year='2007', vol='18', pages='477–485', doi='10.1097/CAD.0b013e32801403f4', source_type='PubMed abstract (no PMC full text)'),
    '40992225': dict(auth='Zhang, S. et al.', title='Efficacy and safety of MRG002 monotherapy in treating patients with locally advanced or metastatic breast cancer with low HER2 expression: a multi-center, non-randomized, open-label phase II clinical trial', journal='Eur. J. Cancer', year='2025', vol='230', pages='115794', doi='10.1016/j.ejca.2025.115794', source_type='PubMed abstract (no PMC full text)'),
    '38772302': dict(auth='Qu, W. et al.', title='A single-arm, multicenter, phase 2 clinical study of recombinant humanized anti-HER2 monoclonal antibody-MMAE conjugate (MRG002) in HER2-positive unresectable locally advanced or metastatic urothelial carcinoma', journal='Eur. J. Cancer', year='2024', vol='205', pages='114096', doi='10.1016/j.ejca.2024.114096', source_type='PubMed abstract (no PMC full text)'),
    '27793850': dict(auth='Liu, J. F. et al.', title='Phase I study of safety and pharmacokinetics of the anti-MUC16 antibody-drug conjugate DMUC5754A in patients with platinum-resistant ovarian cancer or unresectable pancreatic cancer', journal='Ann. Oncol.', year='2016', vol='27', pages='2124–2130', doi='10.1093/annonc/mdw401', source_type='PubMed abstract (no PMC full text)'),
    '41791381': dict(auth='Duan, J. et al.', title='HS-20093, a B7-H3-targeted antibody-drug conjugate in lung cancer: Results from the ARTEMIS-001 phase 1a/b trial', journal='Cancer Cell', year='2026', vol='44', pages='846–857.e3', doi='10.1016/j.ccell.2026.02.006', source_type='PubMed abstract (no PMC full text)'),
}

ref_rows = []
key_to_sid = {}


def get_or_add_ref(key, row_dict):
    if key not in key_to_sid:
        sid = f'S{len(ref_rows) + 1}'
        key_to_sid[key] = sid
        row_dict['supp_ref_id'] = sid
        ref_rows.append(row_dict)
    return key_to_sid[key]


row_ref_ids = []
for _, r in clin.iterrows():
    src = str(r.source_id) if pd.notna(r.source_id) else ''
    if not src:
        sid = get_or_add_ref('NO_CLIN_PUB', dict(
            source_type='No monotherapy safety report',
            identifier='PubMed / openFDA search (2026-11)',
            effective_or_year='2026-11',
            doi_or_set_id='N/A',
            citation='Systematic PubMed and openFDA drug-label search conducted November 2026 identified no peer-reviewed monotherapy clinical safety publication for each of the 12 excluded investigational ADCs citing S17 in Table S1.',
        ))
        row_ref_ids.append(sid)
    elif 'FDA label' in src:
        m_brand = re.search(r'FDA label\s+([A-Za-z0-9-]+)', src)
        brand = m_brand.group(1).upper() if m_brand else src
        meta = label_meta.get(brand, {})
        eff = meta.get('eff_fmt', '')
        if not eff:
            m_eff = re.search(r'effective\s+([0-9-]+)', src)
            eff = m_eff.group(1) if m_eff else ''
        appl = meta.get('appl', '')
        set_id = meta.get('set_id', '')
        appl_txt = f' ({appl})' if appl else ''
        sid_txt = f' DailyMed/openFDA set identifier: {set_id}.' if set_id else ''
        sid = get_or_add_ref(f'LABEL:{brand}', dict(
            source_type='FDA prescribing information',
            identifier=f'{brand}{appl_txt}',
            effective_or_year=f'effective {eff}',
            doi_or_set_id=set_id,
            citation=f'U.S. Food and Drug Administration. {brand} ({r.drug}) prescribing information{appl_txt}; openFDA label effective {eff}.{sid_txt}',
        ))
        row_ref_ids.append(sid)
    else:
        pmids = re.findall(r'PMID\s*(\d+)', src)
        ids_for_row = []
        for pmid in pmids:
            info = pmc_meta.get(pmid) or NON_PMC_PUBMED.get(pmid)
            if info:
                vol_pg = f" **{info['vol']}**" if info.get('vol') else ''
                pg = f", {info['pages']}" if info.get('pages') else ''
                pmcid_txt = f" ({info['pmcid']})" if info.get('pmcid') else ''
                cit = f"{info['auth']} {info['title']}. *{info['journal']}*{vol_pg}{pg} ({info['year']}). PMID {pmid}{pmcid_txt}; doi:{info['doi']}."
                sid = get_or_add_ref(f'PMID:{pmid}', dict(
                    source_type=info['source_type'],
                    identifier=f'PMID {pmid}',
                    effective_or_year=info['year'],
                    doi_or_set_id=info['doi'],
                    citation=cit,
                ))
                ids_for_row.append(sid)
        row_ref_ids.append('; '.join(ids_for_row) if ids_for_row else '')

clin['supp_ref_id'] = row_ref_ids
supp_refs_df = pd.DataFrame(ref_rows)[['supp_ref_id', 'source_type', 'identifier', 'effective_or_year', 'doi_or_set_id', 'citation']]

merge_cols = [
    'drug', 'pclass', 'cens', 'y_pct', 'W_basal', 'W_basal_quiescent', 'W_basal_prolif',
    'W_basal_det', 'W_alltissue', 'HPA_sc_esoph_basal', 'HPA_esoph_DAB', 'HPA_oral_DAB', 'HPA_skin_DAB',
    'GTEx_esoph_mucosa', 'TCGA_normal_oral', 'CPTAC_NAT_protein',
]
rec_cols = ['drug', 'Comparison_Status', 'Adjudicated_By', 'Adjudication_Note']
s1 = clin.merge(reconciled[rec_cols], on='drug', how='left').merge(adc_feat[merge_cols], on='drug', how='left').rename(columns={'W_basal_prolif': 'W_basal_cycling'})
# Place supp_ref_id cleanly right after source_id
cols_s1 = list(s1.columns)
cols_s1.remove('supp_ref_id')
src_idx = cols_s1.index('source_id') + 1
cols_s1.insert(src_idx, 'supp_ref_id')
s1 = s1[cols_s1]

s5_adcs = p3_axes[[
    'drug', 'target_gene', 'tg', 'pclass', 'y_pct', 'cens',
    'abs_oral_basal', 'tum_oscc', 'ratio_oscc', 'tum_puram', 'ratio_puram',
    'ratio_choi_internal', 'ratio_choi_pc01', 'ratio_choi_pc001',
    'ratio_choi_linear_pc01', 'tcga_paired_lfc',
]]

s6_cols = [
    'drug', 'target_gene', 'tg', 'pclass', 'payload_class', 'y_pct', 'cens', 'abs_oral_basal',
    'pred_loto', 'pred_loto_payload_only', 'pred_loto_4payload',
    'sigma_loto', 'sigma_loto_payload', 'log_score_loto', 'log_score_loto_payload',
    'pred_loo', 'pred_payload_only', 'pred_loo_4payload',
    'log_score_loo', 'log_score_loo_payload',
]

p5_sum = pd.read_csv(R + 'part5_calib_summary.tsv', sep='\t')
s3c_protocol = pd.DataFrame([
    dict(
        plan_item='Y1 (Primary clinical endpoint)',
        prespecified_plan='Any-grade stomatitis/oral mucositis % in the largest monotherapy safety population (label > phase 3 > phase 2 > phase 1). Not listed in AE table -> left-censored below the table threshold.',
        implemented_analysis='Applied strictly to openFDA labels + PMC full texts (with PubMed abstracts for two left-censored entries) across 66 ADCs (Wendi Dou pre-linkage, Xiang Yu blinded, consensus adjudication; 44 evaluable: 15 numeric, 29 left-censored). Composite oral/pharyngeal row used for sigvotatug (22%; sensitivity at 11% single PT LR=15.47); abstract >20% lists used for sofituzumab & vandortuzumab (sensitivity excluding n=42, LR=13.40).',
        output_location='Results (Study Design; Absolute Oral Basal Expression); Methods (Blinded Dual Clinical Stomatitis Data Extraction); Table S1, Tables S1a–e, Table S3b',
        status='Implemented as planned (with minor operational rules & S3b checks)',
    ),
    dict(
        plan_item='Y2 & Negative Controls (FAERS)',
        prespecified_plan='FAERS reporting proportion of oral-mucositis PTs (STOMATITIS, MUCOSAL INFLAMMATION, MOUTH ULCERATION, APHTHOUS ULCER, GLOSSITIS, ORAL MUCOSAL ERYTHEMA, ORAL MUCOSAL BLISTERING) among reports where the ADC is the only suspect drug and no other antineoplastic is listed; negative controls peripheral neuropathy and neutropenia PTs.',
        implemented_analysis='Exact FDA FAERS Quarterly ASCII extracts (56 quarters, 2012Q4–2026Q3; 17,994,015 deduplicated cases after removing 99,817 FDA-deleted case IDs and retaining latest CASEVERSION/FDA_DT). Primary cohort restricts strictly to deduplicated sole-suspect (ROLE_COD=PS with zero SS) monotherapy ADC reports excluding all co-reported antineoplastic/anti-HER2 agents (14 FDA-approved ADCs, 12 unique targets, 37,051 cases, 448 oral events; disitamab vedotin had 0 PS reports), evaluated via payload-adjusted exposure-offset NB2 regression (10,000 target permutations), unadjusted Spearman, equal-weight OLS, whole-FAERS ROR/PRR/BCPNN IC025, 2 prespecified + 4 exploratory negative controls, and 6 READUS-PV sensitivity cohorts (Tables S8a–d).',
        output_location='Results (Basal Expression Generalizes to Held-Out Targets and Pharmacovigilance); Methods (Statistical Modeling, Calibration, FAERS, and Surfaceome Re-Ranking); Figure 4C–D, Tables S8a–d',
        status='Implemented as planned (expanded with full FDA deduplication, whole-FAERS ROR/IC025, 4 exploratory controls, and READUS-PV sensitivity suite)',
    ),
    dict(
        plan_item='Primary Predictor (X) & Replications',
        prespecified_plan='Target-gene expression in basal keratinocytes of healthy oral mucosa (GSE164241, healthy buccal + gingival donors only; pseudobulk log2(mean CP10K + 1)). Bispecific -> max of the two genes. CD44v6 -> CD44. Replications: Choi GSE181919 normal oral mucosa keratinocytes; HPA IHC oral mucosa staining level; TCGA-HNSC adjacent-normal oral bulk RNA; CPTAC HNSCC NAT proteome.',
        implemented_analysis='All five computed identically (with legacy gene symbol PVRL4 harmonized to HGNC symbol NECTIN4 across expression matrices). Choi normal epithelium restricted to oral-cavity subsite (choi_epi_NL_OC_mean, n=199 cells, 3 donors; LR=6.34, P=0.0118) with pan-HNSCC normal epithelium (n=810 cells, 9 donors; LR=8.14) retained alongside in Table S2a.',
        output_location='Results (Absolute Oral Basal Expression; Epithelial Resolution Gradient); Methods (Single-Cell, Bulk RNA, Proteomic, Quantitative IHC, and Reference Atlas Processing); Table 1, Table S2a',
        status='Implemented as planned',
    ),
    dict(
        plan_item='Test 1 & Test 2 (Primary & LOTO)',
        prespecified_plan='Test 1: Spearman(X, Y1) across ADCs; censored values set to threshold/2 (sensitivity: threshold, 0, Tobit on logit scale). Test 2: logit(Y1) ~ X + payload class (topo-I / tubulin / DNA-damaging); report X coefficient, likelihood-ratio test, and leave-one-target-out.',
        implemented_analysis='Exact 3-level payload Tobit (LR=13.55, beta=0.994), all 3 Spearman imputations (rho=0.457-0.535), and prespecified leave-one-target-out (implemented with zero-antigen-overlap folds across 44 ADCs including bispecific EGFR;ERBB3 exclusion, plus exploratory distribution-free all-44 C-index, exact uncensored >=10%/15%/20% ROC-AUCs, out-of-fold Tobit log-scores, and 44-fold LOO-ADC in Tables S6a, S6b, S6c).',
        output_location='Results (Absolute Oral Basal Expression; Leave-One-Target-Out Calibration); Methods (Statistical Modeling); Figure 2B, Figure 4A–B, Table 1, Table S2a, Table S3b, Table S6a, Table S6b, Table S6c',
        status='Implemented as planned (expanded with target permutation, cluster bootstrap, & log-scores)',
    ),
    dict(
        plan_item='Test 3 (Within-payload contrasts)',
        prespecified_plan='Within-payload contrasts (DXd, MMAE, DM1/DM4, SN-38).',
        implemented_analysis='Power-feasible subclasses fit directly: Topo-I (n=12, 7 num; LR=7.59), Tubulin (n=25, 4 num; LR=4.71), MMAE (n=11, 3 num; LR=2.84, P=0.092), post hoc Topo-I excl SN-38 (n=10). Too few numeric events for regression in DM1/DM4 (n=7, 1 num), SN-38 (n=2, 1 num), and DXd (n=4, 3 num; reported descriptively).',
        output_location='Results (Absolute Oral Basal Expression); Methods (Prespecified Protocol vs. Exploratory and Post Hoc Analyses); Figure 3D, Table S3b',
        status='Partial due to sample size (insufficient numeric rows in DM1/DM4, SN-38, DXd)',
    ),
    dict(
        plan_item='Test 4 (Tissue specificity scan)',
        prespecified_plan='Repeat test 1 with X from each GTEx tissue; expect squamous mucosa (esophagus mucosa, minor salivary gland, skin) to rank top.',
        implemented_analysis='Reported both prespecified Test 1 Spearman(X, Y1) across all 54 GTEx tissues (11/54 BH-FDR q<0.05 prior to payload adjustment, including non-squamous tissues) and payload-adjusted Tobit LR with BH-FDR across 54 GTEx tissues (isolating 4/54 stratified squamous tissues at q<0.05: Esophagus-Mucosa #1, Skin #2/#4, Vagina #3) in Table S4, plus exploratory 155-cell-type HPA scan.',
        output_location='Results (Epithelial Resolution Gradient, Quantitative Squamous IHC, and Tissue Specificity); Methods (Prespecified Protocol; Statistical Modeling); Figure 3B–C, Table S4',
        status='Implemented (expanded with payload-adjusted Tobit LR + BH-FDR alongside Spearman)',
    ),
    dict(
        plan_item='Test 5 & Part 2 (OSCC window & surfaceome)',
        prespecified_plan='Test 5: Y2 replication of tests 1-2; negative-control outcomes. Part 2: for ADC targets and surfaceome genes, tumour (TCGA oral, CPTAC tumour proteome, malignant cells in GSE103322 / GSE181919) vs normal oral basal keratinocytes; predicted stomatitis liability from the part-1 model; window = tumour excess over predicted liability; DepMap OSCC lines for model selection.',
        implemented_analysis='FAERS Tests 1–2 & controls run on FDA-deduplicated sole-suspect monotherapy quarterly ASCII extracts (2012Q4–2026Q3, n=14 FDA-approved ADCs across 12 unique targets, 37,051 cases); 420-gene OSCC surfaceome (with PVRL4 harmonized to NECTIN4) profiled across Choi/Puram/TCGA/CPTAC/GTEx-4 & 24-organ/GTEx-54/DepMap 24Q4 (1,178 cell lines) using coherent full-cohort predictions & LOTO sensitivity columns.',
        output_location='Results (Leave-One-Target-Out Calibration and Pharmacovigilance Replication; OSCC Surfaceome Coupling and Exploratory Target Re-Ranking); Methods (Surfaceome Curation); Figure 4C–D, Figure 5, Table 3, Table S7',
        status='Deviated: DepMap used for pan-cell-line common-essential screening (1,178 lines) instead of OSCC-line model selection; window computed as tumor − oral basal Δlog₂ alongside 4-criterion re-ranking & 3x3 grids',
    ),
])

s1_export = s1.rename(columns={'tumour_types': 'tumor_types'})
ext_a_export = ext_a.rename(columns={'tumour_types': 'tumor_types'})
reconciled_export = reconciled.rename(columns={'tumour_types': 'tumor_types'})
p4_export = p4.drop(columns=['robust_7of9', 'robust_ext_7of9'], errors='ignore')

readme_rows = [
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S1_Clinical_66ADCs', column_or_title='Table S1: Clinical Stomatitis Cohort (66 ADCs)',
         units='%, N patients, log2(CP10K+1), TPM, OD',
         description=f'Adjudicated monotherapy any-grade stomatitis/oral mucositis dataset across all {len(s1)} clinical-stage and approved ADCs (44 evaluable across 28 targets: 15 numeric, 29 left-censored; 22 excluded), with source hierarchy, exact quotes, linked expression features, and numbered Supplementary Reference IDs (supp_ref_id: S1-S60).'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S1a_Extractor_A_66ADCs', column_or_title='Table S1a: Original Pre-Linkage Clinical Extraction (Wendi Dou)',
         units='%, N patients',
         description=f'Independent pre-linkage extraction by Wendi Dou across all {len(ext_a)} ADCs finalized prior to linking clinical rows to expression matrices.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S1b_Extractor_B_66ADCs', column_or_title='Table S1b: Blinded Second Clinical Extraction (Xiang Yu)',
         units='%, N patients',
         description=f'Independent second extraction by Xiang Yu across all {len(ext_b)} ADCs blinded to gene/protein expression matrices and Wendi Dou entries.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S1c_Three_Dataset_Results', column_or_title='Table S1c: Parallel Modeling Across Three Clinical Dataset Versions',
         units='LR (df=1), P, rho, percentage points',
         description='Side-by-side comparison of primary Tobit, no-TACSTD2 sensitivity, target permutation, and LOTO metrics across Wendi Dou original, Xiang Yu unadjudicated, and Adjudicated Final datasets.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S1d_Reconciled_66ADCs', column_or_title='Table S1d: Dual-Extractor Reconciliation & Adjudication Log (66 ADCs)',
         units='%, N patients',
         description='Side-by-side comparison of Wendi Dou, Xiang Yu, and Adjudicated Final fields with discrepancy flags and source-verified adjudication rationale.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S1e_Dual_Agreement', column_or_title='Table S1e: Inter-Extractor Agreement Statistics',
         units='Cohen kappa, weighted kappa, ICC(2,1)/ICC(2,k), MAE (pp)',
         description='Pre-reconciliation inter-extractor reliability metrics across all 66 ADCs (status, evaluability, source tier, censoring threshold, and shared numeric rates).'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S2a_Features', column_or_title='Table S2a: Univariate Payload-Adjusted Tobit Models (21 Features)',
         units='logit slope beta, LR (df=1), Spearman rho',
         description='Payload-adjusted left-censored logit Tobit regressions across 21 oral, squamous, and reference expression features (matching Main Table 1).'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S2b_Resolution_Ladder', column_or_title='Table S2b: 15-Modality Squamous Resolution Ladder & Pairwise Delta-LR CIs',
         units='LR (df=1), AIC, delta-LR (95% target-cluster bootstrap CI)',
         description='Comparison of 15 cellular resolution and measurement modalities against Williams oral basal keratinocytes with matched-subset B=2,000 target-cluster bootstrap 95% CIs for delta-LR.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S2c_HPA_DAB_Models', column_or_title='Table S2c: Deterministic Ruifrok-Johnston rgb2hed DAB IHC Models',
         units='Optical density (OD), LR (df=1), Holm & permutation P',
         description='Payload-adjusted Tobit models across 9 tissue-by-summary DAB metrics (oral mucosa, esophagus, skin x mean basal, 75th-percentile basal, whole epithelium), Enhanced/Supported antibody restriction, and epithelial-fraction checks.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S2d_HPA_DAB_Genes', column_or_title='Table S2d: Gene-Level HPA Normal Squamous DAB Optical Densities (28 Genes)',
         units='OD, epithelial area fraction (0-1), log2(CP10K+1)',
         description='Per-gene mean/median DAB optical densities, core counts, epithelial area fractions, and antibody reliability ratings across 509 HPA v26.0 normal oral mucosa, esophagus, and skin TMA sections.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S2e_HPA_DAB_ParamGrid', column_or_title='Table S2e: 81-Combination Algorithmic Segmentation Stability Grid',
         units='Quantile, ratio, px, grayscale (0-255), LR (df=1)',
         description='Sensitivity grid across 3x3x3x3=81 segmentation parameter combinations for esophageal and oral basal DAB Tobit models (esophagus: 81/81 P<0.05).'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S3a_Target_Inference', column_or_title='Table S3a: Target-Unit Permutation & Cluster-Bootstrap Inference',
         units='LR (df=1), logit beta, P',
         description='28-target permutation (B=10,000) and target-cluster bootstrap (B=2,000) inference for primary 3-level and post hoc 4-level payload Tobit models.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S3b_Sensitivity_Subgroups', column_or_title='Table S3b: Sensitivity, Drop-One-Target, Payload & Outcome-Heterogeneity Checks',
         units='logit beta (95% profile or cluster-boot CI), LR (df=1), P',
         description='Complete sensitivity suite including all 28 drop-one-target exclusions, payload-family subsets, 7 outcome-heterogeneity checks, imputation/site/scaling checks, and pseudobulk no-TACSTD2 check.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S3c_Protocol_Alignment', column_or_title='Table S3c: Prespecified Protocol-to-Implementation Alignment',
         units='Text audit log',
         description='Line-by-line audit mapping every item in ANALYSIS_PLAN.md to its implemented analysis, deviations/extensions, and manuscript/table locations.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S4_CellType_Tissue_Scan', column_or_title='Table S4: Unbiased HPA (155 Cell Types) & GTEx (54 Bulk Tissues) Scans',
         units='LR (df=1), Spearman rho, BH-FDR q',
         description='Full unbiased specificity scans across 155 HPA v26.0 single-cell types and 54 GTEx v8 bulk tissues (both full 44-ADC cohort and excluding TACSTD2, plus prespecified unadjusted Spearman).'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S5a_Axes_Single', column_or_title='Table S5a: Single-Predictor Head-to-Head Axis Models (Absolute vs. Ratios)',
         units='logit beta, LR (df=1), Spearman rho',
         description='Payload-adjusted Tobit models swapping only the exposure predictor across absolute oral basal expression, OSCC tumor expression, and tumor-to-normal selectivity ratios (Main Table 2).'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S5b_Joint_Nested', column_or_title='Table S5b: Nested Joint Models (Absolute Basal + Selectivity Ratio + Payload)',
         units='logit beta, added LR (df=1), P',
         description='Direct nested likelihood-ratio tests comparing incremental explanatory value of adding selectivity ratio over absolute basal abundance versus adding basal abundance over ratio.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S5c_Ratio_Scales', column_or_title='Table S5c: Selectivity Ratio Pseudocount, Linear & Basal-Subset Sensitivity',
         units='logit beta, LR (df=1), P',
         description='Sensitivity tests of tumor-to-normal ratios across log pseudocounts (1.0, 0.1, 0.01), linear fold ratios, and subsets restricted to basally expressed targets.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S5d_Axes_Per_ADC', column_or_title='Table S5d: Per-ADC Exposure Axes Across All 44 Modeled ADCs',
         units='%, log2(CP10K+1), log2 fold ratio',
         description='Per-ADC clinical stomatitis outcomes and aligned absolute, tumor, and selectivity-ratio exposure values for all 44 evaluable ADCs.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S6a_LOTO_LOO_Calib', column_or_title='Table S6a: Per-ADC LOTO & LOO-ADC Out-of-Sample Predictions (44 ADCs)',
         units='%, logit SD, logit log-score',
         description='Zero-antigen-overlap leave-one-target-out (LOTO) and exploratory 44-fold leave-one-ADC-out (LOO-ADC) predicted stomatitis rates for basal+payload, payload-only, and 4-level payload models.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S6b_LOTO_Exact_Endpoints', column_or_title='Table S6b: Distribution-Free LOTO Concordance & Exact Uncensored Binary Endpoints',
         units='C-index, ROC-AUC, PR-AUC, P (Mann-Whitney & 5,000-perm)',
         description='Distribution-free left-censored C-index across all 451 ordered pairs and exact uncensored binary endpoint discrimination at >=10% (35/44), >=15% (38/44), and >=20% (44/44) stomatitis.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S6c_Calib_Summary', column_or_title='Table S6c: Summary Calibration Metrics (LOTO vs. LOO-ADC)',
         units='C-index, ROC-AUC, PR-AUC, rho, MAE (pp), McNemar P, log-score',
         description='Aggregated out-of-sample calibration metrics comparing basal+payload, payload-only, and 4-level payload models across LOTO and LOO-ADC schemes.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S7_OSCC_420_Surfaceome', column_or_title='Table S7: Complete Re-Ranking of 420 OSCC-Expressed Surfaceome Proteins',
         units='log2(CP10K+1), %, TPM, log2 fold change, Chronos score, /9 grid',
         description='Multi-omic OSCC surfaceome prioritization table across all 420 dual-cohort expressed surface proteins (after PVRL4-to-NECTIN4 harmonization), including predicted stomatitis with 95% CIs/PIs, LOTO sensitivity estimates, GTEx 4-organ/24-organ/54-tissue maxima, 3x3 stability counts, and DepMap 24Q4 essentiality.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S8a_FAERS_Models', column_or_title='Table S8a: FDA FAERS Quarterly ASCII Pharmacovigilance Models (14 FDA-Approved ADCs, 12 Targets)',
         units='log-rate beta, rate ratio, LR (df=1), P (Wald, LR, 10,000-perm), Spearman rho',
         description='Primary payload-adjusted exposure-offset NB2 count regression alongside unadjusted Spearman and equal-weight empirical log-odds OLS sensitivity models across oral mucositis/stomatitis and 6 control toxicities in deduplicated sole-suspect monotherapy cases (2012Q4–2026Q3, N=37,051).'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S8b_FAERS_Drugs_Disprop', column_or_title='Table S8b: Per-ADC Whole-FAERS Disproportionality (ROR, PRR, BCPNN IC025)',
         units='counts, %, ROR (95% CI), PRR, BCPNN IC / IC025',
         description='Per-ADC deduplicated sole-suspect monotherapy case counts, oral mucositis proportions, and whole-FAERS disproportionality metrics (ROR, PRR, BCPNN IC and IC025 against 17,994,015 deduplicated 2012Q4–2026Q3 background cases) across all 14 FDA-approved ADCs.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S8c_FAERS_Sensitivity', column_or_title='Table S8c: READUS-PV FAERS Sensitivity Suite (7 Cohorts)',
         units='log-rate beta, LR (df=1), 10,000-perm P, Spearman rho, pairwise ROR',
         description='Oral mucositis regression and pairwise ROR stability across primary sole-suspect monotherapy, all-PS monotherapy, HCP-only, serious-only, non-HNSCC/radiotherapy indications, secondary demographic/event/drug/PT deduplication, and unexcluded all-PS reports.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Table_S8d_FAERS_READUS_Flow', column_or_title='Table S8d: READUS-PV Case Attrition Funnel, Baseline Demographics, and Pairwise ROR Contrasts',
         units='case counts, median (IQR) years, %, ROR (95% CI)',
         description='Official FDA DELETE and CASEID deduplication attrition funnel (2012Q4–2026Q3), primary sole-suspect monotherapy cohort demographics (age, sex, reporter qualification, seriousness, top countries), and all 5 pairwise oral ROR contrasts.'),
    dict(section='WORKBOOK INDEX', sheet_or_table='Supplementary_References', column_or_title='Supplementary References (S1-S60)',
         units='Citations (BLA/set_id, PMID/PMCID/DOI)',
         description='Numbered clinical safety sources corresponding to supp_ref_id in Table S1 across all 66 ADCs.'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S1 / S1a / S1b / S1d', column_or_title='tier / listed / cens / stom_any / stom_g3 / threshold_pct / y_pct / pclass',
         units='Category, boolean, incidence (%)',
         description='Clinical curation fields: tier (1 = approved or late-phase pivotal ADC; 2 = clinical-stage investigational ADC); listed ("yes" = explicit oral mucosal AE rate reported; "not_listed" = left-censored below AE table cutoff; "no_data" = excluded, no monotherapy table with N>=15); cens (True = left-censored observation, False = exact numeric observation); stom_any / stom_g3 (any-grade and grade >=3 stomatitis/oral mucositis %); threshold_pct (AE table reporting cutoff % for left-censored trials); y_pct (modeled outcome: stom_any if numeric, threshold_pct if left-censored); pclass (3-level payload family: topoI, tubulin, DNA).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S2c_HPA_DAB_Models', column_or_title='holm_p_LR / holm_perm_p / bh_p_LR',
         units='Adjusted two-sided P value',
         description='Multiplicity-adjusted P values across the 9 primary tissue-by-summary DAB models: Holm-adjusted LR P (holm_p_LR), Holm-adjusted target-permutation P (holm_perm_p), and BH-FDR adjusted LR P (bh_p_LR).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S2b / Table_S2c', column_or_title='dLR_vs_Wbasal / dLR_boot_ci_lo / dLR_boot_ci_hi',
         units='LR difference (df=1)',
         description='Delta-LR = LR(Williams oral basal scRNA) minus LR(comparator feature) refitted on the identical non-missing ADC subset; positive values favor oral basal scRNA. CIs from B=2,000 target-cluster bootstrap.'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S2e_HPA_DAB_ParamGrid', column_or_title='q_epi / r_pos / d_px / g_thr / is_default',
         units='Quantile (0-1), distance ratio (0-1), pixels, grayscale (0-255), boolean',
         description='Algorithmic squamous IHC segmentation grid parameters: q_epi (tissue cellularity quantile threshold: 0.50, 0.58, 0.66), r_pos (basal relative stroma-to-lumen distance ratio cutoff: 0.25, 0.30, 0.35), d_px (maximum Euclidean distance from lamina propria stroma in 4x-downsampled pixels: 14, 18, 22), g_thr (maximum RGB grayscale tissue threshold: 228, 236, 244), and is_default (True for prespecified default 0.58 / 0.30 / 18 / 236).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S4_CellType_Tissue_Scan', column_or_title='LR_noTROP2 / q / q_noTROP2 / q_spearman',
         units='LR (df=1), BH-FDR q value',
         description='Tobit LR and Benjamini-Hochberg FDR q values across full 44-ADC cohort (q) and after excluding all 5 TACSTD2 ADCs (LR_noTROP2, q_noTROP2), plus BH-FDR q for prespecified unadjusted Spearman (q_spearman).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S6a_LOTO_LOO_Calib', column_or_title='pred_loto / pred_loto_payload_only / pred_loto_4payload / sigma_loto / log_score_loto',
         units='Predicted incidence (%), logit residual SD, log-likelihood score',
         description='Out-of-target LOTO conditional median predictions for basal+payload (pred_loto), 3-level payload-only (pred_loto_payload_only), and post hoc 4-level payload+basal (pred_loto_4payload); fitted training-fold Tobit logit error SD (sigma_loto, sigma_loto_payload); and out-of-fold Tobit log-predictive density/CDF score (log_score_loto, log_score_loto_payload; higher/less negative = better). Analogous *_loo columns report 44-fold leave-one-ADC-out estimates.'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S6c_Calib_Summary', column_or_title='c_index_all44_* / n_comparable_pairs / auc_ge15_* / prauc_ge15_* / auc_ge10_*',
         units='Concordance / AUC (0-1), pair count, two-sided P value',
         description='Summary discrimination across *_x (basal+payload), *_payload (payload-only), and *_4pay (post hoc 4-level payload+basal): distribution-free left-censored C-index across n_comparable_pairs=451 ordered pairs (5,000-perm p_perm_c_index), exact uncensored >=15% stomatitis ROC-AUC & PR-AUC across 38/44 ADCs (Mann-Whitney p_mw_auc_ge15, 5,000-perm p_perm_auc_ge15), and exact >=10% ROC-AUC across 35/44 ADCs.'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S6c_Calib_Summary', column_or_title='spearman_* / mae_* / cens_below_* / disc_* / p_mcnemar_twosided / mean_logscore_*',
         units='Spearman rho, MAE (percentage points), count (/29), P, logit log-score',
         description='Quantitative calibration on the 15 numeric ADCs (spearman_*, mae_*), count of 29 left-censored ADCs predicted below their trial reporting threshold (cens_below_*), discordant censored rescues/misses (disc_x_only, disc_payload_only; exact two-sided McNemar p_mcnemar_twosided), and mean out-of-fold Tobit log-scores across all 44, 15 numeric, and 29 censored ADCs (mean_logscore_all_*, mean_logscore_num_*, mean_logscore_cens_*).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S7_OSCC_420_Surfaceome', column_or_title='oral_basal / tum_choi / tum_puram / win_choi / win_puram',
         units='log2(CP10K + 1) and delta-log2',
         description='Shifted log2(mean CP10K + 1) expression in healthy oral basal keratinocytes (oral_basal) and OSCC malignant cells (tum_choi, tum_puram); win_choi and win_puram denote tumor minus oral basal delta-log2(CP10K + 1).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S7_OSCC_420_Surfaceome', column_or_title='vital_max_tpm / vital_max_tissue / pass_vital_15tpm',
         units='GTEx v8 median TPM; boolean',
         description='Maximum GTEx v8 median expression across the four-organ core panel (heart left ventricle, liver, kidney cortex, tibial nerve; labeled "four-organ core panel" in the manuscript) and indicator for <=15 TPM.'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S7_OSCC_420_Surfaceome', column_or_title='ext_max_tpm / nonbrain_ext_max_tpm / brain_max_tpm / gtex54_max_tpm',
         units='GTEx v8 median TPM',
         description='Maximum GTEx v8 median TPM across the extended 24-organ panel (ext_max_tpm: 4 core + 13 brain + 7 non-CNS visceral organs), its 11 non-CNS organs (nonbrain_ext_max_tpm), 13 brain regions (brain_max_tpm), and all 54 GTEx tissues (gtex54_max_tpm).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S7_OSCC_420_Surfaceome', column_or_title='pred_stom_pct / mean_ci_lo..hi / pred_pi_lo..hi / loto_pred_stom_pct',
         units='Predicted topoisomerase-I stomatitis (%)',
         description='Coherent full-cohort predicted stomatitis conditional median (pred_stom_pct) with B=2,000 target-cluster bootstrap 95% CI for the mean (mean_ci_lo/hi) and new-ADC 95% prediction interval (pred_pi_lo/hi); separate zero-antigen-overlap LOTO estimates for trained targets (loto_pred_stom_pct, loto_pi_lo/hi).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S7_OSCC_420_Surfaceome', column_or_title='n_pass / grid_pass / grid_pass_ext',
         units='Counts (out of 4 criteria; out of 9 grid cells)',
         description='Number of exploratory re-ranking criteria passed (n_pass, max 4) and 3x3 threshold stability counts (10%/12%/15% stomatitis x 10/15/30 TPM max, with criteria 1, 3, and 4 held fixed so that targets failing any of these score 0/9) under the four-organ core panel (grid_pass) and extended 24-organ panel (grid_pass_ext).'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S7_OSCC_420_Surfaceome', column_or_title='depmap_common_essential / depmap_median_chronos / depmap_frac_lt_neg05',
         units='Boolean, Chronos gene effect, fraction (0-1)',
         description='DepMap Public 24Q4 CRISPR dependency summary across 1,178 cancer cell lines: common-essential flag, median Chronos score, and fraction of cell lines with Chronos < -0.5.'),
    dict(section='COLUMN DICTIONARY', sheet_or_table='Table_S8a / S8b / S8c', column_or_title='beta / rate_ratio_per_log2 / nb_alpha / ols_beta / p_ols_perm / ROR_whole_FAERS / BCPNN_IC025',
         units='Log-rate slope, rate ratio, dispersion, ROR, PRR, BCPNN IC025',
         description='Exposure-offset NB2 count regression slope (beta, per log2 oral basal unit), rate ratio (exp(beta)), overdispersion (nb_alpha), equal-weight empirical log-odds OLS sensitivity (ols_beta, ols_LR, 10,000-resample 12-target permutation p_ols_perm), and whole-FAERS disproportionality (ROR, PRR, and Bayesian Confidence Propagation Neural Network information component lower 95% bound BCPNN_IC025 > 0 threshold).'),
]
readme_df = pd.DataFrame(readme_rows)

out_xlsx = f'{M}Supplementary_Tables.xlsx'
with pd.ExcelWriter(out_xlsx, engine='openpyxl') as writer:
    readme_df.to_excel(writer, sheet_name='README', index=False)
    s1_export.to_excel(writer, sheet_name='Table_S1_Clinical_66ADCs', index=False)
    ext_a_export.to_excel(writer, sheet_name='Table_S1a_Extractor_A_66ADCs', index=False)
    ext_b.to_excel(writer, sheet_name='Table_S1b_Extractor_B_66ADCs', index=False)
    ds3_sbs.to_excel(writer, sheet_name='Table_S1c_Three_Dataset_Results', index=False)
    reconciled_export.to_excel(writer, sheet_name='Table_S1d_Reconciled_66ADCs', index=False)
    dual_ag.to_excel(writer, sheet_name='Table_S1e_Dual_Agreement', index=False)
    feat_comp.to_excel(writer, sheet_name='Table_S2a_Features', index=False)
    res_ladder.to_excel(writer, sheet_name='Table_S2b_Resolution_Ladder', index=False)
    dab_models.to_excel(writer, sheet_name='Table_S2c_HPA_DAB_Models', index=False)
    dab_genes.to_excel(writer, sheet_name='Table_S2d_HPA_DAB_Genes', index=False)
    dab_grid.to_excel(writer, sheet_name='Table_S2e_HPA_DAB_ParamGrid', index=False)
    tgt_inf.to_excel(writer, sheet_name='Table_S3a_Target_Inference', index=False)
    sens.to_excel(writer, sheet_name='Table_S3b_Sensitivity_Subgroups', index=False)
    s3c_protocol.to_excel(writer, sheet_name='Table_S3c_Protocol_Alignment', index=False)
    scan.to_excel(writer, sheet_name='Table_S4_CellType_Tissue_Scan', index=False)
    axes_comp.to_excel(writer, sheet_name='Table_S5a_Axes_Single', index=False)
    joint_comp.to_excel(writer, sheet_name='Table_S5b_Joint_Nested', index=False)
    ratio_sens.to_excel(writer, sheet_name='Table_S5c_Ratio_Scales', index=False)
    s5_adcs.to_excel(writer, sheet_name='Table_S5d_Axes_Per_ADC', index=False)
    p5[s6_cols].to_excel(writer, sheet_name='Table_S6a_LOTO_LOO_Calib', index=False)
    p5_exact.to_excel(writer, sheet_name='Table_S6b_LOTO_Exact_Endpoints', index=False)
    p5_sum.to_excel(writer, sheet_name='Table_S6c_Calib_Summary', index=False)
    p4_export.to_excel(writer, sheet_name='Table_S7_OSCC_420_Surfaceome', index=False)
    faers_sum.to_excel(writer, sheet_name='Table_S8a_FAERS_Models', index=False)
    faers_disprop.to_excel(writer, sheet_name='Table_S8b_FAERS_Drugs_Disprop', index=False)
    faers_sens.to_excel(writer, sheet_name='Table_S8c_FAERS_Sensitivity', index=False)
    faers_readus_full.to_excel(writer, sheet_name='Table_S8d_FAERS_READUS_Flow', index=False)
    supp_refs_df.to_excel(writer, sheet_name='Supplementary_References', index=False)
    writer.book.properties.creator = ''
    writer.book.properties.title = 'Supplementary Tables S1-S8 and Supplementary References (ADC Stomatitis)'

out_refs_md = f'{M}Supplementary_References.md'
with open(out_refs_md, 'w', encoding='utf-8') as f:
    f.write('# Supplementary References\n\n')
    f.write('Numbered clinical safety sources corresponding to the Supplementary Reference column in Supplementary Table S1 (Clinical 66 ADCs), listing U.S. FDA prescribing information effective dates (with application numbers and set identifiers) and peer-reviewed clinical trial publications (PMIDs, PMCIDs, and DOIs).\n\n')
    for _, rr in supp_refs_df.iterrows():
        f.write(f"**{rr.supp_ref_id}.** {rr.citation}\n\n")

print('Wrote', out_xlsx, 'and', out_refs_md, f'({len(supp_refs_df)} numbered entries covering {len(s1)} ADCs)')
