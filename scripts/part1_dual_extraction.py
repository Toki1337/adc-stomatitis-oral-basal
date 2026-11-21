"""Blinded dual-extractor clinical curation, agreement statistics, adjudication, and three-dataset comparison."""

import os
from pathlib import Path
import numpy as np
import openpyxl
from openpyxl.styles import Alignment, Font, PatternFill
from openpyxl.utils import get_column_letter
import pandas as pd
from scipy import stats
from scipy.special import expit

from tobit_lib import PAYLOAD_MAP, design_matrix, exact_mcnemar_twosided, fit_tobit, lr_test

ROOT = str(Path(__file__).resolve().parents[1])
EXT_A_PATH = f'{ROOT}/clinical/adc_stomatitis_extractor_A.tsv'
EXT_B_PATH = f'{ROOT}/clinical/adc_stomatitis_extractor_B.tsv'
FINAL_PATH = f'{ROOT}/clinical/adc_stomatitis_trials.tsv'
REC_PATH = f'{ROOT}/clinical/adc_stomatitis_reconciled.tsv'
XLSX_PATH = f'{ROOT}/clinical/Dual_Blinded_Extractions_AB.xlsx'
WBASAL_PATH = f'{ROOT}/results/williams_kc_stats.tsv.gz'
OUT_DIR = f'{ROOT}/results/part1'
os.makedirs(OUT_DIR, exist_ok=True)


def cohen_kappa(y1, y2, labels=None, weight=None):
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    if labels is None:
        labels = sorted(set(y1) | set(y2))
    k = len(labels)
    idx = {v: i for i, v in enumerate(labels)}
    O = np.zeros((k, k), dtype=float)
    for a, b in zip(y1, y2):
        O[idx[a], idx[b]] += 1.0
    n = O.sum()
    if n == 0:
        return np.nan, np.nan
    O /= n
    E = O.sum(axis=1, keepdims=True) @ O.sum(axis=0, keepdims=True)
    if weight is None:
        W = 1.0 - np.eye(k)
    elif weight == 'linear':
        W = np.abs(np.subtract.outer(np.arange(k), np.arange(k))) / max(1, k - 1)
    elif weight == 'quadratic':
        W = (np.subtract.outer(np.arange(k), np.arange(k)) ** 2) / max(1, (k - 1) ** 2)
    else:
        raise ValueError(weight)
    denom = float((W * E).sum())
    kappa = 1.0 - float((W * O).sum()) / denom if denom > 0 else np.nan
    return float(kappa), float(np.trace(O))


def boot_kappa(y1, y2, labels=None, weight=None, n_boot=2000, seed=1):
    y1 = np.asarray(y1)
    y2 = np.asarray(y2)
    obs_k, obs_acc = cohen_kappa(y1, y2, labels=labels, weight=weight)
    n = len(y1)
    rng = np.random.default_rng(seed)
    ks = []
    for _ in range(n_boot):
        ix = rng.integers(0, n, size=n)
        k_b, _ = cohen_kappa(y1[ix], y2[ix], labels=labels, weight=weight)
        if np.isfinite(k_b):
            ks.append(k_b)
    lo, hi = np.quantile(ks, [0.025, 0.975]) if len(ks) >= 20 else (np.nan, np.nan)
    return obs_k, float(lo), float(hi), obs_acc, int(n)


def icc_2_1_and_2_k(y1, y2):
    Y = np.column_stack([np.asarray(y1, dtype=float), np.asarray(y2, dtype=float)])
    Y = Y[np.isfinite(Y).all(axis=1)]
    n, k = Y.shape
    if n < 3:
        return np.nan, np.nan, n
    gm = Y.mean()
    ss_tot = float(((Y - gm) ** 2).sum())
    ss_rows = float(k * ((Y.mean(axis=1) - gm) ** 2).sum())
    ss_cols = float(n * ((Y.mean(axis=0) - gm) ** 2).sum())
    ss_err = ss_tot - ss_rows - ss_cols
    ms_rows = ss_rows / (n - 1)
    ms_cols = ss_cols / (k - 1)
    ms_err = ss_err / ((n - 1) * (k - 1))
    d1 = ms_rows + (k - 1) * ms_err + k * (ms_cols - ms_err) / n
    dk = ms_rows + (ms_cols - ms_err) / n
    return (
        float((ms_rows - ms_err) / d1) if d1 > 0 else np.nan,
        float((ms_rows - ms_err) / dk) if dk > 0 else np.nan,
        int(n),
    )


def boot_icc(y1, y2, n_boot=2000, seed=1):
    y1 = np.asarray(y1, dtype=float)
    y2 = np.asarray(y2, dtype=float)
    m = np.isfinite(y1) & np.isfinite(y2)
    y1, y2 = y1[m], y2[m]
    icc1, icck, n = icc_2_1_and_2_k(y1, y2)
    rng = np.random.default_rng(seed)
    v1, vk = [], []
    for _ in range(n_boot):
        ix = rng.integers(0, n, size=n)
        i1, ik, _ = icc_2_1_and_2_k(y1[ix], y2[ix])
        if np.isfinite(i1):
            v1.append(i1)
            vk.append(ik)
    lo1, hi1 = np.quantile(v1, [0.025, 0.975]) if len(v1) >= 20 else (np.nan, np.nan)
    return icc1, float(lo1), float(hi1), icck, n


def prepare_model_df(df_in, W):
    """Build exact feature-linked clinical dataframe on the given extraction table."""
    alias = {'NECTIN4': 'PVRL4'}
    d = df_in[df_in.listed.astype(str).str.strip().str.lower().isin(['yes', 'not_listed'])].copy()
    d['stom_any'] = pd.to_numeric(d['stom_any'], errors='coerce')
    d['threshold_pct'] = pd.to_numeric(d['threshold_pct'], errors='coerce')
    d = d[~(d.stom_any.isna() & d.threshold_pct.isna())].copy()
    d['pclass'] = d.payload_class.map(PAYLOAD_MAP)
    d['cens'] = d.stom_any.isna().astype(bool)
    d['y_pct'] = np.where(d.cens, d.threshold_pct, d.stom_any)
    d['tg'] = d.target_gene.str.split(';').str[0]

    def get_wb(gs):
        vs = []
        for g in str(gs).split(';'):
            g2 = g if g in W.index else alias.get(g, g)
            if g2 in W.index:
                vs.append(float(np.log2(W.loc[g2, 'w_basal_mean'] + 1.0)))
        return max(vs) if vs else np.nan

    d['W_basal'] = [get_wb(gs) for gs in d.target_gene]
    return d.dropna(subset=['W_basal', 'pclass']).reset_index(drop=True)


def evaluate_dataset(dataset_label, stage_desc, df_in, W, n_perm=10000, n_boot=2000):
    """Evaluate primary Tobit, target inference, key subgroups, and LOTO calibration on a 66-ADC table."""
    T = prepare_model_df(df_in, W)
    obs = lr_test(T, 'W_basal', payload=True, compute_se=True)
    sp_rho, sp_p = stats.spearmanr(T.W_basal, np.where(T.cens, T.y_pct / 2.0, T.y_pct))

    # Target-level permutation (B=10,000, seed=1)
    y_v, c_v = T.y_pct.values, T.cens.values
    X_base = design_matrix(T)
    tg_vals = T.groupby('tg').W_basal.first()
    keys = tg_vals.index.values
    rng = np.random.default_rng(1)
    perm_lrs = np.zeros(n_perm)
    for i in range(n_perm):
        pmap = dict(zip(keys, rng.permutation(tg_vals.values)))
        Xp = np.column_stack([X_base, T.tg.map(pmap).values])
        rp = fit_tobit(y_v, c_v, Xp)
        perm_lrs[i] = max(0.0, 2.0 * (obs['nll0'] - rp['nll']))
    perm_p = float((1 + (perm_lrs >= obs['LR']).sum()) / (n_perm + 1))

    # Target-cluster bootstrap (B=2,000, seed=1)
    rng_b = np.random.default_rng(1)
    tgs = T.tg.unique()
    tg_idx = {t: np.where(T.tg.values == t)[0] for t in tgs}
    bs_betas = []
    for _ in range(n_boot):
        pick = rng_b.choice(tgs, len(tgs), replace=True)
        D = T.iloc[np.concatenate([tg_idx[t] for t in pick])].reset_index(drop=True)
        if D.cens.all() or D.pclass.nunique() < 3:
            continue
        rb = lr_test(D, 'W_basal', payload=True)
        if rb['converged'] and np.isfinite(rb['beta']) and 0.05 < rb['sigma'] < 10.0:
            bs_betas.append(rb['beta'])
    ci_lo, ci_hi = [float(v) for v in np.quantile(bs_betas, [0.025, 0.975])]

    # Subgroups & 4-level payload
    T_nt = T[T.tg != 'TACSTD2'].reset_index(drop=True)
    r_nt = lr_test(T_nt, 'W_basal', payload=True)
    T_top = T[T.pclass == 'topoI'].reset_index(drop=True)
    r_top = lr_test(T_top, 'W_basal', payload=False)
    T_tub = T[T.pclass == 'tubulin'].reset_index(drop=True)
    r_tub = lr_test(T_tub, 'W_basal', payload=False)
    r_p4 = lr_test(T, 'W_basal', payload=True, p4=True)

    # Zero-antigen-overlap LOTO cross-validation (excluding all ADCs sharing any antigen with test ADC)
    p_loto, p0_loto = np.zeros(len(T)), np.zeros(len(T))
    gene_alias = {'NECTIN4': 'PVRL4', 'CD44v6': 'CD44'}
    row_genes = [set(gene_alias.get(g, g) for g in str(gs).split(';')) for gs in T.target_gene]
    for i in range(len(T)):
        overlap = np.array([len(row_genes[i] & row_genes[j]) > 0 for j in range(len(T))])
        te, tr = T.iloc[[i]], T[~overlap]
        for xcol, arr in [('W_basal', p_loto), (None, p0_loto)]:
            f_tr = fit_tobit(tr.y_pct.values, tr.cens.values, design_matrix(tr, xcol, payload=True))
            cols_te = [
                np.ones(len(te)),
                (te.pclass == 'tubulin').astype(float).values,
                (te.pclass == 'DNA').astype(float).values,
            ]
            if xcol:
                cols_te.append(te[xcol].values)
            arr[i] = float((expit(np.column_stack(cols_te) @ f_tr['beta']) * 100.0)[0])

    num, cen = T[~T.cens], T[T.cens]
    loto_rho, loto_rho_p = stats.spearmanr(p_loto[num.index], num.y_pct)
    loto_mae = float(np.abs(p_loto[num.index] - num.y_pct).mean())
    loto_mae0 = float(np.abs(p0_loto[num.index] - num.y_pct).mean())
    ok_x = p_loto[cen.index] < cen.y_pct.values
    ok_0 = p0_loto[cen.index] < cen.y_pct.values
    mcn_p = exact_mcnemar_twosided(int((ok_x & (~ok_0)).sum()), int(((~ok_x) & ok_0).sum()))

    return dict(
        dataset=dataset_label,
        stage_blinding=stage_desc,
        n_total=len(df_in),
        n_excluded=int((df_in.listed.astype(str).str.strip().str.lower() == 'no_data').sum()),
        n_evaluable=len(T),
        n_numeric=len(num),
        n_censored=len(cen),
        n_targets=int(len(tgs)),
        beta=obs['beta'],
        se=obs['se'],
        boot_ci_lo=ci_lo,
        boot_ci_hi=ci_hi,
        LR=obs['LR'],
        asym_p=obs['p'],
        perm_p=perm_p,
        spearman_rho=float(sp_rho),
        spearman_p=float(sp_p),
        noTROP2_n=len(T_nt),
        noTROP2_LR=r_nt['LR'],
        noTROP2_p=r_nt['p'],
        topoI_n=len(T_top),
        topoI_LR=r_top['LR'],
        topoI_p=r_top['p'],
        tubulin_n=len(T_tub),
        tubulin_LR=r_tub['LR'],
        tubulin_p=r_tub['p'],
        p4_LR=r_p4['LR'],
        p4_p=r_p4['p'],
        loto_rho=float(loto_rho),
        loto_rho_p=float(loto_rho_p),
        loto_mae=loto_mae,
        loto_payload_mae=loto_mae0,
        loto_cens_x=int(ok_x.sum()),
        loto_cens_n=len(cen),
        loto_cens_payload_x=int(ok_0.sum()),
        loto_mcnemar_p=float(mcn_p),
    )


def build_side_by_side_table(ds_rows, name_A, name_B):
    """Format 3-column human-readable comparison table across Extractor A, Extractor B, and Final."""
    cols = [
        f'Original Extraction ({name_A})',
        f'Second Blinded Extraction ({name_B})',
        f'Adjudicated Final ({name_A} + {name_B})',
    ]

    def fmt_p(p):
        return f'{p:.4f}' if p >= 0.0005 else f'{p:.2e}'

    def col_vals(fn):
        return {c: fn(r) for c, r in zip(cols, ds_rows)}

    specs = [
        ('Curation timing and blinding', lambda r: r['stage_blinding']),
        ('Total ADCs curated / Excluded (no_data)', lambda r: f"{r['n_total']} / {r['n_excluded']}"),
        ('Evaluable ADCs (numeric / left-censored)', lambda r: f"{r['n_evaluable']} ({r['n_numeric']} / {r['n_censored']})"),
        ('Unique primary targets', lambda r: f"{r['n_targets']}"),
        ('Primary Tobit slope beta (SE)', lambda r: f"{r['beta']:.3f} ({r['se']:.3f})"),
        ('Target-cluster bootstrap 95% CI (B = 2,000)', lambda r: f"{r['boot_ci_lo']:.2f} to {r['boot_ci_hi']:.2f}"),
        ('Primary payload-adjusted Tobit LR (df = 1)', lambda r: f"{r['LR']:.2f}"),
        ('Likelihood-ratio asymptotic P', lambda r: fmt_p(r['asym_p'])),
        ('Target-level permutation P (B = 10,000)', lambda r: f"{r['perm_p']:.4f}"),
        ('Spearman rho (P) [threshold/2 imputation]', lambda r: f"{r['spearman_rho']:.3f} ({fmt_p(r['spearman_p'])})"),
        ('Excluding TACSTD2: n, LR (P)', lambda r: f"n = {r['noTROP2_n']}, {r['noTROP2_LR']:.2f} ({fmt_p(r['noTROP2_p'])})"),
        ('Topoisomerase-I subgroup: n, LR (P)', lambda r: f"n = {r['topoI_n']}, {r['topoI_LR']:.2f} ({fmt_p(r['topoI_p'])})"),
        ('Tubulin-inhibitor subgroup: n, LR (P)', lambda r: f"n = {r['tubulin_n']}, {r['tubulin_LR']:.2f} ({fmt_p(r['tubulin_p'])})"),
        ('Post hoc 4-level payload Tobit LR', lambda r: f"{r['p4_LR']:.2f}"),
        ('LOTO numeric Spearman rho (P)', lambda r: f"{r['loto_rho']:.3f} ({fmt_p(r['loto_rho_p'])})"),
        ('LOTO numeric MAE vs. payload-only (pp)', lambda r: f"{r['loto_mae']:.1f} vs. {r['loto_payload_mae']:.1f}"),
        ('LOTO censored below cutoff vs. payload-only (McNemar P)', lambda r: f"{r['loto_cens_x']}/{r['loto_cens_n']} vs. {r['loto_cens_payload_x']}/{r['loto_cens_n']} (P = {r['loto_mcnemar_p']:.2f})"),
    ]
    out_rows = []
    for label, fn in specs:
        row = {'Metric': label}
        row.update(col_vals(fn))
        out_rows.append(row)
    return pd.DataFrame(out_rows)


def write_formatted_workbook(sheets, out_path):
    """Write a clean multi-sheet Excel workbook with consistent header styling and auto-widths."""
    wb = openpyxl.Workbook()
    wb.properties.creator = ''
    wb.properties.description = ''
    wb.remove(wb.active)
    hdr_fill = PatternFill(start_color='EAF2F8', end_color='EAF2F8', fill_type='solid')
    hdr_font = Font(name='Calibri', size=10, bold=True, color='112244')
    cell_font = Font(name='Calibri', size=10)

    for sname, df in sheets:
        ws = wb.create_sheet(title=sname[:31])
        ws.freeze_panes = 'A2'
        df_clean = df.copy()
        for col_idx, col_name in enumerate(df_clean.columns, 1):
            c = ws.cell(row=1, column=col_idx, value=str(col_name))
            c.fill = hdr_fill
            c.font = hdr_font
            c.alignment = Alignment(horizontal='left', vertical='center')

        for row_idx, (_, rec) in enumerate(df_clean.iterrows(), 2):
            for col_idx, col_name in enumerate(df_clean.columns, 1):
                val = rec[col_name]
                if pd.isna(val):
                    val_out = ''
                elif isinstance(val, (np.floating, float)):
                    val_out = float(val)
                elif isinstance(val, (np.integer, int)):
                    val_out = int(val)
                else:
                    val_out = str(val)
                c = ws.cell(row=row_idx, column=col_idx, value=val_out)
                c.font = cell_font
                c.alignment = Alignment(vertical='top')

        for col_idx, col_name in enumerate(df_clean.columns, 1):
            lengths = [len(str(col_name))] + [len(str(v)) for v in df_clean[col_name].dropna().astype(str)]
            w = min(max(lengths) + 3, 55)
            ws.column_dimensions[get_column_letter(col_idx)].width = max(w, 11)

    wb.save(out_path)


def main():
    df_A = pd.read_csv(EXT_A_PATH, sep='\t')
    df_B = pd.read_csv(EXT_B_PATH, sep='\t')
    df_F = pd.read_csv(FINAL_PATH, sep='\t')
    W = pd.read_csv(WBASAL_PATH, sep='\t', index_col=0)

    assert len(df_A) == 66 and len(df_B) == 66 and len(df_F) == 66
    name_A = str(df_A['extractor'].iloc[0]).strip()
    ver_A = str(df_A['verified_by'].iloc[0]).strip()
    name_B = str(df_B['extractor'].iloc[0]).strip()
    ver_B = str(df_B['verified_by'].iloc[0]).strip()
    assert name_A and ver_A and name_B and ver_B

    # 1. Build standardized paired dataframe in alphabetical order for exact agreement CIs
    a_sub = df_A[['drug', 'tier', 'target_gene', 'payload_class', 'population_source_type', 'source_id', 'N', 'threshold_pct', 'listed', 'stom_any', 'stom_g3']].rename(columns={
        'population_source_type': 'A_tier', 'source_id': 'A_source', 'N': 'A_N',
        'threshold_pct': 'A_thresh', 'listed': 'A_listed', 'stom_any': 'A_stom_any', 'stom_g3': 'A_stom_g3',
    })
    a_sub['A_tier'] = a_sub['A_tier'].fillna('no_data').astype(str).str.strip().str.lower()

    b_sub = df_B[['drug', 'source_tier', 'source_id', 'N', 'threshold_pct', 'oral_term', 'is_grouped', 'listed', 'stom_any', 'stom_g3', 'trae_or_allcause']].rename(columns={
        'source_tier': 'B_tier', 'source_id': 'B_source', 'N': 'B_N',
        'threshold_pct': 'B_thresh', 'oral_term': 'B_oral_term', 'is_grouped': 'B_grouped',
        'listed': 'B_listed', 'stom_any': 'B_stom_any', 'stom_g3': 'B_stom_g3', 'trae_or_allcause': 'B_ae_type',
    })
    b_sub['B_tier'] = b_sub['B_tier'].fillna('no_data').astype(str).str.strip().str.lower()
    b_sub['B_listed'] = b_sub['B_listed'].fillna('no_data').astype(str).str.strip().str.lower()

    f_sub = df_F[['drug', 'population_source_type', 'source_id', 'N', 'threshold_pct', 'listed', 'stom_any', 'stom_g3']].rename(columns={
        'population_source_type': 'Final_tier', 'source_id': 'Final_source', 'N': 'Final_N',
        'threshold_pct': 'Final_thresh', 'listed': 'Final_listed', 'stom_any': 'Final_stom_any', 'stom_g3': 'Final_stom_g3',
    })
    f_sub['Final_tier'] = f_sub['Final_tier'].fillna('no_data').astype(str).str.strip().str.lower()

    m_alpha = a_sub.merge(b_sub, on='drug', how='inner').merge(f_sub, on='drug', how='inner').sort_values('drug').reset_index(drop=True)
    assert len(m_alpha) == 66

    # Inter-extractor agreement metrics (prior to adjudication)
    ag_rows = []
    k, lo, hi, acc, n = boot_kappa(m_alpha.A_listed, m_alpha.B_listed, labels=['yes', 'not_listed', 'no_data'])
    ag_rows.append(dict(metric='Status (yes / not_listed / no_data)', n=n, value=k, ci_lo=lo, ci_hi=hi, pct_agree=acc, note="Cohen's kappa (unweighted)"))
    ev_a = np.where(m_alpha.A_listed == 'no_data', 'no_data', 'evaluable')
    ev_b = np.where(m_alpha.B_listed == 'no_data', 'no_data', 'evaluable')
    k, lo, hi, acc, n = boot_kappa(ev_a, ev_b, labels=['evaluable', 'no_data'])
    ag_rows.append(dict(metric='Evaluability (evaluable vs no_data)', n=n, value=k, ci_lo=lo, ci_hi=hi, pct_agree=acc, note="Cohen's kappa"))
    k, lo, hi, acc, n = boot_kappa(m_alpha.A_tier, m_alpha.B_tier, labels=['label', 'phase3', 'phase2', 'phase1', 'no_data'])
    ag_rows.append(dict(metric='Source tier (all 66)', n=n, value=k, ci_lo=lo, ci_hi=hi, pct_agree=acc, note="Cohen's kappa (unweighted)"))
    mb = m_alpha[(m_alpha.A_listed != 'no_data') & (m_alpha.B_listed != 'no_data')]
    k, lo, hi, acc, n = boot_kappa(mb.A_tier, mb.B_tier, labels=['label', 'phase3', 'phase2', 'phase1'])
    ag_rows.append(dict(metric='Source tier (both evaluable)', n=n, value=k, ci_lo=lo, ci_hi=hi, pct_agree=acc, note="Cohen's kappa"))
    mnl = m_alpha[(m_alpha.A_listed == 'not_listed') & (m_alpha.B_listed == 'not_listed') & m_alpha.A_thresh.notna() & m_alpha.B_thresh.notna()]
    k_lw, lo_lw, hi_lw, acc_t, n_t = boot_kappa(mnl.A_thresh.astype(int), mnl.B_thresh.astype(int), labels=[5, 10, 15, 20, 25], weight='linear')
    k_qw, lo_qw, hi_qw, _, _ = boot_kappa(mnl.A_thresh.astype(int), mnl.B_thresh.astype(int), labels=[5, 10, 15, 20, 25], weight='quadratic')
    ag_rows.append(dict(metric='Censoring threshold_pct (linear weighted kappa)', n=n_t, value=k_lw, ci_lo=lo_lw, ci_hi=hi_lw, pct_agree=acc_t, note='Ordinal categories {5,10,15,20,25}'))
    ag_rows.append(dict(metric='Censoring threshold_pct (quadratic weighted kappa)', n=n_t, value=k_qw, ci_lo=lo_qw, ci_hi=hi_qw, pct_agree=acc_t, note='Ordinal categories {5,10,15,20,25}'))
    mnum = m_alpha[m_alpha.A_stom_any.notna() & m_alpha.B_stom_any.notna()]
    icc1, lo_i, hi_i, icck, n_num = boot_icc(mnum.A_stom_any, mnum.B_stom_any)
    rho_n, _ = stats.spearmanr(mnum.A_stom_any, mnum.B_stom_any)
    mae_n = float(np.abs(mnum.A_stom_any - mnum.B_stom_any).mean())
    maxd_n = float(np.abs(mnum.A_stom_any - mnum.B_stom_any).max())
    ag_rows.append(dict(metric='Numeric stom_any ICC(2,1)', n=n_num, value=icc1, ci_lo=lo_i, ci_hi=hi_i, pct_agree=np.nan, note=f'ICC(2,k)={icck:.4f}; Spearman={rho_n:.4f}; MAE={mae_n:.2f} pp; max_diff={maxd_n:.2f} pp'))

    ag_df = pd.DataFrame(ag_rows)
    ag_df['extractors'] = f'{name_A} vs. {name_B}'
    ag_df['verified_by'] = f'{ver_A}; {ver_B}'
    ag_df.to_csv(f'{OUT_DIR}/clinical_dual_agreement.tsv', sep='\t', index=False)

    # 2. Build 66-ADC Reconciled Side-by-Side Table in canonical Final order
    notes_map = {
        'rovalpituzumab tesirine': (
            'Second supersedes (tier rule)',
            f'{name_B} supersedes under higher-tier rule (Phase 2 > Phase 1): TRINITY Phase 2 monotherapy (PMID 31506387 / PMC7105795, N=339, TEAEs >=15%, no oral term) outranks {name_A} Phase 1 (PMID 27932068, N=74, <10%).'
        ),
        'trastuzumab rezetecan': (
            'Second supersedes (tier rule)',
            f'{name_B} supersedes under higher-tier rule (Phase 2 > Phase 1): HELEN-015 Phase 2 single-agent monotherapy (PMID 42811016 / PMC13623843, N=66, Table 4 G1-2 >=3%: Stomatitis 6/66 = 9.1%, G3-4 0%) outranks {name_A} Phase 1 HORIZON-X (PMID 41856971, N=396, <20%).'
        ),
        'trastuzumab emtansine': (
            'Original retained (monotherapy purity)',
            f'{name_A} retained under strict single-agent monotherapy rule (same FDA-label tier): EMILIA (N=490, 14%, metastatic breast cancer) is strict single-agent monotherapy, whereas KATHERINE (N=740, 15%, early breast cancer) permitted concurrent adjuvant radiotherapy and/or hormonal therapy per KADCYLA PI.'
        ),
        'trastuzumab deruxtecan': (
            'Original retained (monotherapy purity)',
            f'{name_A} retained under strict single-agent monotherapy rule (same FDA-label tier): DESTINY-Breast06 (N=434, 15%, metastatic breast cancer) is strict single-agent monotherapy, whereas DESTINY-Breast05 (N=806, 10%, early breast cancer) included concurrent adjuvant radiotherapy in 56% of patients per ENHERTU PI.'
        ),
        'disitamab vedotin': (
            'Original retained (tier rule)',
            f'{name_A} retained under higher-tier rule (Phase 3 > Phase 2): RC48-C006 Phase 3 (PMID 42802176 / PMC13616963, N=53, <20%) outranks {name_B} Phase 2 gastric cohort (PMID 34665942, N=125, <10%).'
        ),
        'depatuxizumab mafodotin': (
            'Original retained (tier rule)',
            f'{name_A} retained under higher-tier rule (Phase 2 > Phase 1): INTELLANCE 2 Phase 2 (PMID 31747009 / PMC7229258, N=84, <10%) outranks {name_B} Phase 1 (PMID 29075855, N=66, <25%).'
        ),
        'zilovertamab vedotin': (
            'Original retained (PMC table verified)',
            f'{name_A} retained after verifying Phase 2 PMC12442023.xml Table 4 (PMID 40762544): TRAEs >=10% across monotherapy cohorts totaling N=102 with no oral term.'
        ),
        'patritumab deruxtecan': (
            'Original retained (sample size N)',
            f'{name_A} retained under larger sample-size rule (both Phase 2 monotherapy): HERTHENA-Lung01 (PMID 37689979, N=225, Fig 5 bar chart 12%) > ICARUS-BREAST01 (PMID 40908353, N=99, 10.1%).'
        ),
        'vobramitamab duocarmazine': (
            'Original retained (RP2D cohort rule)',
            f'{name_A} retained under recommended Phase 2 dose (RP2D) cohort rule (both from Phase 1 PMID 42487497): RP2D 3.0 mg/kg cohort (N=121, 14.9%) vs. all-dose escalation+expansion pool (N=143, 13.3%).'
        ),
        'coltuximab ravtansine': (
            'Original retained (sample size N)',
            f'{name_A} retained under larger sample-size rule (both Phase 2 monotherapy): PMID 29748443 / PMC6068033 Table 3 image (N=61 DLBCL, <10%) > PMID 26775883 (N=36 ALL, <10%).'
        ),
        'aprutumab ixadotin': (
            'Concordant (PMC table cutoff verified)',
            f'Both {name_A} and {name_B} extracted left-censored <10%; source verification against PMC6797631 Table 2 header confirmed explicit ">=10% of patients" TEAE cutoff (N=20).'
        ),
        'polatuzumab vedotin': (
            'Concordant no_data (tier normalized)',
            'Agreed no_data (no monotherapy N>=15 any-grade AE table); source tier harmonized to no_data.'
        ),
        'praluzatamab ravtansine': (
            'Concordant no_data (tier normalized)',
            'Agreed no_data (Table 2 Grade>=3 only); source tier harmonized to no_data.'
        ),
        'vadastuximab talirine': (
            'Concordant no_data (tier normalized)',
            'Agreed no_data (abstract only without any-grade AE table); source tier harmonized to no_data.'
        ),
    }

    order_canon = df_F['drug'].tolist()
    m_canon = a_sub.merge(b_sub, on='drug', how='inner').merge(f_sub, on='drug', how='inner').set_index('drug').loc[order_canon].reset_index()

    rec_rows = []
    for _, r in m_canon.iterrows():
        status, note = notes_map.get(r['drug'], ('Concordant', f'Concordant between {name_A} and {name_B}.'))
        rec_rows.append(dict(
            drug=r['drug'],
            tier=r['tier'],
            target_gene=r['target_gene'],
            payload_class=r['payload_class'],
            Extractor_A_name=name_A,
            Extractor_A_tier=r['A_tier'],
            Extractor_A_N=r['A_N'],
            Extractor_A_listed=r['A_listed'],
            Extractor_A_stom_any=r['A_stom_any'],
            Extractor_A_stom_g3=r['A_stom_g3'],
            Extractor_A_thresh=r['A_thresh'],
            Extractor_A_source=r['A_source'],
            Extractor_B_name=name_B,
            Extractor_B_tier=r['B_tier'],
            Extractor_B_N=r['B_N'],
            Extractor_B_listed=r['B_listed'],
            Extractor_B_stom_any=r['B_stom_any'],
            Extractor_B_stom_g3=r['B_stom_g3'],
            Extractor_B_thresh=r['B_thresh'],
            Extractor_B_oral_term=r['B_oral_term'],
            Extractor_B_ae_type=r['B_ae_type'],
            Extractor_B_source=r['B_source'],
            Final_tier=r['Final_tier'],
            Final_N=r['Final_N'],
            Final_listed=r['Final_listed'],
            Final_stom_any=r['Final_stom_any'],
            Final_stom_g3=r['Final_stom_g3'],
            Final_thresh=r['Final_thresh'],
            Final_source=r['Final_source'],
            Comparison_Status=status,
            Adjudicated_By=f'{ver_A}; {ver_B}',
            Adjudication_Note=note,
        ))
    rec_df = pd.DataFrame(rec_rows)
    rec_df.to_csv(REC_PATH, sep='\t', index=False)

    # 3. Three-Dataset Model & LOTO Comparison in canonical 66-ADC order
    df_B_for_model = df_B.merge(df_F[['drug', 'target_gene', 'payload_class']], on='drug', how='left')
    ds_rows = [
        evaluate_dataset(
            f'Original ({name_A})',
            'Pre-linkage extraction; verified before expression linkage',
            df_A,
            W,
        ),
        evaluate_dataset(
            f'Second ({name_B})',
            f'Post-linkage independent extraction; blinded to expression & {name_A}',
            df_B_for_model,
            W,
        ),
        evaluate_dataset(
            f'Adjudicated Final ({name_A} + {name_B})',
            'Source-verified reconciliation under unified hierarchy',
            df_F,
            W,
        ),
    ]
    ds_df = pd.DataFrame(ds_rows)
    ds_df.to_csv(f'{OUT_DIR}/clinical_three_dataset_comparison.tsv', sep='\t', index=False)

    sbs_df = build_side_by_side_table(ds_rows, name_A, name_B)
    sbs_df.to_csv(f'{OUT_DIR}/clinical_three_dataset_side_by_side.tsv', sep='\t', index=False)

    # 4. Standalone multi-sheet Excel workbook for clinical author review
    write_formatted_workbook(
        [
            ('Extractor_A_Original_66ADCs', df_A),
            ('Extractor_B_Blinded_66ADCs', df_B),
            ('Three_Datasets_SideBySide', sbs_df),
            ('Reconciled_66ADCs_SideBySide', rec_df),
            ('Inter_Extractor_Agreement', ag_df),
        ],
        XLSX_PATH,
    )

    print(f'Wrote {REC_PATH} ({len(rec_df)} ADCs)')
    print(f'Wrote {OUT_DIR}/clinical_dual_agreement.tsv ({len(ag_df)} metrics)')
    print(f'Wrote {OUT_DIR}/clinical_three_dataset_comparison.tsv ({len(ds_df)} datasets)')
    print(f'Wrote {OUT_DIR}/clinical_three_dataset_side_by_side.tsv ({len(sbs_df)} rows)')
    print(f'Wrote {XLSX_PATH} (5 sheets)')
    print('\n=== THREE-DATASET SIDE-BY-SIDE SUMMARY ===')
    print(sbs_df.to_string(index=False))


if __name__ == '__main__':
    main()
