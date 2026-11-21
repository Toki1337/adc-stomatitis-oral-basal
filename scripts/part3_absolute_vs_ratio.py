"""Part 3: Head-to-head comparison, nested joint models, and ratio-scale sensitivity."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

from tobit_lib import PAYLOAD_MAP, fit_tobit, design_matrix, lr_test

pd.set_option('display.width', 260)
ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'
T = pd.read_csv(R + 'part1/adc_table_with_features.tsv', sep='\t')
T['pclass'] = T.payload_class.map(PAYLOAD_MAP)
T['cens'] = T.stom_any.isna().astype(bool)
T['y_pct'] = np.where(T.cens, T.threshold_pct, T.stom_any)

W = pd.read_csv(R + 'williams_kc_stats.tsv.gz', sep='\t', index_col=0)
C = pd.read_csv(R + 'choi_group_stats.tsv.gz', sep='\t', index_col=0)
P = pd.read_csv(R + 'puram_group_stats.tsv.gz', sep='\t', index_col=0)
Tc = pd.read_csv(R + 'tcga_oral_stats.tsv.gz', sep='\t', index_col=0)

alias = {'NECTIN4': 'PVRL4'}


def resolve(df, g):
    if g in df.index:
        return g
    g2 = alias.get(g, g)
    return g2 if g2 in df.index else None


def feat_log(df, col, gs, pc=1.0):
    vs = []
    for g in str(gs).split(';'):
        g2 = resolve(df, g)
        if g2 is not None:
            vs.append(np.log2(float(df.loc[g2, col]) + pc))
    return max(vs) if vs else np.nan


def feat_ratio(df_t, col_t, df_n, col_n, gs, pc=1.0, linear=False):
    vs = []
    for g in str(gs).split(';'):
        gt = resolve(df_t, g)
        gn = resolve(df_n, g)
        if gt is not None and gn is not None:
            vt = float(df_t.loc[gt, col_t])
            vn = float(df_n.loc[gn, col_n])
            if linear:
                vs.append((vt + pc) / (vn + pc))
            else:
                vs.append(np.log2(vt + pc) - np.log2(vn + pc))
    return max(vs) if vs else np.nan


def feat_raw(df, col, gs):
    vs = []
    for g in str(gs).split(';'):
        g2 = resolve(df, g)
        if g2 is not None and not np.isnan(df.loc[g2, col]):
            vs.append(float(df.loc[g2, col]))
    return max(vs) if vs else np.nan


T['tg'] = [resolve(W, str(g).split(';')[0]) or str(g).split(';')[0] for g in T.target_gene]
T['abs_oral_basal'] = [feat_log(W, 'w_basal_mean', g) for g in T.target_gene]
T['w_basal_raw'] = [feat_raw(W, 'w_basal_mean', g) for g in T.target_gene]
T['w_basal_det'] = [feat_raw(W, 'w_basal_det', g) for g in T.target_gene]
T['tum_oscc'] = [feat_log(C, 'choi_malig_OC_HPVneg_mean', g) for g in T.target_gene]
T['ratio_oscc'] = [feat_ratio(C, 'choi_malig_OC_HPVneg_mean', W, 'w_basal_mean', g, pc=1.0) for g in T.target_gene]
T['tum_puram'] = [feat_log(P, 'puram_malig_primary_mean', g) for g in T.target_gene]
T['ratio_puram'] = [feat_ratio(P, 'puram_malig_primary_mean', W, 'w_basal_mean', g, pc=1.0) for g in T.target_gene]
T['ratio_choi_internal'] = [feat_ratio(C, 'choi_malig_OC_HPVneg_mean', C, 'choi_epi_NL_OC_mean', g, pc=1.0) for g in T.target_gene]
T['ratio_choi_pc01'] = [feat_ratio(C, 'choi_malig_OC_HPVneg_mean', W, 'w_basal_mean', g, pc=0.1) for g in T.target_gene]
T['ratio_choi_pc001'] = [feat_ratio(C, 'choi_malig_OC_HPVneg_mean', W, 'w_basal_mean', g, pc=0.01) for g in T.target_gene]
T['ratio_choi_linear_pc01'] = [feat_ratio(C, 'choi_malig_OC_HPVneg_mean', W, 'w_basal_mean', g, pc=0.1, linear=True) for g in T.target_gene]
T['tcga_paired_lfc'] = [feat_raw(Tc, 'tcga_paired_lfc', g) for g in T.target_gene]

# 1. Single-predictor axes comparison
print('=== Single-predictor axes across n=44 ADCs (Tobit, payload-adjusted) ===')
axes = [
    ('abs_oral_basal', 'ABSOLUTE oral basal keratinocyte expression'),
    ('ratio_oscc', 'tumour/normal RATIO (Choi OSCC malignant / oral basal)'),
    ('ratio_puram', 'tumour/normal RATIO (Puram OSCC malignant / oral basal)'),
    ('ratio_choi_internal', 'tumour/normal RATIO (Choi internal OC malignant / OC NL epi)'),
    ('ratio_choi_pc01', 'tumour/normal RATIO (Choi / oral basal, pc=0.1)'),
    ('ratio_choi_pc001', 'tumour/normal RATIO (Choi / oral basal, pc=0.01)'),
    ('ratio_choi_linear_pc01', 'tumour/normal LINEAR ratio ((Choi+0.1)/(basal+0.1))'),
    ('tcga_paired_lfc', 'tumour/normal paired LFC (TCGA oral, n=30 pairs)'),
    ('tum_oscc', 'OSCC tumour expression alone (Choi)'),
    ('tum_puram', 'OSCC tumour expression alone (Puram)'),
]
rows = []
for x, lab in axes:
    d = T.dropna(subset=[x])
    r = lr_test(d, x, compute_se=True)
    imp = np.where(d.cens, d.y_pct / 2.0, d.y_pct)
    rho, prho = stats.spearmanr(d[x], imp)
    rows.append(dict(
        col=x, axis=lab, n=len(d), beta=r['beta'], beta_per_SD=r['beta'] * d[x].std(),
        se=r['se'], LR=r['LR'], p=r['p'], spearman=rho, p_spearman=prho,
    ))
comp_df = pd.DataFrame(rows)
print(comp_df.round(4).to_string(index=False))
comp_df.to_csv(R + 'part3_axes_comparison.tsv', sep='\t', index=False)

# 2. Nested joint models: logit(Y) ~ abs + ratio + payload
print('\n=== Nested joint models: logit(Y) ~ abs + ratio + payload ===')
d44 = T.dropna(subset=['abs_oral_basal']).copy()
j_rows = []
for rcol, rlab in [
    ('ratio_oscc', 'Choi OSCC / oral basal (log2+1)'),
    ('ratio_puram', 'Puram OSCC / oral basal (log2+1)'),
    ('tcga_paired_lfc', 'TCGA-HNSC oral paired LFC (n=30 pairs)'),
]:
    dj = d44.dropna(subset=[rcol])
    r_add_rat = lr_test(dj, rcol, base_cols=['abs_oral_basal'])
    r_add_abs = lr_test(dj, 'abs_oral_basal', base_cols=[rcol])
    r_joint = fit_tobit(dj.y_pct.values, dj.cens.values, design_matrix(dj, ['abs_oral_basal', rcol]))
    j_rows.append(dict(
        ratio_col=rcol, ratio_source=rlab, n=len(dj),
        joint_beta_abs=float(r_joint['beta'][-2]),
        joint_beta_ratio=float(r_joint['beta'][-1]),
        added_LR_ratio_over_abs=r_add_rat['LR'],
        p_added_ratio=r_add_rat['p'],
        added_LR_abs_over_ratio=r_add_abs['LR'],
        p_added_abs=r_add_abs['p'],
    ))
joint_df = pd.DataFrame(j_rows)
print(joint_df.round(4).to_string(index=False))
joint_df.to_csv(R + 'part3_joint_models.tsv', sep='\t', index=False)

# 3. Pseudocount and normal-threshold sensitivity
print('\n=== Ratio sensitivity across 5 scales and normal-expressed subsets ===')
ratio_defs = [
    ('abs_oral_basal', 'Absolute oral basal (reference)', 'log2(CP10K+1)'),
    ('ratio_oscc', 'Log2 ratio (pc=1.0)', 'log2((T+1)/(N+1))'),
    ('ratio_choi_pc01', 'Log2 ratio (pc=0.1)', 'log2((T+0.1)/(N+0.1))'),
    ('ratio_choi_pc001', 'Log2 ratio (pc=0.01)', 'log2((T+0.01)/(N+0.01))'),
    ('ratio_choi_linear_pc01', 'Linear ratio (eps=0.1)', '(T+0.1)/(N+0.1)'),
    ('tcga_paired_lfc', 'TCGA oral paired LFC', 'log2 paired RNA'),
]
rs_rows = []
for sub_name, df_sub in [
    ('All evaluable ADCs', d44),
    ('Normal basal det >= 0.05', d44[d44.w_basal_det >= 0.05]),
    ('Normal basal mean >= 0.05 CP10K', d44[d44.w_basal_raw >= 0.05]),
]:
    for col, lab, scale in ratio_defs:
        ds = df_sub.dropna(subset=[col])
        r = lr_test(ds, col)
        imp = np.where(ds.cens, ds.y_pct / 2.0, ds.y_pct)
        rho, prho = stats.spearmanr(ds[col], imp)
        rs_rows.append(dict(
            subset=sub_name, n=len(ds), n_targets=int(ds.tg.nunique()),
            col=col, metric=lab, scale=scale,
            beta=r['beta'], LR=r['LR'], p=r['p'], spearman=rho, p_spearman=prho,
        ))
rs_df = pd.DataFrame(rs_rows)
print(rs_df.round(4).to_string(index=False))
rs_df.to_csv(R + 'part3_ratio_sensitivity.tsv', sep='\t', index=False)

# Rank agreement across 28 targets
show = d44[['drug', 'target_gene', 'tg', 'pclass', 'y_pct', 'cens', 'abs_oral_basal', 'tum_oscc', 'ratio_oscc']].copy()
u = show.drop_duplicates('tg').set_index('tg')
rho_ar, p_ar = stats.spearmanr(u.abs_oral_basal, u.ratio_oscc)
print('\nTACSTD2 rank among %d primary targets: abs=%d/%d | ratio=%d/%d | Spearman(abs, ratio)=%.3f (p=%.4f)' % (
    len(u), int(u.abs_oral_basal.rank(ascending=False)['TACSTD2']), len(u),
    int(u.ratio_oscc.rank(ascending=False)['TACSTD2']), len(u),
    rho_ar, p_ar,
))
T.to_csv(R + 'part3_axes.tsv', sep='\t', index=False)
