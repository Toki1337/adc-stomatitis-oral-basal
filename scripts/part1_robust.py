"""Part 1b: Target-level inference (B=10,000 permutation + cluster bootstrap), outcome heterogeneity & profile-CI sensitivity."""

from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import optimize, stats
from scipy.special import logit

from tobit_lib import (
    PAYLOAD_MAP,
    design_matrix,
    fit_tobit,
    lr_test,
    profile_ci,
)

ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'
T = pd.read_csv(R + 'part1/adc_table_with_features.tsv', sep='\t')
T['pclass'] = T.payload_class.map(PAYLOAD_MAP)
T['cens'] = T.stom_any.isna().astype(bool)
T['y_pct'] = np.where(T.cens, T.threshold_pct, T.stom_any)
T['tg'] = T.target_gene.str.split(';').str[0]
alias = {'NECTIN4': 'PVRL4'}

# 1. Leave-one-target-out LR across core features
for k in ['W_basal', 'W_basal_det', 'HPA_sc_esoph_basal', 'GTEx_esoph_mucosa', 'W_alltissue', 'TCGA_normal_oral']:
    d = T.dropna(subset=[k])
    out = []
    for tg in d.tg.unique():
        dd = d[d.tg != tg]
        r = lr_test(dd, k)
        out.append((tg, int((d.tg == tg).sum()), round(r['LR'], 2), r['p'], round(r['beta'], 3)))
    o = pd.DataFrame(out, columns=['drop', 'n_dropped', 'LR', 'p', 'beta']).sort_values('p', ascending=False)
    print('==', k, 'worst LOTO:', o.head(3).to_dict('records'))

# 2. Exclude TROP2 entirely
for k in ['W_basal', 'W_basal_det', 'HPA_sc_esoph_basal', 'GTEx_esoph_mucosa', 'W_alltissue']:
    d = T.dropna(subset=[k])
    d = d[d.tg != 'TACSTD2']
    r = lr_test(d, k)
    sp_res = stats.spearmanr(d[k], np.where(d.cens, d.y_pct / 2.0, d.y_pct))
    print('noTROP2', k, len(d), int((~d.cens).sum()), [round(r['LR'], 4), round(r['p'], 4), round(r['beta'], 4)],
          'spearman', round(float(sp_res.statistic), 4), round(float(sp_res.pvalue), 4))

# 3. Primary target-level inference (28 targets): B=10,000 target permutation + B=2,000 cluster bootstrap
obs_res = lr_test(T, 'W_basal', compute_se=True)
obs_lr = obs_res['LR']
obs_nll0 = obs_res['nll0']
y_v = T.y_pct.values
c_v = T.cens.values
X_base = design_matrix(T)

rng = np.random.default_rng(1)
tg_series = T.groupby('tg').W_basal.first()
keys = tg_series.index.values
perm_lrs = np.zeros(10000)
for i in range(10000):
    pmap = dict(zip(keys, rng.permutation(tg_series.values)))
    Xp = np.column_stack([X_base, T.tg.map(pmap).values])
    r_p = fit_tobit(y_v, c_v, Xp)
    perm_lrs[i] = max(0.0, 2.0 * (obs_nll0 - r_p['nll']))
perm_p = float((1 + (perm_lrs >= obs_lr).sum()) / (len(perm_lrs) + 1))

rng_bs = np.random.default_rng(1)
tgs = T.tg.unique()
tg_idx = {t: np.where(T.tg.values == t)[0] for t in tgs}
bs_betas = []
for _ in range(2000):
    pick = rng_bs.choice(tgs, len(tgs), replace=True)
    idx = np.concatenate([tg_idx[t] for t in pick])
    D = T.iloc[idx].reset_index(drop=True)
    if D.cens.all() or D.pclass.nunique() < 3:
        continue
    rb = lr_test(D, 'W_basal')
    if rb['converged'] and np.isfinite(rb['beta']) and 0.05 < rb['sigma'] < 10.0:
        bs_betas.append(rb['beta'])
bs_betas = np.array(bs_betas)
bs_med = float(np.median(bs_betas))
bs_lo, bs_hi = [float(v) for v in np.quantile(bs_betas, [0.025, 0.975])]
bs_frac_le0 = float((bs_betas <= 0).mean())

# 3b. Post hoc 4-level payload model: basal LR, 4-vs-3 overall payload LR (df=1), and target-permutation P
obs4_res = lr_test(T, 'W_basal', payload=True, p4=True, compute_se=True)
lr_payload_4v3 = max(0.0, 2.0 * (obs_res['nll1'] - obs4_res['nll1']))
p_payload_4v3 = float(stats.chi2.sf(lr_payload_4v3, 1))
X_base4 = design_matrix(T, payload=True, p4=True)
rng_p4 = np.random.default_rng(1)
perm4_lrs = np.zeros(10000)
for i in range(10000):
    pmap = dict(zip(keys, rng_p4.permutation(tg_series.values)))
    Xp = np.column_stack([X_base4, T.tg.map(pmap).values])
    r_p4 = fit_tobit(y_v, c_v, Xp)
    perm4_lrs[i] = max(0.0, 2.0 * (obs4_res['nll0'] - r_p4['nll']))
perm4_p = float((1 + (perm4_lrs >= obs4_res['LR']).sum()) / (len(perm4_lrs) + 1))

target_inf = pd.DataFrame([dict(
    n_adcs=len(T), n_targets=len(tgs),
    obs_LR=obs_lr, obs_beta=obs_res['beta'],
    perm_B=10000, perm_p=perm_p,
    perm_null_median_LR=float(np.median(perm_lrs)),
    perm_null_q95_LR=float(np.quantile(perm_lrs, 0.95)),
    boot_B=2000, boot_valid=len(bs_betas),
    boot_beta_median=bs_med, boot_beta_ci_lo=bs_lo, boot_beta_ci_hi=bs_hi,
    boot_frac_le0=bs_frac_le0,
    p4_basal_LR=obs4_res['LR'], p4_basal_beta=obs4_res['beta'], p4_basal_se=obs4_res['se'],
    p4_basal_perm_p=perm4_p, p4_vs_p3_LR=lr_payload_4v3, p4_vs_p3_p=p_payload_4v3,
)])
print('\n=== TARGET-LEVEL INFERENCE (28 targets) ===')
print(target_inf.round(5).to_string(index=False))
target_inf.to_csv(R + 'part1/target_inference.tsv', sep='\t', index=False)

# 4. Sensitivity, subgroup, outcome heterogeneity, and basal definition table with Profile Likelihood 95% CIs
heme_drugs = {
    'inotuzumab ozogamicin', 'gemtuzumab ozogamicin', 'loncastuximab tesirine',
    'brentuximab vedotin', 'belantamab mafodotin', 'camidanlumab tesirine',
    'coltuximab ravtansine', 'naratuximab emtansine', 'pivekimab sunirine',
}

# Exact TRAE vs all-cause and single stomatitis PT classification from source evidence_quote
trae_kws = ['trae', 'treatment-related', 'related treatment-emergent', 'common related adverse']
T['is_trae'] = T.evidence_quote.fillna('').str.lower().apply(lambda s: any(kw in s for kw in trae_kws))


def _is_single_term(row):
    if row.cens:
        return True
    eq_l = str(row.evidence_quote).lower()
    is_grp = (
        ('including' in eq_l) or ('includes' in eq_l) or ('grouped' in eq_l)
        or ('composite' in eq_l) or ('oral mucositis' in eq_l)
        or ('mucositis (' in eq_l) or ('oral ulcer' in eq_l)
    )
    return not is_grp


T['is_single_term'] = T.apply(_is_single_term, axis=1)

sens_rows = []


def add_sens(cat, label, d, k='W_basal', payload=True, p4=False, imp_mode='half', weight_mode=None, p_target_perm=np.nan):
    d = d.dropna(subset=[k]).copy()
    if imp_mode == 'half':
        imp = np.where(d.cens, d.y_pct / 2.0, d.y_pct)
    elif imp_mode == 'thresh':
        imp = d.y_pct
    else:
        imp = np.where(d.cens, 0.0, d.y_pct)
    rho, prho = stats.spearmanr(d[k], imp)
    if cat == 'Imputation':
        sens_rows.append(dict(
            category=cat, specification=label, n=len(d), n_num=int((~d.cens).sum()),
            beta=np.nan, se=np.nan, prof_ci_lo=np.nan, prof_ci_hi=np.nan,
            ci_type='Spearman_imputation_only (Tobit invariant)',
            LR=np.nan, p_LR=np.nan, p_target_perm=np.nan,
            spearman=rho, p_spearman=prho, converged=True,
        ))
        return
    w = None
    if weight_mode == 'N':
        w = d.N.values / d.N.values.mean()
    elif weight_mode == 'sqrtN':
        w = np.sqrt(d.N.values) / np.sqrt(d.N.values).mean()
    r = lr_test(d, k, payload=payload, p4=p4, weights=w, compute_se=True)
    if weight_mode is not None:
        # Weighted pseudo-likelihood does not follow chi-squared(1); report target-cluster bootstrap 95% CI
        rng_w = np.random.default_rng(1)
        tgs_w = d.tg.unique()
        tg_idx_w = {t: np.where(d.tg.values == t)[0] for t in tgs_w}
        bs_w = []
        for _ in range(2000):
            pick = rng_w.choice(tgs_w, len(tgs_w), replace=True)
            idx = np.concatenate([tg_idx_w[t] for t in pick])
            Dw = d.iloc[idx].reset_index(drop=True)
            if Dw.cens.all() or Dw.pclass.nunique() < 3:
                continue
            ww = (Dw.N.values / Dw.N.values.mean()) if weight_mode == 'N' else (np.sqrt(Dw.N.values) / np.sqrt(Dw.N.values).mean())
            rb = lr_test(Dw, k, payload=payload, p4=p4, weights=ww)
            if rb['converged'] and np.isfinite(rb['beta']) and 0.05 < rb['sigma'] < 10.0:
                bs_w.append(rb['beta'])
        plo, phi = [float(v) for v in np.quantile(bs_w, [0.025, 0.975])]
        ci_type = 'cluster_bootstrap_95CI'
        lr_val, p_lr_val = np.nan, np.nan
    else:
        _, plo, phi = profile_ci(d, k, payload=payload, p4=p4)
        ci_type = 'profile_likelihood_95CI'
        # Do not expose non-target P < 1e-5 for 4-level payload specification; pair with target-permutation P instead
        lr_val = r['LR']
        p_lr_val = np.nan if p4 else r['p']
    sens_rows.append(dict(
        category=cat, specification=label, n=len(d), n_num=int((~d.cens).sum()),
        beta=r['beta'], se=r['se'],
        prof_ci_lo=plo, prof_ci_hi=phi, ci_type=ci_type,
        LR=lr_val, p_LR=p_lr_val, p_target_perm=p_target_perm,
        spearman=rho, p_spearman=prho, converged=r['converged'],
    ))


add_sens('Primary', 'Full cohort (3-level payload)', T, p_target_perm=perm_p)
add_sens('Payload_model_posthoc', '4-level payload post hoc (split SN-38 vs DXd/other)', T, p4=True, p_target_perm=perm4_p)
add_sens('Target_subset', 'Exclude TACSTD2 (TROP2)', T[T.tg != 'TACSTD2'])
add_sens('Target_subset', 'Exclude TACSTD2 and EGFR;ERBB3', T[~T.target_gene.isin(['TACSTD2', 'EGFR;ERBB3'])])
add_sens('Target_subset', 'Exclude all TACSTD2 and EGFR', T[~T.target_gene.isin(['TACSTD2', 'EGFR', 'EGFR;ERBB3'])])
add_sens('Tumour_subset', 'Solid tumours only', T[~T.drug.isin(heme_drugs)])
add_sens('Tumour_subset', 'Solid tumours excl TACSTD2', T[(~T.drug.isin(heme_drugs)) & (T.tg != 'TACSTD2')])
add_sens('Tumour_subset', 'Phase 2/3/FDA label only (excl Phase 1)', T[T.population_source_type != 'phase1'])
add_sens('Within_payload', 'Topo-I only (all 12)', T[T.pclass == 'topoI'], payload=False)
add_sens('Within_payload', 'Topo-I excl SN-38, post hoc (n=10)', T[T.payload_class.isin(['DXd', 'topo1-other'])], payload=False)
add_sens('Within_payload', 'Tubulin only (n=25)', T[T.pclass == 'tubulin'], payload=False)
add_sens('Within_payload', 'MMAE only (n=11)', T[T.payload_class == 'MMAE'], payload=False)

# Outcome heterogeneity sensitivity (7 checks total)
add_sens('Outcome_heterogeneity', 'All-cause / label AR only (excl TRAE-only)', T[~T.is_trae])
add_sens('Outcome_heterogeneity', 'TRAE-only studies', T[T.is_trae])
add_sens('Outcome_heterogeneity', 'Single stomatitis PT numeric + censored (n=35; 6 num)', T[T.cens | T.is_single_term])
add_sens('Outcome_heterogeneity', 'Exclude small trials (N < 30; excl SHR-A1921, aprutumab)', T[T.N >= 30])
add_sens('Outcome_heterogeneity', 'Sample-size weighted Tobit (weight = N / mean N)', T, weight_mode='N')
add_sens('Outcome_heterogeneity', 'Sample-size weighted Tobit (weight = sqrt N)', T, weight_mode='sqrtN')
add_sens('Outcome_heterogeneity', 'FDA label source only (n=13)', T[T.population_source_type == 'label'])

# Data source & composite-row sensitivity checks
T_sig11 = T.copy()
T_sig11.loc[T_sig11.drug == 'sigvotatug vedotin', 'y_pct'] = 11.0
add_sens('Data_source_sensitivity', 'Sigvotatug vedotin single stomatitis PT (11% vs 22% composite)', T_sig11)
add_sens('Data_source_sensitivity', 'Exclude abstract-only censored rows (sofituzumab, vandortuzumab)', T[~T.drug.isin(['sofituzumab vedotin', 'vandortuzumab vedotin'])])

add_sens('Imputation', 'Censored set to threshold/2', T, imp_mode='half')
add_sens('Imputation', 'Censored set to threshold', T, imp_mode='thresh')
add_sens('Imputation', 'Censored set to 0', T, imp_mode='zero')


# Basal compartment definition sensitivity from williams_epi.h5ad
def read_h5_col(grp, key):
    item = grp[key]
    if isinstance(item, h5py.Group):
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in item['categories'][:]]
        return pd.Categorical.from_codes(item['codes'][:], categories=cats)
    vals = item[:]
    return [v.decode() if isinstance(v, bytes) else v for v in vals]


with h5py.File(R + 'williams_epi.h5ad', 'r') as f:
    obs_grp = f['obs']
    site = np.array(read_h5_col(obs_grp, 'site'))
    sub = np.array(read_h5_col(obs_grp, 'sub')).astype(str)
    var_names = read_h5_col(f['var'], '_index')
    cg = f['layers']['counts']
    X = sp.csr_matrix((cg['data'][:], cg['indices'][:], cg['indptr'][:]), shape=(len(sub), len(var_names)))
    Xn = sp.diags(1e4 / np.asarray(X.sum(1)).ravel()) @ X
vmap = {g: i for i, g in enumerate(var_names)}


def basal_def_feat(mask, pc=1.0):
    vals = []
    for gs in T.target_gene:
        vs = []
        for g in str(gs).split(';'):
            g2 = g if g in vmap else alias.get(g, g)
            if g2 in vmap:
                vs.append(np.log2(float(Xn[mask, vmap[g2]].mean()) + pc))
        vals.append(max(vs) if vs else np.nan)
    return np.array(vals)


is_basal = np.isin(sub, ['0', '1', '13', '5', '6'])
for label, mask, pc in [
    ('Basal buccal mucosa only (1,048 cells)', is_basal & (site == 'buccal'), 1.0),
    ('Basal gingiva only (1,828 cells)', is_basal & (site == 'gingiva'), 1.0),
    ('Basal quiescent only excl cycling 5/6 (2,051 cells)', np.isin(sub, ['0', '1', '13']), 1.0),
    ('Basal cycling clusters 5/6 only (825 cells)', np.isin(sub, ['5', '6']), 1.0),
    ('Basal excl sub 13 (2,807 cells)', np.isin(sub, ['0', '1', '5', '6']), 1.0),
    ('All epithelial incl other (6,157 cells)', np.ones(len(sub), dtype=bool), 1.0),
    ('Basal log2(CP10K + 0.1)', is_basal, 0.1),
]:
    T['_tmp'] = basal_def_feat(mask, pc=pc)
    add_sens('Basal_definition', label, T, k='_tmp')

# 5. Resolution-ladder no-TROP2 control (whole-mucosa pseudobulk excl TACSTD2) & all 28 drop-one-target exclusions
add_sens('Resolution_noTROP2', 'Whole oral mucosa pseudobulk excl TACSTD2 (n=39)', T[T.tg != 'TACSTD2'], k='W_alltissue')
for tg_drop in sorted(T.tg.unique()):
    n_dr = int((T.tg == tg_drop).sum())
    tg_disp = 'NECTIN4' if tg_drop == 'PVRL4' else tg_drop
    add_sens('Drop_one_target', f'Drop {tg_disp} (k={n_dr} ADC{"s" if n_dr > 1 else ""})', T[T.tg != tg_drop])

S_df = pd.DataFrame(sens_rows)
print('\n=== SENSITIVITY SUMMARY (with Profile Likelihood 95% CIs) ===')
print(S_df.round(4).to_string())
S_df.to_csv(R + 'part1/sensitivity_subgroups.tsv', sep='\t', index=False)

# Descriptive small-N DXd note (kept out of hypothesis testing, n < 8)
d_dxd = T[T.payload_class == 'DXd']
print('\nDescriptive DXd-only subset (n=%d, not hypothesis-tested):' % len(d_dxd))
print(d_dxd[['drug', 'tg', 'y_pct', 'cens', 'W_basal']].to_string())
