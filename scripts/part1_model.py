"""Part 1: Link blinded clinical stomatitis data to expression features and build resolution ladder."""

from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy import stats

from tobit_lib import PAYLOAD_MAP, design_matrix, exact_hess_se, fit_tobit, lr_test

ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'
T = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_trials.tsv', sep='\t')
T = T[T.listed.isin(['yes', 'not_listed'])].copy()
T = T[~(T.stom_any.isna() & T.threshold_pct.isna())].copy()
T['pclass'] = T.payload_class.map(PAYLOAD_MAP)
T['cens'] = T.stom_any.isna().astype(bool)
T['y_pct'] = np.where(T.cens, T.threshold_pct, T.stom_any)
T['tg'] = T.target_gene.str.split(';').str[0]
print('ADCs', len(T), 'numeric', int((~T.cens).sum()), 'censored', int(T.cens.sum()), 'targets', T.tg.nunique())

alias = {'NECTIN4': 'PVRL4'}
W = pd.read_csv(R + 'williams_kc_stats.tsv.gz', sep='\t', index_col=0)
C = pd.read_csv(R + 'choi_group_stats.tsv.gz', sep='\t', index_col=0)
G = pd.read_csv(R + 'gtex_median_tpm.tsv.gz', sep='\t', index_col=0)
Tc = pd.read_csv(R + 'tcga_oral_stats.tsv.gz', sep='\t', index_col=0)
Cp = pd.read_csv(R + 'cptac_oral_stats.tsv.gz', sep='\t', index_col=0)
H = pd.read_csv(R + 'hpa_ihc_scored.tsv.gz', sep='\t')
Hd = pd.read_csv(R + 'part1/hpa_ihc_dab_genes.tsv', sep='\t', index_col=0)
Hs = pd.read_csv(R + 'hpa_sc_ncpm.tsv.gz', sep='\t', index_col=0)
hpa_oral = H[H.Tissue == 'Oral mucosa'].groupby('Gene name').score.max()


def val(df, col, g, log=True, pc=1.0):
    g2 = g if g in df.index else alias.get(g, g)
    if g2 not in df.index:
        return np.nan
    v = df.loc[g2, col]
    return np.log2(v + pc) if log else v


def feat(gs, fn):
    vs = [fn(g) for g in str(gs).split(';')]
    vs = [v for v in vs if not np.isnan(v)]
    return max(vs) if vs else np.nan


# Extract quiescent (0,1,13) and cycling (5,6) basal subclusters from williams_epi.h5ad
def _read_h5_col(grp, key):
    item = grp[key]
    if isinstance(item, h5py.Group):
        cats = [c.decode() if isinstance(c, bytes) else str(c) for c in item['categories'][:]]
        return pd.Categorical.from_codes(item['codes'][:], categories=cats)
    vals = item[:]
    return [v.decode() if isinstance(v, bytes) else v for v in vals]


with h5py.File(R + 'williams_epi.h5ad', 'r') as f:
    sub = np.array(_read_h5_col(f['obs'], 'sub')).astype(str)
    var_names = _read_h5_col(f['var'], '_index')
    cg = f['layers']['counts']
    X_epi = sp.csr_matrix((cg['data'][:], cg['indices'][:], cg['indptr'][:]), shape=(len(sub), len(var_names)))
    Xn_epi = sp.diags(1e4 / np.asarray(X_epi.sum(1)).ravel()) @ X_epi
vmap = {g: i for i, g in enumerate(var_names)}


def epi_mask_val(mask, g, pc=1.0):
    g2 = g if g in vmap else alias.get(g, g)
    if g2 not in vmap:
        return np.nan
    return float(np.log2(float(Xn_epi[mask, vmap[g2]].mean()) + pc))


m_basal_q = np.isin(sub, ['0', '1', '13'])
m_basal_p = np.isin(sub, ['5', '6'])

F = {
    'W_basal': lambda g: val(W, 'w_basal_mean', g),
    'W_basal_quiescent': lambda g: epi_mask_val(m_basal_q, g),
    'W_basal_prolif': lambda g: epi_mask_val(m_basal_p, g),
    'W_basal_donormed': lambda g: val(W, 'w_basal_donormed', g),
    'W_basal_det': lambda g: val(W, 'w_basal_det', g, log=False),
    'W_suprabasal': lambda g: val(W, 'w_suprabasal_mean', g),
    'W_superficial': lambda g: val(W, 'w_superficial_mean', g),
    'W_alltissue': lambda g: val(W, 'w_alltissue_mean', g),
    'Choi_normal_epi': lambda g: val(C, 'choi_epi_NL_OC_mean', g),
    'Choi_panHNSCC_normal_epi': lambda g: val(C, 'choi_epi_NL_mean', g),
    'HPA_IHC_oral': lambda g: hpa_oral.get(g, np.nan),
    'HPA_sc_esoph_basal': lambda g: val(Hs, 'esophageal basal cells', g),
    'HPA_sc_skin_basal_kc': lambda g: val(Hs, 'basal keratinocytes', g),
    'GTEx_esoph_mucosa': lambda g: val(G, 'Esophagus - Mucosa', g),
    'GTEx_minor_salivary': lambda g: val(G, 'Minor Salivary Gland', g),
    'GTEx_skin_sun': lambda g: val(G, 'Skin - Sun Exposed (Lower leg)', g),
    'TCGA_normal_oral': lambda g: val(Tc, 'tcga_nor_mean', g, log=False),
    'CPTAC_NAT_protein': lambda g: val(Cp, 'cptac_prot_nor_mean', g, log=False),
    'HPA_esoph_DAB': lambda g: val(Hd, 'esoph_dab_basal', g, log=False),
    'HPA_skin_DAB': lambda g: val(Hd, 'skin_dab_basal', g, log=False),
    'HPA_oral_DAB': lambda g: val(Hd, 'oral_dab_basal', g, log=False),
}

for k, fn in F.items():
    T[k] = [feat(g, fn) for g in T.target_gene]
T.to_csv(R + 'part1/adc_table_with_features.tsv', sep='\t', index=False)

rows = []
for k in F:
    d = T.dropna(subset=[k])
    res_p = lr_test(d, k, payload=True, compute_se=True)
    res_np = lr_test(d, k, payload=False, compute_se=False)
    sd = float(d[k].std())
    imp = np.where(d.cens, d.y_pct / 2.0, d.y_pct)
    rho, prho = stats.spearmanr(d[k], imp)
    rows.append(dict(
        feature=k, n=len(d), n_num=int((~d.cens).sum()),
        beta=res_p['beta'], beta_per_SD=res_p['beta'] * sd, se=res_p['se'],
        LR=res_p['LR'], p=res_p['p'], aic=res_p['aic'], daic=res_p['daic'],
        LR_nopayload=res_np['LR'], p_nopayload=res_np['p'],
        spearman=rho, p_spearman=prho, converged=res_p['converged'],
    ))

res = pd.DataFrame(rows).sort_values('p')
print(res.round(4).to_string())
res.to_csv(R + 'part1/feature_comparison.tsv', sep='\t', index=False)

# Resolution ladder with cluster-bootstrap (B=2,000 by target) delta-LR vs W_basal
ladder_specs = [
    ('HPA_sc_esoph_basal', 'HPA scRNA: Esophageal basal cells', 'Single-cell RNA (nCPM)'),
    ('W_basal_quiescent', 'Oral quiescent basal KC (scRNA)', 'Single-cell RNA (CP10K)'),
    ('W_basal', 'Oral basal KC overall (scRNA)', 'Single-cell RNA (CP10K)'),
    ('HPA_sc_skin_basal_kc', 'HPA scRNA: Skin basal KC', 'Single-cell RNA (nCPM)'),
    ('TCGA_normal_oral', 'TCGA-HNSC adjacent normal oral', 'Bulk RNA-seq'),
    ('GTEx_esoph_mucosa', 'GTEx bulk: Esophagus mucosa', 'Bulk RNA-seq (TPM)'),
    ('W_basal_prolif', 'Oral cycling basal KC (scRNA)', 'Single-cell RNA (CP10K)'),
    ('W_suprabasal', 'Oral suprabasal KC (scRNA)', 'Single-cell RNA (CP10K)'),
    ('HPA_esoph_DAB', 'HPA esophagus IHC (rgb2hed DAB)', 'Immunohistochemistry (OD)'),
    ('CPTAC_NAT_protein', 'CPTAC oral NAT proteome', 'TMT proteomics (log2)'),
    ('W_superficial', 'Oral superficial KC (scRNA)', 'Single-cell RNA (CP10K)'),
    ('W_alltissue', 'Whole oral mucosa pseudobulk', 'Unstratified scRNA (CP10K)'),
    ('HPA_skin_DAB', 'HPA skin IHC (rgb2hed DAB)', 'Immunohistochemistry (OD)'),
    ('HPA_IHC_oral', 'HPA oral mucosa IHC', 'Immunohistochemistry'),
    ('HPA_oral_DAB', 'HPA oral mucosa IHC (rgb2hed DAB)', 'Immunohistochemistry (OD)'),
]

obs_lad = {}
for col, lab, mod in ladder_specs:
    d = T.dropna(subset=[col])
    rc = lr_test(d, col, payload=True)
    rb = lr_test(d, 'W_basal', payload=True)
    obs_lad[col] = (lab, mod, len(d), rc['LR'], rc['p'], rc['aic'], rc['daic'], rb['LR'] - rc['LR'])

rng = np.random.default_rng(1)
tgs = T.tg.unique()
tg_idx = {t: np.where(T.tg.values == t)[0] for t in tgs}
boot_diffs = {col: [] for col, _, _ in ladder_specs}

for _ in range(2000):
    pick = rng.choice(tgs, len(tgs), replace=True)
    idx = np.concatenate([tg_idx[t] for t in pick])
    D = T.iloc[idx].reset_index(drop=True)
    if D.cens.all() or D.pclass.nunique() < 2:
        continue
    r0_all = fit_tobit(D.y_pct.values, D.cens.values, design_matrix(D))
    r1_b_all = fit_tobit(D.y_pct.values, D.cens.values, design_matrix(D, 'W_basal'))
    if not (r0_all['converged'] and r1_b_all['converged']):
        continue
    lr_b_all = max(0.0, 2.0 * (r0_all['nll'] - r1_b_all['nll']))
    for col, _, _ in ladder_specs:
        if col == 'W_basal':
            boot_diffs[col].append(0.0)
            continue
        if not D[col].isna().any():
            r1_c = fit_tobit(D.y_pct.values, D.cens.values, design_matrix(D, col))
            if r1_c['converged']:
                boot_diffs[col].append(lr_b_all - max(0.0, 2.0 * (r0_all['nll'] - r1_c['nll'])))
        else:
            Dc = D.dropna(subset=[col])
            if len(Dc) < 15 or Dc.cens.all() or Dc[col].std() < 1e-8:
                continue
            rc0 = fit_tobit(Dc.y_pct.values, Dc.cens.values, design_matrix(Dc))
            rc1 = fit_tobit(Dc.y_pct.values, Dc.cens.values, design_matrix(Dc, col))
            rb1 = fit_tobit(Dc.y_pct.values, Dc.cens.values, design_matrix(Dc, 'W_basal'))
            if rc0['converged'] and rc1['converged'] and rb1['converged']:
                boot_diffs[col].append(max(0.0, 2.0 * (rc0['nll'] - rb1['nll'])) - max(0.0, 2.0 * (rc0['nll'] - rc1['nll'])))

lad_rows = []
for col, _, _ in ladder_specs:
    lab, mod, n, lr, p, aic, daic, dlr = obs_lad[col]
    arr = np.array(boot_diffs[col])
    lo, hi = (0.0, 0.0) if col == 'W_basal' else np.quantile(arr, [0.025, 0.975])
    excl_zero = bool((col != 'W_basal') and ((lo > 0) or (hi < 0)))
    lad_rows.append(dict(
        feature=col, label=lab, modality=mod, n=n, LR=lr, p_LR=p, aic=aic, daic=daic,
        dLR_basal_minus_src=dlr, dLR_boot_ci_lo=lo, dLR_boot_ci_hi=hi,
        ci_excludes_zero=excl_zero, n_boot=len(arr),
    ))
lad_df = pd.DataFrame(lad_rows)
print('\n=== RESOLUTION LADDER (cluster bootstrap B=2,000 by target) ===')
print(lad_df.round(4).to_string())
lad_df.to_csv(R + 'part1/resolution_ladder.tsv', sep='\t', index=False)
