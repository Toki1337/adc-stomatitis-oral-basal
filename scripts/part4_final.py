"""Part 4: OSCC surfaceome structural coupling, full-cohort & LOTO 95% PIs, PVRL4->NECTIN4 harmonization, and 3x3 grids."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import expit

from tobit_lib import PAYLOAD_MAP, design_matrix, fit_tobit

pd.set_option('display.width', 330)
ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'
surf = set(pd.read_csv(R + 'surfaceome_genes.txt', header=None)[0])
W = pd.read_csv(R + 'williams_kc_stats.tsv.gz', sep='\t', index_col=0)
C = pd.read_csv(R + 'choi_group_stats.tsv.gz', sep='\t', index_col=0)
P = pd.read_csv(R + 'puram_group_stats.tsv.gz', sep='\t', index_col=0)
Tc = pd.read_csv(R + 'tcga_oral_stats.tsv.gz', sep='\t', index_col=0)
Cp = pd.read_csv(R + 'cptac_oral_stats.tsv.gz', sep='\t', index_col=0)
G = pd.read_csv(R + 'gtex_median_tpm.tsv.gz', sep='\t', index_col=0)
Dm = pd.read_csv(R + 'depmap_24q4_summary.tsv', sep='\t').set_index('gene')
T = pd.read_csv(R + 'part3_axes.tsv', sep='\t')

# Harmonize legacy symbol PVRL4 -> NECTIN4 consistently across RNA matrices (matching Part 1 / Part 3 alias rule)
LEGACY_ALIASES = {'PVRL4': 'NECTIN4'}


def harmonize_legacy_index(df):
    new_idx = [LEGACY_ALIASES[g] if (g in LEGACY_ALIASES and LEGACY_ALIASES[g] not in df.index) else g for g in df.index]
    out = df.copy()
    out.index = new_idx
    return out


W = harmonize_legacy_index(W)
C = harmonize_legacy_index(C)
P = harmonize_legacy_index(P)
Tc = harmonize_legacy_index(Tc)

harm_log = pd.DataFrame([
    dict(
        legacy_symbol='PVRL4',
        current_hgnc_symbol='NECTIN4',
        tables_mapped='Williams GSE164241; Choi GSE181919; Puram GSE103322; TCGA-HNSC',
        choi_malig_det=float(C.loc['NECTIN4', 'choi_malig_OC_HPVneg_det']),
        puram_malig_det=float(P.loc['NECTIN4', 'puram_malig_primary_det']),
        passes_25pct_pool=bool(
            C.loc['NECTIN4', 'choi_malig_OC_HPVneg_det'] >= 0.25 and P.loc['NECTIN4', 'puram_malig_primary_det'] >= 0.25
        ),
        passes_50pct_screen=bool(
            C.loc['NECTIN4', 'choi_malig_OC_HPVneg_det'] >= 0.50 and P.loc['NECTIN4', 'puram_malig_primary_det'] >= 0.50
        ),
    )
])
harm_log.to_csv(R + 'part4_symbol_harmonization.tsv', sep='\t', index=False)

d = T.dropna(subset=['abs_oral_basal']).reset_index(drop=True)
Xm_full = design_matrix(d, 'abs_oral_basal')
full_fit = fit_tobit(d.y_pct.values, d.cens.values, Xm_full)
b_full = full_fit['beta']


def pred_fn(xb, beta=b_full):
    return expit(beta[0] + beta[-1] * np.asarray(xb)) * 100.0


# Cluster-bootstrap distribution helper (B=2,000 by target) for full cohort & LOTO folds
def cluster_boot_params(df_sub, seed=1, B=2000):
    rng = np.random.default_rng(seed)
    tgs = df_sub.tg.unique()
    tg_idx = {t: np.where(df_sub.tg.values == t)[0] for t in tgs}
    y_a = df_sub.y_pct.values
    c_a = df_sub.cens.values
    pc_a = df_sub.pclass.values
    X_a = design_matrix(df_sub, 'abs_oral_basal')
    bb_list, ss_list = [], []
    for _ in range(B):
        pick = rng.choice(tgs, len(tgs), replace=True)
        idx = np.concatenate([tg_idx[t] for t in pick])
        if c_a[idx].all() or len(np.unique(pc_a[idx])) < 3:
            continue
        r = fit_tobit(y_a[idx], c_a[idx], X_a[idx])
        if r['converged'] and np.isfinite(r['beta']).all() and 0.05 < r['sigma'] < 10.0:
            bb_list.append(r['beta'])
            ss_list.append(r['sigma'])
    return np.array(bb_list), np.array(ss_list)


genes = sorted(surf & set(W.index) & set(C.index) & set(P.index))
M = pd.DataFrame(index=genes)
M['oral_basal_raw'] = W.loc[genes, 'w_basal_mean']
M['oral_basal'] = np.log2(M.oral_basal_raw + 1.0)
M['oral_basal_det'] = W.loc[genes, 'w_basal_det']
M['pred_stom_pct'] = pred_fn(M.oral_basal)
M['tum_choi_raw'] = C.loc[genes, 'choi_malig_OC_HPVneg_mean']
M['tum_choi'] = np.log2(M.tum_choi_raw + 1.0)
M['tum_choi_det'] = C.loc[genes, 'choi_malig_OC_HPVneg_det']
M['tum_puram_raw'] = P.loc[genes, 'puram_malig_primary_mean']
M['tum_puram'] = np.log2(M.tum_puram_raw + 1.0)
M['tum_puram_det'] = P.loc[genes, 'puram_malig_primary_det']

# Shifted log2 differences (Delta log2(CP10K+1)) and linear ratios (T+eps)/(N+eps)
M['win_choi'] = M.tum_choi - M.oral_basal
M['win_puram'] = M.tum_puram - M.oral_basal
M['lin_fc_choi_raw'] = M.tum_choi_raw / M.oral_basal_raw.replace(0, np.nan)
M['lin_fc_choi_01'] = (M.tum_choi_raw + 0.1) / (M.oral_basal_raw + 0.1)
M['lin_fc_puram_01'] = (M.tum_puram_raw + 0.1) / (M.oral_basal_raw + 0.1)

M['tcga_lfc'] = Tc.reindex(genes)['tcga_paired_lfc']
M['tcga_frac_up'] = Tc.reindex(genes)['tcga_paired_frac_up']
M['cptac_prot_lfc'] = Cp.reindex(genes)['cptac_prot_paired_lfc']
M['cptac_prot_frac_up'] = Cp.reindex(genes)['cptac_prot_paired_frac_up']
M['cptac_prot_n'] = Cp.reindex(genes)['cptac_prot_paired_n']
M['cptac_rna_lfc'] = Cp.reindex(genes)['cptac_rna_paired_lfc']
# Four-organ & extended 24-organ GTEx systemic + CNS panels (excluding squamous skin/esophagus/vagina/salivary mucosa)
VITAL4_COLS = {
    'Heart - Left Ventricle': 'Heart',
    'Liver': 'Liver',
    'Kidney - Cortex': 'Kidney',
    'Nerve - Tibial': 'Nerve',
}
NONBRAIN_OTHER_COLS = {
    'Lung': 'Lung',
    'Stomach': 'Stomach',
    'Colon - Transverse': 'Colon',
    'Small Intestine - Terminal Ileum': 'Small intestine',
    'Pancreas': 'Pancreas',
    'Spleen': 'Spleen',
    'Whole Blood': 'Whole blood',
}
BRAIN_COLS = {
    'Brain - Amygdala': 'Amygdala',
    'Brain - Anterior cingulate cortex (BA24)': 'Anterior cingulate',
    'Brain - Caudate (basal ganglia)': 'Caudate',
    'Brain - Cerebellar Hemisphere': 'Cerebellar hemisphere',
    'Brain - Cerebellum': 'Cerebellum',
    'Brain - Cortex': 'Cerebral cortex',
    'Brain - Frontal Cortex (BA9)': 'Frontal cortex',
    'Brain - Hippocampus': 'Hippocampus',
    'Brain - Hypothalamus': 'Hypothalamus',
    'Brain - Nucleus accumbens (basal ganglia)': 'Nucleus accumbens',
    'Brain - Putamen (basal ganglia)': 'Putamen',
    'Brain - Spinal cord (cervical c-1)': 'Spinal cord',
    'Brain - Substantia nigra': 'Substantia nigra',
}
EXT24_COLS = {**VITAL4_COLS, **NONBRAIN_OTHER_COLS, **BRAIN_COLS}
NONBRAIN_EXT11_COLS = {**VITAL4_COLS, **NONBRAIN_OTHER_COLS}

for nm, col in [
    ('gtex_heart', 'Heart - Left Ventricle'),
    ('gtex_liver', 'Liver'),
    ('gtex_kidney', 'Kidney - Cortex'),
    ('gtex_lung', 'Lung'),
    ('gtex_nerve', 'Nerve - Tibial'),
    ('gtex_esoph', 'Esophagus - Mucosa'),
]:
    M[nm] = G.reindex(genes)[col]

G_sub = G.reindex(genes)
v4_mat = G_sub[list(VITAL4_COLS.keys())]
nb_mat = G_sub[list(NONBRAIN_EXT11_COLS.keys())]
br_mat = G_sub[list(BRAIN_COLS.keys())]
ex_mat = G_sub[list(EXT24_COLS.keys())]

M['vital_max_tpm'] = v4_mat.max(axis=1)
M['vital_max_tissue'] = v4_mat.dropna(how='all').idxmax(axis=1).map(VITAL4_COLS).reindex(genes)
M['nonbrain_ext_max_tpm'] = nb_mat.max(axis=1)
M['nonbrain_ext_max_tissue'] = nb_mat.dropna(how='all').idxmax(axis=1).map(NONBRAIN_EXT11_COLS).reindex(genes)
M['brain_max_tpm'] = br_mat.max(axis=1)
M['brain_max_tissue'] = br_mat.dropna(how='all').idxmax(axis=1).map(BRAIN_COLS).reindex(genes)
M['ext_max_tpm'] = ex_mat.max(axis=1)
M['ext_max_tissue'] = ex_mat.dropna(how='all').idxmax(axis=1).map(EXT24_COLS).reindex(genes)
M['gtex54_max_tpm'] = G_sub.max(axis=1)
M['gtex54_max_tissue'] = G_sub.dropna(how='all').idxmax(axis=1).reindex(genes)

E = M[(M.tum_choi_det >= 0.25) & (M.tum_puram_det >= 0.25)].copy()
print('surfaceome genes measured in all 3 single-cell sources: %d ; expressed in OSCC malignant cells (>=25%% both): %d' % (len(M), len(E)))
print('\n=== structural coupling & window scale comparison ===')
print('Spearman(OSCC malignant, healthy oral basal KC): Choi %.3f  Puram %.3f' % (
    stats.spearmanr(E.tum_choi, E.oral_basal)[0], stats.spearmanr(E.tum_puram, E.oral_basal)[0],
))
print('Delta log2(CP10K+1) >= 1 (Choi): %.1f%% (%d/%d) | dual (Choi>=1 & Puram>=0.5): %.1f%% (%d/%d)' % (
    100 * (E.win_choi >= 1).mean(), int((E.win_choi >= 1).sum()), len(E),
    100 * ((E.win_choi >= 1) & (E.win_puram >= 0.5)).mean(), int(((E.win_choi >= 1) & (E.win_puram >= 0.5)).sum()), len(E),
))
print('Delta log2(CP10K+1) >= 2 (Choi): %.1f%% (%d/%d) | dual (Choi>=2 & Puram>=1): %.1f%% (%d/%d)' % (
    100 * (E.win_choi >= 2).mean(), int((E.win_choi >= 2).sum()), len(E),
    100 * ((E.win_choi >= 2) & (E.win_puram >= 1)).mean(), int(((E.win_choi >= 2) & (E.win_puram >= 1)).sum()), len(E),
))
print('Linear unshifted (T/N >= 2, Choi raw): %.1f%% (%d/%d) | Linear (T+0.1)/(N+0.1) >= 2: Choi %.1f%% (%d/%d), dual %.1f%% (%d/%d)' % (
    100 * ((E.tum_choi_raw / np.maximum(E.oral_basal_raw, 1e-12)) >= 2).mean(),
    int(((E.tum_choi_raw / np.maximum(E.oral_basal_raw, 1e-12)) >= 2).sum()), len(E),
    100 * (E.lin_fc_choi_01 >= 2).mean(), int((E.lin_fc_choi_01 >= 2).sum()), len(E),
    100 * ((E.lin_fc_choi_01 >= 2) & (E.lin_fc_puram_01 >= 2)).mean(),
    int(((E.lin_fc_choi_01 >= 2) & (E.lin_fc_puram_01 >= 2)).sum()), len(E),
))
print('median predicted stomatitis across %d targets: %.1f%% (IQR %.1f-%.1f)' % (
    len(E), E.pred_stom_pct.median(), *np.percentile(E.pred_stom_pct, [25, 75]),
))

# Compute coherent full-cohort 95% CIs/PIs AND separate zero-antigen-overlap LOTO intervals for in-training targets
train_max_basal = float(d.abs_oral_basal.max())
E['extrapolated_above_train_max'] = E.oral_basal > (train_max_basal + 1e-9)
E['fit_source'] = 'full_cohort_44ADC'
print(f'\nComputing coherent full-cohort 95% PIs (train max basal = {train_max_basal:.3f}; extrapolated genes = {int(E.extrapolated_above_train_max.sum())})...')
full_bb, full_ss = cluster_boot_params(d, seed=1, B=2000)

train_gene_sets = [set(str(gs).split(';')) for gs in d.target_gene]
all_train_genes = sorted({g for s in train_gene_sets for g in s} & set(E.index))
loto_cache = {}
for tg in all_train_genes:
    sub_d = d[~np.array([tg in s for s in train_gene_sets])].reset_index(drop=True)
    loto_fit = fit_tobit(sub_d.y_pct.values, sub_d.cens.values, design_matrix(sub_d, 'abs_oral_basal'))
    bb, ss = cluster_boot_params(sub_d, seed=1, B=2000)
    loto_cache[tg] = (loto_fit['beta'], bb, ss)

rng_pi = np.random.default_rng(42)
rng_loto = np.random.default_rng(42)
mean_lo, mean_hi, pi_lo, pi_hi = [], [], [], []
loto_pt, loto_lo, loto_hi = [], [], []
for gene, row in E.iterrows():
    x = float(row.oral_basal)
    lp_f = full_bb[:, 0] + full_bb[:, -1] * x
    m_ci = np.quantile(expit(lp_f) * 100.0, [0.025, 0.975])
    z_f = lp_f + rng_pi.normal(0.0, full_ss)
    p_pi = np.quantile(expit(z_f) * 100.0, [0.025, 0.975])
    mean_lo.append(float(m_ci[0]))
    mean_hi.append(float(m_ci[1]))
    pi_lo.append(float(p_pi[0]))
    pi_hi.append(float(p_pi[1]))

    if gene in loto_cache:
        b_pt, bb_l, ss_l = loto_cache[gene]
        loto_pt.append(float(expit(b_pt[0] + b_pt[-1] * x) * 100.0))
        z_l = bb_l[:, 0] + bb_l[:, -1] * x + rng_loto.normal(0.0, ss_l)
        p_lpi = np.quantile(expit(z_l) * 100.0, [0.025, 0.975])
        loto_lo.append(float(p_lpi[0]))
        loto_hi.append(float(p_lpi[1]))
    else:
        loto_pt.append(np.nan)
        loto_lo.append(np.nan)
        loto_hi.append(np.nan)

E['mean_ci_lo'] = mean_lo
E['mean_ci_hi'] = mean_hi
E['pred_pi_lo'] = pi_lo
E['pred_pi_hi'] = pi_hi
E['loto_pred_stom_pct'] = loto_pt
E['loto_pi_lo'] = loto_lo
E['loto_pi_hi'] = loto_hi

# Four exploratory ranking criteria + 3x3 threshold stability grids (4-organ & extended 24-organ panels)
E['pass_cover'] = (E.tum_choi_det >= 0.5) & (E.tum_puram_det >= 0.5)
E['pass_tox'] = E.pred_stom_pct <= 12.0
E['pass_tcga'] = (E.tcga_lfc > 0) & (E.tcga_frac_up >= 0.6)
E['pass_prot'] = (E.cptac_prot_lfc > 0) & (E.cptac_prot_n >= 15) & (E.cptac_prot_frac_up >= 0.5)
E['n_pass'] = E[['pass_cover', 'pass_tox', 'pass_tcga', 'pass_prot']].sum(axis=1)

grid_count = np.zeros(len(E), dtype=int)
grid_count_ext = np.zeros(len(E), dtype=int)
for stom_cut in [10.0, 12.0, 15.0]:
    for vit_cut in [10.0, 15.0, 30.0]:
        ok = E['pass_cover'] & E['pass_tcga'] & E['pass_prot'] & (E.pred_stom_pct <= stom_cut) & (E.vital_max_tpm <= vit_cut)
        ok_ext = E['pass_cover'] & E['pass_tcga'] & E['pass_prot'] & (E.pred_stom_pct <= stom_cut) & (E.ext_max_tpm <= vit_cut)
        grid_count += ok.astype(int)
        grid_count_ext += ok_ext.astype(int)
E['grid_pass'] = grid_count
E['robust_7of9'] = E.grid_pass >= 7
E['grid_pass_ext'] = grid_count_ext
E['robust_ext_7of9'] = E.grid_pass_ext >= 7
E['cptac_prot_note'] = np.where(
    (E.cptac_prot_lfc > 0) & (E.cptac_prot_lfc < 0.05),
    'marginal / not upregulated at protein level (LFC < 0.05)',
    np.where(E.pass_prot, 'upregulated', 'not upregulated'),
)
E['depmap_common_essential'] = Dm.reindex(E.index)['common_essential']
E['depmap_median_chronos'] = Dm.reindex(E.index)['median_chronos']
E['depmap_frac_lt_neg05'] = Dm.reindex(E.index)['frac_lt_neg05']

win = E[E.n_pass == 4].sort_values('tum_choi', ascending=False)
print('\n=== 4/4 criteria winners: n = %d / %d ; 4-organ >=7/9: %d ; 24-organ ext >=6/9: %d (<=15 TPM: %d) ===' % (
    len(win), len(E), int(E.robust_7of9.sum()), int((E.grid_pass_ext >= 6).sum()), int((win.ext_max_tpm <= 15.0).sum()),
))
show_cols = [
    'oral_basal', 'pred_stom_pct', 'pred_pi_lo', 'pred_pi_hi',
    'tum_choi_det', 'tum_puram_det', 'tcga_lfc', 'cptac_prot_lfc',
    'vital_max_tpm', 'vital_max_tissue', 'ext_max_tpm', 'ext_max_tissue',
    'grid_pass', 'grid_pass_ext', 'cptac_prot_note',
]
print(E[(E.grid_pass >= 1) | (E.n_pass == 4)].sort_values(['grid_pass', 'grid_pass_ext', 'tum_choi'], ascending=[False, False, False])[show_cols].round(2).to_string())

E.to_csv(R + 'part4_final_ranking.tsv', sep='\t')
win.to_csv(R + 'part4_winners.tsv', sep='\t')
