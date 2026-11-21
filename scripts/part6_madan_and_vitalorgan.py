"""Part 6: Benchmark against Madan et al. (2023) dual selectivity+safety nominees and vital-organ/DepMap checks."""

from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'
E = pd.read_csv(R + 'part4_final_ranking.tsv', sep='\t', index_col=0)

madan_meta = {
    'CDH3': 'Madan Dual-pass (Selectivity + Safety)',
    'COL17A1': 'Madan Dual-pass (Selectivity + Safety)',
    'DSC3': 'Madan Dual-pass (Selectivity + Safety)',
    'DSG3': 'Madan Dual-pass (Selectivity + Safety)',
    'PKP3': 'Madan Dual-pass (Intracellular plaque; excl.)',
    'TACSTD2': 'Madan Selectivity-only + Clinical',
}

print('=== Madan et al. (2023) HNSCC targets on the calibrated absolute axis ===')
M_sub = E.reindex(list(madan_meta.keys()))[
    ['oral_basal', 'pred_stom_pct', 'pred_pi_lo', 'pred_pi_hi', 'pass_cover', 'pass_tox', 'pass_tcga', 'pass_prot', 'n_pass']
].copy()
M_sub['madan_status'] = [madan_meta[g] for g in M_sub.index]
print(M_sub.round(2).to_string())

W = E[E.n_pass == 4]
v15 = W[W.vital_max_tpm <= 15].sort_values('tum_choi', ascending=False)
ext15 = W[W.ext_max_tpm <= 15].sort_values('tum_choi', ascending=False)
nb15 = W[W.nonbrain_ext_max_tpm <= 15].sort_values('tum_choi', ascending=False)
rob = W[W.robust_7of9].sort_values('tum_choi', ascending=False)
print('\n4/4 winners: %d | <=15 TPM (4-organ): %d | <=15 TPM (11 non-CNS ext): %d | <=15 TPM (24 ext): %d | 4-organ >=7/9: %d' % (
    len(W), len(v15), len(nb15), len(ext15), len(rob),
))
print('<=15 TPM (4-organ):', ', '.join(v15.index))
print('<=15 TPM (24-organ ext):', ', '.join(f"{g} ({W.loc[g, 'ext_max_tpm']:.2f} {W.loc[g, 'ext_max_tissue']})" for g in ext15.index))
print('4-organ Robust (>=7/9) candidates:', ', '.join(f"{g} (4-org {int(v)}/9, ext {int(W.loc[g, 'grid_pass_ext'])}/9)" for g, v in rob['grid_pass'].items()))

# DepMap 24Q4 common-essential check from repository summary table
depmap_path = f'{ROOT}/results/depmap_24q4_summary.tsv'
dm = pd.read_csv(depmap_path, sep='\t').set_index('gene')
print('\n=== DepMap 24Q4 check on 12 vital-low 4/4 winners ===')
for g in v15.index:
    if g in dm.index:
        r = dm.loc[g]
        print(f"  {g:8s}: common_essential={bool(r['common_essential'])} | median_Chronos={r['median_chronos']:+.3f} | frac<-0.5={r['frac_lt_neg05']:.3f}")
