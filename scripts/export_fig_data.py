"""Export synchronized figure auxiliary datasets in fast, deterministic form."""

import os
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.formula.api as smf
from scipy.special import expit

from tobit_lib import PAYLOAD_MAP, design_matrix, fit_tobit

ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'


# 1. Export epithelial compartment marker & ADC target stats from williams_epi.h5ad (if present)
H5_PATH = R + 'williams_epi.h5ad'
EPI_OUT = R + 'fig_epi_markers.tsv'
if os.path.exists(H5_PATH):
    def read_h5_col(grp, key):
        item = grp[key]
        if isinstance(item, h5py.Group):
            cats = [c.decode() if isinstance(c, bytes) else str(c) for c in item['categories'][:]]
            return pd.Categorical.from_codes(item['codes'][:], categories=cats)
        vals = item[:]
        return [v.decode() if isinstance(v, bytes) else v for v in vals]

    with h5py.File(H5_PATH, 'r') as f:
        sub = np.array(read_h5_col(f['obs'], 'sub')).astype(str)
        var_names = read_h5_col(f['var'], '_index')
        cg = f['layers']['counts']
        X = sp.csr_matrix((cg['data'][:], cg['indices'][:], cg['indptr'][:]), shape=(len(sub), len(var_names)))
        lib = np.asarray(X.sum(1)).ravel()
        Xn = sp.diags(1e4 / lib) @ X

    vmap = {g: i for i, g in enumerate(var_names)}
    comp_map = {'basal': ['0', '1', '13', '5', '6'], 'suprabasal': ['2', '3', '9'], 'superficial': ['4', '8', '11']}
    comp = np.full(len(sub), 'other', dtype=object)
    for k, v in comp_map.items():
        comp[np.isin(sub, v)] = k

    genes_fig1 = [
        ('KRT5', 'Basal marker'), ('KRT14', 'Basal marker'), ('KRT15', 'Basal marker'),
        ('COL17A1', 'Basal marker'), ('TP63', 'Basal marker'), ('MKI67', 'Cycling basal'),
        ('KRT13', 'Differentiation'), ('KRT4', 'Differentiation'), ('SPRR3', 'Superficial'),
        ('CRNN', 'Superficial'), ('TACSTD2', 'ADC target'), ('CD44', 'ADC target'),
        ('GPNMB', 'ADC target'), ('EGFR', 'ADC target'), ('ERBB2', 'ADC target'),
        ('PTK7', 'ADC target'), ('PVRL4', 'ADC target (NECTIN4)'), ('F3', 'ADC target'),
        ('FOLR1', 'ADC target'),
    ]
    rows = []
    for g, grp in genes_fig1:
        idx = vmap[g]
        for c in ['basal', 'suprabasal', 'superficial']:
            m = (comp == c)
            vals = Xn[m, idx].toarray().ravel()
            rows.append(dict(
                gene=g if g != 'PVRL4' else 'NECTIN4',
                gene_group=grp,
                compartment=c,
                mean_log2=float(np.log2(vals.mean() + 1.0)),
                det_frac=float((vals > 0).mean()),
            ))
    pd.DataFrame(rows).to_csv(EPI_OUT, sep='\t', index=False)
elif os.path.exists(EPI_OUT):
    print('williams_epi.h5ad absent; preserving pre-computed', EPI_OUT)
else:
    raise FileNotFoundError(f'Missing both {H5_PATH} and {EPI_OUT}')

# 2. Export fitted Tobit curves for Figure 2B
T = pd.read_csv(R + 'part3_axes.tsv', sep='\t')
d = T.dropna(subset=['abs_oral_basal']).copy()
Xm = design_matrix(d, 'abs_oral_basal')
fit = fit_tobit(d.y_pct.values, d.cens.values, Xm)
b = fit['beta']
xs = np.linspace(0.0, 2.45, 200)
crows = []
for pc, dt, dd in [('topoI', 0.0, 0.0), ('DNA', 0.0, 1.0), ('tubulin', 1.0, 0.0)]:
    pred = expit(b[0] + b[1] * dt + b[2] * dd + b[3] * xs) * 100.0
    for xval, pval in zip(xs, pred):
        crows.append(dict(pclass=pc, abs_oral_basal=xval, pred_pct=pval))
pd.DataFrame(crows).to_csv(R + 'fig_tobit_curves.tsv', sep='\t', index=False)

# 3. Export FAERS summary (exposure-offset Negative-Binomial GLM + 10k target permutations + OLS/Spearman sensitivity) & pairwise RORs
from part1_faers2 import run_faers_analysis
run_faers_analysis()
print('Exported figure auxiliary tables successfully.')
