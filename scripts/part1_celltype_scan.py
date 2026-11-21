"""Part 1c: Unbiased specificity scan across 155 HPA single-cell types and 54 GTEx bulk tissues with BH-FDR."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

from tobit_lib import PAYLOAD_MAP, bh_fdr, lr_test

ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'
T = pd.read_csv(R + 'part1/adc_table_with_features.tsv', sep='\t')
T['pclass'] = T.payload_class.map(PAYLOAD_MAP)
T['cens'] = T.stom_any.isna().astype(bool)
T['y_pct'] = np.where(T.cens, T.threshold_pct, T.stom_any)
T['tg'] = T.target_gene.str.split(';').str[0]
alias = {'NECTIN4': 'PVRL4'}

Hs = pd.read_csv(R + 'hpa_sc_ncpm.tsv.gz', sep='\t', index_col=0)
G = pd.read_csv(R + 'gtex_median_tpm.tsv.gz', sep='\t', index_col=0)


def feat(gs, fn):
    vs = [fn(g) for g in str(gs).split(';')]
    vs = [v for v in vs if not np.isnan(v)]
    return max(vs) if vs else np.nan


def scan(df, label):
    rows = []
    for ct in df.columns:
        fn = lambda g: (
            np.log2(df.loc[g if g in df.index else alias.get(g, g), ct] + 1.0)
            if (g in df.index or alias.get(g, g) in df.index) else np.nan
        )
        T['_x'] = [feat(g, fn) for g in T.target_gene]
        d = T.dropna(subset=['_x'])
        if d._x.std() == 0:
            continue
        r_all = lr_test(d, '_x')
        dn = d[d.tg != 'TACSTD2']
        r_no = lr_test(dn, '_x')
        imp = np.where(d.cens, d.y_pct / 2.0, d.y_pct)
        sp_rho, sp_p = stats.spearmanr(d['_x'], imp)
        rows.append(dict(
            source=label, celltype=ct,
            LR=r_all['LR'], p=r_all['p'],
            LR_noTROP2=r_no['LR'], p_noTROP2=r_no['p'],
            spearman_rho=float(sp_rho), p_spearman=float(sp_p),
        ))
    out = pd.DataFrame(rows)
    out['q'] = bh_fdr(out['p'].values)
    out['q_noTROP2'] = bh_fdr(out['p_noTROP2'].values)
    out['q_spearman'] = bh_fdr(out['p_spearman'].values)
    out['rank_all'] = out.LR.rank(ascending=False)
    out['rank_noTROP2'] = out.LR_noTROP2.rank(ascending=False)
    return out.sort_values('rank_all')


S = pd.concat([scan(Hs, 'HPA_singlecell'), scan(G, 'GTEx')], ignore_index=True)
S.to_csv(R + 'part1/celltype_tissue_scan.tsv', sep='\t', index=False)

for src in ['HPA_singlecell', 'GTEx']:
    s = S[S.source == src]
    print('==', src, 'total:', len(s), '| q < 0.05:', int((s.q < 0.05).sum()))
    print(s[['celltype', 'LR', 'p', 'q', 'LR_noTROP2', 'p_noTROP2', 'q_noTROP2', 'rank_all', 'rank_noTROP2']].head(12).round(5).to_string(index=False))
