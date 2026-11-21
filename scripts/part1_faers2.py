"""Part 1 FAERS pharmacovigilance replication: exposure-offset Negative-Binomial GLM + equal-weight OLS sensitivity (both under 10,000 target permutations) + pairwise ROR."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

from tobit_lib import PAYLOAD_MAP, fit_nb2_offset

ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'


def run_faers_analysis():
    F = pd.read_csv(f'{R}faers_counts.tsv', sep='\t')
    allT = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_trials.tsv', sep='\t').set_index('drug')
    W = pd.read_csv(f'{R}williams_kc_stats.tsv.gz', sep='\t', index_col=0)
    d = F[F.set == 'noco'].set_index('drug').join(allT[['target_gene', 'payload_class']])
    d['pclass'] = d.payload_class.map(PAYLOAD_MAP)
    d = d[d.n_reports >= 50].copy()
    d['W_basal'] = np.log2(W.reindex(d.target_gene.replace({'NECTIN4': 'PVRL4'}))['w_basal_mean'].values + 1.0)
    d['tg'] = d.target_gene

    X0 = np.column_stack([np.ones(len(d)), (d.pclass == 'topoI').astype(float).values, (d.pclass == 'DNA').astype(float).values])
    X1 = np.column_stack([X0, d.W_basal.values])
    Q0, _ = np.linalg.qr(X0)
    log_n = np.log(d.n_reports.values.astype(float))
    tg_s = d.groupby('tg').W_basal.first()
    tgs_f = tg_s.index.values
    w_vec = d.W_basal.values

    frows = []
    for out in ['oral', 'neutropenia', 'neuropathy', 'skin', 'ocular', 'ild', 'diarrhea']:
        k_arr = d['n_' + out].values.astype(float)
        rate = k_arr / d.n_reports.values.astype(float)
        rho_sp, p_sp = stats.spearmanr(w_vec, rate)

        # 1. Primary exposure-offset Negative-Binomial (NB2) count regression + B=10,000 12-target permutations
        r0 = fit_nb2_offset(k_arr, log_n, X0, compute_se=False)
        r1 = fit_nb2_offset(k_arr, log_n, X1, compute_se=True)
        b_nb = float(r1['beta'][-1])
        se_nb = float(r1['se_last'])
        p_wald = float(2.0 * stats.norm.sf(abs(b_nb) / se_nb))
        lr_nb = max(0.0, 2.0 * (r0['nll'] - r1['nll']))
        p_lr = float(stats.chi2.sf(lr_nb, 1))

        warm = np.r_[r0['beta'], 0.0, np.log(1.0 / max(r0['alpha'], 1e-6))]
        rng_nb = np.random.default_rng(1)
        null_lr = np.zeros(10000)
        for b_ix in range(10000):
            pm = dict(zip(tgs_f, rng_nb.permutation(tg_s.values)))
            Xb = np.column_stack([X0, d.tg.map(pm).values])
            rb = fit_nb2_offset(k_arr, log_n, Xb, compute_se=False, warm_theta=warm)
            null_lr[b_ix] = max(0.0, 2.0 * (r0['nll'] - rb['nll']))
        p_perm_target = float((1 + int((null_lr >= lr_nb).sum())) / 10001)

        # 2. Unweighted empirical log-odds OLS sensitivity + identical 10,000 12-target LR/F permutations
        y_ols = np.log((k_arr + 0.5) / (d.n_reports.values - k_arr + 0.5))
        d['y_ols'] = y_ols
        m_ols = smf.ols('y_ols ~ C(pclass) + W_basal', d).fit()
        b_ols = float(m_ols.params['W_basal'])
        se_ols = float(m_ols.bse['W_basal'])
        y_res = y_ols - Q0 @ (Q0.T @ y_ols)
        rss0 = float(y_res @ y_res)
        w_res0 = w_vec - Q0 @ (Q0.T @ w_vec)
        rss1 = float(rss0 - b_ols * b_ols * float(w_res0 @ w_res0))
        ols_lr = max(0.0, float(len(d) * np.log(rss0 / max(rss1, 1e-15))))

        rng_ols = np.random.default_rng(1)
        null_ols_lr = np.zeros(10000)
        for b_ix in range(10000):
            pm = dict(zip(tgs_f, rng_ols.permutation(tg_s.values)))
            wb = d.tg.map(pm).values
            wr = wb - Q0 @ (Q0.T @ wb)
            ww = float(wr @ wr)
            if ww < 1e-12:
                null_ols_lr[b_ix] = 0.0
                continue
            bb = float((wr @ y_res) / ww)
            rss_b = float(rss0 - bb * bb * ww)
            null_ols_lr[b_ix] = max(0.0, float(len(d) * np.log(rss0 / max(rss_b, 1e-15))))
        p_ols_perm = float((1 + int((null_ols_lr >= ols_lr - 1e-12).sum())) / 10001)

        frows.append(dict(
            outcome=out,
            n_adcs=len(d),
            n_targets=len(tgs_f),
            total_events=int(k_arr.sum()),
            beta=b_nb,
            se=se_nb,
            ci_lo=b_nb - 1.96 * se_nb,
            ci_hi=b_nb + 1.96 * se_nb,
            rate_ratio_per_log2=float(np.exp(b_nb)),
            p_wald=p_wald,
            LR=lr_nb,
            p_LR=p_lr,
            p_perm=p_perm_target,
            nb_alpha=float(r1['alpha']),
            spearman_rho=float(rho_sp),
            p_spearman=float(p_sp),
            ols_beta=b_ols,
            ols_se=se_ols,
            ols_LR=ols_lr,
            p_ols=float(m_ols.pvalues['W_basal']),
            p_ols_perm=p_ols_perm,
        ))

    faers_df = pd.DataFrame(frows)
    faers_df.to_csv(f'{R}fig_faers_summary.tsv', sep='\t', index=False)

    # Pairwise ROR contrasts (4 higher-vs-lower-basal + 1 same-target TACSTD2)
    def _calc_ror(k1, n1, k2, n2):
        c1, c2 = n1 - k1, n2 - k2
        r = (k1 / c1) / (k2 / c2)
        se = np.sqrt(1.0 / k1 + 1.0 / c1 + 1.0 / k2 + 1.0 / c2)
        return float(r), float(r * np.exp(-1.96 * se)), float(r * np.exp(1.96 * se))

    pair_specs = [
        ('datopotamab deruxtecan', 'trastuzumab deruxtecan', 'Dato-DXd vs T-DXd\n(both DXd)', 'Higher vs lower oral basal'),
        ('enfortumab vedotin', 'brentuximab vedotin', 'EV vs BV\n(both MMAE)', 'Higher vs lower oral basal'),
        ('enfortumab vedotin', 'polatuzumab vedotin', 'EV vs Pola\n(both MMAE)', 'Higher vs lower oral basal'),
        ('trastuzumab emtansine', 'mirvetuximab soravtansine', 'T-DM1 vs Mirv\n(DM1 vs DM4)', 'Higher vs lower oral basal'),
        ('datopotamab deruxtecan', 'sacituzumab govitecan', 'Dato-DXd vs SG [same target]\n(DXd vs SN-38, TACSTD2)', 'Same target (TACSTD2)'),
    ]
    prows = []
    for d1, d2, label, ptype in pair_specs:
        k1, n1 = int(d.loc[d1, 'n_oral']), int(d.loc[d1, 'n_reports'])
        k2, n2 = int(d.loc[d2, 'n_oral']), int(d.loc[d2, 'n_reports'])
        r_val, lo_val, hi_val = _calc_ror(k1, n1, k2, n2)
        prows.append(dict(
            drug1=d1, drug2=d2, contrast=label, pair_type=ptype,
            basal1=float(d.loc[d1, 'W_basal']), basal2=float(d.loc[d2, 'W_basal']),
            k1=k1, n1=n1, pct1=round(100.0 * k1 / n1, 2),
            k2=k2, n2=n2, pct2=round(100.0 * k2 / n2, 2),
            ror=round(r_val, 2), lo=round(lo_val, 2), hi=round(hi_val, 2),
        ))
    pairs_df = pd.DataFrame(prows)
    pairs_df.to_csv(f'{R}fig_faers_pairs.tsv', sep='\t', index=False)
    return faers_df, pairs_df


if __name__ == '__main__':
    fdf, pdf = run_faers_analysis()
    print(fdf[['outcome', 'beta', 'se', 'p_wald', 'LR', 'p_LR', 'p_perm', 'spearman_rho', 'ols_beta', 'ols_LR', 'p_ols', 'p_ols_perm']].round(4).to_string(index=False))
    print('\nPairwise oral ROR contrasts:')
    print(pdf[['drug1', 'drug2', 'ror', 'lo', 'hi']].round(2).to_string(index=False))
