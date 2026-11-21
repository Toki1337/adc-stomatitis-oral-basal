"""Part 5: Prespecified LOTO calibration (zero-antigen-overlap folds), exact uncensored endpoints, & LOO-ADC."""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from scipy.special import expit, logit

from tobit_lib import design_matrix, exact_mcnemar_twosided, fit_tobit

pd.set_option('display.width', 260)
ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'
T = pd.read_csv(R + 'part3_axes.tsv', sep='\t')
W = pd.read_csv(R + 'williams_kc_stats.tsv.gz', sep='\t', index_col=0)
d = T.dropna(subset=['abs_oral_basal']).reset_index(drop=True)
ALIAS = {'NECTIN4': 'PVRL4', 'CD44v6': 'CD44'}
row_gene_sets = [frozenset(ALIAS.get(g, g) for g in str(gs).split(';')) for gs in d.target_gene]


def predict_split(train_df, test_df, x_col=None, p4=False):
    Xm_tr = design_matrix(train_df, x_col, payload=True, p4=p4)
    fit = fit_tobit(train_df.y_pct.values, train_df.cens.values, Xm_tr)
    cols_te = [np.ones(len(test_df), dtype=float)]
    if (train_df.pclass == 'tubulin').any():
        cols_te.append((test_df.pclass == 'tubulin').astype(float).values)
    if (train_df.pclass == 'DNA').any():
        cols_te.append((test_df.pclass == 'DNA').astype(float).values)
    if p4 and (train_df.payload_class == 'SN-38').any():
        cols_te.append((test_df.payload_class == 'SN-38').astype(float).values)
    if x_col:
        cols_te.append(test_df[x_col].astype(float).values)
    Xt = np.column_stack(cols_te)
    lp = Xt @ fit['beta']
    return expit(lp) * 100.0, lp, float(fit['sigma'])


def row_log_score(y_pct, is_cens, lp, sigma):
    yv = logit(np.clip(float(y_pct) / 100.0, 0.001, 0.999))
    if is_cens:
        return float(stats.norm.logcdf((yv - float(lp)) / float(sigma)))
    return float(stats.norm.logpdf(yv, loc=float(lp), scale=float(sigma)))


def left_censored_c_index(y, cens, score):
    """Exact concordance over all unambiguously ordered pairs without distributional assumptions."""
    n = len(y)
    conc = 0.0
    total = 0
    for i in range(n):
        if cens[i]:
            continue
        for j in range(n):
            if i == j:
                continue
            if ((not cens[j]) and (y[i] > y[j])) or (cens[j] and (y[i] >= y[j])):
                total += 1
                if score[i] > score[j]:
                    conc += 1.0
                elif score[i] == score[j]:
                    conc += 0.5
    return (conc / total if total > 0 else np.nan), total


def exact_roc_pr_auc(y_bin, scores):
    """Exact Mann-Whitney ROC-AUC (with midrank ties), two-sided P, and step-function PR-AUC."""
    y_b = np.asarray(y_bin, dtype=int)
    sc = np.asarray(scores, dtype=float)
    pos = sc[y_b == 1]
    neg = sc[y_b == 0]
    u_res = stats.mannwhitneyu(pos, neg, alternative='two-sided')
    auc = float(u_res.statistic / (len(pos) * len(neg)))
    p_mw = float(u_res.pvalue)
    # Exact step-function average precision (matching sklearn average_precision_score)
    uniq = np.sort(np.unique(sc))[::-1]
    ap = 0.0
    prev_rec = 0.0
    n_pos = float(len(pos))
    for thr in uniq:
        m = sc >= thr
        tp = float(y_b[m].sum())
        prec = tp / float(m.sum())
        rec = tp / n_pos
        ap += (rec - prev_rec) * prec
        prev_rec = rec
    return auc, p_mw, float(ap)


# 1. Primary LOTO (zero-antigen-overlap folds: excluding all ADCs sharing any antigen with the held-out ADC)
p_loto, p0_loto, p4_loto = np.zeros(len(d)), np.zeros(len(d)), np.zeros(len(d))
ll_loto, ll0_loto, sig_loto, sig0_loto = np.zeros(len(d)), np.zeros(len(d)), np.zeros(len(d)), np.zeros(len(d))
for i in range(len(d)):
    overlap = np.array([len(row_gene_sets[i] & row_gene_sets[j]) > 0 for j in range(len(d))])
    te = d.iloc[[i]]
    tr = d[~overlap]
    px, lpx, sx = predict_split(tr, te, 'abs_oral_basal', p4=False)
    p0, lp0, s0 = predict_split(tr, te, None, p4=False)
    p4v, _, _ = predict_split(tr, te, 'abs_oral_basal', p4=True)
    p_loto[i], p0_loto[i], p4_loto[i] = px[0], p0[0], p4v[0]
    sig_loto[i], sig0_loto[i] = sx, s0
    ll_loto[i] = row_log_score(te.y_pct.values[0], bool(te.cens.values[0]), lpx[0], sx)
    ll0_loto[i] = row_log_score(te.y_pct.values[0], bool(te.cens.values[0]), lp0[0], s0)

d['pred_loto'] = p_loto
d['pred_loto_payload_only'] = p0_loto
d['pred_loto_4payload'] = p4_loto
d['sigma_loto'] = sig_loto
d['sigma_loto_payload'] = sig0_loto
d['log_score_loto'] = ll_loto
d['log_score_loto_payload'] = ll0_loto

# 2. Secondary Leave-One-ADC-Out (44 folds, supplement reference)
p_loo, p0_loo, p4_loo = np.zeros(len(d)), np.zeros(len(d)), np.zeros(len(d))
ll_loo, ll0_loo = np.zeros(len(d)), np.zeros(len(d))
for i in range(len(d)):
    tr = d.drop(i)
    te = d.loc[[i]]
    px, lpx, sx = predict_split(tr, te, 'abs_oral_basal', p4=False)
    p0, lp0, s0 = predict_split(tr, te, None, p4=False)
    p4v, _, _ = predict_split(tr, te, 'abs_oral_basal', p4=True)
    p_loo[i], p0_loo[i], p4_loo[i] = px[0], p0[0], p4v[0]
    ll_loo[i] = row_log_score(te.y_pct.values[0], bool(te.cens.values[0]), lpx[0], sx)
    ll0_loo[i] = row_log_score(te.y_pct.values[0], bool(te.cens.values[0]), lp0[0], s0)

d['pred_loo'] = p_loo
d['pred_payload_only'] = p0_loo
d['pred_loo_4payload'] = p4_loo
d['log_score_loo'] = ll_loo
d['log_score_loo_payload'] = ll0_loo

# 3. Exact 28-target LOTO refit permutation test (B=5,000) for distribution-free C-index & exact binary ROC-AUCs
all_genes = sorted({ALIAS.get(g, g) for gs in d.target_gene for g in str(gs).split(';')})
gene_x = np.array([float(np.log2(W.loc[g, 'w_basal_mean'] + 1.0)) for g in all_genes])

# Group identical antigen sets so each distinct LOTO training fold is fitted once per permutation
uniq_gsets = sorted(set(row_gene_sets), key=lambda s: sorted(s))
fold_specs = []
for gs_u in uniq_gsets:
    te_idx = np.where([s == gs_u for s in row_gene_sets])[0]
    tr_idx = np.where([len(s & gs_u) == 0 for s in row_gene_sets])[0]
    fold_specs.append((tr_idx, te_idx))

X_pay_all = design_matrix(d)
y_v = d.y_pct.values
c_v = d.cens.values

c_loto_obs, n_comp_pairs = left_censored_c_index(y_v, c_v, p_loto)
c_loto_0, _ = left_censored_c_index(y_v, c_v, p0_loto)
c_loto_4, _ = left_censored_c_index(y_v, c_v, p4_loto)

m10 = (~d.cens.values) | (d.y_pct.values <= 10.0)
y10 = ((~d.cens.values) & (d.y_pct.values >= 10.0))[m10].astype(int)
auc10_x, pmw10_x, ap10_x = exact_roc_pr_auc(y10, p_loto[m10])
auc10_0, pmw10_0, ap10_0 = exact_roc_pr_auc(y10, p0_loto[m10])

m15 = (~d.cens.values) | (d.y_pct.values <= 15.0)
y15 = ((~d.cens.values) & (d.y_pct.values >= 15.0))[m15].astype(int)
auc15_x, pmw15_x, ap15_x = exact_roc_pr_auc(y15, p_loto[m15])
auc15_0, pmw15_0, ap15_0 = exact_roc_pr_auc(y15, p0_loto[m15])

m20 = (~d.cens.values) | (d.y_pct.values <= 20.0)
y20 = ((~d.cens.values) & (d.y_pct.values >= 20.0))[m20].astype(int)
auc20_x, pmw20_x, ap20_x = exact_roc_pr_auc(y20, p_loto[m20])
auc20_0, pmw20_0, ap20_0 = exact_roc_pr_auc(y20, p0_loto[m20])

# Precompute fast vectorized pair indices for C-index inside permutation loop
pair_i, pair_j = [], []
for ii in range(len(y_v)):
    if c_v[ii]:
        continue
    for jj in range(len(y_v)):
        if ii != jj and (((not c_v[jj]) and (y_v[ii] > y_v[jj])) or (c_v[jj] and (y_v[ii] >= y_v[jj]))):
            pair_i.append(ii)
            pair_j.append(jj)
pair_i = np.array(pair_i, dtype=int)
pair_j = np.array(pair_j, dtype=int)

B_PERM = 5000
rng = np.random.default_rng(1)
c_null = np.zeros(B_PERM)
a10_null = np.zeros(B_PERM)
a15_null = np.zeros(B_PERM)
a20_null = np.zeros(B_PERM)

for b in range(B_PERM):
    perm_gx = dict(zip(all_genes, rng.permutation(gene_x)))
    xp = np.array([max(perm_gx[g] for g in gs_set) for gs_set in row_gene_sets])
    X_all = np.column_stack([X_pay_all, xp])
    preds_b = np.zeros(len(d))
    for tr_idx, te_idx in fold_specs:
        r_fit = fit_tobit(y_v[tr_idx], c_v[tr_idx], X_all[tr_idx])
        preds_b[te_idx] = expit(X_all[te_idx] @ r_fit['beta']) * 100.0
    diff_p = preds_b[pair_i] - preds_b[pair_j]
    c_null[b] = (np.sum(diff_p > 0) + 0.5 * np.sum(diff_p == 0)) / len(pair_i)
    a10_null[b] = float(stats.mannwhitneyu(preds_b[m10][y10 == 1], preds_b[m10][y10 == 0]).statistic / (y10.sum() * (len(y10) - y10.sum())))
    a15_null[b] = float(stats.mannwhitneyu(preds_b[m15][y15 == 1], preds_b[m15][y15 == 0]).statistic / (y15.sum() * (len(y15) - y15.sum())))
    a20_null[b] = float(stats.mannwhitneyu(preds_b[m20][y20 == 1], preds_b[m20][y20 == 0]).statistic / (y20.sum() * (len(y20) - y20.sum())))

p_perm_c = float((1 + int((c_null >= c_loto_obs).sum())) / (B_PERM + 1))
p_perm_a10 = float((1 + int((a10_null >= auc10_x).sum())) / (B_PERM + 1))
p_perm_a15 = float((1 + int((a15_null >= auc15_x).sum())) / (B_PERM + 1))
p_perm_a20 = float((1 + int((a20_null >= auc20_x).sum())) / (B_PERM + 1))

exact_ep_df = pd.DataFrame([
    dict(
        endpoint='All-44-ADC Left-Censored Concordance (C-index)',
        censoring_assumption='None (451 unambiguously ordered pairs)',
        n_evaluable=len(d),
        n_targets=int(d.tg.nunique()),
        n_pos_or_pairs=int(n_comp_pairs),
        n_neg=np.nan,
        loto_metric_x=c_loto_obs,
        loto_metric_payload_only=c_loto_0,
        delta_vs_payload=c_loto_obs - c_loto_0,
        loto_metric_4payload=c_loto_4,
        p_mannwhitney_x=np.nan,
        p_mannwhitney_payload=np.nan,
        p_target_perm_5000=p_perm_c,
        loto_prauc_x=np.nan,
        loto_prauc_payload=np.nan,
    ),
    dict(
        endpoint='Exact Uncensored Stomatitis >= 10% (ROC-AUC)',
        censoring_assumption='None (excludes 9 censored ADCs with >10% table cutoff)',
        n_evaluable=int(m10.sum()),
        n_targets=int(d.loc[m10, 'tg'].nunique()),
        n_pos_or_pairs=int(y10.sum()),
        n_neg=int((1 - y10).sum()),
        loto_metric_x=auc10_x,
        loto_metric_payload_only=auc10_0,
        delta_vs_payload=auc10_x - auc10_0,
        loto_metric_4payload=exact_roc_pr_auc(y10, p4_loto[m10])[0],
        p_mannwhitney_x=pmw10_x,
        p_mannwhitney_payload=pmw10_0,
        p_target_perm_5000=p_perm_a10,
        loto_prauc_x=ap10_x,
        loto_prauc_payload=ap10_0,
    ),
    dict(
        endpoint='Exact Uncensored Stomatitis >= 15% (ROC-AUC)',
        censoring_assumption='None (excludes 6 censored ADCs with 20% table cutoff)',
        n_evaluable=int(m15.sum()),
        n_targets=int(d.loc[m15, 'tg'].nunique()),
        n_pos_or_pairs=int(y15.sum()),
        n_neg=int((1 - y15).sum()),
        loto_metric_x=auc15_x,
        loto_metric_payload_only=auc15_0,
        delta_vs_payload=auc15_x - auc15_0,
        loto_metric_4payload=exact_roc_pr_auc(y15, p4_loto[m15])[0],
        p_mannwhitney_x=pmw15_x,
        p_mannwhitney_payload=pmw15_0,
        p_target_perm_5000=p_perm_a15,
        loto_prauc_x=ap15_x,
        loto_prauc_payload=ap15_0,
    ),
    dict(
        endpoint='Exact Uncensored Stomatitis >= 20% (ROC-AUC)',
        censoring_assumption='None (all 44/44 ADCs have exact binary status)',
        n_evaluable=int(m20.sum()),
        n_targets=int(d.loc[m20, 'tg'].nunique()),
        n_pos_or_pairs=int(y20.sum()),
        n_neg=int((1 - y20).sum()),
        loto_metric_x=auc20_x,
        loto_metric_payload_only=auc20_0,
        delta_vs_payload=auc20_x - auc20_0,
        loto_metric_4payload=exact_roc_pr_auc(y20, p4_loto[m20])[0],
        p_mannwhitney_x=pmw20_x,
        p_mannwhitney_payload=pmw20_0,
        p_target_perm_5000=p_perm_a20,
        loto_prauc_x=ap20_x,
        loto_prauc_payload=ap20_0,
    ),
])
exact_ep_df.to_csv(R + 'part5_loto_exact_endpoints.tsv', sep='\t', index=False)
print('\n=== DISTRIBUTION-FREE / EXACT UNCENSORED LOTO ENDPOINTS (B=5,000 target permutations) ===')
print(exact_ep_df.round(4).to_string(index=False))


def summarize_scheme(df_in, col_x, col_0, col_4, ll_col_x, ll_col_0, label, is_loto=False):
    num = df_in[~df_in.cens]
    cen = df_in[df_in.cens]
    rho_x, p_rho_x = stats.spearmanr(num[col_x], num.y_pct)
    rho_0, p_rho_0 = stats.spearmanr(num[col_0], num.y_pct)
    rho_4, _ = stats.spearmanr(num[col_4], num.y_pct)
    mae_x = float(np.abs(num[col_x] - num.y_pct).mean())
    mae_0 = float(np.abs(num[col_0] - num.y_pct).mean())
    mae_4 = float(np.abs(num[col_4] - num.y_pct).mean())

    ok_x = (cen[col_x] < cen.y_pct).values
    ok_0 = (cen[col_0] < cen.y_pct).values
    ok_4 = (cen[col_4] < cen.y_pct).values
    disc_x_only = int((ok_x & (~ok_0)).sum())
    disc_0_only = int(((~ok_x) & ok_0).sum())
    p_mcnemar_2s = exact_mcnemar_twosided(disc_x_only, disc_0_only)

    ll_all_x = float(df_in[ll_col_x].mean())
    ll_all_0 = float(df_in[ll_col_0].mean())
    ll_num_x = float(num[ll_col_x].mean())
    ll_num_0 = float(num[ll_col_0].mean())
    ll_cen_x = float(cen[ll_col_x].mean())
    ll_cen_0 = float(cen[ll_col_0].mean())

    cx, n_pairs = left_censored_c_index(df_in.y_pct.values, df_in.cens.values, df_in[col_x].values)
    c0, _ = left_censored_c_index(df_in.y_pct.values, df_in.cens.values, df_in[col_0].values)
    c4, _ = left_censored_c_index(df_in.y_pct.values, df_in.cens.values, df_in[col_4].values)
    a15_s, pmw15_s, ap15_s = exact_roc_pr_auc(y15, df_in[col_x].values[m15])
    a15_0_s, _, ap15_0_s = exact_roc_pr_auc(y15, df_in[col_0].values[m15])
    a10_s, pmw10_s, ap10_s = exact_roc_pr_auc(y10, df_in[col_x].values[m10])
    a10_0_s, _, _ = exact_roc_pr_auc(y10, df_in[col_0].values[m10])

    print(f'\n=== {label} ===')
    print(f'All-44 C-index ({n_pairs} pairs): x={cx:.3f} vs payload={c0:.3f} (4-pay={c4:.3f}) | Exact >=15% AUC (n={int(m15.sum())}): {a15_s:.3f} vs {a15_0_s:.3f} (p={pmw15_s:.4f})')
    print(f'Numeric ADCs (n={len(num)}): Spearman x={rho_x:.3f} (p={p_rho_x:.4f}), MAE={mae_x:.1f} pp | payload-only={rho_0:.3f} (p={p_rho_0:.4f}), MAE={mae_0:.1f} pp | 4-payload={rho_4:.3f}, MAE={mae_4:.1f} pp')
    print(f'Censored ADCs (n={len(cen)}): x={int(ok_x.sum())}/{len(cen)} ({100*ok_x.mean():.0f}%) | payload-only={int(ok_0.sum())}/{len(cen)} ({100*ok_0.mean():.0f}%) | 4-payload={int(ok_4.sum())}/{len(cen)} ({100*ok_4.mean():.0f}%)')
    print(f'Log-score (all/num/cens): x=({ll_all_x:.3f}, {ll_num_x:.3f}, {ll_cen_x:.3f}) vs payload=({ll_all_0:.3f}, {ll_num_0:.3f}, {ll_cen_0:.3f})')
    print(f'Discordant (x vs payload-only): {disc_x_only} vs {disc_0_only} | exact two-sided McNemar p={p_mcnemar_2s:.4f}')
    return dict(
        scheme=label, n_numeric=len(num), n_censored=len(cen),
        c_index_all44_x=float(cx), c_index_all44_payload=float(c0), c_index_all44_4pay=float(c4),
        n_comparable_pairs=int(n_pairs),
        p_perm_c_index=(p_perm_c if is_loto else np.nan),
        auc_ge15_x=float(a15_s), auc_ge15_payload=float(a15_0_s),
        prauc_ge15_x=float(ap15_s), prauc_ge15_payload=float(ap15_0_s),
        p_mw_auc_ge15=float(pmw15_s), p_perm_auc_ge15=(p_perm_a15 if is_loto else np.nan),
        auc_ge10_x=float(a10_s), auc_ge10_payload=float(a10_0_s),
        p_mw_auc_ge10=float(pmw10_s), p_perm_auc_ge10=(p_perm_a10 if is_loto else np.nan),
        spearman_x=float(rho_x), p_spearman_x=float(p_rho_x), mae_x=mae_x,
        spearman_payload=float(rho_0), p_spearman_payload=float(p_rho_0), mae_payload=mae_0,
        spearman_4pay=float(rho_4), mae_4pay=mae_4,
        cens_below_x=int(ok_x.sum()), cens_below_payload=int(ok_0.sum()), cens_below_4pay=int(ok_4.sum()),
        disc_x_only=disc_x_only, disc_payload_only=disc_0_only,
        p_mcnemar_twosided=p_mcnemar_2s,
        mean_logscore_all_x=ll_all_x, mean_logscore_all_payload=ll_all_0,
        mean_logscore_num_x=ll_num_x, mean_logscore_num_payload=ll_num_0,
        mean_logscore_cens_x=ll_cen_x, mean_logscore_cens_payload=ll_cen_0,
    )


s_loto = summarize_scheme(d, 'pred_loto', 'pred_loto_payload_only', 'pred_loto_4payload', 'log_score_loto', 'log_score_loto_payload', 'Leave-One-Target-Out (LOTO, prespecified; zero-antigen-overlap folds)', is_loto=True)
s_loo = summarize_scheme(d, 'pred_loo', 'pred_payload_only', 'pred_loo_4payload', 'log_score_loo', 'log_score_loo_payload', 'Leave-One-ADC-Out (LOO-ADC, 44 folds; exploratory)', is_loto=False)

pd.DataFrame([s_loto, s_loo]).to_csv(R + 'part5_calib_summary.tsv', sep='\t', index=False)
d.to_csv(R + 'part5_loo_predictions.tsv', sep='\t', index=False)
d.to_csv(R + 'part5_loto_predictions.tsv', sep='\t', index=False)
