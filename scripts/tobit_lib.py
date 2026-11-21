"""Shared left-censored Tobit regression, inference, bootstrap, and FDR utilities."""

import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.special import digamma, expit, gammaln, log_ndtr, logit

LOG_SQRT_2PI = 0.5 * np.log(2.0 * np.pi)
GAM_MIN, GAM_MAX = -6.0, 4.0

PAYLOAD_MAP = {
    'DXd': 'topoI', 'SN-38': 'topoI', 'topo1-other': 'topoI',
    'MMAE': 'tubulin', 'MMAF': 'tubulin', 'DM1': 'tubulin', 'DM4': 'tubulin',
    'auristatin-other': 'tubulin', 'eribulin': 'tubulin', 'hemiasterlin': 'tubulin',
    'PBD': 'DNA', 'calicheamicin': 'DNA', 'duocarmycin': 'DNA', 'IGN': 'DNA',
}


def design_matrix(df, x_cols=None, payload=True, p4=False):
    """Build design matrix with intercept, active payload indicators, and predictors."""
    cols = [np.ones(len(df), dtype=float)]
    if payload:
        if (df.pclass == 'tubulin').any():
            cols.append((df.pclass == 'tubulin').astype(float).values)
        if (df.pclass == 'DNA').any():
            cols.append((df.pclass == 'DNA').astype(float).values)
        if p4 and ('payload_class' in df.columns) and (df.payload_class == 'SN-38').any():
            cols.append((df.payload_class == 'SN-38').astype(float).values)
    if x_cols:
        if isinstance(x_cols, str):
            x_cols = [x_cols]
        for c in x_cols:
            cols.append(df[c].astype(float).values)
    return np.column_stack(cols)


def fit_tobit(y_pct, cens, Xm, weights=None, gtol=1e-6):
    """Fit left-censored Gaussian Tobit on logit-transformed incidence percentages."""
    yl = logit(np.clip(np.asarray(y_pct, dtype=float) / 100.0, 1e-3, 0.999))
    c = np.asarray(cens, dtype=bool)
    w = np.ones(len(yl), dtype=float) if weights is None else np.asarray(weights, dtype=float)

    Xw = Xm.copy().astype(float)
    means = np.zeros(Xm.shape[1])
    stds = np.ones(Xm.shape[1])
    for j in range(1, Xm.shape[1]):
        means[j] = Xm[:, j].mean()
        sd = Xm[:, j].std()
        stds[j] = sd if sd > 1e-12 else 1.0
        Xw[:, j] = (Xm[:, j] - means[j]) / stds[j]

    def nll_grad(p):
        b = p[:-1]
        gam = p[-1]
        if gam < GAM_MIN:
            gam_c = GAM_MIN
            pen = 100.0 * (gam - GAM_MIN) ** 2
            gpen = 200.0 * (gam - GAM_MIN)
        elif gam > GAM_MAX:
            gam_c = GAM_MAX
            pen = 100.0 * (gam - GAM_MAX) ** 2
            gpen = 200.0 * (gam - GAM_MAX)
        else:
            gam_c = gam
            pen = 0.0
            gpen = 0.0

        s = np.exp(gam_c)
        u = (yl - (Xw @ b)) / s
        lc = log_ndtr(u)
        lp = -0.5 * u * u - LOG_SQRT_2PI
        ll = w * np.where(c, lc, lp - gam_c)
        lam = np.exp(np.clip(lp - lc, -500.0, 500.0))
        g_mu = w * np.where(c, lam / s, -u / s)
        g_gam = (w * np.where(c, lam * u, 1.0 - u * u)).sum() if (GAM_MIN <= gam <= GAM_MAX) else gpen
        return -ll.sum() + pen, np.r_[Xw.T @ g_mu, g_gam]

    # Weighted least-squares starting point
    sw = np.sqrt(w)
    b0 = np.linalg.lstsq(Xw * sw[:, None], yl * sw, rcond=None)[0]
    p0 = np.r_[b0, 0.0]
    r = optimize.minimize(nll_grad, p0, jac=True, method='BFGS', options={'gtol': gtol, 'maxiter': 200})
    if not r.success:
        r2 = optimize.minimize(lambda p: nll_grad(p)[0], r.x, method='Nelder-Mead', options={'maxiter': 400})
        r = optimize.minimize(nll_grad, r2.x, jac=True, method='BFGS', options={'gtol': 1e-4, 'maxiter': 200})

    b_orig = np.zeros(Xm.shape[1])
    for j in range(1, Xm.shape[1]):
        b_orig[j] = r.x[j] / stds[j]
    b_orig[0] = r.x[0] - np.sum(b_orig[1:] * means[1:])
    sigma = float(np.exp(np.clip(r.x[-1], GAM_MIN, GAM_MAX)))

    # Return worked-scale solution and unscaled Hessian helper info too
    return dict(
        nll=float(r.fun),
        beta=b_orig,
        sigma=sigma,
        gam=float(np.clip(r.x[-1], GAM_MIN, GAM_MAX)),
        p_work=r.x,
        stds=stds,
        Xw=Xw,
        yl=yl,
        c=c,
        w=w,
        converged=bool(r.success),
    )


def exact_hess_se(fit_res, eps=1e-5):
    """Central-difference Hessian standard errors transformed back to original covariate scale."""
    Xw = fit_res['Xw']
    yl = fit_res['yl']
    c = fit_res['c']
    w = fit_res['w']
    p_opt = fit_res['p_work']
    stds = fit_res['stds']

    def nll_w(p):
        b = p[:-1]
        gam = np.clip(p[-1], GAM_MIN, GAM_MAX)
        s = np.exp(gam)
        u = (yl - (Xw @ b)) / s
        lc = log_ndtr(u)
        lp = -0.5 * u * u - LOG_SQRT_2PI
        return -float((w * np.where(c, lc, lp - gam)).sum())

    k = len(p_opt)
    H = np.zeros((k, k))
    for i in range(k):
        for j in range(i, k):
            ei = np.zeros(k); ei[i] = eps
            ej = np.zeros(k); ej[j] = eps
            H[i, j] = H[j, i] = (
                nll_w(p_opt + ei + ej) - nll_w(p_opt + ei - ej)
                - nll_w(p_opt - ei + ej) + nll_w(p_opt - ei - ej)
            ) / (4.0 * eps * eps)
    try:
        cov = np.linalg.inv(H)
        se_w = np.sqrt(np.maximum(np.diag(cov), 0.0))[:-1]
        return se_w / stds
    except np.linalg.LinAlgError:
        return np.full(len(stds), np.nan)


def profile_ci(df, x_col, payload=True, p4=False, weights=None, level=0.95):
    """Profile-likelihood confidence interval for the coefficient of x_col."""
    X1 = design_matrix(df, x_col, payload=payload, p4=p4)
    y = df.y_pct.values
    cens = df.cens.values
    res1 = fit_tobit(y, cens, X1, weights=weights)
    f1 = res1['nll']
    bhat = float(res1['beta'][-1])
    crit = stats.chi2.ppf(level, 1)

    X_other = X1[:, :-1]
    xv = X1[:, -1]
    yl = res1['yl']

    def delta_dev(bval):
        # Adjust outcome by fixed offset bval * xv and refit nuisance parameters
        y_adj_pct = expit(yl - bval * xv) * 100.0
        # Exact shifted y on logit scale: reuse direct nll optimization on offset
        return _fit_with_offset(yl - bval * xv, cens, X_other, weights, res1['gam']) - f1

    def obj(bval):
        return 2.0 * delta_dev(bval) - crit

    def find_bound(direction):
        step = 0.25 * max(0.1, abs(bhat))
        x_in = bhat
        x_out = bhat + direction * step
        for _ in range(18):
            val = obj(x_out)
            if np.isfinite(val) and val > 0:
                return float(optimize.brentq(obj, min(x_in, x_out), max(x_in, x_out), xtol=1e-4))
            x_in = x_out
            x_out += direction * step * (1.5 ** _)
            if abs(x_out - bhat) > 25.0:
                break
        return np.nan

    lo = find_bound(-1.0)
    hi = find_bound(+1.0)
    return bhat, lo, hi


def _fit_with_offset(yl_adj, cens, Xm, weights, gam_init=0.0):
    c = np.asarray(cens, dtype=bool)
    w = np.ones(len(yl_adj), dtype=float) if weights is None else np.asarray(weights, dtype=float)
    Xw = Xm.copy().astype(float)
    for j in range(1, Xm.shape[1]):
        sd = Xm[:, j].std()
        Xw[:, j] = (Xm[:, j] - Xm[:, j].mean()) / (sd if sd > 1e-12 else 1.0)

    def nll_grad(p):
        b = p[:-1]
        gam = p[-1]
        if gam < GAM_MIN:
            gam_c = GAM_MIN; pen = 100.0 * (gam - GAM_MIN) ** 2; gpen = 200.0 * (gam - GAM_MIN)
        elif gam > GAM_MAX:
            gam_c = GAM_MAX; pen = 100.0 * (gam - GAM_MAX) ** 2; gpen = 200.0 * (gam - GAM_MAX)
        else:
            gam_c = gam; pen = 0.0; gpen = 0.0
        s = np.exp(gam_c)
        u = (yl_adj - (Xw @ b)) / s
        lc = log_ndtr(u)
        lp = -0.5 * u * u - LOG_SQRT_2PI
        ll = w * np.where(c, lc, lp - gam_c)
        lam = np.exp(np.clip(lp - lc, -500.0, 500.0))
        g_mu = w * np.where(c, lam / s, -u / s)
        g_gam = (w * np.where(c, lam * u, 1.0 - u * u)).sum() if (GAM_MIN <= gam <= GAM_MAX) else gpen
        return -ll.sum() + pen, np.r_[Xw.T @ g_mu, g_gam]

    sw = np.sqrt(w)
    b0 = np.linalg.lstsq(Xw * sw[:, None], yl_adj * sw, rcond=None)[0]
    r = optimize.minimize(nll_grad, np.r_[b0, gam_init], jac=True, method='BFGS', options={'gtol': 1e-6, 'maxiter': 500})
    return float(r.fun)


def lr_test(df, x_col, base_cols=None, payload=True, p4=False, weights=None, compute_se=False):
    """Compare nested Tobit models adding x_col over base_cols + payload."""
    X0 = design_matrix(df, base_cols, payload=payload, p4=p4)
    X1 = design_matrix(df, ([base_cols] if isinstance(base_cols, str) else list(base_cols or [])) + [x_col], payload=payload, p4=p4)
    y = df.y_pct.values
    cens = df.cens.values
    r0 = fit_tobit(y, cens, X0, weights=weights)
    r1 = fit_tobit(y, cens, X1, weights=weights)
    lr = max(0.0, 2.0 * (r0['nll'] - r1['nll']))
    p = float(stats.chi2.sf(lr, 1))
    k0 = X0.shape[1] + 1
    k1 = X1.shape[1] + 1
    aic0 = 2.0 * k0 + 2.0 * r0['nll']
    aic1 = 2.0 * k1 + 2.0 * r1['nll']
    se_x = float(exact_hess_se(r1)[-1]) if compute_se else np.nan
    return dict(
        LR=lr,
        p=p,
        beta=float(r1['beta'][-1]),
        se=se_x,
        sigma=r1['sigma'],
        aic=aic1,
        daic=aic1 - aic0,
        nll0=r0['nll'],
        nll1=r1['nll'],
        beta_full=r1['beta'],
        converged=bool(r0['converged'] and r1['converged']),
    )


def exact_mcnemar_twosided(a_only_correct, b_only_correct):
    """Exact two-sided binomial McNemar test on discordant binary pairs."""
    n_disc = int(a_only_correct) + int(b_only_correct)
    if n_disc == 0:
        return 1.0
    k = min(int(a_only_correct), int(b_only_correct))
    return float(stats.binomtest(k, n_disc, 0.5, alternative='two-sided').pvalue)


def bh_fdr(pvals):
    """Benjamini-Hochberg false discovery rate adjustment."""
    return stats.false_discovery_control(np.asarray(pvals, dtype=float), method='bh')


def fit_nb2_offset(k, log_n, X, compute_se=False, warm_theta=None):
    """Fit exposure-offset Negative-Binomial (NB2) GLM via analytic-gradient L-BFGS-B."""
    k_arr = np.asarray(k, dtype=float)
    ln_arr = np.asarray(log_n, dtype=float)
    Xm = np.asarray(X, dtype=float)
    p = Xm.shape[1]

    def nll_and_grad(theta):
        b = theta[:p]
        log_r = float(np.clip(theta[p], -12.0, 12.0))
        r = np.exp(log_r)
        eta = np.clip(ln_arr + Xm @ b, -20.0, 20.0)
        mu = np.exp(eta)
        ll = (
            gammaln(k_arr + r) - gammaln(r) - gammaln(k_arr + 1.0)
            + r * np.log(r / (r + mu))
            + np.where(k_arr > 0, k_arr * (eta - np.log(r + mu)), 0.0)
        )
        resid = r * (k_arr - mu) / (r + mu)
        g_b = -(Xm.T @ resid)
        dl_dr = digamma(k_arr + r) - digamma(r) + np.log(r / (r + mu)) + (mu - k_arr) / (r + mu)
        g_lr = -float(np.sum(dl_dr) * r)
        return -float(np.sum(ll)), np.append(g_b, g_lr)

    bounds = [(-15.0, 15.0)] * p + [(-10.0, 10.0)]
    if warm_theta is not None:
        best = optimize.minimize(nll_and_grad, warm_theta, jac=True, method='L-BFGS-B', bounds=bounds)
    else:
        p0 = np.zeros(p + 1)
        p0[0] = np.log(max(k_arr.sum(), 0.5) / np.exp(ln_arr).sum())
        best = None
        for init_lr in (0.0, -1.5, 2.0):
            st = p0.copy()
            st[p] = init_lr
            res = optimize.minimize(nll_and_grad, st, jac=True, method='L-BFGS-B', bounds=bounds)
            if best is None or res.fun < best.fun:
                best = res

    se_last = np.nan
    if compute_se:
        eps = 1e-5
        th = best.x
        H = np.zeros((p + 1, p + 1))
        for j in range(p + 1):
            ep = np.zeros(p + 1)
            ep[j] = eps
            _, gp = nll_and_grad(th + ep)
            _, gm = nll_and_grad(th - ep)
            H[:, j] = (gp - gm) / (2.0 * eps)
        H = 0.5 * (H + H.T)
        try:
            cov = np.linalg.inv(H)
            se_last = float(np.sqrt(max(cov[p - 1, p - 1], 0.0)))
        except np.linalg.LinAlgError:
            se_last = np.nan

    return dict(
        nll=float(best.fun),
        beta=best.x[:p],
        alpha=float(np.exp(-np.clip(best.x[p], -12.0, 12.0))),
        se_last=se_last,
        theta=best.x,
        converged=bool(best.success),
    )
