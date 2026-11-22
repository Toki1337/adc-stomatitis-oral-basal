"""Part 1 FDA FAERS Quarterly ASCII (2012Q4-2026Q3) pharmacovigilance replication:
deduplicated sole-suspect monotherapy exposure-offset Negative-Binomial GLM +
equal-weight OLS sensitivity (both under 10,000 target permutations) + unadjusted
Spearman + whole-FAERS ROR/PRR/BCPNN IC025 + pairwise RORs + READUS-PV sensitivity suite.
"""

import os
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf

from tobit_lib import PAYLOAD_MAP, fit_nb2_offset

ROOT = Path(__file__).resolve().parents[1]
R = f"{ROOT}/results/"

COHORT_SPECS = [
    ("sole_suspect_noco", "Primary: FDA-deduplicated sole-suspect (PS, zero SS) monotherapy excluding co-reported antineoplastics"),
    ("ps_noco", "Sensitivity 1: FDA-deduplicated Primary Suspect (PS, SS allowed) excluding co-reported antineoplastics"),
    ("sole_suspect_noco_hcp", "Sensitivity 2: Primary restricted to healthcare professional reporters (MD/PH/OT)"),
    ("sole_suspect_noco_serious", "Sensitivity 3: Primary restricted to serious reports (any FDA OUTC code)"),
    ("sole_suspect_noco_nohnscc", "Sensitivity 4: Primary excluding head/neck/oral-cancer and radiotherapy indications"),
    ("sole_suspect_noco_secdedup", "Sensitivity 5: Primary plus secondary exact demographic/event/drug/PT deduplication"),
    ("ps_all", "Sensitivity 6 (negative-confounding check): All deduplicated PS reports without co-antineoplastic exclusion"),
]


def _calc_pair_ror(k1, n1, k2, n2):
    c1, c2 = n1 - k1, n2 - k2
    if k1 <= 0 or k2 <= 0 or c1 <= 0 or c2 <= 0:
        return np.nan, np.nan, np.nan
    r = (k1 / c1) / (k2 / c2)
    se = np.sqrt(1.0 / k1 + 1.0 / c1 + 1.0 / k2 + 1.0 / c2)
    return float(r), float(r * np.exp(-1.96 * se)), float(r * np.exp(1.96 * se))


def _calc_whole_faers_metrics(a, n_drug, n_ev, N_tot):
    b = n_drug - a
    c = n_ev - a
    d = N_tot - n_drug - n_ev + a
    if a > 0 and b > 0 and c > 0 and d > 0:
        ror = (a * d) / (b * c)
        se = np.sqrt(1.0 / a + 1.0 / b + 1.0 / c + 1.0 / d)
        ror_lo = ror * np.exp(-1.96 * se)
        ror_hi = ror * np.exp(1.96 * se)
        prr = (a / n_drug) / (c / (N_tot - n_drug))
    else:
        ror, ror_lo, ror_hi, prr = np.nan, np.nan, np.nan, np.nan
    e_ij = (n_drug * n_ev) / N_tot
    ic = float(np.log2((a + 0.5) / (e_ij + 0.5)))
    ic025 = float(ic - 3.3 * (a + 0.5) ** (-0.5) - 2.0 * (a + 0.5) ** (-1.5))
    return ror, ror_lo, ror_hi, prr, ic, ic025


def _fit_outcome_models(d, out, n_perm=10000):
    X0 = np.column_stack([
        np.ones(len(d)),
        (d.pclass == "topoI").astype(float).values,
        (d.pclass == "DNA").astype(float).values,
    ])
    X1 = np.column_stack([X0, d.W_basal.values])
    Q0, _ = np.linalg.qr(X0)
    n_arr = d.n_reports.values.astype(float)
    log_n = np.log(n_arr)
    k_arr = d["n_" + out].values.astype(float)
    w_vec = d.W_basal.values

    tg_s = d.groupby("tg").W_basal.first()
    tgs_f = tg_s.index.values
    tg_vals = tg_s.values
    tg_idx = np.array([np.where(tgs_f == t)[0][0] for t in d.tg.values])

    rate = k_arr / n_arr
    rho_sp, p_sp = stats.spearmanr(w_vec, rate)

    # 1. Primary exposure-offset NB2 count regression + B=10,000 target permutations
    r0 = fit_nb2_offset(k_arr, log_n, X0, compute_se=False)
    r1 = fit_nb2_offset(k_arr, log_n, X1, compute_se=True)
    b_nb = float(r1["beta"][-1])
    se_nb = float(r1["se_last"])
    p_wald = float(2.0 * stats.norm.sf(abs(b_nb) / se_nb))
    lr_nb = max(0.0, 2.0 * (r0["nll"] - r1["nll"]))
    p_lr = float(stats.chi2.sf(lr_nb, 1))

    warm = np.r_[r0["beta"], 0.0, np.log(1.0 / max(r0["alpha"], 1e-6))]
    rng_nb = np.random.default_rng(1)
    null_lr = np.zeros(n_perm)
    for b_ix in range(n_perm):
        perm_vals = rng_nb.permutation(tg_vals)
        Xb = np.column_stack([X0, perm_vals[tg_idx]])
        rb = fit_nb2_offset(k_arr, log_n, Xb, compute_se=False, warm_theta=warm)
        null_lr[b_ix] = max(0.0, 2.0 * (r0["nll"] - rb["nll"]))
    p_perm_target = float((1 + int((null_lr >= lr_nb).sum())) / (n_perm + 1))

    # 2. Equal-weight empirical log-odds OLS sensitivity + identical target permutations
    y_ols = np.log((k_arr + 0.5) / (n_arr - k_arr + 0.5))
    df_ols = pd.DataFrame({"y_ols": y_ols, "pclass": d.pclass.values, "W_basal": w_vec})
    m_ols = smf.ols("y_ols ~ C(pclass) + W_basal", df_ols).fit()
    b_ols = float(m_ols.params["W_basal"])
    se_ols = float(m_ols.bse["W_basal"])
    y_res = y_ols - Q0 @ (Q0.T @ y_ols)
    rss0 = float(y_res @ y_res)
    w_res0 = w_vec - Q0 @ (Q0.T @ w_vec)
    rss1 = float(rss0 - b_ols * b_ols * float(w_res0 @ w_res0))
    ols_lr = max(0.0, float(len(d) * np.log(rss0 / max(rss1, 1e-15))))

    rng_ols = np.random.default_rng(1)
    null_ols_lr = np.zeros(n_perm)
    for b_ix in range(n_perm):
        perm_vals = rng_ols.permutation(tg_vals)
        wb = perm_vals[tg_idx]
        wr = wb - Q0 @ (Q0.T @ wb)
        ww = float(wr @ wr)
        if ww < 1e-12:
            continue
        bb = float((wr @ y_res) / ww)
        rss_b = float(rss0 - bb * bb * ww)
        null_ols_lr[b_ix] = max(0.0, float(len(d) * np.log(rss0 / max(rss_b, 1e-15))))
    p_ols_perm = float((1 + int((null_ols_lr >= ols_lr - 1e-12).sum())) / (n_perm + 1))

    return dict(
        outcome=out,
        n_adcs=len(d),
        n_targets=len(tgs_f),
        total_reports=int(n_arr.sum()),
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
        nb_alpha=float(r1["alpha"]),
        spearman_rho=float(rho_sp),
        p_spearman=float(p_sp),
        ols_beta=b_ols,
        ols_se=se_ols,
        ols_LR=ols_lr,
        p_ols=float(m_ols.pvalues["W_basal"]),
        p_ols_perm=p_ols_perm,
    )


def _fit_outcome_task(args):
    d, out, n_perm = args
    return _fit_outcome_models(d, out, n_perm=n_perm)


def _quick_subset_nb2(sub):
    Xs0 = np.column_stack([
        np.ones(len(sub)),
        (sub.pclass == "topoI").astype(float).values,
        (sub.pclass == "DNA").astype(float).values,
    ])
    Xs1 = np.column_stack([Xs0, sub.W_basal.values])
    ln = np.log(sub.n_reports.values.astype(float))
    kk = sub.n_oral.values.astype(float)
    r0 = fit_nb2_offset(kk, ln, Xs0, compute_se=False)
    r1 = fit_nb2_offset(kk, ln, Xs1, compute_se=True)
    lr = max(0.0, 2.0 * (r0["nll"] - r1["nll"]))
    return float(r1["beta"][-1]), float(lr), float(stats.chi2.sf(lr, 1))


def run_faers_analysis():
    F = pd.read_csv(f"{R}faers_counts.tsv", sep="\t")
    bg = pd.read_csv(f"{R}faers_background.tsv", sep="\t").set_index("metric")["value"]
    bg_tot = int(bg["dedup_total_cases"])
    bg_oral = int(bg["dedup_oral"])

    allT = pd.read_csv(f"{ROOT}/clinical/adc_stomatitis_trials.tsv", sep="\t").set_index("drug")
    W = pd.read_csv(f"{R}williams_kc_stats.tsv.gz", sep="\t", index_col=0)

    def _prep(set_name):
        d = F[F.set == set_name].set_index("drug").join(allT[["target_gene", "payload_class"]])
        d["pclass"] = d.payload_class.map(PAYLOAD_MAP)
        d = d[d.n_reports >= 50].copy()
        d["W_basal"] = np.log2(W.reindex(d.target_gene.replace({"NECTIN4": "PVRL4"}))["w_basal_mean"].values + 1.0)
        d["tg"] = d.target_gene
        return d

    # 1. Primary cohort + sensitivity cohorts in parallel (single-threaded BLAS per worker)
    d_prim = _prep("sole_suspect_noco")
    outcomes = ["oral", "neutropenia", "neuropathy", "skin", "ocular", "ild", "diarrhea"]
    prim_tasks = [(d_prim, out, 10000) for out in outcomes]
    sens_cohorts = [(cname, desc, _prep(cname)) for cname, desc in COHORT_SPECS if cname != "sole_suspect_noco"]
    sens_tasks = [(dc, "oral", 10000) for _, _, dc in sens_cohorts]
    all_tasks = prim_tasks + sens_tasks
    with ProcessPoolExecutor(max_workers=min(7, os.cpu_count() or 1)) as ex:
        all_res = list(ex.map(_fit_outcome_task, all_tasks))
    frows = all_res[:len(outcomes)]
    sens_res_map = {"sole_suspect_noco": frows[0]}
    for (cname, _, _), res in zip(sens_cohorts, all_res[len(outcomes):]):
        sens_res_map[cname] = res
    faers_df = pd.DataFrame(frows)
    faers_df.to_csv(f"{R}fig_faers_summary.tsv", sep="\t", index=False)

    # 2. Per-drug whole-FAERS disproportionality table (ROR, PRR, BCPNN IC/IC025)
    drows = []
    for drug, row in d_prim.sort_values("W_basal", ascending=False).iterrows():
        ror, rlo, rhi, prr, ic, ic025 = _calc_whole_faers_metrics(
            int(row["n_oral"]), int(row["n_reports"]), bg_oral, bg_tot
        )
        drows.append(dict(
            drug=drug,
            target_gene=row["target_gene"],
            payload_class=row["payload_class"],
            pclass=row["pclass"],
            W_basal=round(float(row["W_basal"]), 4),
            n_reports=int(row["n_reports"]),
            n_oral=int(row["n_oral"]),
            oral_pct=round(100.0 * row["n_oral"] / row["n_reports"], 2),
            ROR_whole_FAERS=round(ror, 2) if np.isfinite(ror) else np.nan,
            ROR_ci_lo=round(rlo, 2) if np.isfinite(rlo) else np.nan,
            ROR_ci_hi=round(rhi, 2) if np.isfinite(rhi) else np.nan,
            PRR=round(prr, 2) if np.isfinite(prr) else np.nan,
            BCPNN_IC=round(ic, 2),
            BCPNN_IC025=round(ic025, 2),
            signal_IC025_pos="YES" if ic025 > 0 else "NO",
        ))
    disprop_df = pd.DataFrame(drows)
    disprop_df.to_csv(f"{R}faers_per_drug_disprop.tsv", sep="\t", index=False)

    # 3. Pairwise ROR contrasts (4 higher-vs-lower-basal + 1 same-target TACSTD2)
    pair_specs = [
        ("datopotamab deruxtecan", "trastuzumab deruxtecan", "Dato-DXd vs T-DXd\n(both DXd)", "Higher vs lower oral basal"),
        ("enfortumab vedotin", "brentuximab vedotin", "EV vs BV\n(both MMAE)", "Higher vs lower oral basal"),
        ("enfortumab vedotin", "polatuzumab vedotin", "EV vs Pola\n(both MMAE)", "Higher vs lower oral basal"),
        ("trastuzumab emtansine", "mirvetuximab soravtansine", "T-DM1 vs Mirv\n(DM1 vs DM4)", "Higher vs lower oral basal"),
        ("datopotamab deruxtecan", "sacituzumab govitecan", "Dato-DXd vs SG [same target]\n(DXd vs SN-38, TACSTD2)", "Same target (TACSTD2)"),
    ]
    prows = []
    for d1, d2, label, ptype in pair_specs:
        k1, n1 = int(d_prim.loc[d1, "n_oral"]), int(d_prim.loc[d1, "n_reports"])
        k2, n2 = int(d_prim.loc[d2, "n_oral"]), int(d_prim.loc[d2, "n_reports"])
        r_val, lo_val, hi_val = _calc_pair_ror(k1, n1, k2, n2)
        prows.append(dict(
            drug1=d1, drug2=d2, contrast=label, pair_type=ptype,
            basal1=float(d_prim.loc[d1, "W_basal"]), basal2=float(d_prim.loc[d2, "W_basal"]),
            k1=k1, n1=n1, pct1=round(100.0 * k1 / n1, 2),
            k2=k2, n2=n2, pct2=round(100.0 * k2 / n2, 2),
            ror=round(r_val, 2), lo=round(lo_val, 2), hi=round(hi_val, 2),
        ))
    pairs_df = pd.DataFrame(prows)
    pairs_df.to_csv(f"{R}fig_faers_pairs.tsv", sep="\t", index=False)

    # 4. READUS-PV sensitivity suite across all 7 cohorts (oral outcome, 10,000 target permutations)
    srows = []
    for cname, desc in COHORT_SPECS:
        dc = _prep(cname)
        res = sens_res_map[cname]
        r_dt, lo_dt, hi_dt = _calc_pair_ror(
            int(dc.loc["datopotamab deruxtecan", "n_oral"]), int(dc.loc["datopotamab deruxtecan", "n_reports"]),
            int(dc.loc["trastuzumab deruxtecan", "n_oral"]), int(dc.loc["trastuzumab deruxtecan", "n_reports"]),
        )
        r_eb, lo_eb, hi_eb = _calc_pair_ror(
            int(dc.loc["enfortumab vedotin", "n_oral"]), int(dc.loc["enfortumab vedotin", "n_reports"]),
            int(dc.loc["brentuximab vedotin", "n_oral"]), int(dc.loc["brentuximab vedotin", "n_reports"]),
        )
        bd_b, bd_lr, bd_p = _quick_subset_nb2(dc.drop(index=["datopotamab deruxtecan"], errors="ignore"))
        bt_b, bt_lr, bt_p = _quick_subset_nb2(dc[dc.target_gene != "TACSTD2"])
        srows.append(dict(
            cohort=cname,
            description=desc,
            n_adcs=res["n_adcs"],
            n_targets=res["n_targets"],
            total_reports=res["total_reports"],
            oral_events=res["total_events"],
            beta=round(res["beta"], 4),
            se=round(res["se"], 4),
            ci_lo=round(res["ci_lo"], 4),
            ci_hi=round(res["ci_hi"], 4),
            rate_ratio_per_log2=round(res["rate_ratio_per_log2"], 3),
            p_wald=res["p_wald"],
            LR=round(res["LR"], 3),
            p_LR=res["p_LR"],
            p_perm=res["p_perm"],
            spearman_rho=round(res["spearman_rho"], 4),
            p_spearman=res["p_spearman"],
            ols_beta=round(res["ols_beta"], 4),
            ols_LR=round(res["ols_LR"], 3),
            p_ols=res["p_ols"],
            p_ols_perm=res["p_ols_perm"],
            dato_vs_tdxd_ror=f"{r_dt:.2f} ({lo_dt:.2f}–{hi_dt:.2f})",
            ev_vs_bv_ror=f"{r_eb:.2f} ({lo_eb:.2f}–{hi_eb:.2f})",
            drop_dato_beta=round(bd_b, 3),
            drop_dato_p_LR=round(bd_p, 4),
            drop_tacstd2_beta=round(bt_b, 3),
            drop_tacstd2_p_LR=round(bt_p, 4),
        ))
    sens_df = pd.DataFrame(srows)
    sens_df.to_csv(f"{R}faers_sensitivity.tsv", sep="\t", index=False)

    return faers_df, pairs_df, disprop_df, sens_df


if __name__ == "__main__":
    fdf, pdf, ddf, sdf = run_faers_analysis()
    print("=== Primary Sole-Suspect Monotherapy FAERS Summary ===")
    print(fdf[["outcome", "total_events", "beta", "se", "p_wald", "LR", "p_LR", "p_perm", "spearman_rho", "p_spearman", "ols_beta", "ols_LR", "p_ols", "p_ols_perm"]].round(4).to_string(index=False))
    print("\n=== Per-Drug Whole-FAERS Disproportionality ===")
    print(ddf.to_string(index=False))
    print("\n=== Pairwise Oral ROR Contrasts ===")
    print(pdf[["drug1", "drug2", "ror", "lo", "hi"]].to_string(index=False))
    print("\n=== READUS-PV Sensitivity Suite ===")
    print(sdf[["cohort", "n_adcs", "total_reports", "oral_events", "beta", "LR", "p_LR", "p_perm", "spearman_rho", "p_spearman", "p_ols_perm"]].to_string(index=False))
