"""Part 1 HPA IHC color deconvolution: deterministic Ruifrok & Johnston rgb2hed DAB quantification across 509 squamous cores."""

import multiprocessing as mp
import os
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage, stats
from statsmodels.stats.multitest import multipletests

try:
    import pymupdf
except ImportError:
    import fitz as pymupdf

from tobit_lib import PAYLOAD_MAP, design_matrix, fit_tobit, lr_test

ROOT = str(Path(__file__).resolve().parents[1])
R = f'{ROOT}/results'
P1 = f'{R}/part1'
QC_DIR = f'{R}/dab_qc'
os.makedirs(P1, exist_ok=True)
os.makedirs(QC_DIR, exist_ok=True)

# Ruifrok & Johnston (2001) optical density color deconvolution matrix (standard RGB -> H-E-DAB)
HED_FROM_RGB = np.linalg.inv(np.array([
    [0.650, 0.704, 0.286],  # Hematoxylin
    [0.072, 0.990, 0.105],  # Eosin
    [0.268, 0.570, 0.776],  # DAB
], dtype=np.float64)).T

NO_IHC_GENES = ['DLL3', 'F3', 'FOLR1', 'GPNMB', 'ROR1', 'TNFRSF17']
QC_GENES = {'TACSTD2', 'CDH3', 'PTPRZ1', 'COL17A1', 'DSC3', 'ERBB2'}


def rgb_to_hed(rgb_u8):
    rgb = np.maximum(rgb_u8.astype(np.float64) / 255.0, 1e-4)
    od = -np.log10(rgb)
    return np.tensordot(od, HED_FROM_RGB, axes=([2], [1]))


def analyze_core(jpg_path, stride=4):
    """Segment epithelial nuclei-dense compartment and basal margin (dist_stroma <= 18 px, rel_pos <= 0.30)."""
    pix = pymupdf.Pixmap(jpg_path)
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3]
    arr_ds = arr[::stride, ::stride, :]
    hed = rgb_to_hed(arr_ds)
    hem = hed[:, :, 0]
    dab = np.maximum(hed[:, :, 2], 0.0)

    gray = arr_ds.mean(axis=2)
    tissue = gray < 236.0
    bg = ~tissue

    bg_lbl, n_bg = ndimage.label(bg)
    if n_bg > 0:
        sizes = np.bincount(bg_lbl.ravel())
        lumen = (sizes[bg_lbl] > 180) & bg
    else:
        lumen = bg

    cell_score = ndimage.gaussian_filter(hem + 0.35 * dab, sigma=5.0)
    if tissue.sum() < 500:
        return dict(
            ok=False, tissue_frac=0.0, epi_frac=0.0,
            dab_epi_mean=np.nan, dab_basal_mean=np.nan,
            dab_suprabasal_mean=np.nan, dab_basal_p75=np.nan,
        )

    t_vals = cell_score[tissue]
    thr = max(float(np.quantile(t_vals, 0.58)), 0.03)
    epi_raw = tissue & (cell_score >= thr)
    epi_clean = ndimage.binary_opening(epi_raw, iterations=2)
    epi_clean = ndimage.binary_closing(epi_clean, iterations=3)
    lbl, n_c = ndimage.label(epi_clean)
    if n_c > 0:
        sz = np.bincount(lbl.ravel())
        sz[0] = 0
        epi_mask = (sz[lbl] >= 300)
    else:
        epi_mask = epi_clean

    if epi_mask.sum() < 400:
        return dict(
            ok=False,
            tissue_frac=float(tissue.mean()),
            epi_frac=float(epi_mask.mean()),
            dab_epi_mean=float(dab[tissue].mean()),
            dab_basal_mean=float(dab[tissue].mean()),
            dab_suprabasal_mean=float(dab[tissue].mean()),
            dab_basal_p75=float(np.quantile(dab[tissue], 0.75)),
            arr_ds=arr_ds, epi_mask=epi_mask,
            basal_mask=np.zeros_like(epi_mask), supra_mask=np.zeros_like(epi_mask), dab=dab,
        )

    stroma_mask = tissue & (~epi_mask)
    if stroma_mask.sum() < 100 or lumen.sum() < 100:
        return dict(
            ok=False,
            tissue_frac=float(tissue.mean()),
            epi_frac=float(epi_mask.mean()),
            dab_epi_mean=float(dab[epi_mask].mean()),
            dab_basal_mean=float(dab[epi_mask].mean()),
            dab_suprabasal_mean=float(dab[epi_mask].mean()),
            dab_basal_p75=float(np.quantile(dab[epi_mask], 0.75)),
            arr_ds=arr_ds, epi_mask=epi_mask,
            basal_mask=epi_mask, supra_mask=np.zeros_like(epi_mask), dab=dab,
        )

    d_stroma = ndimage.distance_transform_edt(~stroma_mask)
    d_lumen = ndimage.distance_transform_edt(~lumen)
    rel_pos = d_stroma / (d_stroma + d_lumen + 1e-6)

    basal_mask = epi_mask & ((rel_pos <= 0.30) & (d_stroma <= 18.0))
    supra_mask = epi_mask & (~basal_mask) & (rel_pos <= 0.85)
    if basal_mask.sum() < 50:
        basal_mask = epi_mask & (rel_pos <= 0.35)
    if supra_mask.sum() < 50:
        supra_mask = epi_mask & (~basal_mask)

    return dict(
        ok=True,
        tissue_frac=float(tissue.mean()),
        epi_frac=float(epi_mask.mean()),
        dab_epi_mean=float(dab[epi_mask].mean()),
        dab_basal_mean=float(dab[basal_mask].mean()) if basal_mask.sum() > 0 else np.nan,
        dab_suprabasal_mean=float(dab[supra_mask].mean()) if supra_mask.sum() > 0 else np.nan,
        dab_basal_p75=float(np.quantile(dab[basal_mask], 0.75)) if basal_mask.sum() > 0 else np.nan,
        arr_ds=arr_ds, epi_mask=epi_mask, basal_mask=basal_mask, supra_mask=supra_mask, dab=dab,
    )


def save_qc_overlay(res, out_png):
    if 'arr_ds' not in res:
        return
    vis = res['arr_ds'].copy().astype(np.float32)
    b = res['basal_mask']
    s = res['supra_mask']
    vis[b, 0] = 0.45 * vis[b, 0] + 0.55 * 20
    vis[b, 1] = 0.45 * vis[b, 1] + 0.55 * 210
    vis[b, 2] = 0.45 * vis[b, 2] + 0.55 * 40
    vis[s, 0] = 0.55 * vis[s, 0] + 0.45 * 30
    vis[s, 1] = 0.55 * vis[s, 1] + 0.45 * 140
    vis[s, 2] = 0.55 * vis[s, 2] + 0.45 * 235
    out_u8 = np.clip(vis, 0, 255).astype(np.uint8)
    h, w, _ = out_u8.shape
    pix = pymupdf.Pixmap(pymupdf.csRGB, w, h, out_u8.tobytes(), False)
    pix.save(out_png)


def target_perm_test(df, col, n_perm=10000, seed=1):
    d = df.dropna(subset=[col]).reset_index(drop=True)
    obs = lr_test(d, col, payload=True, compute_se=True)
    obs_lr = obs['LR']
    nll0 = obs['nll0']
    y_v = d.y_pct.values
    c_v = d.cens.values
    X_base = design_matrix(d)
    tg_vals = d.groupby('tg')[col].first()
    keys = tg_vals.index.values
    rng = np.random.default_rng(seed)
    perm_lrs = np.zeros(n_perm)
    for i in range(n_perm):
        pmap = dict(zip(keys, rng.permutation(tg_vals.values)))
        Xp = np.column_stack([X_base, d.tg.map(pmap).values])
        rp = fit_tobit(y_v, c_v, Xp)
        perm_lrs[i] = max(0.0, 2.0 * (nll0 - rp['nll']))
    perm_p = float((1 + (perm_lrs >= obs_lr).sum()) / (n_perm + 1))
    return obs, perm_p


def cluster_boot_dlr_multi(df_all, cols, n_boot=2000, seed=1):
    """Target-cluster bootstrap across all 28 primary targets in a single pass (identical RNG draws to part1_model.py)."""
    obs = {}
    for col in cols:
        d_sub = df_all.dropna(subset=[col, 'W_basal'])
        obs[col] = float(lr_test(d_sub, 'W_basal', payload=True)['LR'] - lr_test(d_sub, col, payload=True)['LR'])
    tgs = df_all.tg.unique()
    tg_idx = {t: np.where(df_all.tg.values == t)[0] for t in tgs}
    rng = np.random.default_rng(seed)
    boot_diffs = {col: [] for col in cols}
    for _ in range(n_boot):
        pick = rng.choice(tgs, len(tgs), replace=True)
        idx = np.concatenate([tg_idx[t] for t in pick])
        D = df_all.iloc[idx].reset_index(drop=True)
        if D.cens.all() or D.pclass.nunique() < 2:
            continue
        r0_all = fit_tobit(D.y_pct.values, D.cens.values, design_matrix(D))
        r1_b_all = fit_tobit(D.y_pct.values, D.cens.values, design_matrix(D, 'W_basal'))
        if not (r0_all['converged'] and r1_b_all['converged']):
            continue
        lr_b_all = max(0.0, 2.0 * (r0_all['nll'] - r1_b_all['nll']))
        mask_cache = {}
        for col in cols:
            na_mask = D[col].isna().values
            if not na_mask.any():
                r1_c = fit_tobit(D.y_pct.values, D.cens.values, design_matrix(D, col))
                if r1_c['converged']:
                    boot_diffs[col].append(lr_b_all - max(0.0, 2.0 * (r0_all['nll'] - r1_c['nll'])))
            else:
                m_key = na_mask.tobytes()
                Dc = D.iloc[~na_mask]
                if len(Dc) < 15 or Dc.cens.all() or Dc[col].std() < 1e-8:
                    continue
                if m_key not in mask_cache:
                    rc0 = fit_tobit(Dc.y_pct.values, Dc.cens.values, design_matrix(Dc))
                    rb1 = fit_tobit(Dc.y_pct.values, Dc.cens.values, design_matrix(Dc, 'W_basal'))
                    mask_cache[m_key] = (rc0, rb1)
                else:
                    rc0, rb1 = mask_cache[m_key]
                if rc0['converged'] and rb1['converged']:
                    rc1 = fit_tobit(Dc.y_pct.values, Dc.cens.values, design_matrix(Dc, col))
                    if rc1['converged']:
                        boot_diffs[col].append(max(0.0, 2.0 * (rc0['nll'] - rb1['nll'])) - max(0.0, 2.0 * (rc0['nll'] - rc1['nll'])))
    res_map = {}
    for col in cols:
        arr = np.array(boot_diffs[col])
        lo, hi = np.quantile(arr, [0.025, 0.975]) if len(arr) >= 20 else (np.nan, np.nan)
        res_map[col] = (obs[col], float(lo), float(hi), int(len(arr)))
    return res_map


def _grid_worker(args):
    gene, tg_grp, img_path = args
    pix = pymupdf.Pixmap(img_path)
    arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, pix.n)[:, :, :3]
    arr_ds = arr[::4, ::4, :]
    hed = rgb_to_hed(arr_ds)
    hem = hed[:, :, 0]
    dab = np.maximum(hed[:, :, 2], 0.0)
    gray = arr_ds.mean(axis=2)
    cs = ndimage.gaussian_filter(hem + 0.35 * dab, sigma=5.0)
    out = {}
    for g_thr in [228.0, 236.0, 244.0]:
        tissue = gray < g_thr
        if tissue.sum() < 500:
            for q_epi in [0.50, 0.58, 0.66]:
                for r_pos in [0.25, 0.30, 0.35]:
                    for d_px in [14.0, 18.0, 22.0]:
                        out[(q_epi, r_pos, d_px, g_thr)] = np.nan
            continue
        bg = ~tissue
        bg_lbl, n_bg = ndimage.label(bg)
        if n_bg > 0:
            sizes = np.bincount(bg_lbl.ravel())
            lumen = (sizes[bg_lbl] > 180) & bg
        else:
            lumen = bg
        t_cs = cs[tissue]
        t_dab_mean = float(dab[tissue].mean())
        for q_epi in [0.50, 0.58, 0.66]:
            thr = max(float(np.quantile(t_cs, q_epi)), 0.03)
            epi_raw = tissue & (cs >= thr)
            epi_clean = ndimage.binary_closing(ndimage.binary_opening(epi_raw, iterations=2), iterations=3)
            lbl, n_c = ndimage.label(epi_clean)
            if n_c > 0:
                sz = np.bincount(lbl.ravel())
                sz[0] = 0
                epi_mask = (sz[lbl] >= 300)
            else:
                epi_mask = epi_clean
            if epi_mask.sum() < 400:
                for r_pos in [0.25, 0.30, 0.35]:
                    for d_px in [14.0, 18.0, 22.0]:
                        out[(q_epi, r_pos, d_px, g_thr)] = t_dab_mean
                continue
            stroma_mask = tissue & (~epi_mask)
            if stroma_mask.sum() < 100 or lumen.sum() < 100:
                em = float(dab[epi_mask].mean())
                for r_pos in [0.25, 0.30, 0.35]:
                    for d_px in [14.0, 18.0, 22.0]:
                        out[(q_epi, r_pos, d_px, g_thr)] = em
                continue
            d_stroma = ndimage.distance_transform_edt(~stroma_mask)[epi_mask]
            d_lumen = ndimage.distance_transform_edt(~lumen)[epi_mask]
            rel_pos = d_stroma / (d_stroma + d_lumen + 1e-6)
            d_epi = dab[epi_mask]
            for r_pos in [0.25, 0.30, 0.35]:
                for d_px in [14.0, 18.0, 22.0]:
                    bm = (rel_pos <= r_pos) & (d_stroma <= d_px)
                    if bm.sum() < 50:
                        bm = (rel_pos <= (r_pos + 0.05))
                    out[(q_epi, r_pos, d_px, g_thr)] = float(d_epi[bm].mean()) if bm.sum() > 0 else np.nan
    return gene, tg_grp, out


def run_param_grid(man, T):
    sub_man = man[man.tissue_group.isin(['Esophagus', 'Oral mucosa'])].reset_index(drop=True)
    tasks = [(r['gene'], r['tissue_group'], f"{ROOT}/{r['rel_path']}") for _, r in sub_man.iterrows()]
    with mp.Pool(min(16, os.cpu_count() or 4)) as pool:
        res_all = pool.map(_grid_worker, tasks)
    grid_rows = []
    for q_epi in [0.50, 0.58, 0.66]:
        for r_pos in [0.25, 0.30, 0.35]:
            for d_px in [14.0, 18.0, 22.0]:
                for g_thr in [228.0, 236.0, 244.0]:
                    key = (q_epi, r_pos, d_px, g_thr)
                    es_v, or_v = {}, {}
                    for gene, tg_grp, out in res_all:
                        v = out[key]
                        if np.isfinite(v):
                            (es_v if tg_grp == 'Esophagus' else or_v).setdefault(gene, []).append(v)
                    es_m = {g: np.mean(vs) for g, vs in es_v.items()}
                    or_m = {g: np.mean(vs) for g, vs in or_v.items()}
                    T_sub = T.copy()
                    T_sub['es'] = [max([es_m[g] for g in str(gs).split(';') if g in es_m], default=np.nan) for gs in T_sub.target_gene]
                    T_sub['or'] = [max([or_m[g] for g in str(gs).split(';') if g in or_m], default=np.nan) for gs in T_sub.target_gene]
                    de = T_sub.dropna(subset=['es'])
                    do = T_sub.dropna(subset=['or'])
                    re = lr_test(de, 'es', payload=True, compute_se=True)
                    ro = lr_test(do, 'or', payload=True, compute_se=True)
                    is_def = (q_epi == 0.58) and (r_pos == 0.30) and (d_px == 18.0) and (g_thr == 236.0)
                    grid_rows.append(dict(
                        q_epi=q_epi, r_pos=r_pos, d_px=d_px, g_thr=g_thr, is_default=is_def,
                        esoph_n_adc=len(de), esoph_beta=float(re['beta']), esoph_se=float(re['se']),
                        esoph_LR=float(re['LR']), esoph_p=float(re['p']),
                        oral_n_adc=len(do), oral_beta=float(ro['beta']), oral_se=float(ro['se']),
                        oral_LR=float(ro['LR']), oral_p=float(ro['p']),
                    ))
    grid_df = pd.DataFrame(grid_rows)
    grid_df.to_csv(f'{P1}/hpa_ihc_dab_param_grid.tsv', sep='\t', index=False)
    print('=== 81-COMBINATION SEGMENTATION SENSITIVITY GRID ===')
    print(f"Esophagus LR median={grid_df.esoph_LR.median():.2f} (range {grid_df.esoph_LR.min():.2f}-{grid_df.esoph_LR.max():.2f}), p<0.05: {(grid_df.esoph_p < 0.05).sum()}/{len(grid_df)}, p<0.01: {(grid_df.esoph_p < 0.01).sum()}/{len(grid_df)}")
    print(f"Oral LR median={grid_df.oral_LR.median():.2f} (range {grid_df.oral_LR.min():.2f}-{grid_df.oral_LR.max():.2f}), p<0.05: {(grid_df.oral_p < 0.05).sum()}/{len(grid_df)}")
    return grid_df


def main():
    man = pd.read_csv(f'{ROOT}/dl/hpa_ihc/manifest_509.tsv', sep='\t')
    rows = []
    qc_done = set()
    for _, r in man.iterrows():
        img_path = f"{ROOT}/{r['rel_path']}"
        res = analyze_core(img_path)
        rows.append(dict(
            slide_id=r['slide_id'],
            gene=r['gene'],
            antibody_id=r['antibody_id'],
            hpa_reliability=r['hpa_reliability'],
            tissue=r['tissue'],
            tissue_group=r['tissue_group'],
            hpa_staining=r['hpa_staining'],
            seg_ok=res['ok'],
            tissue_frac=res['tissue_frac'],
            epi_frac=res['epi_frac'],
            dab_epi_mean=res['dab_epi_mean'],
            dab_basal_mean=res['dab_basal_mean'],
            dab_suprabasal_mean=res['dab_suprabasal_mean'],
            dab_basal_p75=res['dab_basal_p75'],
        ))
        tag = f"{r['gene']}_{r['tissue_group']}"
        if res['ok'] and (r['gene'] in QC_GENES) and (r['tissue_group'] in ['Oral mucosa', 'Esophagus']) and (tag not in qc_done):
            qc_done.add(tag)
            save_qc_overlay(res, f"{QC_DIR}/{r['slide_id']}_{r['gene']}_{r['tissue_group'].replace(' ', '_')}.png")

    slides = pd.DataFrame(rows)
    slides.to_csv(f'{P1}/hpa_ihc_dab_slides.tsv', sep='\t', index=False)
    print(f'Quantified {len(slides)} HPA squamous IHC cores via rgb2hed DAB deconvolution (seg_ok={int(slides.seg_ok.sum())}/{len(slides)}).')

    # Gene-level aggregation across 34 clinical/benchmark genes
    W = pd.read_csv(f'{R}/williams_kc_stats.tsv.gz', sep='\t', index_col=0)
    H = pd.read_csv(f'{R}/hpa_ihc_scored.tsv.gz', sep='\t')
    hpa_oral_sum = H[H.Tissue == 'Oral mucosa'].groupby('Gene name').score.max()
    alias = {'NECTIN4': 'PVRL4'}
    full_34 = sorted(set(slides.gene.unique()) | set(NO_IHC_GENES))

    def _agg(sub_df, prefix):
        g = sub_df.groupby('gene').agg(
            n=('slide_id', 'count'),
            epi_frac_median=('epi_frac', 'median'),
            epi_frac_mean=('epi_frac', 'mean'),
            dab_basal=('dab_basal_mean', 'mean'),
            dab_suprabasal=('dab_suprabasal_mean', 'mean'),
            dab_epi=('dab_epi_mean', 'mean'),
            dab_p75=('dab_basal_p75', 'mean'),
        )
        g.columns = [f'{prefix}_{c}' for c in g.columns]
        return g

    slides_es = slides[slides.hpa_reliability.isin(['enhanced', 'supported'])]
    slides_oral = slides[slides.tissue_group == 'Oral mucosa']
    oral_epige15 = slides_oral[slides_oral.epi_frac >= 0.15].groupby('gene')['dab_basal_mean'].mean().rename('oral_dab_epige15')

    def _area_w_mean(sub_g):
        w = np.maximum(sub_g['epi_frac'].values, 1e-6)
        v = sub_g['dab_basal_mean'].values
        m = np.isfinite(v) & np.isfinite(w)
        return float(np.average(v[m], weights=w[m])) if m.any() else np.nan

    oral_area_wt = slides_oral.groupby('gene')[['epi_frac', 'dab_basal_mean']].apply(_area_w_mean).rename('oral_dab_area_wt')

    gdf = pd.DataFrame(index=full_34)
    gdf.index.name = 'gene'
    for part in [
        _agg(slides_oral, 'oral'),
        _agg(slides[slides.tissue_group == 'Esophagus'], 'esoph'),
        _agg(slides[slides.tissue_group == 'Skin'], 'skin'),
        _agg(slides_es[slides_es.tissue_group == 'Esophagus'], 'es_esoph'),
        _agg(slides_es[slides_es.tissue_group == 'Oral mucosa'], 'es_oral'),
        oral_epige15,
        oral_area_wt,
    ]:
        gdf = gdf.join(part, how='left')

    rel_map = slides.groupby('gene').hpa_reliability.apply(lambda s: ','.join(sorted(set(s))))
    gdf['hpa_reliability'] = gdf.index.map(rel_map).fillna('no_hpa_ihc')
    gdf['W_basal'] = [float(np.log2(W.loc[alias.get(g, g), 'w_basal_mean'] + 1.0)) if alias.get(g, g) in W.index else np.nan for g in gdf.index]
    gdf['HPA_IHC_oral_summary'] = [hpa_oral_sum.get(g, np.nan) for g in gdf.index]
    gdf['esoph_dab_rank_asc'] = gdf['esoph_dab_basal'].rank(method='min', ascending=True)
    gdf['low_basal_rna_lt005'] = (gdf['W_basal'] < 0.05) & gdf['esoph_dab_basal'].notna()
    gdf.to_csv(f'{P1}/hpa_ihc_dab_genes.tsv', sep='\t')

    # Tissue coverage comparison: Oral mucosa vs Esophagus vs Skin (slide-level & gene-level paired Wilcoxon across 28 genes)
    both_ef = gdf[['oral_epi_frac_median', 'esoph_epi_frac_median']].dropna()
    wx_stat, wx_p = stats.wilcoxon(both_ef['esoph_epi_frac_median'], both_ef['oral_epi_frac_median'])
    oral_ef = slides_oral['epi_frac'].dropna()
    esoph_ef = slides[slides.tissue_group == 'Esophagus']['epi_frac'].dropna()
    _, mw_p = stats.mannwhitneyu(esoph_ef, oral_ef)

    prefix_map = {'Oral mucosa': 'oral', 'Esophagus': 'esoph', 'Skin': 'skin'}
    cov_rows = []
    for tg_name in ['Oral mucosa', 'Esophagus', 'Skin']:
        sub = slides[slides.tissue_group == tg_name]
        pref = prefix_map[tg_name]
        cov_rows.append(dict(
            tissue_group=tg_name,
            n_cores=len(sub),
            n_genes=int(sub.gene.nunique()),
            tissue_frac_median=float(sub.tissue_frac.median()),
            epi_frac_median=float(sub.epi_frac.median()),
            epi_frac_iqr_lo=float(sub.epi_frac.quantile(0.25)),
            epi_frac_iqr_hi=float(sub.epi_frac.quantile(0.75)),
            epi_frac_mean=float(sub.epi_frac.mean()),
            gene_epifrac_median=float(gdf[f'{pref}_epi_frac_median'].median()),
            dab_basal_median=float(sub.dab_basal_mean.median()),
        ))
    cov_df = pd.DataFrame(cov_rows)
    cov_df['mw_p_esoph_vs_oral_epifrac'] = float(mw_p)
    cov_df['wilcoxon_stat_gene_epifrac'] = float(wx_stat)
    cov_df['wilcoxon_p_esoph_vs_oral_epifrac'] = float(wx_p)
    cov_df.to_csv(f'{P1}/hpa_ihc_dab_coverage.tsv', sep='\t', index=False)
    print(f'Median epi_frac (slides): Oral={oral_ef.median():.3f} (IQR {oral_ef.quantile(0.25):.3f}-{oral_ef.quantile(0.75):.3f}) vs Esophagus={esoph_ef.median():.3f} (IQR {esoph_ef.quantile(0.25):.3f}-{esoph_ef.quantile(0.75):.3f})')
    print(f'Gene-level paired Wilcoxon (n={len(both_ef)}): Oral={both_ef.oral_epi_frac_median.median():.3f} vs Esophagus={both_ef.esoph_epi_frac_median.median():.3f}, stat={wx_stat:.1f}, p={wx_p:.2e}')

    # Clinical Tobit models on adjudicated cohort (all 9 tissue x summary combos + sensitivity)
    T = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_trials.tsv', sep='\t')
    T = T[T.listed.isin(['yes', 'not_listed']) & ~(T.stom_any.isna() & T.threshold_pct.isna())].copy()
    T['pclass'] = T.payload_class.map(PAYLOAD_MAP)
    T['cens'] = T.stom_any.isna().astype(bool)
    T['y_pct'] = np.where(T.cens, T.threshold_pct, T.stom_any)
    T['tg'] = T.target_gene.str.split(';').str[0]
    T['W_basal'] = [
        max([float(np.log2(W.loc[alias.get(g, g), 'w_basal_mean'] + 1.0)) for g in str(gs).split(';') if alias.get(g, g) in W.index])
        for gs in T.target_gene
    ]
    T['HPA_IHC_oral_summary'] = [
        (lambda vs: max(vs) if vs else np.nan)([float(hpa_oral_sum.get(g, np.nan)) for g in str(gs).split(';') if pd.notna(hpa_oral_sum.get(g, np.nan))])
        for gs in T.target_gene
    ]

    def feat_from_gdf(gs, col):
        vs = [gdf.loc[g, col] for g in str(gs).split(';') if g in gdf.index]
        vs = [float(v) for v in vs if pd.notna(v)]
        return max(vs) if vs else np.nan

    specs = [
        ('W_basal_matched', 'scRNA oral basal KC reference (matched HPA IHC subset)', None, False, False),
        ('HPA_IHC_oral_summary', 'HPA oral mucosa IHC summary score (ordinal 0-3)', None, True, False),
        ('esoph_dab_basal', 'rgb2hed DAB: HPA esophagus basal mean OD [primary]', 'esoph_dab_basal', True, True),
        ('esoph_dab_p75', 'rgb2hed DAB: HPA esophagus basal 75th percentile OD', 'esoph_dab_p75', True, True),
        ('esoph_dab_epi', 'rgb2hed DAB: HPA esophagus whole-epithelium mean OD', 'esoph_dab_epi', True, True),
        ('oral_dab_basal', 'rgb2hed DAB: HPA oral mucosa basal mean OD', 'oral_dab_basal', True, True),
        ('oral_dab_p75', 'rgb2hed DAB: HPA oral mucosa basal 75th percentile OD', 'oral_dab_p75', True, True),
        ('oral_dab_epi', 'rgb2hed DAB: HPA oral mucosa whole-epithelium mean OD', 'oral_dab_epi', True, True),
        ('skin_dab_basal', 'rgb2hed DAB: HPA skin basal mean OD', 'skin_dab_basal', True, True),
        ('skin_dab_p75', 'rgb2hed DAB: HPA skin basal 75th percentile OD', 'skin_dab_p75', True, True),
        ('skin_dab_epi', 'rgb2hed DAB: HPA skin whole-epithelium mean OD', 'skin_dab_epi', True, True),
        ('es_esoph_dab_basal', 'rgb2hed DAB: HPA esophagus basal mean OD (Enhanced/Supported Abs)', 'es_esoph_dab_basal', True, False),
        ('oral_dab_epige15', 'rgb2hed DAB: HPA oral mucosa basal mean OD (epi_frac >= 0.15 cores)', 'oral_dab_epige15', True, False),
        ('oral_dab_area_wt', 'rgb2hed DAB: HPA oral mucosa basal mean OD (epi_frac area-weighted)', 'oral_dab_area_wt', True, False),
    ]

    for col, _, gcol, _, _ in specs:
        if gcol is not None:
            T[col] = [feat_from_gdf(gs, gcol) for gs in T.target_gene]

    boot_cols = [col for col, _, _, do_perm, _ in specs if do_perm]
    dlr_boot_map = cluster_boot_dlr_multi(T, boot_cols, n_boot=2000, seed=1)
    matched_mask = T['HPA_IHC_oral_summary'].notna()
    mod_rows = []
    for col, desc, gcol, do_perm, in_9fam in specs:
        if col == 'W_basal_matched':
            d = T[matched_mask].dropna(subset=['W_basal']).reset_index(drop=True)
            xcol = 'W_basal'
        else:
            d = T.dropna(subset=[col]).reset_index(drop=True)
            xcol = col
        if do_perm:
            res, perm_p = target_perm_test(d, xcol, n_perm=10000, seed=1)
            dlr, dlo, dhi, _ = dlr_boot_map[xcol]
        else:
            res = lr_test(d, xcol, payload=True, compute_se=True)
            perm_p = np.nan
            dlr = lr_test(d, 'W_basal', payload=True)['LR'] - res['LR']
            dlo, dhi = np.nan, np.nan
        imp = np.where(d.cens, d.y_pct / 2.0, d.y_pct)
        rho_clin, _ = stats.spearmanr(d[xcol], imp)
        rho_wb_adc, _ = stats.spearmanr(d[xcol], d['W_basal'])
        if gcol is not None:
            gg = gdf[[gcol, 'W_basal']].dropna()
            rho_wb_gene, p_wb_gene = stats.spearmanr(gg[gcol], gg['W_basal'])
        elif col == 'HPA_IHC_oral_summary':
            gg = T[['tg', 'HPA_IHC_oral_summary', 'W_basal']].dropna().drop_duplicates('tg')
            rho_wb_gene, p_wb_gene = stats.spearmanr(gg['HPA_IHC_oral_summary'], gg['W_basal'])
        else:
            rho_wb_gene, p_wb_gene = 1.0, 0.0

        mod_rows.append(dict(
            feature=col, description=desc, is_9combo_family=in_9fam,
            primary_dab_metric=(col == 'esoph_dab_basal'),
            n_adc=len(d), n_num=int((~d.cens).sum()), n_targets=int(d.tg.nunique()),
            rho_vs_Wbasal_genes=float(rho_wb_gene), p_vs_Wbasal_genes=float(p_wb_gene),
            rho_vs_Wbasal_adcs=float(rho_wb_adc), rho_vs_clin=float(rho_clin),
            beta=float(res['beta']), beta_per_SD=float(res['beta'] * d[xcol].std()),
            se=float(res['se']), LR=float(res['LR']), p_LR=float(res['p']),
            perm_p=perm_p, holm_p_LR=np.nan, bh_p_LR=np.nan, holm_perm_p=np.nan,
            dLR_vs_Wbasal=float(dlr), dLR_ci_lo=dlo, dLR_ci_hi=dhi,
        ))

    mod_df = pd.DataFrame(mod_rows)
    fam_mask = mod_df['is_9combo_family'].values
    mod_df.loc[fam_mask, 'holm_p_LR'] = multipletests(mod_df.loc[fam_mask, 'p_LR'].values, method='holm')[1]
    mod_df.loc[fam_mask, 'bh_p_LR'] = multipletests(mod_df.loc[fam_mask, 'p_LR'].values, method='fdr_bh')[1]
    mod_df.loc[fam_mask, 'holm_perm_p'] = multipletests(mod_df.loc[fam_mask, 'perm_p'].values, method='holm')[1]
    mod_df.to_csv(f'{P1}/hpa_ihc_dab_models.tsv', sep='\t', index=False)
    print('\n=== HPA SQUAMOUS IHC rgb2hed DAB TOBIT MODELS ===')
    print(mod_df[['feature', 'n_adc', 'n_targets', 'LR', 'p_LR', 'holm_p_LR', 'perm_p', 'holm_perm_p', 'dLR_vs_Wbasal']].round(4).to_string(index=False))

    run_param_grid(man, T)


if __name__ == '__main__':
    main()
