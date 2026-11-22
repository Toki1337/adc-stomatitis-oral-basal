"""Aggregate all manuscript statistics, CIs, PIs, and table metrics into results/numbers.json."""

import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.special import expit
from tobit_lib import design_matrix, fit_tobit

ROOT = Path(__file__).resolve().parents[1]
R = f'{ROOT}/results/'

clin = pd.read_csv(f'{ROOT}/clinical/adc_stomatitis_trials.tsv', sep='\t')
dual_ag = pd.read_csv(R + 'part1/clinical_dual_agreement.tsv', sep='\t').set_index('metric')
ds3 = pd.read_csv(R + 'part1/clinical_three_dataset_comparison.tsv', sep='\t')
dab_cov = pd.read_csv(R + 'part1/hpa_ihc_dab_coverage.tsv', sep='\t').set_index('tissue_group')
dab_mod = pd.read_csv(R + 'part1/hpa_ihc_dab_models.tsv', sep='\t').set_index('feature')
dab_gene = pd.read_csv(R + 'part1/hpa_ihc_dab_genes.tsv', sep='\t', index_col=0)
dab_grid = pd.read_csv(R + 'part1/hpa_ihc_dab_param_grid.tsv', sep='\t')
adc_feat = pd.read_csv(R + 'part1/adc_table_with_features.tsv', sep='\t')
feat_comp = pd.read_csv(R + 'part1/feature_comparison.tsv', sep='\t').set_index('feature')
res_ladder = pd.read_csv(R + 'part1/resolution_ladder.tsv', sep='\t').set_index('feature')
tgt_inf = pd.read_csv(R + 'part1/target_inference.tsv', sep='\t').iloc[0].to_dict()
sens = pd.read_csv(R + 'part1/sensitivity_subgroups.tsv', sep='\t').set_index('specification')
scan = pd.read_csv(R + 'part1/celltype_tissue_scan.tsv', sep='\t')
axes_comp = pd.read_csv(R + 'part3_axes_comparison.tsv', sep='\t').set_index('col')
joint_comp = pd.read_csv(R + 'part3_joint_models.tsv', sep='\t').set_index('ratio_col')
ratio_sens = pd.read_csv(R + 'part3_ratio_sensitivity.tsv', sep='\t')
calib_sum = pd.read_csv(R + 'part5_calib_summary.tsv', sep='\t')
loto_exact = pd.read_csv(R + 'part5_loto_exact_endpoints.tsv', sep='\t')
loto_df = pd.read_csv(R + 'part5_loto_predictions.tsv', sep='\t')
faers = pd.read_csv(R + 'fig_faers_summary.tsv', sep='\t').set_index('outcome')
faers_pairs = pd.read_csv(R + 'fig_faers_pairs.tsv', sep='\t')
faers_disprop = pd.read_csv(R + 'faers_per_drug_disprop.tsv', sep='\t')
faers_sens = pd.read_csv(R + 'faers_sensitivity.tsv', sep='\t')
faers_bg = pd.read_csv(R + 'faers_background.tsv', sep='\t').set_index('metric')['value'].to_dict()
faers_readus = pd.read_csv(R + 'faers_readus_baseline.tsv', sep='\t')
p4 = pd.read_csv(R + 'part4_final_ranking.tsv', sep='\t', index_col=0)
surf = pd.read_csv(R + 'surfaceome_genes.txt', header=None)[0].tolist()

eval_clin = clin[clin.listed.isin(['yes', 'not_listed']) & ~(clin.stom_any.isna() & clin.threshold_pct.isna())]

hpa_scan = scan[scan.source == 'HPA_singlecell'].sort_values('rank_all')
gtex_scan = scan[scan.source == 'GTEx'].sort_values('rank_all')

win49 = p4[p4.n_pass == 4]
win12 = win49[win49.vital_max_tpm <= 15.0].sort_values('tum_choi', ascending=False)
win_ext15 = win49[win49.ext_max_tpm <= 15.0].sort_values('tum_choi', ascending=False)
win_nb15 = win49[win49.nonbrain_ext_max_tpm <= 15.0].sort_values('tum_choi', ascending=False)
win_rob = win49[win49.robust_7of9].sort_values('tum_choi', ascending=False)
win_rob6 = win49[win49.grid_pass >= 6].sort_values('tum_choi', ascending=False)
win_ext6 = win49[win49.grid_pass_ext >= 6].sort_values('tum_choi', ascending=False)

# Target-level correlation between abs_oral_basal and ratio_oscc across 28 primary targets
p3_all = pd.read_csv(R + 'part3_axes.tsv', sep='\t').dropna(subset=['abs_oral_basal'])
main_fit = fit_tobit(p3_all.y_pct.values, p3_all.cens.values, design_matrix(p3_all, 'abs_oral_basal'))
b_main = main_fit['beta']
u28 = p3_all.drop_duplicates('tg').set_index('tg')
from scipy import stats
rho_abs_ratio, rho_abs_ratio_p = stats.spearmanr(u28.abs_oral_basal, u28.ratio_oscc)
rho_abs_ratio = float(rho_abs_ratio)
rho_abs_ratio_p = float(rho_abs_ratio_p)

nums = {
    'cohort': {
        'n_total_adcs': int(len(clin)),
        'n_unique_targets_all': int(clin.target_gene.nunique()),
        'n_evaluable_adcs': int(len(eval_clin)),
        'n_numeric_adcs': int((eval_clin.listed == 'yes').sum()),
        'n_censored_adcs': int((eval_clin.listed == 'not_listed').sum()),
        'n_excluded_adcs': int((clin.listed == 'no_data').sum()),
        'n_unique_targets_eval': int(adc_feat.target_gene.str.split(';').str[0].nunique()),
        'total_monotherapy_patients': int(eval_clin.N.sum()),
        'median_trial_N': int(eval_clin.N.median()),
        'min_trial_N': int(eval_clin.N.min()),
        'max_trial_N': int(eval_clin.N.max()),
        'min_numeric_stom_pct': float(eval_clin.stom_any.dropna().min()),
        'max_numeric_stom_pct': float(eval_clin.stom_any.dropna().max()),
    },
    'clinical_dual_agreement': dual_ag.reset_index().to_dict(orient='records'),
    'clinical_three_datasets': ds3.to_dict(orient='records'),
    'hpa_dab': {
        'n_sections_total': int(dab_cov['n_cores'].sum()),
        'coverage_by_tissue': dab_cov.reset_index().to_dict(orient='records'),
        'gene_epifrac_oral_median_pct': round(float(100 * dab_cov.loc['Oral mucosa', 'gene_epifrac_median']), 1),
        'gene_epifrac_esoph_median_pct': round(float(100 * dab_cov.loc['Esophagus', 'gene_epifrac_median']), 1),
        'wilcoxon_gene_epifrac_p': float(dab_cov.loc['Oral mucosa', 'wilcoxon_p_esoph_vs_oral_epifrac']),
        'ptprz1_esoph_dab': round(float(dab_gene.loc['PTPRZ1', 'esoph_dab_basal']), 3),
        'ptprz1_esoph_rank': int(dab_gene.loc['PTPRZ1', 'esoph_dab_rank_asc']),
        'low_rna_genes_n': int(dab_gene['low_basal_rna_lt005'].sum()),
        'low_rna_esoph_dab_median': round(float(dab_gene.loc[dab_gene['low_basal_rna_lt005'], 'esoph_dab_basal'].median()), 3),
        'param_grid_81': {
            'n_grid': int(len(dab_grid)),
            'esoph_n_p05': int((dab_grid['esoph_p'] < 0.05).sum()),
            'esoph_pct_p05': round(float(100 * (dab_grid['esoph_p'] < 0.05).mean()), 1),
            'esoph_n_p01': int((dab_grid['esoph_p'] < 0.01).sum()),
            'esoph_median_LR': round(float(dab_grid['esoph_LR'].median()), 2),
            'esoph_min_LR': round(float(dab_grid['esoph_LR'].min()), 2),
            'esoph_max_LR': round(float(dab_grid['esoph_LR'].max()), 2),
            'oral_n_p05': int((dab_grid['oral_p'] < 0.05).sum()),
            'oral_median_LR': round(float(dab_grid['oral_LR'].median()), 2),
            'oral_min_LR': round(float(dab_grid['oral_LR'].min()), 2),
            'oral_max_LR': round(float(dab_grid['oral_LR'].max()), 2),
        },
        'models': dab_mod.reset_index().to_dict(orient='records'),
        'key_genes': dab_gene.reindex([
            'TACSTD2', 'EGFR', 'ERBB2', 'ERBB3', 'CD276', 'NECTIN4', 'MET', 'ITGB6',
            'CDH3', 'PTK7', 'FAT2', 'PTPRZ1', 'COL17A1', 'DSC3', 'DSG3',
        ]).reset_index().rename(columns={'index': 'gene'}).to_dict(orient='records'),
    },
    'primary_model': {
        'W_basal_LR': round(float(feat_comp.loc['W_basal', 'LR']), 2),
        'W_basal_p': float(feat_comp.loc['W_basal', 'p']),
        'W_basal_beta': round(float(feat_comp.loc['W_basal', 'beta']), 3),
        'W_basal_se': round(float(feat_comp.loc['W_basal', 'se']), 3),
        'W_basal_beta_per_SD': round(float(feat_comp.loc['W_basal', 'beta_per_SD']), 3),
        'W_basal_spearman': round(float(feat_comp.loc['W_basal', 'spearman']), 3),
        'W_basal_p_spearman': float(feat_comp.loc['W_basal', 'p_spearman']),
        'W_basal_det_LR': round(float(feat_comp.loc['W_basal_det', 'LR']), 2),
        'W_basal_det_p': float(feat_comp.loc['W_basal_det', 'p']),
        'W_basal_donormed_LR': round(float(feat_comp.loc['W_basal_donormed', 'LR']), 2),
        'W_basal_quiescent_LR': round(float(feat_comp.loc['W_basal_quiescent', 'LR']), 2),
        'W_basal_prolif_LR': round(float(feat_comp.loc['W_basal_prolif', 'LR']), 2),
        'Choi_normal_epi_LR': round(float(feat_comp.loc['Choi_normal_epi', 'LR']), 2),
        'Choi_normal_epi_p': round(float(feat_comp.loc['Choi_normal_epi', 'p']), 4),
        'Choi_normal_epi_spearman': round(float(feat_comp.loc['Choi_normal_epi', 'spearman']), 3),
        'Choi_panHNSCC_LR': round(float(feat_comp.loc['Choi_panHNSCC_normal_epi', 'LR']), 2),
        'intercept_b0': round(float(b_main[0]), 3),
        'floor_x0_topoI_pct': round(float(expit(b_main[0]) * 100.0), 1),
        'floor_x0_tubulin_pct': round(float(expit(b_main[0] + b_main[1]) * 100.0), 1),
        'floor_x0_dna_pct': round(float(expit(b_main[0] + b_main[2]) * 100.0), 1),
    },
    'target_inference': {k: (float(v) if isinstance(v, (float, np.floating)) else int(v)) for k, v in tgt_inf.items()},
    'posthoc_4payload': {
        'basal_LR': round(float(tgt_inf['p4_basal_LR']), 2),
        'basal_beta': round(float(tgt_inf['p4_basal_beta']), 2),
        'basal_se': round(float(tgt_inf['p4_basal_se']), 2),
        'basal_perm_p': round(float(tgt_inf['p4_basal_perm_p']), 4),
        'payload_4v3_LR': round(float(tgt_inf['p4_vs_p3_LR']), 2),
        'payload_4v3_p': round(float(tgt_inf['p4_vs_p3_p']), 6),
    },
    'weighted_tobit': {
        'n_weight_beta': round(float(sens.loc['Sample-size weighted Tobit (weight = N / mean N)', 'beta']), 2),
        'n_weight_boot_ci': [
            round(float(sens.loc['Sample-size weighted Tobit (weight = N / mean N)', 'prof_ci_lo']), 2),
            round(float(sens.loc['Sample-size weighted Tobit (weight = N / mean N)', 'prof_ci_hi']), 2),
        ],
        'sqrtn_weight_beta': round(float(sens.loc['Sample-size weighted Tobit (weight = sqrt N)', 'beta']), 2),
        'sqrtn_weight_boot_ci': [
            round(float(sens.loc['Sample-size weighted Tobit (weight = sqrt N)', 'prof_ci_lo']), 2),
            round(float(sens.loc['Sample-size weighted Tobit (weight = sqrt N)', 'prof_ci_hi']), 2),
        ],
    },
    'n_outcome_heterogeneity': int((sens.reset_index().category == 'Outcome_heterogeneity').sum()),
    'drop_one_target': {
        'n_targets_tested': int((sens.reset_index().category == 'Drop_one_target').sum()),
        'min_LR': round(float(sens[sens.category == 'Drop_one_target']['LR'].min()), 2),
        'max_p_LR': round(float(sens[sens.category == 'Drop_one_target']['p_LR'].max()), 4),
        'worst_spec': str(sens[sens.category == 'Drop_one_target']['LR'].idxmin()),
        'alltissue_noTROP2_LR': round(float(sens.loc['Whole oral mucosa pseudobulk excl TACSTD2 (n=39)', 'LR']), 2),
        'alltissue_noTROP2_p': round(float(sens.loc['Whole oral mucosa pseudobulk excl TACSTD2 (n=39)', 'p_LR']), 3),
    },
    'resolution_ladder': res_ladder.reset_index().to_dict(orient='records'),
    'sensitivity_subgroups': sens.reset_index().to_dict(orient='records'),
    'specificity_scan': {
        'hpa_n_total': int(len(hpa_scan)),
        'hpa_n_q05': int((hpa_scan.q < 0.05).sum()),
        'hpa_n_q05_noTROP2': int((hpa_scan.q_noTROP2 < 0.05).sum()),
        'hpa_min_q_noTROP2': round(float(hpa_scan.q_noTROP2.min()), 4),
        'hpa_top5': hpa_scan[['celltype', 'LR', 'p', 'q', 'LR_noTROP2', 'q_noTROP2', 'rank_all', 'rank_noTROP2']].head(7).to_dict(orient='records'),
        'gtex_n_total': int(len(gtex_scan)),
        'gtex_n_q05': int((gtex_scan.q < 0.05).sum()),
        'gtex_n_q05_noTROP2': int((gtex_scan.q_noTROP2 < 0.05).sum()),
        'gtex_min_q_noTROP2': round(float(gtex_scan.q_noTROP2.min()), 4),
        'gtex_top6': gtex_scan[['celltype', 'LR', 'p', 'q', 'LR_noTROP2', 'q_noTROP2', 'rank_all', 'rank_noTROP2']].head(6).to_dict(orient='records'),
    },
    'axes_comparison': axes_comp.reset_index().to_dict(orient='records'),
    'joint_models': joint_comp.reset_index().to_dict(orient='records'),
    'ratio_sensitivity': ratio_sens.to_dict(orient='records'),
    'orthogonality_28targets': {
        'spearman_abs_vs_ratio': round(rho_abs_ratio, 3),
        'spearman_abs_vs_ratio_p': round(rho_abs_ratio_p, 4),
        'TACSTD2_abs_rank': int(u28.abs_oral_basal.rank(ascending=False)['TACSTD2']),
        'TACSTD2_ratio_rank': int(u28.ratio_oscc.rank(ascending=False)['TACSTD2']),
    },
    'calibration': calib_sum.to_dict(orient='records'),
    'loto_exact_endpoints': loto_exact.to_dict(orient='records'),
    'loto_numeric_predictions': loto_df[~loto_df.cens][
        ['drug', 'target_gene', 'y_pct', 'pred_loto']
    ].to_dict(orient='records'),
    'loto_censored_misses': loto_df[loto_df.cens & (loto_df.pred_loto >= loto_df.y_pct)][
        ['drug', 'tg', 'pclass', 'y_pct', 'pred_loto', 'pred_loto_4payload']
    ].to_dict(orient='records'),
    'faers': faers.reset_index().to_dict(orient='records'),
    'faers_pairs': faers_pairs.to_dict(orient='records'),
    'faers_per_drug_disprop': faers_disprop.to_dict(orient='records'),
    'faers_sensitivity': faers_sens.to_dict(orient='records'),
    'faers_background': faers_bg,
    'faers_readus_baseline': faers_readus.to_dict(orient='records'),
    'surfaceome': {
        'n_surfaceome': int(len(surf)),
        'n_oscc_expressed': int(len(p4)),
        'spearman_choi_basal': round(float(stats.spearmanr(p4.tum_choi, p4.oral_basal)[0]), 3),
        'spearman_puram_basal': round(float(stats.spearmanr(p4.tum_puram, p4.oral_basal)[0]), 3),
        'n_dlog2_ge1_choi': int((p4.win_choi >= 1).sum()),
        'pct_dlog2_ge1_choi': round(float(100 * (p4.win_choi >= 1).mean()), 1),
        'n_dlog2_ge1_dual': int(((p4.win_choi >= 1) & (p4.win_puram >= 0.5)).sum()),
        'pct_dlog2_ge1_dual': round(float(100 * ((p4.win_choi >= 1) & (p4.win_puram >= 0.5)).mean()), 1),
        'n_dlog2_ge2_choi': int((p4.win_choi >= 2).sum()),
        'pct_dlog2_ge2_choi': round(float(100 * (p4.win_choi >= 2).mean()), 1),
        'n_lin_ge2_choi_raw': int(((p4.tum_choi_raw / np.maximum(p4.oral_basal_raw, 1e-12)) >= 2).sum()),
        'pct_lin_ge2_choi_raw': round(float(100 * ((p4.tum_choi_raw / np.maximum(p4.oral_basal_raw, 1e-12)) >= 2).mean()), 1),
        'n_lin_ge2_choi_01': int((p4.lin_fc_choi_01 >= 2).sum()),
        'pct_lin_ge2_choi_01': round(float(100 * (p4.lin_fc_choi_01 >= 2).mean()), 1),
        'n_lin_ge2_dual_01': int(((p4.lin_fc_choi_01 >= 2) & (p4.lin_fc_puram_01 >= 2)).sum()),
        'pct_lin_ge2_dual_01': round(float(100 * ((p4.lin_fc_choi_01 >= 2) & (p4.lin_fc_puram_01 >= 2)).mean()), 1),
        'median_pred_stom_pct': round(float(p4.pred_stom_pct.median()), 1),
        'iqr_pred_stom_pct': [round(float(v), 1) for v in np.percentile(p4.pred_stom_pct, [25, 75])],
        'n_extrapolated_train_max': int(p4.extrapolated_above_train_max.sum()),
        'n_extrapolated_winners': int(win49.extrapolated_above_train_max.sum()),
        'n_winners_4of4': int(len(win49)),
        'n_winners_vital15': int(len(win12)),
        'winners_vital15_genes': win12.index.tolist(),
        'n_winners_nonbrain_ext15': int(len(win_nb15)),
        'winners_nonbrain_ext15_genes': win_nb15.index.tolist(),
        'n_winners_ext15': int(len(win_ext15)),
        'winners_ext15_genes': win_ext15.index.tolist(),
        'n_robust_7of9': int(len(win_rob)),
        'robust_7of9_genes': win_rob.index.tolist(),
        'n_robust_6of9': int(len(win_rob6)),
        'robust_6of9_genes': win_rob6.index.tolist(),
        'n_ext_grid_6of9': int(len(win_ext6)),
        'ext_grid_6of9_genes': win_ext6.index.tolist(),
    },
    'table3_targets': p4.reindex([
        'CDH3', 'PTPRZ1', 'PTK7', 'FAT2', 'FZD6', 'DSG2', 'PTPRK', 'ITGA2',
        'CDCP1', 'IL1RAP', 'IGSF3', 'ATP2B1', 'TNFRSF12A', 'F3',
        'TACSTD2', 'COL17A1', 'DSC3', 'DSG3', 'EGFR', 'CD44', 'GPNMB',
        'ERBB2', 'ERBB3', 'CD276', 'MET', 'ITGB6', 'AXL', 'NECTIN4',
    ]).reset_index().rename(columns={'index': 'gene'}).to_dict(orient='records'),
}

with open(R + 'numbers.json', 'w', encoding='utf-8') as f:
    json.dump(nums, f, indent=2)
print('Wrote', R + 'numbers.json')
