suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggrepel)
  library(cowplot)
  library(scales)
})

args_cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args_cmd, value = TRUE)
ROOT <- if (length(file_arg) > 0) normalizePath(file.path(dirname(sub("^--file=", "", file_arg[1])), "..")) else getwd()
R <- file.path(ROOT, "results", "")
OUT <- file.path(ROOT, "figures", "")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 7.5) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title.position = "plot",
      plot.title = element_text(face = "bold", size = base_size + 0.6, hjust = 0, margin = margin(l = 12, b = 2.5)),
      plot.subtitle = element_text(size = base_size - 1.0, color = "grey30", hjust = 0, margin = margin(l = 12, b = 4)),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "grey15", size = base_size - 0.8),
      legend.title = element_text(face = "bold", size = base_size - 0.6),
      legend.text = element_text(size = base_size - 0.9),
      legend.key.size = unit(0.28, "cm"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.5),
      plot.margin = margin(4, 7, 4, 5)
    )
}

pal_pclass <- c("topoI" = "#D55E00", "DNA" = "#E69F00", "tubulin" = "#0072B2")
pal_pclass_lab <- c("topoI" = "Topo-I inhibitor", "DNA" = "DNA-damaging", "tubulin" = "Tubulin inhibitor")

save_fig <- function(plt, name, w, h) {
  ggsave(paste0(OUT, name, ".pdf"), plt, width = w, height = h, device = cairo_pdf)
  ggsave(paste0(OUT, name, ".png"), plt, width = w, height = h, dpi = 300, bg = "white")
  message("Saved ", name)
}

fmt_p <- function(p) {
  sapply(p, function(v) {
    if (is.na(v)) return("NA")
    if (v < 0.001) {
      ex <- floor(log10(v))
      mant <- v / (10^ex)
      sup_digits <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴",
                      "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹")
      ex_str <- paste0("⁻", paste0(sup_digits[strsplit(as.character(abs(ex)), "")[[1]]], collapse = ""))
      sprintf("%.2f × 10%s", mant, ex_str)
    } else if (v < 0.01) {
      sprintf("%.4f", v)
    } else {
      sprintf("%.3f", v)
    }
  })
}

short_adc <- function(drug, tg) {
  s <- dplyr::recode(drug,
    "trastuzumab deruxtecan" = "T-DXd",
    "trastuzumab emtansine" = "T-DM1",
    "trastuzumab rezetecan" = "T-rez",
    "datopotamab deruxtecan" = "Dato-DXd",
    "sacituzumab govitecan" = "SG",
    "sacituzumab tirumotecan" = "Sac-TMT",
    "patritumab deruxtecan" = "HER3-DXd",
    "izalontamab brengitecan" = "Iza-bren",
    "vobramitamab duocarmazine" = "Vobra-duo",
    "glembatumumab vedotin" = "Glemba-V",
    "gemtuzumab ozogamicin" = "GO",
    "inotuzumab ozogamicin" = "InO",
    "tisotumab vedotin" = "Tisotumab-V",
    "sigvotatug vedotin" = "Sigvotatug-V",
    "becotatug vedotin" = "Becotatug-V",
    "cofetuzumab pelidotin" = "Cofetuzumab-P",
    "ifinatamab deruxtecan" = "I-DXd",
    .default = gsub(" .*", "", drug)
  )
  paste0(s, " (", tg, ")")
}

# ==============================================================================
# FIGURE 1: Study architecture, oral mucosal compartments, & clinical landscape
# ==============================================================================
clin_all <- read_tsv(file.path(ROOT, "clinical", "adc_stomatitis_trials.tsv"), show_col_types = FALSE)
p3 <- read_tsv(paste0(R, "part3_axes.tsv"), show_col_types = FALSE)
epi_mk <- read_tsv(paste0(R, "fig_epi_markers.tsv"), show_col_types = FALSE)

n_num_all <- sum(clin_all$listed == "yes", na.rm = TRUE)
n_cen_all <- sum(clin_all$listed == "not_listed", na.rm = TRUE)
n_nod_all <- sum(!(clin_all$listed %in% c("yes", "not_listed")), na.rm = TRUE)
n_eval_all <- n_num_all + n_cen_all
lbl_num_all <- sprintf("Numeric rate (n=%d)", n_num_all)
lbl_cen_all <- sprintf("Left-censored <threshold (n=%d)", n_cen_all)
lbl_nod_all <- sprintf("No monotherapy data (n=%d)", n_nod_all)

clin_sum <- clin_all %>%
  mutate(
    status = case_when(
      listed == "yes" ~ lbl_num_all,
      listed == "not_listed" ~ lbl_cen_all,
      TRUE ~ lbl_nod_all
    ),
    status = factor(status, levels = c(lbl_num_all, lbl_cen_all, lbl_nod_all)),
    src = case_when(
      population_source_type == "label" ~ "FDA label",
      population_source_type == "phase3" ~ "Phase 3",
      population_source_type == "phase2" ~ "Phase 2",
      population_source_type == "phase1" ~ "Phase 1",
      TRUE ~ "Unavailable"
    ),
    src = factor(src, levels = c("FDA label", "Phase 3", "Phase 2", "Phase 1", "Unavailable"))
  )

clin_totals <- clin_sum %>% count(status)
p1a <- clin_sum %>%
  count(status, src) %>%
  ggplot(aes(x = status, y = n, fill = src)) +
  geom_col(width = 0.60, color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(n >= 2, n, "")), position = position_stack(vjust = 0.5), size = 2.25, color = "white", fontface = "bold") +
  geom_text(data = clin_totals, aes(x = status, y = n, label = paste0("n=", n)), inherit.aes = FALSE, vjust = -0.45, size = 2.45, fontface = "bold") +
  scale_fill_manual(values = c("FDA label" = "#1B4F72", "Phase 3" = "#2874A6", "Phase 2" = "#5DADE2", "Phase 1" = "#AED6F1", "Unavailable" = "grey75"),
                    name = NULL, na.translate = FALSE) +
  guides(fill = guide_legend(nrow = 2)) +
  scale_y_continuous(limits = c(0, 35), expand = expansion(mult = c(0, 0.05))) +
  scale_x_discrete(labels = c(sprintf("Numeric\n(n=%d)", n_num_all), sprintf("Left-censored\n(n=%d)", n_cen_all), sprintf("Excluded\n(n=%d)", n_nod_all))) +
  labs(title = sprintf("Blinded dual-curated cohort (%d ADCs)", nrow(clin_all)),
       subtitle = sprintf("%d monotherapy ADCs enter Tobit regression", n_eval_all),
       x = NULL, y = "Number of ADCs") +
  theme_pub(base_size = 7.4) +
  theme(legend.position = "top",
        legend.justification = "left",
        legend.margin = margin(0, 0, -3, 0),
        legend.key.size = unit(0.24, "cm"),
        legend.text = element_text(size = 6.2))

gene_order <- c("KRT5", "KRT14", "KRT15", "COL17A1", "TP63", "MKI67",
                "KRT13", "KRT4", "SPRR3", "CRNN",
                "TACSTD2", "CD44", "GPNMB", "EGFR", "ERBB2", "PTK7", "NECTIN4", "F3", "FOLR1")
epi_mk <- epi_mk %>%
  mutate(
    gene = factor(gene, levels = rev(gene_order)),
    compartment = factor(compartment, levels = c("basal", "suprabasal", "superficial"),
                         labels = c("Basal\n(n=2,876)", "Suprabasal\n(n=1,700)", "Superficial\n(n=968)"))
  )

p1b <- ggplot(epi_mk, aes(x = compartment, y = gene, size = det_frac * 100, color = mean_log2)) +
  geom_point() +
  geom_hline(yintercept = 9.5, linetype = "dashed", color = "grey60", linewidth = 0.4) +
  scale_color_gradientn(colors = c("#F7F7F7", "#FDB863", "#E66101", "#7F0000"), name = "Mean\nlog₂(CP10K+1)",
                        guide = guide_colorbar(barheight = unit(1.5, "cm"), barwidth = unit(0.25, "cm"))) +
  scale_size_continuous(range = c(0.8, 4.0), name = "Detected\ncells (%)", breaks = c(10, 40, 70, 95)) +
  labs(title = "Oral mucosal epithelial stratification",
       subtitle = "Healthy buccal & gingival scRNA-seq (21 donors)",
       x = NULL, y = NULL) +
  theme_pub(base_size = 7.4) +
  theme(axis.text.y = element_text(face = "italic", size = 6.4),
        axis.text.x = element_text(size = 6.4),
        legend.title = element_text(size = 6.3),
        legend.text = element_text(size = 6.1),
        legend.spacing.y = unit(0.12, "cm"),
        panel.grid.major = element_line(color = "grey94", linewidth = 0.25))

p1c <- ggplot(p3, aes(x = reorder(paste0(drug, "  [", target_gene, "]"), ifelse(cens, y_pct - 0.1, y_pct)), y = y_pct, color = pclass, shape = cens)) +
  geom_segment(aes(xend = reorder(paste0(drug, "  [", target_gene, "]"), ifelse(cens, y_pct - 0.1, y_pct)), y = 0, yend = y_pct),
               alpha = 0.35, linewidth = 0.45) +
  geom_point(size = 1.75) +
  coord_flip() +
  scale_color_manual(values = pal_pclass, labels = pal_pclass_lab, name = "Payload class") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), labels = c(`FALSE` = "Observed %", `TRUE` = "Censored (< threshold)"), name = "Endpoint type") +
  scale_y_continuous(limits = c(0, 70), expand = expansion(mult = c(0.01, 0.04))) +
  labs(title = "Clinical stomatitis incidence across 44 ADCs",
       subtitle = "Triangles denote upper reporting thresholds (not listed in AE table)",
       x = NULL, y = "Any-grade stomatitis / oral mucositis (%)") +
  theme_pub(base_size = 7.4) +
  theme(axis.text.y = element_text(size = 5.9),
        legend.position = c(0.73, 0.22),
        legend.title = element_text(size = 6.3),
        legend.text = element_text(size = 6.0),
        legend.key.size = unit(0.24, "cm"),
        legend.spacing.y = unit(0.08, "cm"),
        legend.box.background = element_rect(fill = "white", color = "grey85", linewidth = 0.3))

fig1_left <- plot_grid(p1a, NULL, p1b, ncol = 1, rel_heights = c(0.84, 0.035, 1.32), labels = c("A", "", "B"), label_size = 10)
fig1 <- plot_grid(fig1_left, p1c, ncol = 2, rel_widths = c(0.95, 1.15), labels = c("", "C"), label_size = 10)
save_fig(fig1, "Fig1_cohort_and_compartments", 7.09, 7.8)


# ==============================================================================
# FIGURE 2: Head-to-head — Absolute basal expression vs Tumor/Normal ratio
# ==============================================================================
a3 <- read_tsv(paste0(R, "part3_axes_comparison.tsv"), show_col_types = FALSE)
curves <- read_tsv(paste0(R, "fig_tobit_curves.tsv"), show_col_types = FALSE)
tgt_inf <- read_tsv(paste0(R, "part1/target_inference.tsv"), show_col_types = FALSE)
a3_abs <- a3 %>% filter(col == "abs_oral_basal")
a3_rat <- a3 %>% filter(col == "ratio_oscc")

a3_sub <- a3 %>%
  filter(col %in% c("abs_oral_basal", "tum_oscc", "ratio_puram", "ratio_oscc", "ratio_choi_internal", "ratio_choi_pc01", "tcga_paired_lfc")) %>%
  mutate(
    short_lab = recode(col,
      "abs_oral_basal" = "Absolute oral basal\n(Williams scRNA)",
      "tum_oscc" = "OSCC tumor alone\n(Choi malignant)",
      "ratio_puram" = "Ratio: Puram OSCC /\noral basal (log₂+1)",
      "ratio_oscc" = "Ratio: Choi OSCC /\noral basal (log₂+1)",
      "ratio_choi_internal" = "Ratio: Choi OSCC /\nChoi normal epi",
      "ratio_choi_pc01" = "Ratio: Choi / basal\n(pseudocount 0.1)",
      "tcga_paired_lfc" = "TCGA oral paired LFC\n(30 tumor/normal)"
    ),
    grp = ifelse(col == "abs_oral_basal", "Absolute", ifelse(col == "tum_oscc", "Tumor only", "Selectivity ratio")),
    short_lab = factor(short_lab, levels = rev(short_lab))
  )

p2a <- ggplot(a3_sub, aes(x = short_lab, y = LR, fill = grp)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_label(aes(label = sprintf("LR=%.1f (P=%s, β=%+.2f)", LR, fmt_p(p), beta)),
             hjust = -0.03, size = 2.0, color = "grey10", fill = "white",
             linewidth = 0, label.padding = unit(0.05, "lines"), show.legend = FALSE) +
  coord_flip() +
  scale_fill_manual(values = c("Absolute" = "#0072B2", "Tumor only" = "#56B4E9", "Selectivity ratio" = "grey65"), name = NULL) +
  scale_x_discrete(expand = expansion(add = c(0.55, 0.85))) +
  scale_y_continuous(limits = c(0, 32.5), breaks = c(0, 5, 10, 15, 20), expand = expansion(mult = c(0, 0.02))) +
  annotate("text", x = 7.42, y = 4.2, label = "P=0.05", size = 2.0, color = "grey40", hjust = 0, vjust = 0) +
  labs(title = "Head-to-head Tobit likelihood ratio (n=44 ADCs)",
       subtitle = "Same outcome, drugs, and payload covariates; only X swapped",
       x = NULL, y = "Likelihood-ratio statistic (df = 1)") +
  theme_pub(base_size = 7.2) +
  theme(axis.text.y = element_text(size = 6.0),
        legend.position = "top",
        legend.justification = "left",
        legend.key.size = unit(0.24, "cm"),
        legend.margin = margin(b = -4))

b_drugs <- c("datopotamab deruxtecan", "SHR-A1921", "sacituzumab tirumotecan",
             "sacituzumab govitecan", "izalontamab brengitecan",
             "glembatumumab vedotin", "gemtuzumab ozogamicin",
             "sigvotatug vedotin", "trastuzumab deruxtecan",
             "trastuzumab rezetecan", "tisotumab vedotin", "ESG401")
p3_2b <- p3 %>% mutate(
  lbl = ifelse(drug %in% b_drugs, short_adc(drug, target_gene), ""),
  nx = case_when(
    drug == "trastuzumab deruxtecan" ~ 0.25,
    drug == "trastuzumab rezetecan" ~ 0.42,
    drug == "tisotumab vedotin" ~ 0.10,
    drug == "glembatumumab vedotin" ~ 0.18,
    drug == "ESG401" ~ -0.32,
    drug == "sacituzumab govitecan" ~ -0.05,
    TRUE ~ 0
  ),
  ny = case_when(
    drug == "trastuzumab deruxtecan" ~ 12.5,
    drug == "trastuzumab rezetecan" ~ -8.5,
    drug == "tisotumab vedotin" ~ -11.5,
    drug == "glembatumumab vedotin" ~ -12.5,
    drug == "ESG401" ~ -6.0,
    drug == "sacituzumab govitecan" ~ -11.5,
    TRUE ~ 0
  )
)

p2b <- ggplot(p3_2b, aes(x = abs_oral_basal, y = y_pct, color = pclass)) +
  geom_line(data = curves, aes(x = abs_oral_basal, y = pred_pct, color = pclass), linewidth = 0.7, alpha = 0.85) +
  geom_point(aes(shape = cens), size = 1.8, alpha = 0.9) +
  geom_text_repel(aes(label = lbl),
                  nudge_x = p3_2b$nx, nudge_y = p3_2b$ny,
                  size = 1.98, max.overlaps = Inf, box.padding = 0.34, point.padding = 0.22,
                  force = 12, force_pull = 0.35, max.time = 5, max.iter = 40000, seed = 101,
                  bg.color = "white", bg.r = 0.16,
                  min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(values = pal_pclass, labels = pal_pclass_lab, name = "Payload") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), labels = c(`FALSE` = "Observed", `TRUE` = "Censored (<)"), name = "Status") +
  scale_y_continuous(limits = c(-6.5, 76), breaks = c(0, 20, 40, 60)) +
  labs(title = "Absolute oral basal expression predicts stomatitis",
       subtitle = sprintf("LR = %.2f (target-perm P = %.4f); Spearman ρ = %.3f (P = %s)",
                          a3_abs$LR, tgt_inf$perm_p[1], a3_abs$spearman, fmt_p(a3_abs$p_spearman)),
       x = "Oral basal keratinocyte expression, log₂(CP10K + 1)",
       y = "Stomatitis incidence / threshold (%)") +
  theme_pub(base_size = 7.2) +
  theme(legend.position = c(0.32, 0.81),
        legend.box = "horizontal",
        legend.title = element_text(size = 6.1),
        legend.text = element_text(size = 5.8),
        legend.key.size = unit(0.22, "cm"),
        legend.spacing.x = unit(0.08, "cm"),
        legend.background = element_rect(fill = alpha("white", 0.85), color = NA))

c_drugs <- c("datopotamab deruxtecan", "SHR-A1921", "sacituzumab govitecan",
             "izalontamab brengitecan", "sigvotatug vedotin", "tisotumab vedotin",
             "cofetuzumab pelidotin", "gemtuzumab ozogamicin")
p3_2c <- p3 %>% mutate(
  lbl = ifelse(drug %in% c_drugs, short_adc(drug, target_gene), ""),
  nx = case_when(
    drug == "izalontamab brengitecan" ~ 0.34,
    drug == "gemtuzumab ozogamicin" ~ -0.24,
    drug == "sigvotatug vedotin" ~ -0.02,
    TRUE ~ 0
  ),
  ny = case_when(
    drug == "izalontamab brengitecan" ~ 3.5,
    drug == "gemtuzumab ozogamicin" ~ 8.5,
    drug == "sigvotatug vedotin" ~ 10.5,
    TRUE ~ 0
  )
)

p2c <- ggplot(p3_2c, aes(x = ratio_oscc, y = y_pct, color = pclass)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(shape = cens), size = 1.8, alpha = 0.9) +
  geom_smooth(aes(group = 1), method = "lm", se = FALSE, color = "grey45", linetype = "dashed", linewidth = 0.55) +
  geom_text_repel(aes(label = lbl),
                  nudge_x = p3_2c$nx, nudge_y = p3_2c$ny,
                  size = 1.98, max.overlaps = Inf, box.padding = 0.32, point.padding = 0.22,
                  force = 10, force_pull = 0.35, max.time = 5, max.iter = 40000, seed = 42,
                  bg.color = "white", bg.r = 0.16,
                  min.segment.length = 0, segment.size = 0.3, segment.alpha = 0.6, show.legend = FALSE) +
  scale_color_manual(values = pal_pclass, guide = "none") +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 2), guide = "none") +
  scale_y_continuous(limits = c(-3, 74), breaks = c(0, 20, 40, 60)) +
  labs(title = "Tumor / normal ratio fails protective prediction",
       subtitle = sprintf("Choi ratio LR = %.2f (P = %s); non-protective positive slope", a3_rat$LR, fmt_p(a3_rat$p)),
       x = "Tumor / normal log₂ ratio (Choi OSCC − oral basal)",
       y = "Stomatitis incidence / threshold (%)") +
  theme_pub(base_size = 7.2)

u28 <- p3 %>%
  distinct(tg, .keep_all = TRUE) %>%
  mutate(
    tg_disp = ifelse(tg == "PVRL4", "NECTIN4", tg),
    highlight = case_when(
      tg == "TACSTD2" ~ "TROP2 (worst abs 1/28,\nselectivity rank 2/28)",
      tg %in% c("F3", "PTK7", "ITGB6", "CD276") ~ "Lower-basal selective targets\n(0.07–0.36 log₂ units)",
      TRUE ~ "Other ADC targets"
    ),
    lbl_2d = ifelse(abs_oral_basal > 0.3 | abs(ratio_oscc) > 0.2 | tg %in% c("F3", "ITGB6", "CD276"), tg_disp, "")
  )
rho_orth <- cor(u28$abs_oral_basal, u28$ratio_oscc, method = "spearman")

p2d <- ggplot(u28, aes(x = abs_oral_basal, y = ratio_oscc)) +
  geom_hline(yintercept = 0, linetype = "dotted", color = "grey55") +
  geom_point(aes(color = highlight), size = 2.0) +
  geom_text_repel(aes(label = lbl_2d, color = highlight), size = 2.1, max.overlaps = Inf,
                  box.padding = 0.25, point.padding = 0.18, seed = 42,
                  bg.color = "white", bg.r = 0.15, show.legend = FALSE) +
  scale_color_manual(values = c("TROP2 (worst abs 1/28,\nselectivity rank 2/28)" = "#D55E00",
                                "Lower-basal selective targets\n(0.07–0.36 log₂ units)" = "#009E73",
                                "Other ADC targets" = "grey40"), name = NULL) +
  labs(title = "Uncorrelated axes across 28 primary ADC targets",
       subtitle = sprintf("Spearman ρ(absolute axis, ratio axis) = %.3f", rho_orth),
       x = "Absolute oral basal expression, log₂(CP10K + 1)",
       y = "Tumor / normal log₂ ratio (Choi OSCC − oral basal)") +
  theme_pub(base_size = 7.2) +
  theme(legend.position = c(0.66, 0.84),
        legend.text = element_text(size = 5.9),
        legend.key.size = unit(0.24, "cm"))

fig2_top <- plot_grid(p2a, p2b, ncol = 2, labels = c("A", "B"), label_size = 10)
fig2_bot <- plot_grid(p2c, p2d, ncol = 2, labels = c("C", "D"), label_size = 10)
fig2 <- plot_grid(fig2_top, NULL, fig2_bot, ncol = 1, rel_heights = c(1, 0.035, 1))
save_fig(fig2, "Fig2_absolute_vs_ratio", 7.09, 6.7)


# ==============================================================================
# FIGURE 3: Resolution ladder, FDR-adjusted scans, profile-CI robustness & nested joint test
# ==============================================================================
res_lad <- read_tsv(paste0(R, "part1/resolution_ladder.tsv"), show_col_types = FALSE)
scan_df <- read_tsv(paste0(R, "part1/celltype_tissue_scan.tsv"), show_col_types = FALSE)
sens_df <- read_tsv(paste0(R, "part1/sensitivity_subgroups.tsv"), show_col_types = FALSE)
joint_df <- read_tsv(paste0(R, "part3_joint_models.tsv"), show_col_types = FALSE)

# Panel A: Squamous vs non-squamous / IHC resolution ladder (sorted by LR)
res_p <- res_lad %>%
  arrange(desc(LR)) %>%
  mutate(
    mod_grp = case_when(
      grepl("Single-cell|Unstratified", modality) ~ "Single-cell RNA",
      grepl("Bulk", modality) ~ "Bulk RNA",
      TRUE ~ "Protein (NAT / IHC)"
    ),
    label = factor(label, levels = rev(label))
  )

p3a <- ggplot(res_p, aes(x = label, y = LR, fill = mod_grp)) +
  geom_col(width = 0.66) +
  geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.1f", LR)), hjust = -0.15, size = 2.0) +
  coord_flip() +
  scale_fill_manual(values = c("Single-cell RNA" = "#0072B2", "Bulk RNA" = "#56B4E9", "Protein (NAT / IHC)" = "#009E73"), name = NULL) +
  scale_y_continuous(limits = c(0, 22.0), expand = expansion(mult = c(0, 0.03))) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(title = "Squamous resolution & modality ladder",
       subtitle = "Squamous RNA, NAT proteome & esoph. DAB predict",
       x = NULL, y = "Tobit LR statistic (df = 1)") +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 6.1),
        axis.text.y = element_text(size = 5.9),
        legend.position = "bottom",
        legend.justification = "left",
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 6.0),
        legend.margin = margin(t = -4))

# Panel B: HPA 155 single-cell scan (all q < 0.05 bars shown)
hpa_all <- scan_df %>% filter(source == "HPA_singlecell")
n_hpa_sig <- sum(hpa_all$q < 0.05, na.rm = TRUE)
n_hpa_sig_no <- sum(hpa_all$q_noTROP2 < 0.05, na.rm = TRUE)
min_hpa_q_no <- min(hpa_all$q_noTROP2, na.rm = TRUE)
sq_ocular_cts <- c(
  "esophageal basal cells",
  "suprabasal keratinocytes",
  "esophageal suprabasal cells",
  "ocular epithelial cells",
  "basal keratinocytes",
  "esophageal apical cells"
)
hpa_top <- hpa_all %>%
  arrange(desc(LR)) %>%
  slice_head(n = max(14, n_hpa_sig)) %>%
  mutate(
    grp_ct = case_when(
      celltype %in% sq_ocular_cts ~ "Stratified squamous / ocular (q < 0.05)",
      celltype == "urothelial cells" ~ "Urothelial epithelium (#2, q < 0.05)",
      TRUE ~ "Other epithelial (q < 0.05)"
    ),
    celltype = factor(celltype, levels = rev(celltype))
  )

p3b <- ggplot(hpa_top, aes(x = celltype, y = LR, fill = grp_ct)) +
  geom_col(width = 0.68) +
  geom_point(aes(y = LR_noTROP2), shape = 18, size = 1.65, color = "black") +
  coord_flip() +
  scale_fill_manual(values = c(
    "Stratified squamous / ocular (q < 0.05)" = "#D55E00",
    "Urothelial epithelium (#2, q < 0.05)" = "#0072B2",
    "Other epithelial (q < 0.05)" = "grey65"
  ), name = NULL) +
  scale_y_continuous(limits = c(0, 20.5)) +
  guides(fill = guide_legend(ncol = 1)) +
  labs(title = sprintf("Unbiased scan: %d HPA cell types", nrow(hpa_all)),
       subtitle = sprintf("%d pass q < 0.05; no TROP2: %d pass (min q = %.3f)", n_hpa_sig, n_hpa_sig_no, min_hpa_q_no),
       x = NULL, y = "Tobit LR statistic") +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 6.1),
        axis.text.y = element_text(size = 5.9),
        legend.position = "bottom",
        legend.justification = "left",
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 5.9),
        legend.margin = margin(t = -4))

# Panel C: GTEx 54-tissue scan (BH-FDR q < 0.05)
gtex_all <- scan_df %>% filter(source == "GTEx")
n_gtex_sig <- sum(gtex_all$q < 0.05, na.rm = TRUE)
min_gtex_q_no <- min(gtex_all$q_noTROP2, na.rm = TRUE)
sig_lbl_gtex <- sprintf("BH-FDR q < 0.05 (n=%d/%d)", n_gtex_sig, nrow(gtex_all))
gtex_top <- gtex_all %>%
  arrange(desc(LR)) %>%
  slice_head(n = 12) %>%
  mutate(
    fdr_sig = ifelse(q < 0.05, sig_lbl_gtex, "q ≥ 0.05"),
    celltype = factor(celltype, levels = rev(celltype))
  )

p3c <- ggplot(gtex_top, aes(x = celltype, y = LR, fill = fdr_sig)) +
  geom_col(width = 0.68) +
  geom_point(aes(y = LR_noTROP2), shape = 18, size = 1.65, color = "black") +
  geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", color = "grey45", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = setNames(c("#D55E00", "grey70"), c(sig_lbl_gtex, "q ≥ 0.05")), name = NULL) +
  scale_y_continuous(limits = c(0, 14.5)) +
  guides(fill = guide_legend(nrow = 2)) +
  labs(title = sprintf("Unbiased scan: %d GTEx bulk tissues", nrow(gtex_all)),
       subtitle = sprintf("%d pass q < 0.05; no TROP2: esoph. #1 (q = %.2f)", n_gtex_sig, min_gtex_q_no),
       x = NULL, y = "Tobit LR statistic") +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 6.1),
        axis.text.y = element_text(size = 6.0),
        legend.position = "bottom",
        legend.justification = "left",
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 6.0),
        legend.margin = margin(t = -4))

# Panel D: Profile-likelihood & target-cluster bootstrap 95% CI sensitivity & outcome heterogeneity
drop_specs <- c(
  "Censored set to threshold", "Censored set to 0",
  "Basal excl sub 13 (2,807 cells)",
  "All epithelial incl other (6,157 cells)"
)
sens_keep <- sens_df %>%
  filter(!specification %in% drop_specs,
         !category %in% c("Imputation", "Data_source_sensitivity", "Drop_one_target", "Resolution_noTROP2")) %>%
  mutate(
    spec_disp = specification %>%
      gsub("tumours", "tumors", .) %>%
      gsub("Sample-size weighted Tobit \\(weight = N / mean N\\)", "Weighted Tobit (w = N / mean N)", .) %>%
      gsub("Sample-size weighted Tobit \\(weight = sqrt N\\)", "Weighted Tobit (w = √N)", .) %>%
      gsub("Single stomatitis PT numeric \\+ censored \\(n=35; 6 num\\)", "Single stomatitis PT (6 num)", .) %>%
      gsub("Single stomatitis PT numeric \\+ censored", "Single stomatitis PT + censored", .) %>%
      gsub("4-level payload post hoc \\(split SN-38 vs DXd/other\\)", "4-level payload (split SN-38)", .) %>%
      gsub("Exclude small trials \\(N < 30; excl SHR-A1921, aprutumab\\)", "Trials with N ≥ 30", .) %>%
      gsub("log2\\(CP10K \\+ 0\\.1\\)", "log₂(CP10K + 0.1)", .) %>%
      gsub(" \\((n=\\d+|all \\d+)\\)", "", .),
    spec_lab = case_when(
      ci_type == "cluster_bootstrap_95CI" ~ sprintf("%s (n=%d, boot CI)", spec_disp, n),
      !is.na(p_target_perm) ~ sprintf("%s (n=%d, perm P=%.4f)", spec_disp, n, p_target_perm),
      TRUE ~ sprintf("%s (n=%d, P=%s)", spec_disp, n, fmt_p(p_LR))
    ),
    spec_lab = factor(spec_lab, levels = rev(spec_lab))
  )

p3d <- ggplot(sens_keep, aes(x = spec_lab, y = beta, color = category)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbar(aes(ymin = prof_ci_lo, ymax = prof_ci_hi), width = 0.28, linewidth = 0.48) +
  geom_point(size = 1.45) +
  coord_flip(ylim = c(-0.4, 5.0)) +
  scale_color_brewer(palette = "Dark2", guide = "none") +
  labs(title = "Sensitivity & outcome heterogeneity",
       subtitle = "Tobit slope (95% profile CI; cluster-boot CI for weighted)",
       x = NULL, y = "Logit slope β (95% CI)") +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 6.0),
        axis.text.y = element_text(size = 5.65),
        plot.margin = margin(4, 6, 4, 4))

# Panel E: Joint nested models — added LR of ratio vs added LR of absolute
joint_long <- joint_df %>%
  select(ratio_source, added_LR_ratio_over_abs, p_added_ratio, added_LR_abs_over_ratio, p_added_abs) %>%
  mutate(
    src_short = recode(ratio_source,
      "Choi OSCC / oral basal (log2+1)" = "Choi OSCC /\noral basal",
      "Puram OSCC / oral basal (log2+1)" = "Puram OSCC /\noral basal",
      "TCGA-HNSC oral paired LFC (n=30 pairs)" = "TCGA oral\npaired LFC"
    )
  ) %>%
  pivot_longer(
    cols = c(added_LR_ratio_over_abs, added_LR_abs_over_ratio),
    names_to = "term", values_to = "added_LR"
  ) %>%
  mutate(
    p_val = ifelse(term == "added_LR_ratio_over_abs", p_added_ratio, p_added_abs),
    term_lab = ifelse(term == "added_LR_abs_over_ratio",
                      "Add Absolute over Ratio", "Add Ratio over Absolute"),
    src_short = factor(src_short, levels = c("Choi OSCC /\noral basal", "Puram OSCC /\noral basal", "TCGA oral\npaired LFC")),
    term_lab = factor(term_lab, levels = c("Add Absolute over Ratio", "Add Ratio over Absolute"))
  )

p3e <- ggplot(joint_long, aes(x = src_short, y = added_LR, fill = term_lab)) +
  geom_col(position = position_dodge(width = 0.70), width = 0.60) +
  geom_hline(yintercept = qchisq(0.95, 1), linetype = "dashed", color = "grey45", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f\n(P=%s)", added_LR, fmt_p(p_val))),
            position = position_dodge(width = 0.70), vjust = -0.22, size = 2.15) +
  scale_fill_manual(values = c("Add Absolute over Ratio" = "#0072B2", "Add Ratio over Absolute" = "grey65"), name = NULL) +
  scale_y_continuous(limits = c(0, 17.5), expand = expansion(mult = c(0, 0.04))) +
  labs(title = "Nested joint models: abs + ratio + payload",
       subtitle = sprintf("Added Tobit LR (df=1): ratio adds ≤%.2f; absolute adds %.1f–%.1f",
                          max(joint_df$added_LR_ratio_over_abs), min(joint_df$added_LR_abs_over_ratio), max(joint_df$added_LR_abs_over_ratio)),
       x = NULL, y = "Added LR statistic (df = 1)") +
  theme_pub(base_size = 7.3) +
  theme(legend.position = "top",
        legend.key.size = unit(0.24, "cm"),
        legend.text = element_text(size = 6.3, margin = margin(r = 8)),
        legend.margin = margin(b = -4),
        plot.margin = margin(4, 65, 4, 65))

fig3_r1 <- plot_grid(p3a, p3b, ncol = 2, align = "h", axis = "tb", labels = c("A", "B"), label_size = 10)
fig3_r2 <- plot_grid(p3c, p3d, ncol = 2, rel_widths = c(0.82, 1.18), labels = c("C", "D"), label_size = 10)
fig3 <- plot_grid(fig3_r1, NULL, fig3_r2, NULL, p3e, ncol = 1, rel_heights = c(0.325, 0.02, 0.405, 0.02, 0.23), labels = c("", "", "", "", "E"), label_size = 10)
save_fig(fig3, "Fig3_resolution_and_specificity", 7.09, 9.2)


# ==============================================================================
# FIGURE 4: Out-of-sample LOTO calibration & FAERS pharmacovigilance
# ==============================================================================
p5 <- read_tsv(paste0(R, "part5_loto_predictions.tsv"), show_col_types = FALSE)
calib_sum <- read_tsv(paste0(R, "part5_calib_summary.tsv"), show_col_types = FALSE)
loto_row <- calib_sum %>% slice(1)
faers <- read_tsv(paste0(R, "fig_faers_summary.tsv"), show_col_types = FALSE)

num5 <- p5 %>%
  filter(!cens) %>%
  mutate(
    lbl = short_adc(drug, target_gene),
    lx = case_when(
      drug == "SHR-A1921"                 ~ 49.0,
      drug == "datopotamab deruxtecan"    ~ 72.5,
      drug == "sacituzumab tirumotecan"   ~ 72.5,
      drug == "sacituzumab govitecan"     ~ 76.5,
      drug == "izalontamab brengitecan"   ~ 20.5,
      drug == "sigvotatug vedotin"        ~ 17.5,
      drug == "gemtuzumab ozogamicin"     ~ 23.0,
      drug == "glembatumumab vedotin"     ~ 58.0,
      drug == "trastuzumab deruxtecan"    ~ 47.0,
      drug == "becotatug vedotin"         ~ 49.0,
      drug == "vobramitamab duocarmazine" ~ 47.0,
      drug == "inotuzumab ozogamicin"     ~ 40.5,
      drug == "patritumab deruxtecan"     ~ 49.0,
      drug == "trastuzumab rezetecan"     ~ 45.0,
      drug == "trastuzumab emtansine"     ~  8.0,
      TRUE ~ pred_loto
    ),
    ly = case_when(
      drug == "SHR-A1921"                 ~ 76.5,
      drug == "datopotamab deruxtecan"    ~ 56.0,
      drug == "sacituzumab tirumotecan"   ~ 41.5,
      drug == "sacituzumab govitecan"     ~ 26.5,
      drug == "izalontamab brengitecan"   ~ 54.5,
      drug == "sigvotatug vedotin"        ~ 41.5,
      drug == "gemtuzumab ozogamicin"     ~ 32.5,
      drug == "glembatumumab vedotin"     ~ 33.5,
      drug == "trastuzumab deruxtecan"    ~ 26.5,
      drug == "becotatug vedotin"         ~ 19.5,
      drug == "vobramitamab duocarmazine" ~ 12.5,
      drug == "inotuzumab ozogamicin"     ~  6.0,
      drug == "patritumab deruxtecan"     ~  0.0,
      drug == "trastuzumab rezetecan"     ~ -5.5,
      drug == "trastuzumab emtansine"     ~ -5.5,
      TRUE ~ y_pct
    )
  ) %>%
  arrange(case_when(
    drug == "becotatug vedotin"         ~ 1L,
    drug == "vobramitamab duocarmazine" ~ 2L,
    drug == "trastuzumab deruxtecan"    ~ 3L,
    TRUE ~ 0L
  ))
cen5 <- p5 %>% filter(cens)

# Panel A: Primary LOTO numeric ADCs
p4a <- ggplot(num5, aes(x = pred_loto, y = y_pct, color = pclass)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_point(size = 1.95) +
  geom_text_repel(aes(label = lbl),
                  nudge_x = num5$lx - num5$pred_loto,
                  nudge_y = num5$ly - num5$y_pct,
                  size = 1.92, max.overlaps = Inf, box.padding = 0.20, point.padding = 0.22,
                  force = 0, force_pull = 0, max.iter = 0,
                  bg.color = "white", bg.r = 0.14,
                  min.segment.length = 0, segment.size = 0.30, segment.alpha = 0.65,
                  show.legend = FALSE) +
  scale_color_manual(values = pal_pclass, labels = pal_pclass_lab, name = NULL) +
  scale_x_continuous(breaks = c(0, 25, 50, 75)) +
  scale_y_continuous(breaks = c(0, 25, 50, 75)) +
  coord_equal(xlim = c(-6, 94), ylim = c(-7, 90)) +
  labs(title = sprintf("Leave-one-target-out: numeric ADCs (n=%d)", loto_row$n_numeric),
       subtitle = sprintf("Spearman ρ = %.3f (P=%.3f); all-44 C-index = %.3f vs %.3f (P=%.4f)",
                          loto_row$spearman_x, loto_row$p_spearman_x,
                          loto_row$c_index_all44_x, loto_row$c_index_all44_payload, loto_row$p_perm_c_index),
       x = "Out-of-target predicted stomatitis (%)",
       y = "Observed clinical stomatitis (%)") +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 5.95),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 6.0),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-3, 0, -4, 0))

# Panel B: Primary LOTO censored ADCs
n_misses_cen <- sum(cen5$pred_loto >= cen5$y_pct)
miss_lbl <- sprintf("Above threshold (%d misses)", n_misses_cen)
cen5_p <- cen5 %>%
  mutate(
    rescued = (pred_loto_payload_only >= y_pct) & (pred_loto < y_pct),
    missed = pred_loto >= y_pct,
    grp_status = case_when(
      missed ~ miss_lbl,
      rescued ~ "Rescued by target expression",
      TRUE ~ "Correctly below threshold"
    )
  )

p4b <- ggplot(cen5_p, aes(x = reorder(paste0(drug, " (<", as.integer(y_pct), "%)"), pred_loto))) +
  geom_segment(aes(xend = reorder(paste0(drug, " (<", as.integer(y_pct), "%)"), pred_loto),
                   y = pred_loto, yend = y_pct), color = "grey75", linewidth = 0.45) +
  geom_point(aes(y = y_pct), shape = 4, size = 1.65, color = "grey30", stroke = 0.75) +
  geom_point(aes(y = pred_loto_payload_only), shape = 1, size = 1.55, color = "grey55") +
  geom_point(aes(y = pred_loto, color = grp_status), size = 1.8) +
  coord_flip() +
  scale_color_manual(values = setNames(c("#0072B2", "#009E73", "#D55E00"),
                                       c("Correctly below threshold", "Rescued by target expression", miss_lbl)),
                     name = "LOTO placement") +
  scale_y_continuous(limits = c(0, 88)) +
  labs(title = sprintf("LOTO censored ADCs (n=%d) & exact ≥15%% AUC", loto_row$n_censored),
       subtitle = sprintf("Exact ≥15%% AUC = %.3f vs %.3f (P=%.4f, n=38); %d/%d < cutoff",
                          loto_row$auc_ge15_x, loto_row$auc_ge15_payload, loto_row$p_mw_auc_ge15,
                          loto_row$cens_below_x, loto_row$n_censored),
       x = NULL, y = "Predicted stomatitis (%) vs trial threshold (×)") +
  theme_pub(base_size = 7.0) +
  theme(plot.subtitle = element_text(size = 5.95),
        axis.text.y = element_text(size = 5.65),
        legend.position = c(0.68, 0.25),
        legend.title = element_text(size = 6.0),
        legend.text = element_text(size = 5.7),
        legend.key.size = unit(0.22, "cm"))

# Panel C: FAERS negative & positive controls
faers_p <- faers %>%
  mutate(
    outcome_lab = recode(outcome,
      "oral" = "Oral mucositis /\nstomatitis",
      "skin" = "Cutaneous / rash\n(expl. ctrl)",
      "diarrhea" = "Diarrhea\n(expl. ctrl)",
      "neutropenia" = "Neutropenia\n(presp. neg. ctrl)",
      "neuropathy" = "Peripheral neuropathy\n(presp. neg. ctrl)",
      "ild" = "ILD / pneumonitis\n(expl. ctrl)",
      "ocular" = "Ocular toxicity\n(expl. ctrl)"
    ),
    is_oral = outcome == "oral",
    outcome_lab = factor(outcome_lab, levels = rev(c(
      "Oral mucositis /\nstomatitis", "Cutaneous / rash\n(expl. ctrl)", "Ocular toxicity\n(expl. ctrl)",
      "Neutropenia\n(presp. neg. ctrl)", "Diarrhea\n(expl. ctrl)",
      "Peripheral neuropathy\n(presp. neg. ctrl)", "ILD / pneumonitis\n(expl. ctrl)"
    )))
  )

p4c <- ggplot(faers_p, aes(x = outcome_lab, y = beta, color = is_oral)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_label(aes(label = ifelse(is_oral, sprintf("LR P=%.3f (perm P=%.4f)", p_LR, p_perm), sprintf("LR P=%.2f (perm P=%.2f)", p_LR, p_perm))),
             vjust = -0.68, size = 1.98, color = "grey20", fill = "white", linewidth = 0, label.padding = unit(0.06, "lines")) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi), width = 0.22, linewidth = 0.55) +
  geom_point(size = 2.1) +
  scale_x_discrete(expand = expansion(add = c(0.55, 0.90))) +
  coord_flip() +
  scale_color_manual(values = c(`TRUE` = "#D55E00", `FALSE` = "grey45"), guide = "none") +
  labs(title = "FAERS: oral mucositis replicates specifically",
       subtitle = "NB2 log-rate slope per log₂ oral basal unit (n=14 ADCs)",
       x = NULL, y = "Log-rate slope per log₂ oral basal unit (95% CI)") +
  theme_pub(base_size = 7.1) +
  theme(axis.text.y = element_text(size = 6.0, lineheight = 0.92))

# Panel D: FAERS pairwise oral ROR — 4 higher-vs-lower basal contrasts vs 1 same-target TACSTD2 contrast
ror_df <- read_tsv(paste0(R, "fig_faers_pairs.tsv"), show_col_types = FALSE) %>%
  mutate(contrast = factor(contrast, levels = rev(contrast)))

p4d <- ggplot(ror_df, aes(x = contrast, y = ror, color = pair_type)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = 1.5, linetype = "dotted", color = "grey55") +
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.22, linewidth = 0.55) +
  geom_point(size = 2.1) +
  geom_text(aes(y = ifelse(ror < 10, 1.32, ror),
                hjust = ifelse(ror < 10, 0, 0.5),
                label = sprintf("%.2f× (%.2f–%.2f)", ror, lo, hi)),
            vjust = -1.1, size = 2.0, show.legend = FALSE) +
  coord_flip() +
  scale_color_manual(values = c("Higher vs lower oral basal" = "#0072B2", "Same target (TACSTD2)" = "#D55E00"), name = NULL) +
  scale_y_log10(limits = c(0.45, 78), breaks = c(1, 3, 10, 30)) +
  labs(title = "FAERS pairwise oral reporting odds ratios",
       subtitle = "Between-target basal contrasts (blue) vs same-target TROP2 (orange)",
       x = NULL, y = "Oral reporting odds ratio (log scale, 95% CI)") +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 5.9),
        axis.text.y = element_text(size = 6.1, lineheight = 0.92),
        legend.position = "top",
        legend.justification = "left",
        legend.key.size = unit(0.22, "cm"),
        legend.text = element_text(size = 5.95),
        legend.margin = margin(b = -4))

fig4_top <- plot_grid(p4a, p4b, ncol = 2, rel_widths = c(0.98, 1.02), labels = c("A", "B"), label_size = 10)
fig4_bot <- plot_grid(p4c, p4d, ncol = 2, rel_widths = c(1.02, 0.98), labels = c("C", "D"), label_size = 10)
fig4 <- plot_grid(fig4_top, NULL, fig4_bot, ncol = 1, rel_heights = c(0.54, 0.03, 0.43))
save_fig(fig4, "Fig4_calibration_and_faers", 7.09, 7.2)


# ==============================================================================
# FIGURE 5: OSCC surfaceome coupling, re-ranking, Madan head-to-head, & winners
# ==============================================================================
e4 <- read_tsv(paste0(R, "part4_final_ranking.tsv"), show_col_types = FALSE) %>%
  rename(gene = 1)

# Panel A: Structural coupling across 420 OSCC surfaceome targets
n_win_4 <- sum(e4$n_pass == 4, na.rm = TRUE)
n_win_vlow <- sum(e4$n_pass == 4 & e4$gtex_heart <= 15 & e4$gtex_liver <= 15 & e4$gtex_kidney <= 15 & e4$gtex_nerve <= 15, na.rm = TRUE)
n_win_other <- n_win_4 - n_win_vlow
lbl_vlow <- sprintf("4/4 winner + 4-organ ≤15 TPM (n=%d)", n_win_vlow)
lbl_woth <- sprintf("4/4 winner (n=%d other)", n_win_other)
rho_choi_s <- cor(e4$tum_choi, e4$oral_basal, method = "spearman")
rho_puram_s <- cor(e4$tum_puram, e4$oral_basal, method = "spearman")
pct_d1_s <- 100 * mean((e4$tum_choi - e4$oral_basal) >= 1.0)

e4_p <- e4 %>%
  mutate(
    win_choi = tum_choi - oral_basal,
    grp = case_when(
      n_pass == 4 & gtex_heart <= 15 & gtex_liver <= 15 & gtex_kidney <= 15 & gtex_nerve <= 15 ~ lbl_vlow,
      n_pass == 4 ~ lbl_woth,
      gene %in% c("TACSTD2", "EGFR", "CD44", "GPNMB", "ERBB2") ~ "Pursued target (failed ≥1)",
      TRUE ~ "Other OSCC surfaceome"
    )
  )

lbl_5a <- e4_p %>%
  filter(gene %in% c("TACSTD2", "CD44", "GPNMB", "EGFR", "ERBB2", "TNFRSF12A", "F3", "CDH3", "PTK7", "FAT2", "ITGA3", "SLC7A8"))

p5a <- ggplot(e4_p, aes(x = oral_basal, y = tum_choi, color = grp)) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "grey65", linewidth = 0.4) +
  geom_abline(slope = 1, intercept = 1, linetype = "dashed", color = "#009E73", linewidth = 0.45) +
  geom_abline(slope = 1, intercept = 2, linetype = "dotted", color = "#0072B2", linewidth = 0.45) +
  geom_point(aes(alpha = grp == "Other OSCC surfaceome"), size = 1.4) +
  geom_text_repel(data = lbl_5a, aes(label = gene), size = 2.0, fontface = "bold", max.overlaps = 30,
                  box.padding = 0.34, point.padding = 0.20, force = 8, seed = 42,
                  bg.color = "white", bg.r = 0.14, min.segment.length = 0, show.legend = FALSE) +
  scale_alpha_manual(values = c(`TRUE` = 0.45, `FALSE` = 0.95), guide = "none") +
  scale_color_manual(values = setNames(
    c("#009E73", "#56B4E9", "#D55E00", "grey70"),
    c(lbl_vlow, lbl_woth, "Pursued target (failed ≥1)", "Other OSCC surfaceome")
  ), name = NULL) +
  annotate("text", x = 3.35, y = 4.72, label = "Δlog₂(CP10K+1) ≥ 1", color = "#009E73", size = 1.95, angle = 37) +
  annotate("text", x = 3.15, y = 5.45, label = "Δlog₂(CP10K+1) ≥ 2", color = "#0072B2", size = 1.95, angle = 37) +
  labs(title = sprintf("OSCC surfaceome structural coupling (n=%d)", nrow(e4)),
       subtitle = sprintf("Spearman ρ = %.3f (Choi), %.3f (Puram); %.1f%% reach Δlog₂ ≥ 1", rho_choi_s, rho_puram_s, pct_d1_s),
       x = "Healthy oral basal KC, log₂(CP10K + 1)",
       y = "OSCC malignant cells (Choi), log₂(CP10K + 1)") +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 5.95),
        legend.position = c(0.37, 0.87),
        legend.text = element_text(size = 5.75),
        legend.key.size = unit(0.20, "cm"),
        legend.background = element_rect(fill = NA, color = NA))

# Panel B: Clinically pursued + Madan 2023 nominees
bench_genes <- c("TACSTD2", "COL17A1", "DSC3", "DSG3", "CD44", "GPNMB", "EGFR", "ERBB2", "ERBB3", "AXL", "CD276", "MET", "ITGB6", "CDH3", "PTK7", "F3", "FAT2")
bench <- e4 %>%
  filter(gene %in% bench_genes) %>%
  mutate(
    source_tag = case_when(
      gene %in% c("COL17A1", "DSC3", "DSG3") ~ paste0(gene, " [Madan dual]"),
      gene == "CDH3" ~ "CDH3 [Madan dual + 4/4]",
      gene == "TACSTD2" ~ "TACSTD2 [Clin + Madan sel.]",
      gene == "FAT2" ~ "FAT2 [4/4; prot +0.02]",
      TRUE ~ gene
    ),
    source_tag = factor(source_tag, levels = rev(c(
      "TACSTD2 [Clin + Madan sel.]", "COL17A1 [Madan dual]", "DSC3 [Madan dual]", "DSG3 [Madan dual]",
      "CD44", "GPNMB", "EGFR", "ERBB2", "ERBB3", "AXL", "CD276", "MET", "ITGB6",
      "PTK7", "F3", "CDH3 [Madan dual + 4/4]", "FAT2 [4/4; prot +0.02]"
    )))
  ) %>%
  select(source_tag, pred_stom_pct, pass_cover, pass_tox, pass_tcga, pass_prot) %>%
  pivot_longer(cols = starts_with("pass_"), names_to = "crit", values_to = "passed") %>%
  mutate(crit = recode(crit,
    "pass_cover" = "Tumor\ncoverage\n(≥50%)",
    "pass_tox" = "Pred.\nstomatitis\n(≤12%)",
    "pass_tcga" = "TCGA\npaired\nRNA up",
    "pass_prot" = "CPTAC\npaired\nprot. up"
  ),
  crit = factor(crit, levels = c("Tumor\ncoverage\n(≥50%)", "Pred.\nstomatitis\n(≤12%)", "TCGA\npaired\nRNA up", "CPTAC\npaired\nprot. up")))

p5b <- ggplot(bench, aes(x = crit, y = source_tag, fill = passed)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = ifelse(passed, "PASS", "FAIL")), size = 2.0, fontface = "bold", color = "white") +
  scale_fill_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#D55E00"), guide = "none") +
  labs(title = "Clinically pursued & Madan 2023 targets",
       subtitle = "Exploratory 4 criteria; PTK7, F3, CDH3 & FAT2 pass all four",
       x = NULL, y = NULL) +
  theme_pub(base_size = 7.1) +
  theme(plot.subtitle = element_text(size = 5.95),
        axis.text.y = element_text(size = 5.9),
        axis.text.x = element_text(size = 5.8, lineheight = 0.90))

# Panel C: 12 4-organ <=15 TPM winners with 95% PIs, 24-organ extended max TPM & 3x3 stability scores
top_show <- c("CDH3", "PTPRZ1", "PTK7", "FAT2", "FZD6", "DSG2", "PTPRK", "ITGA2", "CDCP1", "IL1RAP", "IGSF3", "ATP2B1", "TNFRSF12A", "F3", "TACSTD2")
prof <- e4 %>%
  filter(gene %in% top_show) %>%
  mutate(
    gene_lbl = case_when(
      gene == "FAT2" ~ sprintf("FAT2 [%d/9 | ext %d/9, prot +0.02]", grid_pass, grid_pass_ext),
      gene %in% top_show[1:12] ~ sprintf("%s [%d/9 | ext %d/9]", gene, grid_pass, grid_pass_ext),
      TRUE ~ gene
    ),
    max_vital = pmax(gtex_heart, gtex_liver, gtex_kidney, gtex_nerve)
  )
n_rob9 <- sum(prof$grid_pass == 9 & prof$gene %in% top_show[1:12])
n_rob6 <- sum(prof$grid_pass >= 6 & prof$gene %in% top_show[1:12])
n_ext15 <- sum(prof$ext_max_tpm <= 15 & prof$gene %in% top_show[1:12])
ord_lbl <- prof %>% arrange(match(gene, top_show)) %>% pull(gene_lbl)
prof <- prof %>% mutate(gene_lbl = factor(gene_lbl, levels = rev(ord_lbl)))

p5c <- ggplot(prof, aes(x = pred_stom_pct, y = gene_lbl)) +
  geom_vline(xintercept = 12, linetype = "dashed", color = "#D55E00", linewidth = 0.45) +
  geom_errorbar(aes(xmin = pred_pi_lo, xmax = pred_pi_hi), width = 0.25, color = "grey65", linewidth = 0.45) +
  geom_point(aes(size = tum_choi_det * 100, color = tcga_lfc)) +
  geom_text(aes(x = pred_pi_hi, label = sprintf("  %.1f%% | 4-org %.1f (ext24 %.1f %s) TPM", pred_stom_pct, max_vital, ext_max_tpm, ext_max_tissue)), hjust = 0, size = 2.0) +
  scale_color_gradient2(low = "#2166AC", mid = "grey80", high = "#B2182B", midpoint = 0, name = "TCGA\npaired LFC",
                        guide = guide_colorbar(barheight = unit(1.3, "cm"), barwidth = unit(0.24, "cm"))) +
  scale_size_continuous(range = c(1.6, 4.0), name = "Choi malignant\ndetection (%)") +
  scale_x_continuous(limits = c(0, 178), breaks = c(0, 25, 50, 75, 100)) +
  labs(title = "Re-ranked candidates vs TROP2 (point estimate, 95% PI & 4- vs 24-organ max TPM)",
       subtitle = sprintf("Overlapping 95%% PIs preclude fine ranking; 4-org: %d/12 (9/9), %d/12 (≥6/9); 24-org ≤15 TPM: %d/12", n_rob9, n_rob6, n_ext15),
       x = "Predicted topo-I stomatitis (%, 95% prediction interval)", y = NULL) +
  theme_pub(base_size = 7.2) +
  theme(axis.text.y = element_text(size = 6.1),
        legend.position = c(0.89, 0.45),
        legend.title = element_text(size = 6.1),
        legend.text = element_text(size = 5.8),
        legend.spacing.y = unit(0.10, "cm"))

fig5_top <- plot_grid(p5a, p5b, ncol = 2, rel_widths = c(1.03, 0.97), labels = c("A", "B"), label_size = 10)
fig5 <- plot_grid(fig5_top, NULL, p5c, ncol = 1, rel_heights = c(0.47, 0.03, 0.50), labels = c("", "", "C"), label_size = 10)
save_fig(fig5, "Fig5_oscc_surfaceome_reranking", 7.09, 8.1)
