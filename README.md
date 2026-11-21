# Absolute Oral Basal Keratinocyte Target Expression and ADC Stomatitis

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24429695.svg)](https://doi.org/10.5281/zenodo.24429695)

Reproducible analysis code, dual-curated clinical safety tables, processed summary matrices, figures, and supplementary data archives for:

**Absolute Target Expression in Normal Oral Basal Keratinocytes Predicts Antibody-Drug Conjugate Stomatitis and Re-Ranks Targets for Oral Cavity Carcinoma**

**Authors:** Wendi Dou†*, Xiang Yu†, Shangmin Zhang† (†Equal contribution; *Corresponding author: dou_v@student.sechenov.ru)  
**Affiliation:** I.M. Sechenov First Moscow State Medical University (Sechenov University), 8/2 Trubetskaya St., Moscow, 119991, Russia

## Repository Structure

- `ANALYSIS_PLAN.md` — Prespecified analytical protocol finalized prior to linking clinical stomatitis outcomes to gene or protein expression.
- `clinical/adc_stomatitis_extractor_A.tsv` — Original pre-linkage clinical safety extraction by Wendi Dou across all 66 ADCs (`Table S1a`).
- `clinical/adc_stomatitis_extractor_B.tsv` — Independent second extraction by Xiang Yu, blinded to gene/protein expression data and to Wendi Dou's entries, across all 66 ADCs (`Table S1b`).
- `clinical/adc_stomatitis_trials.tsv` — Final adjudicated monotherapy stomatitis table across 66 ADCs (44 evaluable across 28 primary targets: 15 numeric, 29 left-censored).
- `clinical/adc_stomatitis_reconciled.tsv` — Side-by-side comparison of Wendi Dou, Xiang Yu, and Adjudicated Final curations with adjudication notes across all 66 ADCs (`Table S1d`).
- `clinical/Dual_Blinded_Extractions_AB.xlsx` — Multi-sheet workbook containing Wendi Dou's extraction, Xiang Yu's blinded extraction, three-dataset parallel modeling comparison (`Table S1c`), 66-ADC side-by-side reconciliation (`Table S1d`), and inter-extractor agreement metrics (`Table S1e`).
- `dl/hpa_ihc/manifest_509.tsv` — Image manifest and Human Protein Atlas (HPA) URLs for the 509 primary normal squamous IHC sections (`189` oral mucosa, `166` esophagus, `154` skin across 28 genes).
- `scripts/` — End-to-end Python and R analytical, figure-generation, table-building, and archive-verification pipeline (`Makefile` orchestrates execution).
- `results/` — Processed summary tables, Tobit model outputs, Ruifrok & Johnston `rgb2hed` DAB outputs (`hpa_ihc_dab_models.tsv`, `hpa_ihc_dab_genes.tsv`, `hpa_ihc_dab_param_grid.tsv`), leave-one-target-out predictions, surfaceome prioritization table, `dab_qc/` overlay examples, and machine-readable `numbers.json`.
- `figures/` — Publication PDF and 300-dpi PNG figures (`Fig1`–`Fig5`).
- `manuscript/` — Multi-sheet `Supplementary_Tables.xlsx` (Tables S1–S8, README dictionary, and Supplementary References S1–S60) and `Supplementary_References.md`. *(The manuscript text is submitted separately for peer review; all analytical code, clinical data, and statistical outputs here correspond identically to the submitted study.)*

## Large Binary Data Accessions

Two large binary inputs exceed standard Git file-size limits and are omitted from the lightweight repository; each can be regenerated deterministically from public sources:

1. **Processed oral epithelial AnnData (`results/williams_epi.h5ad`, ~159 MB):**
   Download raw count matrices for healthy buccal mucosa and gingiva donors from NCBI GEO `GSE164241` (Williams et al., *Cell* 2021), then run `scripts/sc_williams.py`, `scripts/sc_williams_epi.py`, and `scripts/williams_stats.py`. Pre-computed compartment summary matrices (`results/williams_kc_stats.tsv.gz`) are included directly in `results/`.
2. **Raw HPA normal squamous IHC JPEGs (`dl/hpa_ihc/raw_jpg/`, 509 images):**
   Fetch directly from the HPA image URLs listed in `dl/hpa_ihc/manifest_509.tsv`. Pre-computed slide-level and gene-level `rgb2hed` optical density tables (`results/part1/hpa_ihc_dab_slides.tsv`, `hpa_ihc_dab_genes.tsv`, `hpa_ihc_dab_models.tsv`) are included in `results/`.

### Public Multi-Omic Accessions

NCBI GEO `GSE181919` (Choi et al.) and `GSE103322` (Puram et al.), TCGA-HNSC via UCSC Xena / NCI GDC, CPTAC-HNSCC oral cavity via PDC000221 / LinkedOmics, GTEx v8, HPA v26.0 Single Cell Type & IHC tables, UniProtKB/Swiss-Prot reviewed human proteome, DepMap 24Q4 CRISPR Chronos, openFDA Drug Label & FAERS APIs, and PubMed Central full texts.

## Reproduction

All scripts resolve paths dynamically from the repository root with zero hardcoded local directories. Install Python dependencies (`pip install -r requirements.txt`) and R packages listed in `R_sessionInfo.txt`, then follow either reproduction track:

1. **Track 1 — Reproduce all modeling, figures, supplementary tables, and archive verification from included processed summary matrices (~25 min, including 5,000 LOTO permutations in `part5` and 7 × 10,000 FAERS target permutations via `export_fig_data.py` in `numbers`):**
   ```bash
   make part3 part5 part4 numbers figures supp verify
   ```
   *(Notice `results/` includes all processed transcriptomic, proteomic, GTEx, HPA, FAERS, and `rgb2hed` DAB summary matrices, plus `results/surfaceome_genes.txt`, `results/fig_epi_markers.tsv`, and `results/supp_refs_metadata.json`, which `surfaceome_build.py`, `export_fig_data.py`, and `build_supp_tables.py` automatically preserve when raw UniProt / `.h5ad` / `dl/` intermediates are omitted.)*

2. **Track 2 — Full end-to-end rebuild from raw downloads (additionally runs `part1`):**
   Fetch `GSE164241` (`results/williams_epi.h5ad`), raw UniProt TSV (`dl/uniprot/human_all.tsv.gz`), and the 509 HPA JPEGs listed in `dl/hpa_ihc/manifest_509.tsv`, install optional upstream dependencies (`scanpy`, `anndata`, `harmonypy`), and run:
   ```bash
   make all
   ```
