# Locked analysis plan (written 2026-11-17, before the clinical stomatitis table was returned and before any expression-outcome linkage)

Unit of analysis: one ADC (monotherapy). Clinical extraction is done blind to expression (rules in the extraction protocol).

Outcomes
- Y1 (primary): any-grade stomatitis/oral mucositis % in the largest monotherapy safety population (label > phase 3 > phase 2 > phase 1). Not listed in AE table -> left-censored below the table threshold.
- Y2 (secondary, independent source): FAERS reporting proportion of oral-mucositis PTs (STOMATITIS, MUCOSAL INFLAMMATION, MOUTH ULCERATION, APHTHOUS ULCER, GLOSSITIS, ORAL MUCOSAL ERYTHEMA, ORAL MUCOSAL BLISTERING) among reports where the ADC is the only suspect drug and no other antineoplastic is listed.
- Negative-control outcomes (FAERS, same report set): peripheral neuropathy and neutropenia PTs (payload-driven, not expected to track oral-mucosa target expression).

Primary predictor (X): target-gene expression in basal keratinocytes of healthy oral mucosa (GSE164241, healthy buccal + gingival donors only; pseudobulk log2(mean CP10K + 1)). Bispecific -> max of the two genes. CD44v6 -> CD44.
Replication predictors: Choi GSE181919 normal oral mucosa keratinocytes; HPA IHC oral mucosa staining level; TCGA-HNSC adjacent-normal oral bulk RNA; CPTAC HNSCC NAT proteome.

Tests
1. Spearman(X, Y1) across ADCs; censored values set to threshold/2 (sensitivity: threshold, 0, Tobit on logit scale).
2. Model: logit(Y1) ~ X + payload class (topo-I / tubulin / DNA-damaging); report X coefficient, likelihood-ratio test, and leave-one-target-out.
3. Within-payload contrasts (DXd, MMAE, DM1/DM4, SN-38).
4. Tissue specificity: repeat test 1 with X from each GTEx tissue; expect squamous mucosa (esophagus mucosa, minor salivary gland, skin) to rank top.
5. Y2 replication of tests 1-2; negative-control outcomes.

Part 2 (OSCC window): for ADC targets and surfaceome genes, tumour (TCGA oral, CPTAC tumour proteome, malignant cells in GSE103322 / GSE181919) vs normal oral basal keratinocytes; predicted stomatitis liability from the part-1 model; window = tumour excess over predicted liability; DepMap OSCC lines for model selection.
