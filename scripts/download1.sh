set -x
ROOT_DIR="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT_DIR/dl"
mkdir -p GSE181919 GSE103322 GSE164241 hpa gtex tcga
G=https://ftp.ncbi.nlm.nih.gov/geo/series
wget -q -c -O GSE181919/GSE181919_UMI_counts.txt.gz $G/GSE181nnn/GSE181919/suppl/GSE181919_UMI_counts.txt.gz
wget -q -c -O GSE181919/GSE181919_Barcode_metadata.txt.gz $G/GSE181nnn/GSE181919/suppl/GSE181919_Barcode_metadata.txt.gz
wget -q -c -O GSE103322/GSE103322_HNSCC_all_data.txt.gz $G/GSE103nnn/GSE103322/suppl/GSE103322_HNSCC_all_data.txt.gz
HPA=https://v26.proteinatlas.org/download/tsv
wget -q -c -O hpa/normal_ihc_data.tsv.zip $HPA/normal_ihc_data.tsv.zip
wget -q -c -O hpa/rna_single_cell_type.tsv.zip $HPA/rna_single_cell_type.tsv.zip
wget -q -c -O hpa/rna_tissue_consensus.tsv.zip $HPA/rna_tissue_consensus.tsv.zip
wget -q -c -O hpa/proteinatlas.tsv.zip $HPA/proteinatlas.tsv.zip
wget -q -c -O gtex/gene_median_tpm.gct.gz https://storage.googleapis.com/adult-gtex/bulk-gex/v8/rna-seq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz
wget -q -c -O tcga/HiSeqV2.gz https://tcga.xenahubs.net/download/TCGA.HNSC.sampleMap/HiSeqV2.gz
wget -q -c -O tcga/HNSC_clinicalMatrix https://tcga.xenahubs.net/download/TCGA.HNSC.sampleMap/HNSC_clinicalMatrix
wget -q -c -O GSE164241/GSE164241_RAW.tar $G/GSE164nnn/GSE164241/suppl/GSE164241_RAW.tar
wget -q -c -O GSE164241/GSE164241_submission_update.xlsx $G/GSE164nnn/GSE164241/suppl/GSE164241_submission_update.xlsx
echo DONE
