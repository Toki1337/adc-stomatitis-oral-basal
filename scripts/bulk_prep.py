from pathlib import Path
import pandas as pd, numpy as np
ROOT = Path(__file__).resolve().parents[1]
R=f'{ROOT}/results/'; D=f'{ROOT}/dl/'
# TCGA
E=pd.read_csv(D+'tcga/HiSeqV2.gz',sep='\t',index_col=0)
C=pd.read_csv(D+'tcga/HNSC_clinicalMatrix',sep='\t',index_col=0)
oral={'Alveolar Ridge','Buccal Mucosa','Floor of mouth','Hard Palate','Lip','Oral Cavity','Oral Tongue'}
C=C.reindex(E.columns)
isoral=C.anatomic_neoplasm_subdivision.isin(oral)
st=C.sample_type
tum=E.loc[:,(isoral&(st=='Primary Tumor')).values]; nor=E.loc[:,(isoral&(st=='Solid Tissue Normal')).values]
print('TCGA oral tumour',tum.shape[1],'normal',nor.shape[1])
# paired
pid=lambda s:s[:12]
pairs=[(t,n) for n in nor.columns for t in tum.columns if pid(t)==pid(n)]
print('pairs',len(pairs))
T=pd.DataFrame({'tcga_tum_mean':tum.mean(1),'tcga_tum_q25':tum.quantile(0.25,axis=1),'tcga_nor_mean':nor.mean(1)})
d=np.array([tum[t].values-nor[n].values for t,n in pairs])
T['tcga_paired_lfc']=d.mean(0); T['tcga_paired_frac_up']=(d>0).mean(0)
T.to_csv(R+'tcga_oral_stats.tsv.gz',sep='\t')
# CPTAC (oral cavity)
cli=pd.read_csv(D+'cptac/HS_CPTAC_HNSCC_CLI.tsi',sep='\t',index_col=0).iloc[1:]
oc=set(cli.index[cli.tumor_site_curated=='Oral cavity'])
out={}
for kind in ['Proteomics_TMT_Gene_level','RNAseq_RSEM_UQ_log2']:
    t=pd.read_csv(D+f'cptac/HS_CPTAC_HNSCC_{kind}_Tumor.cct',sep='\t',index_col=0)
    n=pd.read_csv(D+f'cptac/HS_CPTAC_HNSCC_{kind}_Normal.cct',sep='\t',index_col=0)
    t=t[[c for c in t.columns if c in oc]]; n=n[[c for c in n.columns if c in oc]]
    pp=[c for c in t.columns if c in n.columns]
    print(kind,'oral T',t.shape[1],'N',n.shape[1],'paired',len(pp))
    k='prot' if kind.startswith('Prot') else 'rna'
    out[f'cptac_{k}_tum_mean']=t.mean(1); out[f'cptac_{k}_nor_mean']=n.mean(1)
    dd=t[pp].values-n[pp].values
    out[f'cptac_{k}_paired_lfc']=pd.Series(np.nanmean(dd,1),index=t.index)
    out[f'cptac_{k}_paired_n']=pd.Series(np.sum(~np.isnan(dd),1),index=t.index)
    out[f'cptac_{k}_paired_frac_up']=pd.Series(np.nanmean(np.where(np.isnan(dd),np.nan,dd>0),1),index=t.index)
pd.DataFrame(out).to_csv(R+'cptac_oral_stats.tsv.gz',sep='\t')
# GTEx medians
G=pd.read_csv(D+'gtex/gene_median_tpm.gct.gz',sep='\t',skiprows=2,index_col=1).drop(columns=['Name'])
G=G.groupby(level=0).max(); G.to_csv(R+'gtex_median_tpm.tsv.gz',sep='\t')
# HPA IHC
H=pd.read_csv(D+'hpa/normal_ihc_data.tsv',sep='\t')
lv={'Not detected':0,'Low':1,'Medium':2,'High':3}
H['score']=H.Level.map(lv)
H.to_csv(R+'hpa_ihc_scored.tsv.gz',sep='\t',index=False)
print(H[H.Tissue=='Oral mucosa']['Cell type'].value_counts())
S=pd.read_csv(D+'hpa/rna_single_cell_type.tsv',sep='\t').pivot_table(index='Gene name',columns='Cell type',values='nCPM',aggfunc='max')
S.to_csv(R+'hpa_sc_ncpm.tsv.gz',sep='\t')
print('done')
