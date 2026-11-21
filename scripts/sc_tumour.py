from pathlib import Path
import pandas as pd, numpy as np, scipy.sparse as sp, gzip
ROOT = Path(__file__).resolve().parents[1]
R=f'{ROOT}/results/'
# ---- Choi GSE181919
meta=pd.read_csv(f'{ROOT}/dl/GSE181919/GSE181919_Barcode_metadata.txt.gz',sep='\t',index_col=0)
rows=[];genes=[]
it=pd.read_csv(f'{ROOT}/dl/GSE181919/GSE181919_UMI_counts.txt.gz',sep='\t',index_col=0,chunksize=2000)
for ch in it:
    genes+=list(ch.index); rows.append(sp.csr_matrix(ch.values.astype(np.float32)))
    cols=list(ch.columns)
X=sp.vstack(rows).T.tocsr()   # cells x genes
cells=[c.replace('.','-') for c in cols]
meta=meta.reindex(cells)
print('Choi',X.shape, meta['cell.type'].isna().sum())
lib=np.asarray(X.sum(1)).ravel(); Xn=sp.diags(1e4/lib)@X
def grp_stats(Xn,mask,name):
    sub=Xn[mask]
    mean=np.asarray(sub.mean(0)).ravel(); det=np.asarray((sub>0).mean(0)).ravel()
    return pd.DataFrame({name+'_mean':mean,name+'_det':det},index=genes)
groups={
 'choi_malig_OC_HPVneg':(meta['cell.type']=='Malignant.cells')&(meta.subsite=='OC')&(meta.hpv=='HPV-')&(meta['tissue.type']=='CA'),
 'choi_malig_all_CA':(meta['cell.type']=='Malignant.cells')&(meta['tissue.type']=='CA'),
 'choi_epi_NL':(meta['cell.type']=='Epithelial.cells')&(meta['tissue.type']=='NL'),
 'choi_epi_NL_OC':(meta['cell.type']=='Epithelial.cells')&(meta['tissue.type']=='NL')&(meta.subsite=='OC'),
 'choi_fibro':(meta['cell.type']=='Fibroblasts'),
}
out=[]
for k,m in groups.items():
    m=m.fillna(False).values; print(k,m.sum(), meta.loc[m,'sample.id'].nunique())
    out.append(grp_stats(Xn,m,k))
# per-patient pseudobulk for malignant OC HPV-neg (for homogeneity across patients)
C=pd.concat(out,axis=1); C=C.groupby(level=0).mean()
C.to_csv(R+'choi_group_stats.tsv.gz',sep='\t')
mm=groups['choi_malig_OC_HPVneg'].fillna(False).values
pp={}
for s in meta.loc[mm,'sample.id'].unique():
    m2=mm&(meta['sample.id']==s).values
    if m2.sum()>=20: pp[s]=np.asarray(Xn[m2].mean(0)).ravel()
pd.DataFrame(pp,index=genes).groupby(level=0).mean().to_csv(R+'choi_malig_OC_perpatient.tsv.gz',sep='\t')
# ---- Puram GSE103322
P=pd.read_csv(f'{ROOT}/dl/GSE103322/GSE103322_HNSCC_all_data.txt.gz',sep='\t',index_col=0,low_memory=False)
ann=P.iloc[:5]; E=P.iloc[5:].astype(np.float32)
E.index=[g.strip("'") for g in E.index]
ln=ann.loc['Lymph node'].astype(int); ca=ann.loc['classified  as cancer cell'].astype(int); nct=ann.loc['non-cancer cell type']
tpm=10*(2**E-1)
lin=tpm/100.0   # TPM/100 ~ CP10K scale
import re
pat=pd.Series([re.sub(r'^HN(\d+)$', r'HNSCC\1', c.split('_')[0]) for c in E.columns],index=E.columns)
res={}
for name,mask in {'puram_malig_primary':(ca==1)&(ln==0),'puram_malig_LN':(ca==1)&(ln==1),'puram_fibro':(nct=='Fibroblast')}.items():
    sub=lin.loc[:,mask.values]; res[name+'_mean']=sub.mean(1); res[name+'_det']=(sub>0).mean(1)
    pats_valid=pat[mask.values]
    print(name,mask.sum(),pats_valid[pats_valid!='HNSCC'].nunique())
pd.DataFrame(res).groupby(level=0).mean().to_csv(R+'puram_group_stats.tsv.gz',sep='\t')
mm=((ca==1)&(ln==0)).values
pp={p:lin.loc[:,mm&(pat==p).values].mean(1) for p in pat[mm].unique() if (mm&(pat==p).values).sum()>=20}
pd.DataFrame(pp).groupby(level=0).mean().to_csv(R+'puram_malig_perpatient.tsv.gz',sep='\t')
print('done')
