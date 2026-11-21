from pathlib import Path
import scanpy as sc, pandas as pd, numpy as np, scipy.sparse as sp
ROOT = Path(__file__).resolve().parents[1]
E=sc.read_h5ad(f'{ROOT}/results/williams_epi.h5ad')
A=sc.read_h5ad(f'{ROOT}/results/williams_healthy.h5ad',backed='r')
comp={'basal':['0','1','13','5','6'],'suprabasal':['2','3','9'],'superficial':['4','8','11']}
E.obs['comp']='other'
for k,v in comp.items(): E.obs.loc[E.obs['sub'].isin(v),'comp']=k
print(E.obs.comp.value_counts())
X=E.layers['counts']; lib=np.asarray(X.sum(1)).ravel(); Xn=sp.diags(1e4/lib)@X
genes=E.var_names
res={}
for k in ['basal','suprabasal','superficial']:
    m=(E.obs.comp==k).values
    res[f'w_{k}_mean']=np.asarray(Xn[m].mean(0)).ravel(); res[f'w_{k}_det']=np.asarray((X[m]>0).mean(0)).ravel()
    # per-donor pseudobulk median
    pdm=[]
    for s in E.obs['sample'].unique():
        mm=m&(E.obs['sample']==s).values
        if mm.sum()>=20: pdm.append(np.asarray(Xn[mm].mean(0)).ravel())
    res[f'w_{k}_donormed']=np.median(np.array(pdm),0); res[f'w_{k}_ndonor']=len(pdm)
m=E.obs.comp.isin(['basal','suprabasal','superficial']).values
res['w_kc_mean']=np.asarray(Xn[m].mean(0)).ravel()
# non-epithelial whole-tissue pseudobulk (all healthy cells) for bulk-like comparison
Xa=A.layers['counts'][:] if 'counts' in A.layers else A.X[:]
Xa=sp.csr_matrix(Xa); la=np.asarray(Xa.sum(1)).ravel()
res['w_alltissue_mean']=np.asarray((sp.diags(1e4/la)@Xa).mean(0)).ravel()
R=pd.DataFrame(res,index=genes).groupby(level=0).mean()
R.to_csv(f'{ROOT}/results/williams_kc_stats.tsv.gz',sep='\t')
print(R.loc[['TACSTD2','ERBB2','ERBB3','EGFR','NECTIN4','F3','CD276','ITGB6','MET','FOLR1','CD19','CD33','TNFRSF8','SDC1','CD44','GPNMB','CEACAM5','PTK7','KRT15']].round(2).iloc[:,[0,1,3,5,7,9,10]].to_string())
