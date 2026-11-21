from pathlib import Path
import scanpy as sc, pandas as pd, numpy as np, scipy.sparse as sp
ROOT = Path(__file__).resolve().parents[1]
A=sc.read_h5ad(f'{ROOT}/results/williams_healthy.h5ad')
E=A[A.obs.leiden.isin(['5','28','29'])].copy()
print('epithelial cells',E.n_obs, E.obs['sample'].nunique())
B=E.copy()
sc.pp.highly_variable_genes(B,n_top_genes=2000,batch_key='sample')
B=B[:,B.var.highly_variable].copy(); sc.pp.scale(B,max_value=10); sc.tl.pca(B,n_comps=30)
import harmonypy
ho=harmonypy.run_harmony(B.obsm['X_pca'],B.obs,'sample',random_state=0)
Z=ho.Z_corr; Z=Z.T if Z.shape[0]!=B.n_obs else Z
B.obsm['X_pca_harmony']=Z
sc.pp.neighbors(B,use_rep='X_pca_harmony'); sc.tl.leiden(B,resolution=0.6,key_added='sub')
E.obs['sub']=B.obs['sub']
mk=['KRT5','KRT14','KRT15','COL17A1','TP63','ITGA6','ITGB4','LAMB3','MKI67','TOP2A','KRT6A','KRT13','KRT4','SPRR3','CRNN','KRT10','KRT1','CNFN','PTPRC','COL1A1','PECAM1']
M=sc.get.obs_df(E,keys=mk+['sub'])
print(M.groupby('sub').mean().round(2).to_string()); print(E.obs['sub'].value_counts().sort_index().to_string())
E.write(f'{ROOT}/results/williams_epi.h5ad')
