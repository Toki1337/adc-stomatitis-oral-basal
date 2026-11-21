from pathlib import Path
import scanpy as sc, anndata as ad, pandas as pd, numpy as np, glob, os, scipy.sparse as sp
ROOT = Path(__file__).resolve().parents[1]
D = f'{ROOT}/dl/GSE164241/raw'
adatas=[]
for m in sorted(glob.glob(D+'/*_matrix.mtx.gz')):
    pre=m.replace('_matrix.mtx.gz','')
    s=os.path.basename(pre).split('_',1)[1]
    if not (s.startswith('BM') or s.startswith('GM')): continue
    feat=pre+'_features.tsv.gz' if os.path.exists(pre+'_features.tsv.gz') else pre+'_genes.tsv.gz'
    X=sc.read_mtx(m).X.T.tocsr()
    bc=pd.read_csv(pre+'_barcodes.tsv.gz',header=None)[0].values
    fe=pd.read_csv(feat,header=None,sep='\t')
    a=ad.AnnData(X=X,obs=pd.DataFrame(index=[s+'_'+b for b in bc]),var=pd.DataFrame(index=fe[1].values))
    a.var_names_make_unique()
    a.obs['sample']=s; a.obs['site']='buccal' if s.startswith('BM') else 'gingiva'
    adatas.append(a); print(s,a.shape,flush=True)
A=ad.concat(adatas,join='inner'); A.obs_names_make_unique()
print(A)
sc.pp.filter_cells(A,min_genes=300)
A.var['mt']=A.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(A,qc_vars=['mt'],inplace=True,percent_top=None,log1p=False)
A=A[(A.obs.pct_counts_mt<20)&(A.obs.n_genes_by_counts<7000)].copy()
A.layers['counts']=A.X.copy()
sc.pp.normalize_total(A,target_sum=1e4); sc.pp.log1p(A)
sc.pp.highly_variable_genes(A,n_top_genes=3000,batch_key='sample')
B=A[:,A.var.highly_variable].copy(); sc.pp.scale(B,max_value=10); sc.tl.pca(B,n_comps=40)
try:
    import harmonypy
    sc.external.pp.harmony_integrate(B,'sample'); rep='X_pca_harmony'
except Exception as e:
    print('no harmony',e); rep='X_pca'
sc.pp.neighbors(B,use_rep=rep); sc.tl.leiden(B,resolution=1.0,key_added='leiden')
A.obs['leiden']=B.obs['leiden']; A.obsm['X_umap']=None if False else B.obsm.get('X_umap')
A.write(f'{ROOT}/results/williams_healthy.h5ad')
markers=['KRT5','KRT14','KRT15','COL17A1','TP63','MKI67','TOP2A','KRT13','KRT4','KRT1','KRT10','SPRR3','CRNN','KRT76','PTPRC','CD3E','MS4A1','CD68','LYZ','COL1A1','DCN','PECAM1','VWF','ACTA2','RGS5','TPSAB1','PMEL','MLANA','PLP1','JCHAIN','CPA3']
mk=[m for m in markers if m in A.var_names]
df=pd.DataFrame({g:np.asarray(A[:,g].X.mean(axis=0)).ravel() for g in []})
M=sc.get.obs_df(A,keys=mk+['leiden'])
print(M.groupby('leiden').mean().round(2).to_string())
print(A.obs.leiden.value_counts().sort_index().to_string())
