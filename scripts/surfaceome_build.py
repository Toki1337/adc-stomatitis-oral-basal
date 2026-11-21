# Surfaceome definition used by part4_final.py: reviewed UniProt human entries with a TM helix or GPI anchor,
# plus a signal peptide or cell-membrane annotation, excluding intracellular-organelle-only and primarily
# intracellular-compartment proteins (first subcellular location clause in ER/Golgi/Mitochondrion/Nucleus/Peroxisome/Cytoplasm).
# Writes results/surfaceome_genes.txt (3,152 genes; all 38 clinical ADC targets captured).
import os
import sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_UNIPROT = f'{ROOT}/dl/uniprot/human_all.tsv.gz'
OUT_SURF = f'{ROOT}/results/surfaceome_genes.txt'

if not os.path.exists(RAW_UNIPROT):
    if os.path.exists(OUT_SURF):
        n_existing = len(pd.read_csv(OUT_SURF, header=None))
        print(f'Raw UniProt file absent; preserving pre-computed surfaceome_genes.txt ({n_existing} genes).')
        sys.exit(0)
    raise FileNotFoundError(f'Missing both {RAW_UNIPROT} and {OUT_SURF}')

U = pd.read_csv(RAW_UNIPROT, sep='\t')
U.columns = ['acc', 'gene', 'scl', 'tm', 'lipid', 'kw']
scl = U.scl.fillna('')
kw = U.kw.fillna('')
tmf = U.tm.fillna('')
lip = U.lipid.fillna('')

has_tm = tmf.str.contains('TRANSMEM')
gpi = lip.str.contains('GPI-anchor') | scl.str.contains('GPI-anchor')
sig = kw.str.contains('Signal')
cellmem = scl.str.contains('Cell membrane') | kw.str.contains('Cell membrane')
intra_only = scl.str.contains('Mitochondrion|Endoplasmic reticulum|Golgi|Nucleus|Peroxisome') & ~cellmem & ~scl.str.contains('Membrane;')

def primary_intracellular(s):
    if not isinstance(s, str) or not s:
        return False
    first_clause = s.replace('SUBCELLULAR LOCATION: ', '').split('.')[0].split(';')[0]
    return any(o in first_clause for o in ['Endoplasmic reticulum', 'Golgi', 'Mitochondrion', 'Nucleus', 'Peroxisome', 'Cytoplasm'])

prim_intra = scl.map(primary_intracellular)
surf_mask = (has_tm | gpi) & (sig | cellmem) & ~intra_only & ~prim_intra & U.gene.notna()
surf = sorted(set(U.loc[surf_mask, 'gene']))
print('surfaceome:', len(surf))

adc = ['ALCAM','AXL','CD19','CD22','CD276','CD33','CD37','CD44','CD79B','CDH6','CEACAM5','DLL3','EGFR','ERBB2','ERBB3','F3','FGFR2','FOLR1','GPNMB','IL2RA','IL3RA','ITGB6','MET','MSLN','MUC16','NCAM1','NECTIN4','PTK7','ROR1','ROR2','SDC1','SLC34A2','SLC39A6','STEAP1','TACSTD2','TNFRSF17','TNFRSF8','VTCN1']
miss = [g for g in adc if g not in surf]
print('ADC targets captured: %d/%d  missing: %s' % (len(adc) - len(miss), len(adc), miss))
pd.Series(surf).to_csv(OUT_SURF, index=False, header=False)
