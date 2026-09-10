# Run in the SAME Colab session. Still no motifs scored.
#
# Both sweeps now point at one thing. RUNX2 detection tracks sequencing depth:
#   gene-floor sweep : 17.84% at MIN_GENES=0  ->  21.16% at 1000
#   mito sweep       : 20.25% at MAX_MT=5     ->  14.17% at 100  (strict set is deeper)
# COL1A1 barely moves in either. That is the signature of a lowly expressed gene whose
# DETECTION is depth limited, which RUNX2 is.
#
# This matters far beyond a percentage. The preregistered primary correlates per nucleus
# RUNX motif deviation against RUNX2 expression. If RUNX2 detection tracks RNA depth AND
# motif deviation tracks ATAC depth, a correlation appears with no biology in it at all.
#
# This cell measures the confound before the statistic is fixed.

import numpy as np, pandas as pd
from scipy.stats import spearmanr

MT_CUT = 20      # the plateau start from the mito sweep
fibx = None

def vec(a, g):
    v = a[:, g].X
    return np.asarray(v.todense()).ravel() if hasattr(v,'todense') else np.asarray(v).ravel()

# rebuild the fibroblast set at the loosened cut, reusing the labeling already defined
import scanpy as sc
PANELS = {
 'fibroblast':['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','FBLN1'],
 'keratinocyte':['KRT14','KRT5','KRT1','KRT10','KRT6A'],
 'endothelial':['PECAM1','VWF','CDH5','EGFL7'],
 'pericyte_smc':['RGS5','MYH11','NOTCH3','KCNJ8'],
 'myeloid':['LYZ','CD68','ITGAM','CSF1R','AIF1'],
 't_nk':['CD3D','CD3E','IL7R','CD2','TRAC'],
 'b_plasma':['MS4A1','CD79A','JCHAIN','MZB1'],
 'mast':['TPSAB1','TPSB2','CPA3','MS4A2'],
 'melanocyte':['PMEL','MLANA','TYRP1','DCT'],
 'adnexal':['KRT7','AQP5','SCGB2A2','DCD'],
}
b = rna.copy()
b.var['mt'] = b.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(b, qc_vars=['mt'], inplace=True, percent_top=None, log1p=False)
sel = (b.obs.n_genes_by_counts >= 500) & (b.obs.pct_counts_mt <= MT_CUT)
a = b[sel].copy(); raw = a.copy()
sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
for n, gs in PANELS.items():
    sc.tl.score_genes(a, [g for g in gs if g in a.var_names], score_name=f's_{n}')
S = a.obs[[f's_{n}' for n in PANELS]].to_numpy(); o = np.argsort(-S, axis=1)
lab = np.array(list(PANELS))[o[:,0]]
lab[(S[np.arange(len(S)),o[:,0]] - S[np.arange(len(S)),o[:,1]]) < 0.05] = 'ambiguous'
fibx = raw[lab == 'fibroblast'].copy()
fib_atac = atac[fibx.obs_names].copy()

rna_depth  = np.asarray(fibx.X.sum(axis=1)).ravel()
atac_depth = np.asarray(fib_atac.X.sum(axis=1)).ravel()
print(f'MAX_MT={MT_CUT}: {fibx.n_obs:,} fibroblast nuclei')
print(f'RNA  counts per nucleus  median {np.median(rna_depth):,.0f}')
print(f'ATAC counts per nucleus  median {np.median(atac_depth):,.0f}')
print(f'RNA vs ATAC depth, Spearman rho = {spearmanr(rna_depth, atac_depth).statistic:+.3f}')
print()

print('=' * 74)
print('DOES DETECTION TRACK DEPTH?  rho of expression vs RNA depth, per gene')
print('=' * 74)
rows = []
for g in ['RUNX1','RUNX2','RUNX3','COL1A1','SMAD3','EGR1']:
    if g not in fibx.var_names:
        rows.append((g, np.nan, np.nan, np.nan)); continue
    v = vec(fibx, g)
    rows.append((g, round(100*float((v>0).mean()),2),
                 round(float(spearmanr(v, rna_depth).statistic), 3),
                 round(float(spearmanr((v>0).astype(int), rna_depth).statistic), 3)))
print(pd.DataFrame(rows, columns=['gene','pct nonzero','rho expr~depth','rho detected~depth']
                   ).to_string(index=False))
print()
print('COL1A1 is the reference. It is detected in nearly every nucleus, so its rho with')
print('depth should be small. A much larger rho for RUNX2 means RUNX2 zeros are dropout,')
print('not absence, and depth must be controlled in the primary.')
print()
print('=' * 74)
print('WHAT THE CONTROL COSTS')
print('=' * 74)
print('Partial Spearman controlling for log RNA depth and log ATAC depth removes the')
print('shared-depth path. It cannot remove a real biological association, because motif')
print('deviation z is already background matched on GC and accessibility. If the RUNX2')
print('association survives the partial, it is not a depth artifact. If it does not')
print('survive, it never was one.')
