# Run in the SAME Colab session. Still no motifs scored.
#
# QC discarded 88.2% of barcodes (6,836 -> 810). The margin rule was innocent:
# only 5.1% ambiguous and only 5 would-be fibroblasts lost to it.
#
# The danger is specific. If MIN_GENES=500 is keeping only the highest-complexity
# nuclei, then every detection fraction measured downstream, RUNX2's 20.25% above all,
# was measured in a biased slice. High-complexity nuclei detect more genes by
# construction. This cell asks whether that headline number survives the threshold.

import numpy as np, pandas as pd, scanpy as sc

PANELS = {
 'fibroblast'  : ['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','FBLN1'],
 'keratinocyte': ['KRT14','KRT5','KRT1','KRT10','KRT6A'],
 'endothelial' : ['PECAM1','VWF','CDH5','EGFL7'],
 'pericyte_smc': ['RGS5','MYH11','NOTCH3','KCNJ8'],
 'myeloid'     : ['LYZ','CD68','ITGAM','CSF1R','AIF1'],
 't_nk'        : ['CD3D','CD3E','IL7R','CD2','TRAC'],
 'b_plasma'    : ['MS4A1','CD79A','JCHAIN','MZB1'],
 'mast'        : ['TPSAB1','TPSB2','CPA3','MS4A2'],
 'melanocyte'  : ['PMEL','MLANA','TYRP1','DCT'],
 'adnexal'     : ['KRT7','AQP5','SCGB2A2','DCD'],
}
MARGIN = 0.05

base = rna.copy()
base.var['mt'] = base.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(base, qc_vars=['mt'], inplace=True,
                           percent_top=None, log1p=False)

print('=' * 70)
print('A. WHICH CRITERION DID THE CUTTING?')
print('=' * 70)
ng, mt = base.obs.n_genes_by_counts, base.obs.pct_counts_mt
print(f'MT- genes present in var : {int(base.var["mt"].sum())}')
print(f'fails genes>=500 only    : {int(((ng < 500) & (mt <= 5)).sum()):,}')
print(f'fails mt<=5% only        : {int(((ng >= 500) & (mt > 5)).sum()):,}')
print(f'fails both               : {int(((ng < 500) & (mt > 5)).sum()):,}')
print(f'passes both              : {int(((ng >= 500) & (mt <= 5)).sum()):,}')
print()
print('genes per nucleus, all barcodes in the filtered h5:')
print(ng.describe(percentiles=[.1,.25,.5,.75,.9]).to_string())
print()
print('If the median is well under 500, the threshold is wrong for this assay,')
print('not the data. snMultiome RNA is shallow by construction.')

def label_at(min_genes):
    a = base[(base.obs.n_genes_by_counts >= min_genes) & (base.obs.pct_counts_mt <= 5)].copy()
    if a.n_obs < 50:
        return a, None
    raw = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    for n, gs in PANELS.items():
        p = [g for g in gs if g in a.var_names]
        sc.tl.score_genes(a, p, score_name=f's_{n}')
    S = a.obs[[f's_{n}' for n in PANELS]].to_numpy()
    o = np.argsort(-S, axis=1)
    top, sec = S[np.arange(len(S)), o[:,0]], S[np.arange(len(S)), o[:,1]]
    lab = np.array(list(PANELS))[o[:,0]]
    lab[(top - sec) < MARGIN] = 'ambiguous'
    return raw[lab == 'fibroblast'].copy(), int((lab == 'fibroblast').sum())

def frac(a, g):
    if g not in a.var_names or a.n_obs == 0: return np.nan
    v = a[:, g].X
    v = np.asarray(v.todense()).ravel() if hasattr(v, 'todense') else np.asarray(v).ravel()
    return float((v > 0).mean())

print()
print('=' * 70)
print('B. IS THE RUNX2 NUMBER AN ARTIFACT OF THE THRESHOLD?')
print('=' * 70)
rows = []
for mg in [0, 100, 200, 300, 500, 1000]:
    fibx, n = label_at(mg)
    if n is None:
        rows.append((mg, 0, 0, np.nan, np.nan, np.nan, np.nan)); continue
    rows.append((mg, int(((base.obs.n_genes_by_counts >= mg) &
                          (base.obs.pct_counts_mt <= 5)).sum()), n,
                 frac(fibx,'RUNX1'), frac(fibx,'RUNX2'),
                 frac(fibx,'RUNX3'), frac(fibx,'COL1A1')))
df = pd.DataFrame(rows, columns=['MIN_GENES','nuclei_kept','fibroblasts',
                                 'RUNX1','RUNX2','RUNX3','COL1A1'])
for c in ['RUNX1','RUNX2','RUNX3','COL1A1']:
    df[c] = (df[c]*100).round(2)
print(df.to_string(index=False))

print()
print('=' * 70)
print('C. READING')
print('=' * 70)
r5  = df.loc[df.MIN_GENES == 500, 'RUNX2'].iloc[0]
r2  = df.loc[df.MIN_GENES == 200, 'RUNX2'].iloc[0]
print(f'RUNX2 at MIN_GENES=500: {r5}%')
print(f'RUNX2 at MIN_GENES=200: {r2}%')
if not np.isnan(r2):
    if r2 < 5:
        print('\nRUNX2 falls below the 5% floor once the complexity bias is relaxed.')
        print('The 20.25% was a property of the threshold, not of the tissue.')
    elif r5 - r2 > 8:
        print('\nStrongly threshold dependent. Any RUNX2 claim has to state the cut it')
        print('was measured at, and the cut has to be justified independently.')
    else:
        print('\nStable across thresholds. The detection rate is a property of the data.')
print()
print('No threshold is changed by this cell. It reports; the protocol decides.')
