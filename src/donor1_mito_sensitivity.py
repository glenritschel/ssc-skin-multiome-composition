# Run in the SAME Colab session. Still no motifs scored.
#
# CORRECTION. The previous sweep varied MIN_GENES while holding pct_mt <= 5 FIXED.
# It therefore swept the parameter that was not doing the damage. Look at the row:
# MIN_GENES=0 kept 868 nuclei, not 6,836. The mito filter was still on the whole time.
#
# Actual attribution: mito removed 5,968 of 6,836 barcodes (87.3%). Genes removed 585
# (8.6%). Median genes per nucleus is 1,851, so this is not a shallow assay and the
# 500 floor was near harmless.
#
# This cell sweeps the guilty parameter, with an ambient-contamination control.

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
MARGIN, MIN_GENES = 0.05, 500

# AMBIENT SENTINELS. None of these belong in a fibroblast nucleus. If their detection
# rate inside the fibroblast set climbs as the mito cut is loosened, the extra nuclei
# are carrying soup, and any gain in RUNX2 is soup too, not biology.
SENTINELS = ['KRT14', 'LYZ', 'PTPRC', 'HBB']

base = rna.copy()
base.var['mt'] = base.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(base, qc_vars=['mt'], inplace=True,
                           percent_top=None, log1p=False)

def frac(a, g):
    if g not in a.var_names or a.n_obs == 0: return np.nan
    v = a[:, g].X
    v = np.asarray(v.todense()).ravel() if hasattr(v,'todense') else np.asarray(v).ravel()
    return float((v > 0).mean())

def run(max_mt):
    a = base[(base.obs.n_genes_by_counts >= MIN_GENES) &
             (base.obs.pct_counts_mt <= max_mt)].copy()
    if a.n_obs < 50: return None
    raw = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    for n, gs in PANELS.items():
        sc.tl.score_genes(a, [g for g in gs if g in a.var_names], score_name=f's_{n}')
    S = a.obs[[f's_{n}' for n in PANELS]].to_numpy()
    o = np.argsort(-S, axis=1)
    lab = np.array(list(PANELS))[o[:,0]]
    lab[(S[np.arange(len(S)),o[:,0]] - S[np.arange(len(S)),o[:,1]]) < MARGIN] = 'ambiguous'
    fibx = raw[lab == 'fibroblast'].copy()
    r = {'MAX_MT': max_mt, 'kept': a.n_obs, 'fibro': int((lab=='fibroblast').sum()),
         'med_mt': round(float(np.median(raw.obs.pct_counts_mt)), 1)}
    for g in ['RUNX1','RUNX2','RUNX3','COL1A1']:
        r[g] = round(100*frac(fibx,g), 2)
    for g in SENTINELS:
        r['amb_'+g] = round(100*frac(fibx,g), 2)
    return r

rows = [x for x in (run(m) for m in [5, 10, 20, 30, 50, 100]) if x]
df = pd.DataFrame(rows)
print('=' * 78)
print('SWEEPING THE MITOCHONDRIAL CUT (MIN_GENES held at 500)')
print('=' * 78)
print(df[['MAX_MT','kept','fibro','med_mt','RUNX1','RUNX2','RUNX3','COL1A1']].to_string(index=False))
print()
print('AMBIENT SENTINELS inside the fibroblast set (percent of nuclei nonzero).')
print('These do not belong in a fibroblast. Rising = soup.')
print(df[['MAX_MT','fibro'] + ['amb_'+g for g in SENTINELS]].to_string(index=False))

print()
print('=' * 78)
print('READING')
print('=' * 78)
lo, hi = df.iloc[0], df.iloc[-1]
d_runx2 = hi.RUNX2 - lo.RUNX2
d_amb   = np.nanmean([hi['amb_'+g] - lo['amb_'+g] for g in SENTINELS])
print(f'fibroblasts    {int(lo.fibro):>5} -> {int(hi.fibro):>5}')
print(f'RUNX2 detected {lo.RUNX2:>5.2f}% -> {hi.RUNX2:>5.2f}%   (delta {d_runx2:+.2f})')
print(f'ambient mean   {np.nanmean([lo["amb_"+g] for g in SENTINELS]):>5.2f}% -> '
      f'{np.nanmean([hi["amb_"+g] for g in SENTINELS]):>5.2f}%   (delta {d_amb:+.2f})')
print()
if d_amb > 5 and d_runx2 > 0:
    print('Ambient sentinels climb with the cut. The extra nuclei carry soup, and any')
    print('RUNX2 gain rides on it. Keep the strict cut and accept the smaller n.')
elif d_amb <= 5:
    print('Ambient sentinels are flat. The nuclei excluded by mt<=5% are not obviously')
    print('contaminated, and 87.3% of the data was being discarded by a threshold')
    print('imported from scRNA convention. Loosening it is defensible, and it would be')
    print('logged as a deviation with this table as the reason.')
print()
print('Nothing is changed here. This is evidence for a decision, not the decision.')
