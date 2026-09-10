# Run AFTER the current run finishes. Cheap: no motif scan, no hg38, no deviations.
# Loads each donor's RNA, reruns ONLY the labeling, and prints the full distribution.
#
# Why: fibroblast share ran 0.2% to 33.3% across seven donors of the same tissue under
# the same protocol. SSC4 gave 5 fibroblasts out of 2,760 nuclei. That is not biology,
# and the pipeline never logged where the other nuclei went, so it cannot be diagnosed
# from the run log. Same omission as the motif names.

import os, numpy as np, pandas as pd, scanpy as sc, urllib.request
sc.settings.verbosity = 0
WORK = '/content/drive/MyDrive/RR/runx_paralog'

SAMPLES = {'SSC1':'GSM9338143','SSC2':'GSM9338145','SSC3':'GSM9338147','SSC4':'GSM9338149',
 'SSC5':'GSM9338151','SSC6':'GSM9338153','SSC7':'GSM9338155','SSC8':'GSM9338157',
 'SSC9':'GSM9338159','SSC10':'GSM9338161','HC1':'GSM9338163','HC2':'GSM9338165',
 'HC3':'GSM9338167','HC4':'GSM9338169'}
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
rows, det_rows = [], []
for s, gsm in SAMPLES.items():
    fn = f'{gsm}_{s}_filtered_feature_bc_matrix.h5'
    p = f'{WORK}/{fn}'
    if not os.path.exists(p):
        try:
            urllib.request.urlretrieve(
                f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9338nnn/{gsm}/suppl/{fn}', p)
        except Exception as e:
            print(s, 'download failed', e); continue
    full = sc.read_10x_h5(p, gex_only=False); full.var_names_make_unique()
    rna = full[:, full.var.feature_types == 'Gene Expression'].copy()
    del full
    a = rna.copy()
    a.var['mt'] = a.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], inplace=True, percent_top=None, log1p=False)
    raw_n = a.n_obs
    med_g, med_mt = float(a.obs.n_genes_by_counts.median()), float(a.obs.pct_counts_mt.median())
    a = a[(a.obs.n_genes_by_counts >= 500) & (a.obs.pct_counts_mt <= 20)].copy()
    if a.n_obs < 50:
        print(f'{s}: only {a.n_obs} nuclei survive QC'); continue
    # how well is the fibroblast panel even measured in this sample?
    det_rows.append(dict(donor=s, **{g: float((np.asarray(
        (a[:, g].X > 0).sum(axis=0)).ravel()[0]) / a.n_obs) if g in a.var_names else np.nan
        for g in PANELS['fibroblast']}))
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    for n, gs in PANELS.items():
        sc.tl.score_genes(a, [g for g in gs if g in a.var_names], score_name=f's_{n}')
    S = a.obs[[f's_{n}' for n in PANELS]].to_numpy(); o = np.argsort(-S, axis=1)
    lab = np.array(list(PANELS))[o[:, 0]]
    lab[(S[np.arange(len(S)), o[:,0]] - S[np.arange(len(S)), o[:,1]]) < 0.05] = 'ambiguous'
    vc = pd.Series(lab).value_counts()
    rows.append(dict(donor=s, raw=raw_n, qc=a.n_obs, med_genes=med_g, med_mt=med_mt,
                     **{k: int(vc.get(k, 0)) for k in list(PANELS) + ['ambiguous']}))
    del a, rna

df = pd.DataFrame(rows).set_index('donor')
print('=' * 100)
print('WHERE THE NUCLEI GO, per donor')
print('=' * 100)
print(df.to_string())
print()
print('fibroblast share of QC-passing nuclei:')
print((df.fibroblast / df.qc).map('{:.1%}'.format).to_string())
print()
print('=' * 100)
print('IS THE FIBROBLAST PANEL EVEN DETECTED?  fraction of QC nuclei with a nonzero count')
print('=' * 100)
print(pd.DataFrame(det_rows).set_index('donor').map(lambda v: f'{v:.1%}').to_string())
print()
print('If COL1A1 is detected in a normal fraction of nuclei but few are LABELED')
print('fibroblast, the argmax rule is losing them to a competing panel, and the panel')
print('that wins tells you which one. If COL1A1 detection itself collapses, the sample')
print('is different and the labeling is reporting that correctly.')
