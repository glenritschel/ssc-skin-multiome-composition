"""
RUNX paralog attribution in dcSSc skin fibroblasts, GSE312129.

Implements RR_PREREG_runx_paralog_v2_2026-09-09.md as amended by
RR_DEVIATION_2_2026-09-09.md and its same-day amendment, and
RR_LABELING_PROTOCOL_2026-09-09.md.

Runs as a separate process. Nothing in a notebook kernel can override the
controls below. One donor at a time, checkpointed, so a Colab timeout costs
one donor and not the run.

    python runx_paralog_pipeline.py --work /content/drive/MyDrive/RR/runx_paralog \
                                    --cache /content/drive/MyDrive/RR/_cache
"""
import argparse, gzip, hashlib, json, os, sys, urllib.request
import numpy as np, pandas as pd, scipy.sparse as sp
from scipy.stats import spearmanr, rankdata, norm

# ---------------------------------------------------------------- FIXED
SAMPLES = {                                   # protocol section 7
 'SSC1':'GSM9338143','SSC2':'GSM9338145','SSC3':'GSM9338147','SSC4':'GSM9338149',
 'SSC5':'GSM9338151','SSC6':'GSM9338153','SSC7':'GSM9338155','SSC8':'GSM9338157',
 'SSC9':'GSM9338159','SSC10':'GSM9338161','HC1':'GSM9338163','HC2':'GSM9338165',
 'HC3':'GSM9338167','HC4':'GSM9338169',
}
GROUP = {s: ('dcSSc' if s.startswith('SSC') else 'HC') for s in SAMPLES}

PARALOGS   = ['RUNX1', 'RUNX2', 'RUNX3']      # preregistration section 2
RUNX_MOTIF = ['RUNX1', 'RUNX2', 'RUNX3']      # JASPAR names, averaged into one family score
# positive controls taken from the paper's OWN reported enrichments, not chosen by us
POS_CTRL   = {'SMAD3':'SMAD3', 'EGR1':'EGR1', 'CREB1':'CREB1', 'ELK4':'ELK4'}
NEG_CTRL   = {'HNF4A':'HNF4A', 'POU5F1':'POU5F1'}

PANELS = {
 'fibroblast'  : ['COL1A1','COL1A2','COL3A1','DCN','LUM','PDGFRA','FBLN1'],
 'keratinocyte': ['KRT14','KRT5','KRT1','KRT10','KRT6A'],
 'endothelial' : ['PECAM1','VWF','CDH5','EGFL7'],
 'pericyte_smc': ['RGS5','MYH11','NOTCH3','KCNJ8'],   # ACTA2/TAGLN deliberately absent
 'myeloid'     : ['LYZ','CD68','ITGAM','CSF1R','AIF1'],
 't_nk'        : ['CD3D','CD3E','IL7R','CD2','TRAC'],
 'b_plasma'    : ['MS4A1','CD79A','JCHAIN','MZB1'],
 'mast'        : ['TPSAB1','TPSB2','CPA3','MS4A2'],
 'melanocyte'  : ['PMEL','MLANA','TYRP1','DCT'],
 'adnexal'     : ['KRT7','AQP5','SCGB2A2','DCD'],
}
MIN_GENES   = 500      # protocol s3
MAX_MT      = 20       # DEVIATION 2 change 1
MARGIN      = 0.05     # protocol s3
MIN_FIBRO   = 200      # protocol s4
DETECT_FLOOR= 0.05     # protocol s5
N_BG        = 50
SEED        = 0
NULL_PANEL  = 80       # DEVIATION 3: empirical null, seeded random draw
NULL_MIN_DETECT = 0.05 # a panel motif's own TF must be detected this often
# ----------------------------------------------------------------------


def log(*a): print(*a, flush=True)


def fetch(gsm, sample, work):
    fn  = f'{gsm}_{sample}_filtered_feature_bc_matrix.h5'
    dst = os.path.join(work, fn)
    if not os.path.exists(dst):
        url = f'https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9338nnn/{gsm}/suppl/{fn}'
        log(f'    downloading {fn}')
        urllib.request.urlretrieve(url, dst)
    if os.path.getsize(dst) < 1_000_000:
        raise RuntimeError(f'{fn} is implausibly small')
    return dst


def split_modalities(path):
    import scanpy as sc
    full = sc.read_10x_h5(path, gex_only=False)
    full.var_names_make_unique()
    if 'feature_types' not in full.var:
        raise RuntimeError('no feature_types; not a multiome h5')
    ft = full.var['feature_types']
    if 'Peaks' not in set(ft) or 'Gene Expression' not in set(ft):
        raise RuntimeError(f'missing a modality: {sorted(set(ft))}')
    rna  = full[:, ft == 'Gene Expression'].copy()
    atac = full[:, ft == 'Peaks'].copy()
    assert (rna.obs_names == atac.obs_names).all()
    return rna, atac


def label_fibroblasts(rna):
    """Protocol sections 2 and 3. Per nucleus, rule based. No cluster is inspected."""
    import scanpy as sc
    a = rna.copy()
    a.var['mt'] = a.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(a, qc_vars=['mt'], inplace=True,
                               percent_top=None, log1p=False)
    keep = (a.obs.n_genes_by_counts >= MIN_GENES) & (a.obs.pct_counts_mt <= MAX_MT)
    a = a[keep].copy()
    raw = a.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    for n, gs in PANELS.items():
        present = [g for g in gs if g in a.var_names]
        if not present:
            raise RuntimeError(f'panel {n} has no genes present')
        sc.tl.score_genes(a, present, score_name=f's_{n}')
    S = a.obs[[f's_{n}' for n in PANELS]].to_numpy()
    o = np.argsort(-S, axis=1)
    lab = np.array(list(PANELS))[o[:, 0]]
    lab[(S[np.arange(len(S)), o[:,0]] - S[np.arange(len(S)), o[:,1]]) < MARGIN] = 'ambiguous'
    return raw, lab


def parse_peaks(names):
    rows = []
    for p in names:
        parts = str(p).replace(':', '-').rsplit('-', 2)
        if len(parts) != 3:
            raise RuntimeError(f'cannot parse peak {p!r}')
        rows.append((parts[0], int(parts[1]), int(parts[2])))
    return pd.DataFrame(rows, columns=['chrom', 'start', 'end'])


def motif_matrix(peaks, genome, motifs, keep_names, require_human=True):
    """peaks x motifs boolean hit matrix, per-peak GC, and full provenance.

    Two failure modes the first run exposed:
      - JASPAR CORE holds MORE THAN ONE matrix under a single name (HNF4A appeared
        twice). Keying anything by name silently drops one.
      - JASPAR names carry legacy capitalization ('Runx1'). Case-insensitive selection
        matched it; exact-match lookup downstream then missed it, so the RUNX1 motif
        was scanned and then excluded from the family average.
    Provenance is returned so neither can happen silently again.
    """
    import MOODS.scan, MOODS.tools
    want = {k.upper() for k in keep_names}
    sel  = [m for m in motifs if m.name.upper() in want]

    prov = []
    for m in sel:
        sp = getattr(m, 'species', None)
        human = (sp is None) or ('9606' in str(sp))
        prov.append(dict(name=m.name, key=m.name.upper(), matrix_id=m.matrix_id,
                         species=str(sp), used=bool(human or not require_human)))
    if require_human:
        drop = [p for p in prov if not p['used']]
        for p in drop:
            log(f"    DROPPED non-human matrix {p['name']} {p['matrix_id']} species={p['species']}")
        sel = [m for m, p in zip(sel, prov) if p['used']]

    missing = want - {m.name.upper() for m in sel}
    if missing:
        log(f'    WARNING no usable matrix for: {sorted(missing)}')

    mats, thr, keys = [], [], []
    for m in sel:
        pwm = MOODS.tools.log_odds([list(m.counts[b]) for b in 'ACGT'], [0.25]*4, 1)
        mats.append(pwm); keys.append(m.name.upper())
        thr.append(MOODS.tools.threshold_from_p(pwm, [0.25]*4, 1e-4))
    M  = np.zeros((len(peaks), len(keys)), dtype=bool)
    gc = np.empty(len(peaks), dtype=float)
    for i, (c, s, e) in enumerate(peaks.itertuples(index=False)):
        try:
            seq = str(genome[c][s:e]).upper()
        except Exception:
            seq = ''
        if not seq:
            gc[i] = 0.5; continue
        gc[i] = (seq.count('G') + seq.count('C')) / len(seq)
        for j, h in enumerate(MOODS.scan.scan_dna(seq, mats, [0.25]*4, thr)):
            if len(h): M[i, j] = True
    return M, gc, keys, prov


def deviations(X, M, gc, n_bg=N_BG, seed=SEED, chunk=20):
    """chromVAR-style deviation z. Background peaks matched on GC and accessibility.

    Rewritten sparse and vectorized. The previous version drew one background peak at a
    time in a Python loop, which is fine for ten motifs and does not finish for eighty.
    """
    rng = np.random.default_rng(seed)
    X = sp.csr_matrix(X)
    cell_tot = np.asarray(X.sum(axis=1)).ravel()
    peak_tot = np.asarray(X.sum(axis=0)).ravel()
    grand = peak_tot.sum()
    if grand == 0: raise RuntimeError('empty ATAC matrix')
    pfrac = peak_tot / grand
    n_peaks, n_mot = M.shape

    def qbin(v, n):
        q = np.quantile(v, np.linspace(0, 1, n + 1)[1:-1])
        return np.digitize(v, q)
    key   = qbin(gc, 20) * 20 + qbin(np.log1p(peak_tot), 20)
    members = {b: np.flatnonzero(key == b) for b in np.unique(key)}

    def dev_from(indptr_cols):
        """indptr_cols: list of index arrays, one per column."""
        rows = np.concatenate(indptr_cols) if indptr_cols else np.array([], int)
        cols = np.concatenate([np.full(len(c), j) for j, c in enumerate(indptr_cols)]) \
               if indptr_cols else np.array([], int)
        B = sp.csc_matrix((np.ones(len(rows)), (rows, cols)),
                          shape=(n_peaks, len(indptr_cols)))
        obs = np.asarray((X @ B).todense())
        exp = np.outer(cell_tot, np.asarray(pfrac @ B).ravel())
        with np.errstate(divide='ignore', invalid='ignore'):
            return np.where(exp > 0, (obs - exp) / exp, 0.0)

    hits = [np.flatnonzero(M[:, j]) for j in range(n_mot)]
    obsd = dev_from(hits)
    out  = np.zeros_like(obsd)

    for s in range(0, n_mot, chunk):
        e = min(s + chunk, n_mot)
        cols = []
        for j in range(s, e):
            h = hits[j]
            if not len(h):
                cols.extend([np.array([], int)] * n_bg); continue
            hb = key[h]
            for _ in range(n_bg):
                pick = np.empty(len(h), int)
                for b in np.unique(hb):
                    m = hb == b
                    pool = members[b]
                    pick[m] = pool[rng.integers(0, len(pool), m.sum())]
                cols.append(pick)
        bgd = dev_from(cols).reshape(X.shape[0], e - s, n_bg)
        mu, sd = bgd.mean(2), bgd.std(2)
        out[:, s:e] = np.where(sd > 0, (obsd[:, s:e] - mu) / sd, 0.0)
        del bgd
    return out


def partial_spearman(x, y, Z):
    """Spearman of x,y controlling for columns of Z. Rank all, regress out, Pearson."""
    r  = lambda v: rankdata(v)
    Zr = np.column_stack([r(Z[:, i]) for i in range(Z.shape[1])] + [np.ones(len(x))])
    def resid(v):
        beta, *_ = np.linalg.lstsq(Zr, r(v), rcond=None)
        return r(v) - Zr @ beta
    a, b = resid(x), resid(y)
    if a.std() == 0 or b.std() == 0: return np.nan
    return float(np.corrcoef(a, b)[0, 1])


def fisher_combine(rows):
    """Inverse variance on z, w = k - 3, which is what the preregistration's own
    reasoning implies once per donor cell counts differ."""
    rows = [r for r in rows if not np.isnan(r['rho']) and r['k'] > 4]
    if not rows: return dict(rho=np.nan, p=np.nan, n_donors=0)
    z = np.array([np.arctanh(np.clip(r['rho'], -0.999999, 0.999999)) for r in rows])
    w = np.array([r['k'] - 3 for r in rows], dtype=float)
    zbar = (w * z).sum() / w.sum()
    se   = 1 / np.sqrt(w.sum())
    return dict(rho=float(np.tanh(zbar)), p=float(2 * norm.sf(abs(zbar / se))),
                n_donors=len(rows))


def bh(p):
    p = np.asarray(p, float); ok = ~np.isnan(p)
    q = np.full_like(p, np.nan); v = p[ok]; n = len(v)
    if n == 0: return q
    o = np.argsort(v); adj = np.empty(n)
    run = 1.0
    for i in range(n - 1, -1, -1):
        run = min(run, v[o[i]] * n / (i + 1)); adj[o[i]] = run
    q[ok] = adj
    return q


def run_donor(sample, gsm, work, genome, motifs):
    ck = os.path.join(work, f'ck_{sample}.json')
    if os.path.exists(ck):
        log(f'  {sample}: checkpoint found, skipping')
        return json.load(open(ck))

    log(f'  {sample} ({GROUP[sample]})')
    rna, atac = split_modalities(fetch(gsm, sample, work))
    raw, lab = label_fibroblasts(rna)
    n_fib = int((lab == 'fibroblast').sum())
    rec = dict(sample=sample, group=GROUP[sample], n_nuclei=int(raw.n_obs), n_fibro=n_fib,
               ambiguous=float((lab == 'ambiguous').mean()))
    log(f'    {raw.n_obs:,} nuclei after QC, {n_fib:,} fibroblasts, '
        f'{rec["ambiguous"]:.1%} ambiguous')

    if n_fib < MIN_FIBRO:
        rec.update(excluded=f'fibroblasts {n_fib} < {MIN_FIBRO} floor')
        json.dump(rec, open(ck, 'w'), indent=1); return rec

    fib  = raw[lab == 'fibroblast'].copy()
    fatac = atac[fib.obs_names].copy()

    # DEVIATION 3: empirical null. Every human CORE motif whose OWN transcription factor
    # is detected in >= NULL_MIN_DETECT of this donor's fibroblasts is a candidate; a
    # seeded random NULL_PANEL of them is scored by the same code on the same nuclei.
    def _detect(g):
        if g not in fib.var_names: return 0.0
        v = fib[:, g].X
        v = np.asarray(v.todense()).ravel() if hasattr(v, 'todense') else np.asarray(v).ravel()
        return float((v > 0).mean())

    fixed = {n.upper() for n in RUNX_MOTIF} | {n.upper() for n in POS_CTRL} | \
            {n.upper() for n in NEG_CTRL}
    cand = sorted({m.name.upper() for m in motifs
                   if (getattr(m, 'species', None) is None or '9606' in str(m.species))
                   and m.name.upper() not in fixed
                   and m.name.upper() in set(fib.var_names)
                   and _detect(m.name.upper()) >= NULL_MIN_DETECT})
    panel = list(np.random.default_rng(SEED).permutation(cand)[:NULL_PANEL])
    log(f'    null panel: {len(cand)} eligible motifs, {len(panel)} drawn (seed {SEED})')
    rec['null_panel_names'] = panel

    pk = parse_peaks(fatac.var_names)
    M, gc, keys, prov = motif_matrix(pk, genome, motifs,
                                     RUNX_MOTIF + list(POS_CTRL) + list(NEG_CTRL) + panel)
    log(f'    {M.shape[0]:,} peaks')
    for p in prov:
        if p['key'] in fixed:
            log(f"      {p['name']:<8} {p['matrix_id']:<12} species={p['species']:<28} "
                f"{'USED' if p['used'] else 'dropped'}")
    rec['motif_provenance'] = prov
    D = deviations(fatac.X, M, gc)
    # duplicates under one name are AVERAGED, never silently dropped
    dev = {}
    for k in dict.fromkeys(keys):
        cols = [i for i, kk in enumerate(keys) if kk == k]
        if len(cols) > 1:
            log(f'    {k}: {len(cols)} matrices, averaged')
        dev[k] = D[:, cols].mean(axis=1)

    runx_cols = [dev[n.upper()] for n in RUNX_MOTIF if n.upper() in dev]
    log(f'    RUNX family motifs in the average: '
        f'{[n for n in RUNX_MOTIF if n.upper() in dev]}')
    if not runx_cols:
        raise RuntimeError('no RUNX motif found in JASPAR CORE')
    runx_dev = np.mean(np.column_stack(runx_cols), axis=1)
    if len(runx_cols) > 1:
        cc = np.corrcoef(np.column_stack(runx_cols).T)
        rec['runx_motif_intercorr'] = float(cc[np.triu_indices(len(runx_cols), 1)].min())
        log(f"    RUNX motif deviations intercorrelate at "
            f"{rec['runx_motif_intercorr']:+.3f} (averaging assumption)")

    import scanpy as sc
    ex = fib.copy(); sc.pp.normalize_total(ex, target_sum=1e4); sc.pp.log1p(ex)
    def expr(g):
        if g not in ex.var_names: return None
        v = ex[:, g].X
        return np.asarray(v.todense()).ravel() if hasattr(v, 'todense') else np.asarray(v).ravel()

    rna_depth  = np.log1p(np.asarray(fib.X.sum(axis=1)).ravel())
    atac_depth = np.log1p(np.asarray(fatac.X.sum(axis=1)).ravel())
    Z = np.column_stack([rna_depth, atac_depth])

    # DEVIATION 2 amendment: leg two of the confound, measured
    rec['dev_vs_atac_depth'] = float(spearmanr(runx_dev, atac_depth).statistic)

    out = {}
    for g in PARALOGS:                                   # primary
        v = expr(g)
        if v is None: out[g] = dict(rho=np.nan, note='gene absent'); continue
        d = float((v > 0).mean())
        if d < DETECT_FLOOR:
            out[g] = dict(rho=np.nan, detect=d, note='below detectability floor'); continue
        out[g] = dict(rho=partial_spearman(runx_dev, v, Z), detect=d,
                      rho_uncontrolled=float(spearmanr(runx_dev, v).statistic),
                      k=int(fib.n_obs))
    for tf, mot in {**POS_CTRL, **NEG_CTRL}.items():     # controls: own motif vs own RNA
        v = expr(tf)
        if v is None or mot.upper() not in dev:
            out[tf] = dict(rho=np.nan, note='motif or gene absent'); continue
        d = float((v > 0).mean())
        out[tf] = dict(rho=partial_spearman(dev[mot.upper()], v, Z), detect=d,
                       k=int(fib.n_obs), kind='pos' if tf in POS_CTRL else 'neg')
    # the empirical null, same statistic, same nuclei, same code
    null = {}
    for g in panel:
        if g not in dev: continue
        v = expr(g)
        if v is None: continue
        null[g] = partial_spearman(dev[g], v, Z)
    nv = np.array([x for x in null.values() if not np.isnan(x)])
    rec['null_panel'] = {k: (None if np.isnan(v) else float(v)) for k, v in null.items()}
    if len(nv) >= 20:
        rec['null_median'] = float(np.median(nv))
        rec['null_iqr']    = [float(np.percentile(nv, 25)), float(np.percentile(nv, 75))]
        for g in PARALOGS:
            r = out.get(g, {}).get('rho')
            if r is not None and not (isinstance(r, float) and np.isnan(r)):
                out[g]['null_pct'] = float((nv < r).mean() * 100)
        log(f'    empirical null over {len(nv)} motifs: median {np.median(nv):+.4f}, '
            f'IQR [{np.percentile(nv,25):+.4f}, {np.percentile(nv,75):+.4f}], '
            f'max {nv.max():+.4f}')
        for g in PARALOGS:
            pc = out.get(g, {}).get('null_pct')
            if pc is not None:
                log(f'      {g} sits at the {pc:.1f}th percentile of that null')
        if np.percentile(nv, 75) - np.percentile(nv, 25) < 0.01:
            log('    WARNING the null has almost no spread. The pipeline may not detect')
            log('    this class of relationship at all. See DEVIATION 3.')
    else:
        log(f'    null panel produced only {len(nv)} usable motifs, below 20. '
            'No percentile computed.')

    # DEVIATION 4. The null and the RUNX numbers were not built the same way. Each null
    # motif uses ONE deviation vector against its own TF. The RUNX numbers use the MEAN
    # of two motif deviations, which suppresses noise and inflates rho relative to the
    # null. Add the matched statistic: RUNX2's own single motif against RUNX2 mRNA.
    # No equivalent exists for RUNX1, because JASPAR2024 CORE has no human RUNX1 matrix.
    if 'RUNX2' in dev:
        v = expr('RUNX2')
        if v is not None and float((v > 0).mean()) >= DETECT_FLOOR:
            r_single = partial_spearman(dev['RUNX2'], v, Z)
            out['RUNX2']['rho_own_motif_only'] = r_single
            if len(nv) >= 20:
                out['RUNX2']['null_pct_matched_construction'] = float((nv < r_single) * 1.0).mean() * 100 \
                    if False else float((nv < r_single).mean() * 100)
                log(f'    RUNX2 own motif only, matched to how the null is built: '
                    f'{r_single:+.4f}, {out["RUNX2"]["null_pct_matched_construction"]:.1f}th percentile')
                log('      This is the comparable number. The family average is not.')

    # detection matched null: RUNX2 is seen in 14.6% of nuclei; most panel TFs are not.
    det_map = {g: _detect(g) for g in null}
    if len(nv) >= 20 and 'RUNX2' in out and out['RUNX2'].get('detect'):
        d0 = out['RUNX2']['detect']
        near = sorted(det_map, key=lambda g: abs(det_map[g] - d0))[:20]
        sub = np.array([null[g] for g in near if not np.isnan(null[g])])
        if len(sub) >= 10:
            for stat, lbl in ((out['RUNX2'].get('rho_own_motif_only'), 'own motif'),
                              (out['RUNX2'].get('rho'), 'family avg')):
                if stat is not None and not np.isnan(stat):
                    log(f'    RUNX2 {lbl} vs {len(sub)} detection matched null motifs '
                        f'(detection {min(det_map[g] for g in near):.1%} to '
                        f'{max(det_map[g] for g in near):.1%}): '
                        f'{(sub < stat).mean()*100:.0f}th percentile')
            out['RUNX2']['null_pct_detection_matched'] = float(
                (sub < out['RUNX2'].get('rho_own_motif_only', np.nan)).mean() * 100)

    # DEVIATION 5. The RUNX2 and RUNX3 motif deviations intercorrelate at only ~0.18
    # despite PWM similarity of 0.822. They are not two noisy readings of one quantity,
    # so averaging them mixes two signals rather than denoising one. The family average
    # is demoted to a reported sensitivity. The primary grid uses ONE motif at a time as
    # the predictor, which is both a fair head to head across paralogs (identical
    # predictor) and matched to how every null motif is built.
    #
    # Chosen knowing only the intercorrelation, not knowing any single motif result for
    # RUNX1, which has no human matrix and can only appear on the expression side.
    grid = {}
    for mot in [m for m in ['RUNX2', 'RUNX3'] if m in dev]:
        for g in PARALOGS:
            v = expr(g)
            if v is None: continue
            d = float((v > 0).mean())
            if d < DETECT_FLOOR:
                grid[f'{mot}_motif~{g}_rna'] = dict(rho=None, detect=d,
                                                    note='below detectability floor')
                continue
            r = partial_spearman(dev[mot], v, Z)
            cell = dict(rho=r, detect=d)
            if len(nv) >= 20 and not np.isnan(r):
                cell['null_pct'] = float((nv < r).mean() * 100)
            grid[f'{mot}_motif~{g}_rna'] = cell
    rec['grid'] = grid
    log('    ' + '-' * 60)
    log('    DEVIATION 5 grid, one motif as predictor, matched to the null construction')
    log(f'    {"predictor":<12} {"mRNA":<8} {"rho":>9} {"null pct":>10}  note')
    for k, c in grid.items():
        mot, g = k.split('_motif~')[0], k.split('~')[1].replace('_rna', '')
        r, pc = c.get('rho'), c.get('null_pct')
        log(f'    {mot:<12} {g:<8} {"nan" if r is None else format(r,"+.4f"):>9} '
            f'{"" if pc is None else format(pc,".0f")+"th":>10}  {c.get("note","")}')

    rec['results'] = out
    json.dump(rec, open(ck, 'w'), indent=1)
    log(f"    rho(RUNX deviation, ATAC depth) = {rec['dev_vs_atac_depth']:+.3f}   "
        f"(leg two of the depth confound)")
    log('    ' + '-' * 60)
    log(f'    {"target":<8} {"rho":>8} {"uncontrolled":>13} {"detect":>8}  note')
    for g in PARALOGS + list(POS_CTRL) + list(NEG_CTRL):
        r = out.get(g, {})
        rho = r.get('rho'); unc = r.get('rho_uncontrolled'); det = r.get('detect')
        pc = r.get('null_pct')
        log(f'    {g:<8} {"nan" if rho is None or (isinstance(rho,float) and np.isnan(rho)) else format(rho,"+.4f"):>8} '
            f'{"" if unc is None else format(unc,"+.4f"):>13} '
            f'{"" if det is None else format(det,".1%"):>8}  '
            f'{"" if pc is None else format(pc,".0f")+"th pct":>10}  {r.get("note","")}')
    return rec


def main():
    # Guard: this file is a SCRIPT. Pasting it into a notebook cell hands argparse the
    # kernel's own argv and produces an unreadable SystemExit inside IPython. Same class
    # of failure as running the test in-kernel: say so instead of erroring obscurely.
    if 'ipykernel' in sys.argv[0] or 'colab_kernel' in sys.argv[0]:
        raise RuntimeError(
            'This is a script, not a notebook cell. Running it in-kernel would also let '
            'anything defined in the kernel override the preregistered controls, which is '
            'the whole reason it is a separate process.\n\n'
            'Save it to Drive and run:\n'
            "  !python {WORK}/runx_paralog_pipeline.py --work {WORK} --cache {CACHE} 2>&1 | tee {WORK}/run.log\n"
            'Add --only SSC1 to do one donor first.')
    ap = argparse.ArgumentParser()
    ap.add_argument('--work', required=True)
    ap.add_argument('--cache', required=True)
    ap.add_argument('--only', default=None, help='comma separated sample names')
    a = ap.parse_args()
    os.makedirs(a.work, exist_ok=True); os.makedirs(a.cache, exist_ok=True)

    log('script sha256:', hashlib.sha256(open(__file__, 'rb').read()).hexdigest()[:16])
    log(f'FIXED: MIN_GENES={MIN_GENES} MAX_MT={MAX_MT} MARGIN={MARGIN} '
        f'MIN_FIBRO={MIN_FIBRO} DETECT_FLOOR={DETECT_FLOOR} N_BG={N_BG} SEED={SEED}')
    log(f'PARALOGS={PARALOGS}  POS_CTRL={list(POS_CTRL)}  NEG_CTRL={list(NEG_CTRL)}')

    from pyfaidx import Fasta
    from pyjaspar import jaspardb
    fa = os.path.join(a.cache, 'hg38.fa')
    if not os.path.exists(fa):
        raise FileNotFoundError(f'{fa} missing. Download hg38 into the cache first.')
    genome = Fasta(fa)
    motifs = jaspardb(release='JASPAR2024').fetch_motifs(collection='CORE',
                                                         tax_group=['vertebrates'])
    log(f'genome contigs {len(genome.keys())}, JASPAR motifs {len(motifs)}\n')

    want = a.only.split(',') if a.only else list(SAMPLES)
    recs = [run_donor(s, SAMPLES[s], a.work, genome, motifs) for s in want]

    log('\n' + '=' * 72)
    log('COMBINED, dcSSc donors only, Fisher z, inverse variance w = k - 3')
    log('=' * 72)
    ssc = [r for r in recs if r['group'] == 'dcSSc' and 'results' in r]
    if len(ssc) < 6:
        log(f'STOP: {len(ssc)} dcSSc donors survived the fibroblast floor. '
            'Protocol section 4 requires at least 6. Report this, do not interpret.')
        return

    names = PARALOGS + list(POS_CTRL) + list(NEG_CTRL)
    rows = []
    for g in names:
        rows.append(dict(target=g, **fisher_combine(
            [dict(rho=r['results'][g].get('rho', np.nan),
                  k=r['results'][g].get('k', 0)) for r in ssc if g in r['results']])))
    df = pd.DataFrame(rows)
    df['q'] = bh(df.p.values)
    log(df.to_string(index=False))

    pcts = {g: [r['results'][g]['null_pct'] for r in ssc
                if g in r['results'] and 'null_pct' in r['results'][g]] for g in PARALOGS}
    log('\nDEVIATION 3, position in the empirical null of within-cell motif to own')
    log('expression correlation, one panel per donor:')
    for g, v in pcts.items():
        if v:
            log(f'  {g:<7} median {np.median(v):5.1f}th percentile   per donor '
                f'{[round(x) for x in v]}')
        else:
            log(f'  {g:<7} no percentile available')
    meds = [r['null_median'] for r in ssc if 'null_median' in r]
    if meds:
        log(f'  null median across donors: {np.median(meds):+.4f}')
    if pcts.get('RUNX2') and pcts.get('RUNX1'):
        d = np.median(pcts['RUNX2']) - np.median(pcts['RUNX1'])
        log(f'  RUNX2 minus RUNX1, median percentile gap: {d:+.1f}')
        log('  A gap near zero means the two paralogs are not separable here, which is')
        log('  the preregistered "uninformative" outcome, not evidence for either.')

    log('\nleg two of the depth confound, rho(RUNX deviation, ATAC depth) per donor:')
    legs = [r['dev_vs_atac_depth'] for r in ssc if 'dev_vs_atac_depth' in r]
    log(f'  median {np.median(legs):+.3f}   range {min(legs):+.3f} to {max(legs):+.3f}')
    if abs(np.median(legs)) < 0.05:
        log('  Below 0.05. The depth covariate was unnecessary. State that in the writeup,')
        log('  and report that the uncontrolled result is materially identical.')
    elif abs(np.median(legs)) > 0.20:
        log('  Above 0.20. The covariate was load bearing. State that in the writeup.')

    pos_ok = df[(df.target.isin(POS_CTRL)) & (df.q < 0.05)].shape[0]
    log(f'\npositive controls surviving q<0.05: {pos_ok} of {len(POS_CTRL)}')
    if pos_ok == 0:
        log('PIPELINE FAILURE. No positive control moved. Do not interpret RUNX.')
        log('Preregistration section 5, row one. Report as pipeline failure and stop.')
        return
    neg_bad = df[(df.target.isin(NEG_CTRL)) & (df.q < 0.05)].shape[0]
    if neg_bad:
        log(f'WARNING: {neg_bad} negative control(s) survived. The null is not behaving.')

    log('\nREMINDER, preregistration section 6. GSE312129 is early, treatment naive,')
    log('DIFFUSE cutaneous, lesional forearm. There is no calcinosis in this cohort.')
    log('Nothing above is a statement about a calcified lesion.')
    pd.DataFrame(rows).to_csv(os.path.join(a.work, 'combined_results.csv'), index=False)
    log(f'\nwrote {os.path.join(a.work, "combined_results.csv")}')


if __name__ == '__main__':
    main()
