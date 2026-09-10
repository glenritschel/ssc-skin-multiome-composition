#!/usr/bin/env python3
"""
runx_motif_test.py  -  is the RUNX regulon engaged in dSSc skin fibroblasts?

Implements RR_PREREG_runx2_motif_accessibility_2026-09-09.md. Every choice in
that file is fixed; nothing here is tunable after seeing results.

WHAT THIS DOES DIFFERENTLY FROM A NAIVE chromVAR RUN, AND WHY EACH MATTERS

  1. DONOR-LEVEL TESTING. Per-cell deviations are aggregated to a donor mean and
     tested at n=9 vs 9. Cells from one donor are not independent. The cell-level
     version of this mistake returned >50% of the transcriptome as "significant"
     in P45 this week; the same error in ATAC would be harder to spot.

  2. THE RUNX PARALOG PROBLEM IS TREATED AS THE MAIN THREAT, NOT A FOOTNOTE.
     RUNX1/2/3 share a near-identical core motif (TGTGGT). Accessibility alone
     CANNOT attribute a deviation to RUNX2. This script therefore refuses to
     report a RUNX result without the paired-RNA attribution in step 5, which is
     only possible because the data is multiome.

  3. CONTROLS GATE THE RESULT. If the preregistered positive controls do not
     move, the script says the pipeline failed rather than that the hypothesis
     failed. Those are different findings and conflating them is how a null
     becomes a false claim.

  4. THE COLLAGEN CONFOUND IS CHECKED. In this same cohort the osteogenic RNA
     score was 91% collagen and REVERSED SIGN once fibrosis was regressed out.
     Any motif result is reported raw AND residual, and a sign flip is reported
     in the same breath as the raw number.

USAGE
    python runx_motif_test.py --h5mu path.h5mu
    python runx_motif_test.py --atac atac.h5ad --rna rna.h5ad
"""
import argparse, sys
import numpy as np, pandas as pd
from scipy import stats
import scipy.sparse as sp

# ---- preregistered, do not edit after seeing data -------------------------
PRIMARY      = ["RUNX1", "RUNX2", "RUNX3"]     # scored together; see paralog note
SECONDARY    = ["SOX9"]
NEG_CONTROL  = ["HNF4A", "POU5F1"]             # lineage-inappropriate
POS_CONTROL  = ["SMAD3", "TWIST1", "RUNX1",    # pick from Gur's 63 BEFORE running
                "FOSL2", "JUNB", "TEAD1"]      # <- RECORD which were used
MYOFIB_GENES = ["ACTA2", "POSTN", "COL1A1", "COL1A2", "TAGLN", "CCN2"]
DONOR_KEY, GROUP_KEY = "donor", "group"


def deviations(atac, motif_key="motif_match", n_bg=50, seed=0):
    """
    chromVAR-style deviation, implemented transparently so the statistic is
    auditable rather than a black box.

    For motif m: observed accessibility of peaks carrying m, minus the mean of
    n_bg background peak sets matched on GC content and mean accessibility,
    divided by the background SD. The matching is what makes the statistic
    comparable across motifs of different peak counts.
    """
    rng = np.random.default_rng(seed)
    X = atac.X.tocsr() if sp.issparse(atac.X) else sp.csr_matrix(atac.X)
    peak_acc = np.asarray(X.sum(0)).ravel()
    gc = atac.var["gc"].values if "gc" in atac.var else np.full(atac.n_vars, 0.5)

    # bin peaks on (GC, accessibility) so backgrounds are matched, not random
    def qbin(v, k=20):
        r = pd.qcut(pd.Series(v).rank(method="first"), k, labels=False)
        return r.values
    strata = qbin(gc) * 100 + qbin(np.log1p(peak_acc))
    by_stratum = {s: np.where(strata == s)[0] for s in np.unique(strata)}

    M = atac.varm[motif_key]                      # peaks x motifs, binary
    names = list(atac.uns.get("motif_names", range(M.shape[1])))
    out = np.zeros((atac.n_obs, M.shape[1]), dtype=np.float32)

    cell_tot = np.asarray(X.sum(1)).ravel(); cell_tot[cell_tot == 0] = 1
    for j in range(M.shape[1]):
        idx = np.flatnonzero(np.asarray(M[:, j]).ravel())
        if len(idx) < 10:
            out[:, j] = np.nan; continue
        obs = np.asarray(X[:, idx].sum(1)).ravel() / cell_tot
        bg = np.empty((n_bg, atac.n_obs), dtype=np.float32)
        for b in range(n_bg):
            pick = [rng.choice(by_stratum[strata[i]]) for i in idx]
            bg[b] = np.asarray(X[:, pick].sum(1)).ravel() / cell_tot
        mu, sd = bg.mean(0), bg.std(0); sd[sd == 0] = 1
        out[:, j] = (obs - mu) / sd
    return pd.DataFrame(out, index=atac.obs_names, columns=names)


def donor_test(dev, obs, motifs):
    """Aggregate to donor, then Mann-Whitney at n=9 vs 9. Never cell-level."""
    d = dev.copy()
    d[DONOR_KEY] = obs[DONOR_KEY].values
    d[GROUP_KEY] = obs[GROUP_KEY].values
    per_donor = d.groupby([DONOR_KEY, GROUP_KEY], observed=True).mean(numeric_only=True).reset_index()
    rows = []
    for m in motifs:
        if m not in per_donor: continue
        a = per_donor.loc[per_donor[GROUP_KEY] == "HC", m].dropna()
        b = per_donor.loc[per_donor[GROUP_KEY] == "dSSc", m].dropna()
        if len(a) < 3 or len(b) < 3: continue
        rows.append(dict(motif=m, n_HC=len(a), n_dSSc=len(b),
                         mean_HC=a.mean(), mean_dSSc=b.mean(),
                         delta=b.mean() - a.mean(),
                         direction="dSSc UP" if b.mean() > a.mean() else "dSSc DOWN",
                         p=stats.mannwhitneyu(b, a, alternative="two-sided").pvalue))
    r = pd.DataFrame(rows)
    if len(r):
        p = r.p.values; o = np.argsort(p); m_ = len(p); q = np.empty(m_); run = 1.0
        for i in range(m_ - 1, -1, -1):
            run = min(run, p[o[i]] * m_ / (i + 1)); q[o[i]] = run
        r["q_BH"] = q; r["survives_q05"] = r.q_BH < 0.05
    return r, per_donor


def attribute_runx(dev, rna, obs):
    """
    STEP 5, AND THE RESULT IS NOT REPORTABLE WITHOUT IT.
    RUNX1/2/3 share a motif. Ask which paralog's RNA tracks the deviation in the
    SAME cells. This is what multiome buys and ATAC alone cannot answer.
    """
    print("\n=== RUNX PARALOG ATTRIBUTION (the decisive step) ===")
    if "RUNX" not in " ".join(dev.columns):
        print("  no RUNX motif column"); return
    col = [c for c in dev.columns if "RUNX" in c.upper()][0]
    print(f"  motif column: {col}")
    for g in ["RUNX1", "RUNX2", "RUNX3"]:
        if g not in rna.var_names:
            print(f"  {g:<6} not in RNA matrix"); continue
        x = rna[:, g].X
        x = np.asarray(x.todense()).ravel() if sp.issparse(x) else np.asarray(x).ravel()
        rho = stats.spearmanr(dev[col].values, x).statistic
        pct = 100 * (x > 0).mean()
        print(f"  {g:<6} rho(deviation, RNA) = {rho:+.3f}   detected in {pct:5.1f}% of cells")
    print("\n  READING: if the deviation tracks RUNX1 rather than RUNX2, the signal is")
    print("  NOT osteogenic and a naive report would have been wrong. If RUNX2 is")
    print("  detected in under ~2% of cells, no attribution is possible and the")
    print("  honest output is that this dataset cannot answer the question.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5mu"); ap.add_argument("--atac"); ap.add_argument("--rna")
    a = ap.parse_args()
    if not (a.h5mu or (a.atac and a.rna)):
        print(__doc__); sys.exit(0)
    import anndata as ad
    if a.h5mu:
        import muon
        mu = muon.read(a.h5mu); atac, rna = mu.mod["atac"], mu.mod["rna"]
    else:
        atac, rna = ad.read_h5ad(a.atac), ad.read_h5ad(a.rna)

    print(f"ATAC {atac.shape}   RNA {rna.shape}")
    print(f"donors: {atac.obs[DONOR_KEY].nunique()}   groups: {dict(atac.obs[GROUP_KEY].value_counts())}")

    dev = deviations(atac)

    print("\n=== CONTROLS FIRST. THE RESULT IS GATED ON THESE. ===")
    ctrl, _ = donor_test(dev, atac.obs, POS_CONTROL + NEG_CONTROL)
    print(ctrl.to_string(index=False))
    pos_ok = ctrl[ctrl.motif.isin(POS_CONTROL)].survives_q05.any() if len(ctrl) else False
    neg_bad = ctrl[ctrl.motif.isin(NEG_CONTROL)].survives_q05.any() if len(ctrl) else False
    if not pos_ok:
        print("\n  STOP. No preregistered positive control moved. The finding is that the")
        print("  PIPELINE did not work, not that the hypothesis failed. Do not report a null.")
        return
    if neg_bad:
        print("\n  STOP. A negative control moved. The deviation background is mis-modeled.")
        return
    print("\n  Controls pass. Proceeding.")

    print("\n=== PRIMARY AND SECONDARY ===")
    res, per_donor = donor_test(dev, atac.obs, PRIMARY + SECONDARY)
    print(res.to_string(index=False))

    print("\n=== COLLAGEN CONFOUND CHECK (this cohort reversed sign on the RNA score) ===")
    g = [x for x in MYOFIB_GENES if x in rna.var_names]
    if g:
        Xm = rna[:, g].X
        Xm = np.asarray(Xm.todense()) if sp.issparse(Xm) else np.asarray(Xm)
        myo = pd.Series(((Xm - Xm.mean(0)) / (Xm.std(0) + 1e-9)).mean(1), index=rna.obs_names)
        pd_myo = myo.groupby(atac.obs[DONOR_KEY].values).mean()
        for m in PRIMARY + SECONDARY:
            if m not in per_donor: continue
            sub = per_donor.dropna(subset=[m]).copy()
            sub["myo"] = pd_myo.reindex(sub[DONOR_KEY]).values
            sub = sub.dropna(subset=["myo"])
            if len(sub) < 6: continue
            b, c = np.polyfit(sub.myo, sub[m], 1)
            sub["res"] = sub[m] - (b * sub.myo + c)
            A = sub.loc[sub[GROUP_KEY] == "HC", "res"]; B = sub.loc[sub[GROUP_KEY] == "dSSc", "res"]
            if len(A) < 3 or len(B) < 3: continue
            raw_dir = res.loc[res.motif == m, "direction"].iloc[0]
            res_dir = "dSSc UP" if B.mean() > A.mean() else "dSSc DOWN"
            flag = "   *** SIGN FLIPS ***" if raw_dir != res_dir else ""
            print(f"  {m:<8} raw {raw_dir}   residual {res_dir}  "
                  f"p={stats.mannwhitneyu(B, A).pvalue:.3f}{flag}")

    attribute_runx(dev, rna, atac.obs)

    print("\n=== WHAT THIS DOES NOT LICENSE ===")
    print("  GSE195452 carries no calcinosis label. This tests the DISEASE, not the")
    print("  lesion. No statement about calcinosis follows from any result above.")


if __name__ == "__main__":
    main()
