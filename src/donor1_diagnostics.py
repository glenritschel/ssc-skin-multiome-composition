# Run as one Colab cell in the SAME session as the probe notebook.
# Answers three questions the probe printed past. No motifs are scored.
import numpy as np, pandas as pd

def col(a, g):
    v = a[:, g].X
    return np.asarray(v.todense()).ravel() if hasattr(v, 'todense') else np.asarray(v).ravel()

print('=' * 66)
print('1. WHERE DID THE NUCLEI GO?  237 fibroblasts is low for whole skin.')
print('=' * 66)
print('nuclei in filtered h5      :', f'{full.n_obs:,}')
print('after QC (>=500 genes, <=5% mt):', f'{q.n_obs:,}')
print()
print(q.obs['label'].value_counts().to_string())
print()
amb = (q.obs['label'] == 'ambiguous').mean()
print(f'ambiguous fraction: {amb:.1%}   (protocol section 4 fails above 40%)')
print()
print('margin distribution for nuclei called ambiguous, to see how close they were:')
S = q.obs[[c for c in q.obs.columns if c.startswith('score_')]].to_numpy()
o = np.argsort(-S, axis=1)
marg = S[np.arange(len(S)), o[:, 0]] - S[np.arange(len(S)), o[:, 1]]
am = marg[(q.obs['label'] == 'ambiguous').to_numpy()]
if len(am):
    print(pd.Series(am).describe().to_string())
    print(f'  would-be fibroblast among ambiguous: '
          f'{(np.array(list(q.obs.columns[[c.startswith("score_") for c in q.obs.columns]]))[o[:,0]][(q.obs["label"]=="ambiguous").to_numpy()] == "score_fibroblast").sum()}')

print()
print('=' * 66)
print('2. RUNX1 AND EGR1 BOTH HIT 46.84%.  COINCIDENCE, OR MISALIGNMENT?')
print('=' * 66)
r1, eg = col(fib, 'RUNX1') > 0, col(fib, 'EGR1') > 0
print(f'RUNX1 nonzero nuclei: {r1.sum()}')
print(f'EGR1  nonzero nuclei: {eg.sum()}')
print(f'overlap             : {(r1 & eg).sum()}')
print(f'RUNX1 only          : {(r1 & ~eg).sum()}')
print(f'EGR1 only           : {(~r1 & eg).sum()}')
if (r1 == eg).all():
    raise RuntimeError(
        'STOP: RUNX1 and EGR1 are nonzero in EXACTLY the same nuclei. Two independent '
        'genes do not do that. Suspect a var indexing bug or a duplicated column. '
        'Do not proceed until this is explained.')
print('\nnot identical sets. Equal counts is a coincidence at these numbers.')

print()
print('=' * 66)
print('3. HOW MUCH DOES THIS ONE DONOR ACTUALLY CARRY?')
print('=' * 66)
n2 = int((col(fib, 'RUNX2') > 0).sum())
print(f'RUNX2 positive fibroblast nuclei in SSC1: {n2} of {fib.n_obs}')
print('Per donor power at this n is modest. The power comes from combining the ten')
print('SSc donors by Fisher z, not from any single donor. SSC1 alone decides nothing.')
