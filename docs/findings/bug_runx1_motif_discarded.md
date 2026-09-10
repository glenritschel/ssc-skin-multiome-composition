# BUG: the RUNX1 motif was scanned and then silently discarded from the family score

**Found** 2026-09-09 from the SSC1 run log
**Severity** Invalidates the SSC1 numbers. No result was reported externally.

## What the log said

```
motifs scored: ['SMAD3','RUNX2','CREB1','EGR1','ELK4','HNF4A','HNF4A','POU5F1','Runx1','RUNX3']
RUNX1 0.0740546922950557  RUNX2 0.08949657127925355  RUNX3 nan
```

Two defects are visible in that one line.

## Defect 1: RUNX1 was scanned, then excluded from the thing it was scanned for

Selection was case insensitive:

```python
sel = [m for m in motifs if m.name.upper() in {k.upper() for k in keep_names}]
```

so JASPAR's legacy `Runx1` matched. Lookup afterward was case **sensitive**:

```python
dev = {n: D[:, i] for i, n in enumerate(mnames)}      # key 'Runx1'
runx_cols = [dev[n] for n in RUNX_MOTIF if n in dev]  # asks for 'RUNX1'
```

`'RUNX1' in dev` is False. **The RUNX family deviation score was the average of the
RUNX2 and RUNX3 motifs only.** The RUNX1 motif was computed and thrown away.

The `RUNX1 0.074` in the log is not the RUNX1 motif. It is RUNX1 **expression**
correlated against a RUNX-family deviation that excluded the RUNX1 motif.

## Defect 2: duplicate names, one silently dropped

`HNF4A` appears twice. JASPAR CORE holds more than one matrix under a single name.
A dict keyed by name keeps whichever came last. Only a negative control here, but the
same construct would have silently dropped a duplicate RUNX matrix.

## Defect 3: species was never checked

`tax_group=['vertebrates']` returns human and mouse matrices. Lowercase `Runx1` is
JASPAR's legacy naming and does not by itself say which species. **The pipeline had no
check.** Answering a human paralog question with a mouse position weight matrix would
have produced a number that looked entirely normal.

Run `jaspar_inspect.py` to see matrix IDs and species for all nine names.

## Fixes

1. `motif_matrix` returns full provenance, name, matrix ID, species, used or dropped,
   and the per donor log prints every row. No motif enters a score unnamed.
2. Non-human matrices are dropped by default, with the drop logged by matrix ID.
3. Duplicate names are **averaged**, never dropped, with the count logged.
4. All keys are uppercase, so selection and lookup cannot disagree again.
5. The per donor log now prints all four positive controls, both negative controls, the
   uncontrolled rho beside the partial, and `rho(RUNX deviation, ATAC depth)`. Those
   were computed in the first run and never shown, which is how a wrong family score sat
   unnoticed behind two plausible looking numbers.

## Delete before rerunning

`ck_SSC1.json`. It holds results from the broken family score and the pipeline skips any
donor with a checkpoint.

## What survives from the first run

The infrastructure. 5,761 nuclei after QC, 1,296 fibroblasts at 4.4 percent ambiguous,
62,956 peaks, the scan completing, and the six donor floor stopping the combine at one.
Those are real. Every rho is void.

## The pattern

This is the third time in this thread that a **string match** has failed silently: the
gate zero hint list, the barred term check earlier this week, and now motif names. The
common fix each time was to print what matched rather than that something matched.
