# DELTA: gate zero fired, GSE195452 has no chromatin, and the RUNX motif result is already in print

**Lane** DISCOVERY (SCI)
**Clock** 2026-09-09, 13:30 EDT / 17:30 UTC
**Status** v1 preregistration withdrawn, v2 written, run not yet attempted

---

## What happened

Gate zero in `RUNX_motif_prep.ipynb` raised on cell 3. GSE195452 supplementary files
are exactly two objects plus a manifest:

```
GSE195452_Cell_metadata_v26_anno.txt.gz
GSE195452_RAW.tar
filelist.txt
```

No ATAC. No multiome. No peaks. No fragments. The gate did its job and stopped the
run before hg38 was downloaded, before MOODS was installed, and before any motif was
scored.

## A defect in my own gate

The gate listed the FTP directory and nothing else. GEO also publishes
`filelist.txt`, which enumerates the contents of `GSE*_RAW.tar`. A series whose
per-sample ATAC files sit inside that tar would have shown as NONE.

In this case the verdict is still correct, Gur 2022 is scRNA-seq and CITE-seq, but
the gate could have produced a false negative on a different series and I would have
believed it. `gate_zero_v2.py` reads `filelist.txt`.

## What replaces it

**GSE312129**, deposited 2025-12-12, published in JCI Insight, PMID 41411065:
*Multimodal analyses of early, untreated systemic sclerosis skin identify a
proinflammatory vascular niche of macrophage-fibroblast signaling.*

Paired single-nuclei multiome, snRNA plus snATAC **in the same nuclei**. Ten
treatment-naive dcSSc, four age and sex matched healthy controls, lesional forearm
punch biopsy. GSE312932 is the spatial arm of the same study.

Also found and deprioritized: **GSE99702** (Liu 2020, Nature Communications), bulk
ATAC on flow-sorted skin cell types including fibroblasts, SSc affected forearm and
unaffected back versus healthy. No RNA in the same cells, so it cannot attribute a
paralog. Retained only as a possible independent replication in a different assay.

## The finding that changes the study

The GSE312129 authors ran chromVAR themselves and report, verbatim:

> "Transcription factors downstream of the Hippo and PI3K pathways, such as RUNX1,
> CREB, and EGR1 binding, were enriched in SSc fibroblasts."

**A RUNX motif enrichment in SSc fibroblasts is already published.** That was the
primary hypothesis of the v1 preregistration. Running it would have replicated a known
result.

The authors attributed the signal to **RUNX1**. The paper does not state what
evidence separated RUNX1 from RUNX2 or RUNX3. RUNX2 and RUNX3 are not mentioned
anywhere in it. RUNX1, RUNX2 and RUNX3 share a near-identical TGTGGT core, so the
motif alone cannot make that call.

That gap is the only thing here still open, and it is the thing paired multiome can
actually address.

## Why the group comparison is close to worthless here

At 10 versus 4, the exact two-sided Mann-Whitney p under perfect separation is
**0.00200**. That is the floor. Anything short of every dcSSc donor outranking every
control returns nothing, and a null from that test carries no information.

So v2 makes the primary a **within-nucleus** question instead: across dcSSc fibroblast
nuclei, does RUNX motif deviation track RUNX1, RUNX2, or RUNX3 expression. That uses
thousands of nuclei rather than fourteen donors, and it is the question the pairing
exists to answer.

## The honest expected outcome

Most likely: **RUNX2 is not detectable in these nuclei at all.** RUNX2 is lowly
expressed outside bone and snRNA-seq is sparse. That outcome is preregistered as
"attribution impossible," not as a null about biology.

## The ceiling, unchanged and lower than before

GSE312129 is early, treatment-naive, **diffuse cutaneous**, lesional forearm.
Calcinosis is a late, limited-cutaneous-predominant feature. There is no calcinosis in
this cohort and none in any public single cell or spatial dataset, which is the
deposited Zenodo finding.

Nothing this test returns is a statement about a calcified lesion. It cannot adjudicate
Bao against Nasrallah.

## Files

- `gate_zero_v2.py` --- corrected gate, reads filelist.txt, checks all four candidates
- `RR_PREREG_runx_paralog_v2_2026-09-09.md` --- superseding preregistration
- `RR_PREREG_runx2_motif_accessibility_2026-09-09.md` --- WITHDRAWN, retain for record
- `runx_motif_test.py` --- needs rewrite for the within-nucleus primary, not yet done

## Next

Run `gate_zero_v2.py` in Colab. It is one cell and costs nothing. It reports what is
actually deposited under GSE312129, and the file format it reports determines the
loading code. Nothing further gets built until it prints.

## Sources

- JCI Insight, PMID 41411065, https://insight.jci.org/articles/view/198954
- GSE312129, https://www.omicsdi.org/dataset/geo/GSE312129
- GSE312932, https://www.omicsdi.org/dataset/geo/GSE312932
- Liu 2020, https://www.nature.com/articles/s41467-020-19702-z (GSE99702)
