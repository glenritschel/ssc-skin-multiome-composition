# DEVIATION 2 to the RUNX paralog preregistration

**Logged** 2026-09-09, 15:52 EDT / 19:52 UTC
**Against** `RR_PREREG_runx_paralog_v2_2026-09-09.md` section 7, and
`RR_LABELING_PROTOCOL_2026-09-09.md` section 3
**Status of the run at logging time** No motif has been scored. No RUNX association has
been computed. Nothing about the answer is known.

---

## Change 1: the mitochondrial cut moves from 5 percent to 20 percent

**What the protocol said.** Section 3, item 1: nuclei above 5 percent mitochondrial
counts are dropped, "since this is single nuclei data and high mitochondrial content
indicates ambient contamination or damage."

**What it did.** On SSC1 it removed 5,968 of 6,836 barcodes, 87.3 percent. The gene
floor removed 585, 8.6 percent. I attributed the loss to the gene floor in writing
before measuring it, and that attribution was wrong.

**The evidence table.**

| MAX_MT | nuclei kept | fibroblasts | median mt pct | RUNX2 pct | COL1A1 pct |
|---|---|---|---|---|---|
| 5 | 810 | 237 | 3.1 | 20.25 | 95.78 |
| 10 | 2,752 | 796 | 6.6 | 17.46 | 97.24 |
| 20 | 5,761 | 1,296 | 10.3 | 14.58 | 97.53 |
| 30 | 6,137 | 1,362 | 10.8 | 14.39 | 97.36 |
| 50 | 6,235 | 1,352 | 10.9 | 14.35 | 97.56 |
| 100 | 6,251 | 1,362 | 10.9 | 14.17 | 97.36 |

**Why loosening is not a contamination risk here.** Ambient RNA can only add counts.
If the admitted nuclei were carrying soup, RUNX2 detection would rise. It falls, 20.25
to 14.17. The 1,125 newly admitted fibroblasts carry RUNX2 at 12.9 percent, below the
strict set, which is the opposite of the ambient signature.

**The ambient sentinels, reported honestly.** KRT14 +4.38, LYZ +0.51, PTPRC -1.72,
HBB +4.50. My cell printed "flat" off a mean of +1.92, which averaged across opposite
signs and is a weak summary. Two sentinels rise by four to five points. That is a real
if modest increase in soup, and it is accepted because the direction of the RUNX2
change rules out soup as the driver of the quantity that matters.

**A separate finding, not a reason for this change.** HBB sits at 20.68 percent inside
the fibroblast set at the *strictest* cut. One in five called fibroblast nuclei carries
hemoglobin transcript before any loosening. That is a baseline ambient floor in this
biopsy and it is independent of the mitochondrial question. It is recorded here because
it bears on any absolute expression claim from this dataset.

**Value chosen, and why not by eye.** The table plateaus from 20 upward: 1,296 to 1,362
fibroblasts and RUNX2 14.17 to 14.58 across MAX_MT 20 through 100. The rule is the most
conservative point on the plateau, which is 20. The result is insensitive across the
plateau, so the choice carries little weight either way.

---

## Change 2: depth enters the primary statistic as a covariate

**This is the consequential one.**

Both sweeps say the same thing. RUNX2 detection tracks sequencing depth.

- Gene-floor sweep: 17.84 percent at MIN_GENES 0, rising to 21.16 percent at 1,000.
- Mito sweep: 20.25 percent in the strict, deeper set, falling to 14.17 percent when
  shallower nuclei are admitted.

COL1A1 moves 89.96 to 97.56 across both, a much smaller relative change from a much
higher base. That is the signature of a lowly expressed gene whose zeros are dropout
rather than absence. RUNX2 is lowly expressed outside bone, so this is expected.

**The confound it creates.** The preregistered primary correlates per nucleus RUNX
motif deviation against RUNX2 expression, within donor. If RUNX2 detection tracks RNA
depth, and motif deviation tracks ATAC depth, and RNA and ATAC depth are correlated
across nuclei in a multiome library, then a positive correlation appears with no
biology in it whatsoever. The preregistration did not anticipate this and would have
reported that artifact as an attribution.

**The change.** The primary becomes a **partial Spearman**, controlling per nucleus for
log1p RNA total counts and log1p ATAC total counts. Everything else is unchanged:
still within donor, still combined across donors by Fisher z, still two sided, still BH
across the three paralogs plus controls.

**Why this cannot manufacture a result.** Partialling can only remove the shared-depth
path. Motif deviation z is already background matched on GC and accessibility, so the
ATAC side is partly protected already. If the RUNX2 association survives the partial it
is not a depth artifact. If it does not survive, it never was an association.

**The diagnostic that has to run first.** `depth_confound_check.py` reports rho of
expression against depth for each paralog with COL1A1 as the reference. If RUNX2's rho
is not materially larger than COL1A1's, this change is unnecessary and will be
withdrawn rather than kept for tidiness.

---

## What is not changed

The marker panels. The 0.05 margin rule, which cost 5 fibroblasts out of 237 and is not
where the action was. The 500 gene floor, which removed 8.6 percent and is close to
harmless at a median of 1,851 genes per nucleus. The 5 percent paralog detectability
floor. Every outcome in preregistration section 5, including "attribution impossible."
The ceiling in section 6.

RUNX2 clears the detectability floor at every mitochondrial cut tested, 14.17 to 20.25
percent against a floor of 5. That conclusion does not depend on any choice made here.

---

## Errors logged against myself in reaching this

1. Named MIN_GENES as the cause of an 88 percent loss caused by the mitochondrial cut.
2. Wrote that "snMultiome RNA is shallow by construction" as justification. Median here
   is 1,851 genes per nucleus. The claim was invented to fit the suspicion.
3. Wrote a sensitivity sweep that held the guilty parameter fixed, then read
   "stable across thresholds" off it. The tell was in the output, MIN_GENES=0 keeping
   868 nuclei rather than 6,836, and I did not read it.
4. Printed a "flat" verdict on ambient sentinels from a mean that averaged across
   opposite signs.
