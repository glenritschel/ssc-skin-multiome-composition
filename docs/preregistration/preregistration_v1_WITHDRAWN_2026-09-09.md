# `-> SCI` PREREGISTRATION: does the RUNX motif open in dSSc skin fibroblasts? Fixed before any data is touched.

**Cut Wednesday 2026-09-09, container clock read before titling [Rule 10(d)].**
**Everything below is fixed in advance. No data has been downloaded, opened or inspected.**

---

# 1. WHY THIS TEST EXISTS

**The osteogenic hypothesis is a claim about a transcription factor.** Both lesional datasets measured
RNA and neither measured TF activity.

| study | what it measured | result |
|:--|:--|:--|
| Bao 2026, lupus calcinosis | fibroblast **RNA** | osteogenic and ossification programs up |
| Nasrallah 2025, rat calciphylaxis | **RNA** | *"Runx2 or Sox9 was not differentially expressed"* |

**TF transcript abundance is a poor proxy for TF activity.** A factor can be abundant and idle, or
scarce and busy. **The direct readout is whether the motif is accessible at its targets.**

**And the data exists.** Gur et al., Cell 2022, verbatim: *"we performed single-cell genome-wide
multiome profiling of fresh fibroblasts including ScAFs and other stromal cells derived from the
affected skin of 9 dSSc patients ... and 9 healthy subjects, simultaneously recovering the mRNA and
accessible chromatin regions (by ATAC-seq) from the same individual cells."* **6,710 healthy and 2,417
dSSc QC-positive cells with both modalities.**

---

# 2. ⛔ GATE ZERO, AND IT MAY STOP THIS BEFORE IT STARTS

**The paper's data availability statement, read verbatim, says only:** *"Single cell RNA-seq data that
support the findings of this study is publicly available for download in the NCBI Gene Expression
Omnibus (GEO) with accession number GEO: GSE195452."*

> ⛔ **IT DOES NOT MENTION THE ATAC. The multiome chromatin data may not be deposited under that
> accession, or at all.** Check the GSE195452 sample list for ATAC or multiome samples **before
> anything else in this file is run.**

**If the ATAC is absent, the fallback is GSE312129 and GSE312932** (JCI Insight 2025, early untreated
dcSSc, snRNA plus snATAC plus Visium plus Xenium, 10 SSc and 4 HC). **Smaller, and not paired in the
same cells, which costs section 5.**

---

# 3. THE PREDICTIONS, FIXED NOW

**Primary, two-sided, and it is a real fork.**

> **H1 (Bao-consistent):** RUNX motif accessibility deviation is **HIGHER** in dSSc fibroblasts than
> in healthy fibroblasts.
> **H0 (Nasrallah-consistent):** no difference.
> **H2 (the outcome nobody predicts and which the collagen result makes possible):** **LOWER** in dSSc.

**My prediction, on record before the run: H0. No significant difference.** Reasons, stated so they can
be held against me: the rat found RUNX2 flat at an actual calcified lesion; the GSE195452 osteogenic
RNA score reversed sign once collagen was regressed out; and this cohort has no calcinosis, so any
osteogenic program present would be a disease-wide feature rather than a lesion feature.

**Secondary.** SOX9 motif, same three-way fork, same prediction.

---

# 4. ⛔ THE LIMITATION THAT MOST INVALIDATES NAIVE VERSIONS OF THIS TEST

**RUNX1, RUNX2 and RUNX3 bind a nearly identical core motif (TGTGGT / ACCACA).** chromVAR and every
motif-deviation method **cannot tell them apart from accessibility alone.**

> ⛔ **A POSITIVE RUNX DEVIATION IS NOT EVIDENCE OF RUNX2. It is evidence that something in the RUNX
> family is active, and RUNX1 is the one with broad expression outside bone.**

**Any result reported without addressing this is uninterpretable, and this is the single most likely
way for this analysis to produce a false positive that reads as a discovery.**

---

# 5. ⭐ THE MITIGATION IS THE REASON TO USE MULTIOME AND NOT ATAC ALONE

**Paired RNA and ATAC in the same cell answers the question motif analysis cannot.**

**In cells with high RUNX motif deviation, which RUNX gene is actually transcribed?** Correlate
per-cell RUNX deviation against per-cell RUNX1, RUNX2 and RUNX3 RNA **from the same cells**.

| outcome | reading |
|:--|:--|
| deviation tracks **RUNX2** RNA | the osteogenic reading survives |
| deviation tracks **RUNX1** RNA | ⛔ the signal is not osteogenic and the naive result would have been wrong |
| tracks neither | motif accessibility is not TF-driven here; report and stop |

**This is not a robustness check bolted on at the end. It is the test.** A separate ATAC experiment
could not do it.

---

# 6. CONTROLS, FIXED IN ADVANCE

**POSITIVE CONTROL, and the analysis is void without it.** Gur reports **63 significantly
overrepresented motifs** distinguishing dSSc stromal cells. **At least one must reproduce**, or the
pipeline is broken rather than the hypothesis. **Pick the controls from their list before looking at
RUNX**, and record which were picked.

**NEGATIVE CONTROL.** A lineage-inappropriate motif with no expected dermal fibroblast role.
Preregistered choice: **HNF4A** and **POU5F1**. Either moving means the deviation background is
mis-modeled.

**CONFOUND CONTROL, and it is this session's own lesson.** The osteogenic RNA score in this same cohort
was 91 percent collagen and reversed direction once fibrosis was regressed out. **Report RUNX deviation
raw AND residual after regressing out a myofibroblast axis** (ACTA2, POSTN, COL1A1 module score or the
matched motif family). **If the sign flips, say so in the same sentence as the raw result.**

---

# 7. THE UNIT OF ANALYSIS, AND IT IS NOT NEGOTIABLE

**Per-cell deviations are computed, then aggregated to a DONOR MEAN, then tested at n = 9 versus 9 with
Mann-Whitney.**

> ⛔ **NOT cell-level.** 9,127 cells from 18 donors are not 9,127 independent observations. **This is
> the exact error found in P45 this week, where a cell-level Wilcoxon over 62,639 cells from 13 donors
> returned more than half the transcriptome as significant.** Applying it to a new modality would
> reproduce the defect in a place nobody is looking for it.

**Power, stated honestly in advance: 9 versus 9 is small.** A null result is weak evidence of absence.
**This test can falsify H1 convincingly only if the positive controls land while RUNX does not.**

**Multiple testing: BH across every motif tested, and the q value is reported beside every p.**

---

# 8. WHAT A RESULT WOULD AND WOULD NOT LICENSE

**WOULD:** a statement about whether the RUNX regulon is engaged in **diffuse SSc skin fibroblasts**.
Third independent line beside the collagen reversal and the rat result.

**WOULD NOT:** any statement about calcinosis. **GSE195452 carries no calcinosis label.** This tests
the disease, not the lesion, and the gap note's central claim is untouched either way.

**⚠ AND IT DOES NOT ENTER THE v2 PAPER.** v2 is deposited. If this produces something, it is a
separate result, and the honest framing is that a preregistered test in disease-wide skin cannot
adjudicate a lesion-level disagreement. **It can only make one side more or less plausible.**

---

# 9. RELEASE

**Nothing external until the run is complete and the controls have passed.** A preregistration is not a
result. **If the positive controls fail, the correct output is that the pipeline did not work, and that
gets written down the same as any other outcome.**
