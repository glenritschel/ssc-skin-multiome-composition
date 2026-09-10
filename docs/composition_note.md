---
title: "Systemic sclerosis skin yields six times fewer fibroblast nuclei than healthy skin in the only public paired multiome of SSc skin, and that asymmetry sits underneath every between group claim drawn from it"
author:
  - Glen Charles Ritschel (Ritschel Research, ORCID 0009-0006-5341-8840)
  - Claude (Anthropic)
date: "10 September 2026"
---

## Summary

GSE312129 is the single nuclei multiome arm of a 2025 study of early, treatment naive
diffuse cutaneous systemic sclerosis skin, ten patients and four healthy controls,
lesional forearm punch biopsy (Ashton et al., JCI Insight, PMID 41411065).

Reprocessing all fourteen samples from the deposited filtered matrices with one fixed
labeling rule finds that **fibroblast nuclei are recovered from SSc skin at roughly one
sixth the rate of healthy skin**: median 6.1 percent of quality controlled nuclei in
dcSSc against 38.3 percent in controls, exact Mann-Whitney p = 0.014. The direction is
opposite to the disease. SSc skin is the fibrotic tissue.

The effect is corroborated independently by transcript detection rather than by cell
labels alone. COL1A1 is detected in a median 31.8 percent of dcSSc nuclei against 78.6
percent of control nuclei (p = 0.054). In the two extreme samples, SSC4 and SSC6,
COL1A1 falls to 7.9 and 7.4 percent and those libraries are 85 and 82 percent
keratinocyte. The fibroblasts are not mislabeled. They are not in the tube.

That composition asymmetry is a confound for any dcSSc versus healthy comparison built
on this dataset, including the transcription factor motif enrichments the source study
reports in fibroblasts. Such comparisons contrast a small, recovery selected minority of
SSc fibroblasts against a large and much less selected majority of control fibroblasts.

Dissociation and preparation bias in single cell and single nuclei workflows is
established in general (Denisenko et al., Genome Biology 2020;21:130). What is recorded
here is a specific, quantified, and **differentially distributed** instance of it in a
dataset now entering use as a reference for SSc skin.

## Figure

![Fibroblast recovery across the fourteen GSE312129 donors. Left, absolute fibroblast nuclei on a log scale, with the 200 nucleus analysis floor that was fixed on 9 September 2026 before any donor beyond the first was loaded. Right, fibroblast share against COL1A1 detection, showing that the label and the transcript agree.](composition_fig.png)

## Method

Per sample `filtered_feature_bc_matrix.h5` was taken from the per GSM supplementary
path. Gene expression features were separated from peak features. Nuclei with fewer
than 500 detected genes or more than 20 percent mitochondrial counts were dropped.

Cell identity was assigned per nucleus, never per cluster, by `scanpy.tl.score_genes`
over ten fixed marker panels with matched control gene sets, taking the highest scoring
panel and labeling a nucleus ambiguous when the top two scores differed by less than
0.05. No UMAP was inspected and no cluster was chosen by eye. The panels, thresholds and
rule were fixed in writing before any donor beyond the first was processed.

ACTA2 and TAGLN were deliberately excluded from the pericyte and smooth muscle panel,
because myofibroblasts express both and their inclusion would move exactly the nuclei of
interest into the wrong category. The pericyte panel is RGS5, MYH11, NOTCH3, KCNJ8. The
fibroblast panel is COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA, FBLN1.

The 20 percent mitochondrial threshold replaced an initial 5 percent after the stricter
cut was found to discard 87.3 percent of barcodes in the first donor; the change is
logged with its evidence table in the working record and does not affect the group
comparison reported here, which is computed under the single stated threshold for all
fourteen samples.

## Results

| donor | QC nuclei | fibroblast nuclei | share | COL1A1 detected |
|---|---|---|---|---|
| SSC1 | 5,761 | 1,296 | 22.5% | 64.7% |
| SSC2 | 3,444 | 274 | 8.0% | 49.6% |
| SSC3 | 2,269 | 139 | 6.1% | 28.2% |
| SSC4 | 2,760 | 5 | 0.2% | 7.9% |
| SSC5 | 3,820 | 572 | 15.0% | 35.7% |
| SSC6 | 1,880 | 12 | 0.6% | 7.4% |
| SSC7 | 267 | 89 | 33.3% | 68.5% |
| SSC8 | 2,801 | 170 | 6.1% | 19.7% |
| SSC9 | 6,835 | 419 | 6.1% | 35.4% |
| SSC10 | 5,000 | 239 | 4.8% | 21.5% |
| HC1 | 4,662 | 3,435 | 73.7% | 98.2% |
| HC2 | 3,993 | 2,061 | 51.6% | 88.6% |
| HC3 | 5,684 | 1,423 | 25.0% | 68.5% |
| HC4 | 5,923 | 1,027 | 17.3% | 29.3% |

Median share 6.1 percent dcSSc against 38.3 percent healthy control, exact
Mann-Whitney p = 0.014. Median COL1A1 detection 31.8 against 78.6 percent, p = 0.054.

All four controls exceed 1,000 fibroblast nuclei. Five of ten dcSSc samples fall below
200. SSC7 has the third highest fibroblast share in the study at 33.3 percent while
yielding only 89 fibroblast nuclei, because it produced 267 nuclei in total; share and
absolute count must be read together.

## Interpretation, and its limits

The likely mechanism is mechanical rather than biological. Fibrotic dermis is dense and
collagen cross linked and releases nuclei poorly, while the epidermis above it releases
them normally, so the recovered fraction is enriched for keratinocytes in proportion to
how fibrotic the dermis is. This note does not test that mechanism. It reports the
asymmetry and its size.

Two alternatives are not excluded by these data. Sampling depth within the biopsy could
vary systematically between arms for reasons unrelated to fibrosis. And a genuine
reduction in dermal fibroblast density in lesional skin, rather than a recovery
artifact, would produce the same table. Distinguishing these requires paired histology
or nuclei counts from the same biopsies, which are not deposited.

What does not depend on the mechanism is the consequence. Whatever produced it, the
fibroblast populations compared between arms in this dataset were not sampled
equivalently, and a between group difference in fibroblast state cannot be separated
from a between group difference in which fibroblasts were recovered.

## What was being attempted when this was found

This came out of a preregistered attempt to determine which RUNX paralog accompanies the
RUNX family motif accessibility that the source study reports as enriched in SSc
fibroblasts and attributes to RUNX1. That analysis specified a floor of six donors. Five
cleared it. The preregistered analysis was not run, and no RUNX result is reported here
or anywhere else.

The full record, including four logged deviations, one bug in which the RUNX1 motif was
scanned and silently discarded, and the finding that JASPAR2024 CORE contains no human
RUNX1 matrix, is retained.

## Data and code

GSE312129, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE312129

The labeling audit that produces the table and figure above is a single script operating
only on the deposited filtered matrices. No fragment files, no reference genome and no
motif database are required to reproduce it.

## Conflict of interest

GCR is an independent researcher who has filed provisional patent applications relating
to cell-state-based patient stratification in fibrotic disease, including systemic
sclerosis. **Those filings claim host cell-state measurement. This note reports a sample
recovery artifact in a public dataset and claims no marker, no method, no compound and no
indication; nothing reported here is claimed in them, and nothing here is offered in
support of any filing.**

The authors have no relationship with the authors of the study whose dataset is
reanalyzed here, and no access to it beyond the public deposit.

GCR's research interest in systemic sclerosis is personally motivated by a family member
with the disease.

## References

Ashton et al. Multimodal analyses of early, untreated systemic sclerosis skin identify a
proinflammatory vascular niche of macrophage-fibroblast signaling. JCI Insight, 2025.
PMID 41411065.

Denisenko E, Guo BB, Jones M, et al. Systematic assessment of tissue dissociation and
storage biases in single-cell and single-nucleus RNA-seq workflows. Genome Biology
2020;21:130. doi:10.1186/s13059-020-02048-6. Adult mouse kidney; reports that "cell type
composition differences are observed between single-cell and single-nucleus RNA
sequencing libraries" and that sensitive populations are "underrepresented" after warm
dissociation.
