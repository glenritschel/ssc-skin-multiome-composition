# PREREGISTRATION v2: which RUNX paralog drives the RUNX motif signal in dcSSc skin fibroblasts

**Written** 2026-09-09, 13:30 EDT / 17:30 UTC, before any access to GSE312129.
**Supersedes** `RR_PREREG_runx2_motif_accessibility_2026-09-09.md`, which is withdrawn.
**Author** Glen Ritschel, Ritschel Research.

---

## 0. Why v1 is withdrawn, stated plainly

v1 preregistered the question "is a RUNX motif more accessible in dcSSc skin
fibroblasts than in healthy control fibroblasts," to be run on GSE195452.

Two things killed it.

**GSE195452 has no chromatin data.** Gate zero fired. The deposited supplementary
files are a cell metadata table and a RAW tar. There is no ATAC, no multiome, no
peaks, no fragments.

**The question v1 asked has already been answered in print.** The dataset that does
carry paired chromatin, GSE312129, comes from a published study whose authors ran
chromVAR themselves and report, verbatim: "Transcription factors downstream of the
Hippo and PI3K pathways, such as RUNX1, CREB, and EGR1 binding, were enriched in
SSc fibroblasts." (JCI Insight, PMID 41411065.)

So a RUNX motif is already known to be enriched. Running v1 would have replicated a
published result and called it a finding.

**This document was written after reading that sentence.** That is disqualifying for
any hypothesis about whether a RUNX motif is enriched, and that hypothesis is
therefore demoted here to a pipeline check with a known expected answer. It is not
disqualifying for the question below, which the paper did not ask.

---

## 1. Prior knowledge on the record, before data

Everything known to the author before writing, so that nothing below can later be
presented as having been discovered:

1. A RUNX motif is enriched in SSc fibroblast chromatin in GSE312129. Published.
2. The authors attributed it to **RUNX1**. The paper does not state what evidence
   distinguished RUNX1 from RUNX2 or RUNX3. RUNX2 and RUNX3 are not mentioned.
3. RUNX1, RUNX2 and RUNX3 share a near-identical core binding motif, TGTGGT.
   Accessibility at a RUNX motif cannot by itself say which paralog occupies it.
4. The cohort is 10 treatment-naive dcSSc and 4 healthy controls, lesional forearm.
5. **No calcinosis.** These are early diffuse cutaneous cases. Calcinosis is a late,
   limited-cutaneous-predominant feature. The lesion of interest is not in this
   dataset and is not in any public single cell dataset.

---

## 2. The open question

Not whether a RUNX motif is enriched. **Which paralog is expressed in the cells where
that motif is accessible.**

This matters to the calcinosis program for one reason. RUNX2 is the master osteogenic
transcription factor. RUNX1 and RUNX3 are not. A RUNX motif signal attributed to
RUNX1 is a fibrosis result. The same signal attributed to RUNX2 would be the first
chromatin-level support for the osteogenic hypothesis in SSc skin. The published paper
picked RUNX1 without stating why, and the motif cannot make that call.

---

## 3. Design, and why it is not a group comparison

The obvious design, dcSSc versus HC donor-level motif deviation, is close to useless
here and this is arithmetic, not opinion.

At n = 10 versus 4, the exact two-sided Mann-Whitney p value under **perfect
separation**, every dcSSc donor ranked above every control, is **0.00200**. That is
the smallest p the test can produce. Any imperfection at all and it produces nothing.
The group test therefore has exactly one detectable outcome, and a null result from it
carries no information whatsoever.

So the group comparison is **not the primary**. It is reported, and it is reported with
that sentence attached.

The primary is a **within-cell** question, which does not depend on having many donors:

> Among dcSSc fibroblast nuclei, does per-nucleus RUNX motif deviation track
> RUNX1 expression, RUNX2 expression, RUNX3 expression, or none of them?

GSE312129 measures RNA and ATAC **in the same nuclei**. That pairing is the entire
reason this dataset and not another. Thousands of fibroblast nuclei, not fourteen
donors, is the unit that answers it.

---

## 4. Fixed before data

**Motifs scored.** JASPAR2024 CORE vertebrates. RUNX1 (MA0002), RUNX2 (MA0511),
RUNX3 (MA0684), RUNX1::CBFB if present.

**Positive controls, taken from the paper's own reported enrichments so they are not
chosen by me after the fact:** SMAD3, EGR1, CREB1, and one ETS family member (ELK4
if present, else the first ETS CORE motif alphabetically). These must move. If they
do not, the pipeline is broken and nothing else in the run is interpretable.

**Negative controls:** HNF4A, POU5F1. Lineage inappropriate. These must not move.

**Cells.** Fibroblast nuclei only, by the authors' own cell type labels as deposited.
If no label column is deposited, the run stops and the labeling problem is solved
first, in writing, before any motif is scored. This is the P34 rule.

**Primary statistic.** Per donor, Spearman correlation across that donor's fibroblast
nuclei between RUNX motif deviation z and each of RUNX1, RUNX2, RUNX3 normalized
expression. Ten donor level correlations per paralog. Combined by Fisher z with equal
weights, since at fixed per-donor cell counts the Fisher z variance is 1/(k-3) and
differs across donors only through k. Two-sided. BH across the three paralogs plus
controls.

**Secondary.** The dcSSc versus HC donor level group test, exact Mann-Whitney,
reported with the 0.00200 floor stated in the same line as the p value.

---

## 5. Outcomes fixed in advance

| Result | Reading |
|---|---|
| Positive controls do not move | Pipeline failure. Report as pipeline failure. Do not interpret RUNX. Stop. |
| RUNX2 not detected in fibroblast nuclei at all | **Attribution impossible.** Report as such. This is the expected outcome and it is not a null result about biology, it is a limit of snRNA sparsity. RUNX2 is lowly expressed outside bone. |
| RUNX1 correlates, RUNX2 does not | Supports the published attribution. The osteogenic reading loses its only available chromatin foothold. Report it. |
| RUNX2 correlates, RUNX1 does not | The published attribution is questionable. Still not evidence about calcinosis, because there is no calcinosis in this cohort. Would justify a targeted follow up, nothing more. |
| Both correlate | Uninformative. Report as uninformative. Do not pick the one that suits the hypothesis. |
| Neither correlates, controls fine | The motif signal is not explained by paralog expression in these cells. Report it. |

---

## 6. The ceiling, restated so it cannot be forgotten later

GSE312129 is **early, treatment-naive, diffuse cutaneous, lesional forearm skin.**

There is no calcinosis in it. There is no calcinosis in any public single cell or
spatial dataset, which is the finding already deposited at Zenodo. Nothing this test
returns is a statement about a calcified lesion. At absolute best it makes one reading
of the osteogenic hypothesis marginally more or less plausible in a tissue that is not
the tissue in question.

It cannot adjudicate Bao against Nasrallah. It is not designed to, and it will not be
described as having done so.

---

## 7. Deviations

Any deviation from the above gets appended here with its date and its reason, before
the result is written up. If that section is empty at the end of the run, no deviation
occurred.

*(none as of writing)*
