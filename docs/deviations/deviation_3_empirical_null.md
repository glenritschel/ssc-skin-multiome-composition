# DEVIATION 3: the positive controls test a different quantity than the primary, and the negative controls test nothing

**Logged** 2026-09-09, from the first valid SSC1 run.
**Against** `RR_PREREG_runx_paralog_v2_2026-09-09.md` section 4.
**Status** One donor. No combined result exists. The six donor floor has not been met.

---

## The SSC1 table

| target | partial rho | uncontrolled | detection | nominal p |
|---|---|---|---|---|
| RUNX1 | +0.0697 | +0.0365 | 51.1% | 0.012 |
| RUNX2 | +0.0866 | +0.0739 | 14.6% | 0.002 |
| RUNX3 | nan | | 3.9% | below floor |
| SMAD3 **(pos)** | **-0.0160** | | 33.0% | 0.565 |
| EGR1 **(pos)** | **+0.0007** | | 62.3% | 0.980 |
| CREB1 **(pos)** | +0.0583 | | 24.2% | 0.036 |
| ELK4 **(pos)** | **+0.0185** | | 21.0% | 0.506 |
| HNF4A (neg) | -0.0074 | | **0.2%** | 0.790 |
| POU5F1 (neg) | +0.0283 | | **0.6%** | 0.309 |

Three of four positive controls are flat or negative. The RUNX values sit above every
one of them. Read naively that is a strong result. It is not, and the reason is mine.

---

## Error 1: the positive controls measure the wrong thing

Preregistration section 4 says the positive controls are "taken from the paper's own
reported enrichments so they are not chosen by me after the fact." I treated provenance
as the only thing that mattered and never asked whether the quantity matched.

The paper reported motifs **differentially enriched in SSc fibroblasts versus healthy
controls.** That is a between-group contrast.

The primary here is a **within-cell correlation** between a motif's deviation and its own
transcription factor's mRNA, across nuclei of a single donor.

Those are different quantities. A motif can be strongly differentially accessible between
groups while its per-nucleus deviation has no relationship to its own transcript.
Transcription factor motif accessibility and transcription factor mRNA are known to
correspond weakly at single-cell resolution, because protein abundance, nuclear
localization and cofactor availability all sit between them.

So SMAD3 at -0.016 and EGR1 at +0.0007 are not necessarily evidence that the pipeline is
broken. They may be evidence that I picked controls for a quantity the analysis does not
compute. **I cannot currently tell those two possibilities apart**, and that is the whole
problem: the positive control exists precisely to distinguish them.

Worse, I did not notice this by reading my own preregistration. I noticed it because the
controls failed. That ordering is bad and is recorded as such.

## Error 2: the negative controls are not controls

HNF4A is detected in **0.2 percent** of fibroblast nuclei. POU5F1 in **0.6 percent**.
Both sit far below the 5 percent floor that excludes RUNX3 from the analysis entirely.

Correlating a motif deviation against a gene that is essentially never observed produces
noise around zero by construction, whatever the pipeline is doing. It would return the
same answer if the deviation calculation were replaced with random numbers.

"Lineage inappropriate" was the criterion in the preregistration, and lineage
inappropriate genes are not expressed, which is exactly what makes them useless as a
control for this particular statistic. The criterion defeated its own purpose.

---

## The fix: an empirical null, not a hand-picked control set

Swapping in positive controls until some of them work would be the worst available
response. Instead the design changes so that the control is the data itself.

**Change.** Score, per donor, every human JASPAR CORE motif whose own transcription
factor gene is detected in at least 5 percent of that donor's fibroblast nuclei, capped
at a seeded random panel of 80 to keep the motif scan tractable. Compute the same partial
Spearman for every one of them.

That panel is the **empirical null distribution** of within-cell motif-to-own-expression
correlation in this tissue, computed by the same code on the same nuclei.

**The primary question becomes two things, both fixed now:**

1. Where does the RUNX family deviation's correlation with RUNX2 expression fall in that
   distribution? Reported as a percentile, combined across donors.
2. Is RUNX2's correlation higher than RUNX1's, against the spread of the null?

**What this buys.** If the whole panel is flat, the pipeline cannot detect this class of
relationship and nothing about RUNX is interpretable, which is the honest reading and the
one my broken controls could not deliver. If the panel has real structure and RUNX2 sits
in its upper tail while RUNX1 does not, that is a result with a null behind it rather
than a p value against zero.

**What it costs.** Eight times the motif scanning per donor. The deviation routine is
being rewritten to be sparse and vectorized, since the current background sampling loops
in Python and would not finish at panel scale.

---

## Recorded, per the Deviation 2 amendment

`rho(RUNX deviation, ATAC depth) = -0.033`, **below the 0.05 threshold the amendment
fixed in advance.** The amendment therefore requires this sentence in any writeup: the
depth covariate was unnecessary, and the uncontrolled result is materially the same.

Note the partial rho is **larger** than the uncontrolled for both paralogs, +0.0697
against +0.0365 and +0.0866 against +0.0739. Depth was suppressing the association
slightly, not inflating it, consistent with leg two being negative. My concern in
Deviation 2 was that depth would manufacture a correlation. It did the opposite.

---

## Not changed

Everything in the labeling protocol. The 5 percent detectability floor, which correctly
excluded RUNX3 at 3.9 percent. The six donor floor, which correctly refused to combine at
one. The ceiling in preregistration section 6: no calcinosis in this cohort, and nothing
here is a statement about a calcified lesion.

The two RUNX numbers in the table above are **one donor of ten** and are not a result.
