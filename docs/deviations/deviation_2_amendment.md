# AMENDMENT to Deviation 2, same day, before any motif is scored

**Logged** 2026-09-09, 16:05 EDT / 20:05 UTC

## The withdrawal criterion fired

Deviation 2 change 2 was written with this test attached, verbatim:

> "If RUNX2's rho is not materially larger than COL1A1's, this change is unnecessary and
> will be withdrawn rather than kept for tidiness."

Measured, on 1,296 SSC1 fibroblast nuclei at MAX_MT 20:

| gene | pct nonzero | rho(detected ~ RNA depth) |
|---|---|---|
| RUNX1 | 51.08 | +0.493 |
| **RUNX2** | **14.58** | **+0.174** |
| RUNX3 | 3.86 | -0.001 |
| **COL1A1 (reference)** | **97.53** | **+0.158** |
| SMAD3 | 33.02 | +0.350 |
| EGR1 | 62.35 | +0.362 |

RUNX2 sits at 1.10 times the reference. That is not materially larger.

**The stated rationale for the covariate is withdrawn.** RUNX2 zeros in this sample are
not predominantly dropout. My reasoning across two sweeps, that RUNX2 detection tracks
depth and its zeros are therefore dropout, was wrong. The between-cut movement of the
RUNX2 percentage came from the composition of the nuclei admitted, not from depth
sensitivity in RUNX2 itself.

## The covariate is nonetheless retained, on a different and stated ground

The diagnostic found a confound the criterion was not written to catch.

Depth sensitivity is **differential across the three paralogs being compared**. RUNX1
is at +0.493, three times the reference. RUNX2 is at +0.174, at the reference. The gap
is +0.319.

The primary is not a test of one gene. It is a **comparison** among RUNX1, RUNX2 and
RUNX3 to decide which one accompanies the motif signal. When one arm of a comparison is
three times more depth sensitive than another, and RNA and ATAC depth are coupled at
rho +0.401 across these nuclei, the comparison is biased toward RUNX1 before any
biology enters.

That matters here specifically, because RUNX1 is the paralog the published paper
attributed the signal to. An uncontrolled analysis would tend to confirm the published
attribution for a reason that has nothing to do with transcription factors.

## Honesty about what this move is

Writing a withdrawal test, watching it fire, and keeping the change anyway is the exact
maneuver that hollows out a preregistration. It is defensible here only because all of
the following hold, and if any of them failed the covariate would be dropped:

1. No motif has been scored. No RUNX association exists yet, in any form.
2. The new ground is stated explicitly and is different from the withdrawn one. It is
   not the old rationale rephrased.
3. The change makes the analysis **more** conservative toward the hypothesis of
   interest. It penalizes RUNX1, the incumbent answer, not RUNX2.
4. It is logged before the fact rather than reconstructed after.

## Magnitude, and the leg that is still unmeasured

The artifact needs two legs. Leg one is expression against RNA depth, now measured and
differential. Leg two is motif deviation against ATAC depth, and it is **unmeasured**,
because no motif has been scored.

chromVAR-style deviation z is background matched on GC and accessibility precisely to
kill leg two. If it does its job, leg two is near zero and the whole confound is moot.

Induced spurious correlation at plausible values of leg two:

| rho(deviation, ATAC depth) | spurious RUNX1 | spurious RUNX2 | gap |
|---|---|---|---|
| 0.1 | 0.020 | 0.007 | 0.013 |
| 0.2 | 0.040 | 0.014 | 0.026 |
| 0.3 | 0.059 | 0.021 | 0.038 |
| 0.5 | 0.099 | 0.035 | 0.064 |

Small in absolute terms. Not small relative to the effect sizes a 1,296 nucleus
correlation will resolve.

**Required at motif time:** report rho(deviation, ATAC depth) alongside the primary. If
it is below 0.05, state in the writeup that the depth covariate was unnecessary and
that the uncontrolled result is materially identical. If it is above 0.2, the covariate
was load bearing and that gets stated too.

## Also recorded

RUNX3 at 3.86 percent nonzero and rho -0.001 is below the 5 percent detectability floor
in SSC1 and is excluded for this donor under preregistration section 5, unchanged.

COL1A1 shows rho 0.631 for expression against depth but only 0.158 for detection
against depth. Depth moves its level, not whether it is seen, which is what a gene
detected in 97.53 percent of nuclei should do. The diagnostic behaves correctly.
