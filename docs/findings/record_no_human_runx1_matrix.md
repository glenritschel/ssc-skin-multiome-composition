# RECORD: JASPAR2024 CORE has no human RUNX1 matrix, and the paralog motifs are not near identical

**Logged** 2026-09-09, before any valid rho exists. The SSC1 numbers are void per the
bug report; nothing has been rerun.

## What the source actually contains

| name | matrix ID | species | length |
|---|---|---|---|
| Runx1 | MA0002.3 | **10090, Mus musculus** | 9 |
| RUNX2 | MA0511.2 | 9606 | 9 |
| RUNX3 | MA0684.3 | 9606 | 8 |
| SMAD3 | MA0795.1 | 9606 | 10 |
| EGR1 | MA0162.5 | 9606 | 10 |
| CREB1 | MA0018.5 | 9606 | 8 |
| ELK4 | MA0076.3 | 9606 | 9 |
| HNF4A | MA0114.5 | 9606 | 9 |
| HNF4A | MA1494.2 | 9606 | 14 |
| POU5F1 | MA1115.2 | 9606 | 7 |

**There is no human RUNX1 position weight matrix in the preregistered source.** The only
RUNX1 entry is mouse. That is a property of JASPAR2024 CORE vertebrates, not a coding
error, and the preregistration named that source, so it stands.

## Measured, not asserted

Best aligned column correlation between the position frequency matrices, over all
ungapped offsets and both strands:

| pair | similarity |
|---|---|
| Runx1 (mouse) vs RUNX2 (human) | **0.816** |
| Runx1 (mouse) vs RUNX3 (human) | 0.777 |
| RUNX2 (human) vs RUNX3 (human) | **0.822** |

Two things follow, and both matter.

**Dropping the mouse matrix costs nothing.** Mouse Runx1 sits 0.816 from human RUNX2.
Human RUNX3 sits 0.822 from human RUNX2. The mouse matrix is no more distant from the
human ones than the human ones are from each other. Excluding it removes almost no
information, so `require_human=True` stays as the cleaner default, now on a measurement
rather than a preference.

**"Near identical" was an overstatement of mine.** I have written that phrase repeatedly
in this thread about the RUNX core. At 0.78 to 0.82 the matrices are highly similar,
not near identical. The correct statement is that they are far too similar for
accessibility alone to attribute a paralog, which is what the preregistration's
reasoning actually needs and which these numbers support.

## Consequence for the family score

The RUNX family deviation is the average of **RUNX2 (MA0511.2) and RUNX3 (MA0684.3)**,
both human. That is what the writeup must say. It is not "RUNX1, RUNX2 and RUNX3
motifs," and any figure legend claiming three matrices would be wrong.

This does not weaken the design, because paralog discrimination was never coming from
the motif side. It comes from the paired snRNA in the same nuclei. The motif side only
needs to locate RUNX-family accessible chromatin, and two human matrices at 0.822 to
each other do that.

The pipeline reports `runx_motif_intercorr`, the minimum pairwise correlation of the
deviation vectors actually used. If that comes back low, the two matrices are not
measuring the same thing in this data and the averaging assumption fails. It is the
check on the paragraph above.

## HNF4A

Two human matrices, MA0114.5 at length 9 and MA1494.2 at length 14. The patched pipeline
averages them and logs the count. Previously one was silently discarded.
