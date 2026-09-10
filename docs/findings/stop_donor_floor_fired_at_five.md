# STOP: five dcSSc donors clear the fibroblast floor, the protocol requires six

**Logged** 2026-09-10, after the label audit across all fourteen donors.
**Rule invoked** `RR_LABELING_PROTOCOL_2026-09-09.md` section 4: "Donors surviving the
floor: fewer than 6, the study stops."

---

## The count

Clearing the 200 fibroblast floor: **SSC1 (1,296), SSC2 (274), SSC5 (572), SSC9 (419),
SSC10 (239)**. Five of ten.

Failing: SSC3 (139), SSC4 (5), SSC6 (12), SSC7 (89), SSC8 (170).

All four healthy controls clear it comfortably: HC1 3,435, HC2 2,061, HC3 1,423,
HC4 1,027.

The floor was fixed on 2026-09-09 before any donor beyond SSC1 was loaded. It fires at
five. **The preregistered dcSSc analysis does not run.**

---

## The labeling is not broken. The tissue recovery is.

That was the live hypothesis and the audit rules it out.

The donors that fail are **keratinocyte dominated**, not mislabeled. SSC3 is 1,958
keratinocyte of 2,269 nuclei, 86 percent. SSC4 is 2,359 of 2,760, 85 percent. SSC6 is
1,551 of 1,880, 82 percent. In those same samples the fibroblast panel is genuinely
absent, not merely outvoted: COL1A1 is detected in 7.9 percent of SSC4 nuclei and
7.4 percent of SSC6, against 64.7 percent in SSC1.

Both halves of the diagnostic agree. Few fibroblast labels **and** no fibroblast
transcript. The argmax rule reported what was there.

A 4 mm punch biopsy contains epidermis and dermis. The ratio of epidermal to dermal
nuclei recovered varies enormously between preparations, and that is what this table is
measuring.

---

## The finding that outlives this study

**Fibroblast recovery differs systematically between the groups, in the direction
opposite to the disease.**

| | fibroblast share of QC nuclei | COL1A1 detection |
|---|---|---|
| dcSSc, median of 10 | **6.1%** | **31.8%** |
| HC, median of 4 | **38.3%** | **78.6%** |
| exact Mann-Whitney | p = 0.014 | p = 0.054 |

SSc skin is the fibrotic tissue. It yields **six times fewer** fibroblast nuclei than
healthy skin in this dataset.

The straightforward explanation is dissociation, not biology: fibrotic dermis is dense
and collagen cross linked, so it releases nuclei poorly, while the epidermis above it
releases them normally. The recovered fraction is therefore not a random sample of the
tissue, and it is differentially non random between the arms.

**This bears on the published analysis, not only on ours.** Any dcSSc versus HC
comparison of fibroblasts in GSE312129 is comparing a small, dissociation-selected
minority of SSc fibroblasts against a large and much less selected majority of HC
fibroblasts. That is a composition confound sitting underneath every between-group
statement drawn from this dataset, including the chromVAR enrichments in the paper that
started this whole line of work.

It also retroactively explains the secondary that was already demoted for power. It was
never going to be interpretable.

---

## What the five surviving donors showed, recorded and not interpreted

RUNX2 motif against RUNX2 mRNA, the null-matched statistic, for the three scored before
the floor was known to fail:

| donor | rho | null percentile |
|---|---|---|
| SSC1 | +0.0797 | 95th |
| SSC2 | -0.0307 | 30th |
| SSC5 | +0.1098 | 99th |

Two strongly positive, one negative. **These are not combined and are not a result.**
They are recorded so that nobody, including me, can later present a subset of them as
one.

---

## The three ways forward, and what each costs

**A. Stop.** The preregistered question does not have the donors. Write it up as a
negative methodological result: the paralog attribution cannot be made in the only
public paired multiome of SSc skin, because SSc fibroblast recovery is too poor in over
half the donors. The composition confound above is the publishable finding.

**B. Lower the floor.** Dropping to 150 admits SSC8 at 170 and reaches six. This is a
post hoc threshold change made **after** seeing that five is one short, which is exactly
the move I flagged in the Deviation 2 amendment as the one that hollows out a
preregistration. I do not recommend it and would not defend it.

**C. Write a new preregistration for a different question.** The paralog question is
molecular, not disease specific: does RUNX2 mRNA track RUNX2 motif accessibility in
**skin fibroblasts**. That admits all nine donors clearing the floor, five dcSSc and
four HC, with disease as a covariate or a stratum. It is better powered, it is a fair
question, and the HC donors become an independent replication set rather than a
comparison arm, which sidesteps the composition confound entirely.

It is a **different question** from the one preregistered on 2026-09-09, it needs its
own document written before the combine is run, and the writeup must say plainly that
the dcSSc-only version stopped on its floor first.

---

## Unchanged

Nothing here is a statement about calcinosis. GSE312129 is early, treatment naive,
diffuse cutaneous, lesional forearm. There is no calcinosis in it, and none in any
public single cell or spatial dataset, which remains the deposited Zenodo finding.
