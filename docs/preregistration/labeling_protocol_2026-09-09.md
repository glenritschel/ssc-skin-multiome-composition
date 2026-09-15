# LABELING PROTOCOL for GSE312129, fixed before any nucleus is scored

**Addendum to** `RR_PREREG_runx_paralog_v2_2026-09-09.md`
**Written** 2026-09-09, 14:05 EDT / 18:05 UTC
**Trigger** The preregistration's own stop condition fired.

---

## 0. Why this document exists

Section 4 of the preregistration says:

> "Fibroblast nuclei only, by the authors' own cell type labels as deposited. If no
> label column is deposited, the run stops and the labeling problem is solved first,
> in writing, before any motif is scored. This is the P34 rule."

**No labels are deposited.** GSE312129 contains raw and filtered matrices, ATAC
fragments, and spatial images. There is no annotation table. Compare GSE195452, which
does deposit `GSE195452_Cell_metadata_v26_anno.txt.gz`. GSE312129 deposits nothing
equivalent.

So the run stops here, and this is the writing.

This is logged as deviation 1 in preregistration section 7: the labeling method was
not specified in advance, because the preregistration assumed labels would be
deposited. It is specified now, before any data is loaded.

---

## 1. The failure this is designed to prevent

P34 shipped `TOP_CLUSTERS = ['2', '3']` into an executed run. Clusters 2 and 3 turned
out to be simply the two largest, 5,990 and 5,132 nuclei, 11,122 exactly. Nobody chose
them for a biological reason. The published labels were loaded and never used.

The generalized error is **picking cell populations by looking at them.** Any procedure
where a human inspects a UMAP and names the fibroblast cluster reproduces it.

So no clusters are picked here. Assignment is per nucleus and rule based.

---

## 2. Marker panels, fixed

Every panel below is final. Nothing is added or removed after seeing data.

| Population | Genes |
|---|---|
| **fibroblast** | COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRA, FBLN1 |
| keratinocyte | KRT14, KRT5, KRT1, KRT10, KRT6A |
| endothelial | PECAM1, VWF, CDH5, EGFL7 |
| pericyte / smooth muscle | RGS5, MYH11, NOTCH3, KCNJ8 |
| myeloid | LYZ, CD68, ITGAM, CSF1R, AIF1 |
| T and NK | CD3D, CD3E, IL7R, CD2, TRAC |
| B and plasma | MS4A1, CD79A, JCHAIN, MZB1 |
| mast | TPSAB1, TPSB2, CPA3, MS4A2 |
| melanocyte | PMEL, MLANA, TYRP1, DCT |
| adnexal and sweat gland | KRT7, AQP5, SCGB2A2, DCD |

**ACTA2 and TAGLN are deliberately excluded from the pericyte panel.** Myofibroblasts
express both. Including them would push exactly the nuclei this study cares about into
the pericyte bin. RGS5, MYH11, NOTCH3 and KCNJ8 do not have that collision.

---

## 3. Assignment rule

Per donor, independently. No cross-donor integration, for reasons in section 5.

1. QC filter, fixed: nuclei with fewer than 500 detected genes are dropped. Nuclei with
   more than 5 percent mitochondrial counts are dropped, since this is single **nuclei**
   data and high mitochondrial content indicates ambient contamination or damage.
2. Normalize to 10,000 counts per nucleus, log1p.
3. `scanpy.tl.score_genes` for each of the ten panels, with scanpy's matched control
   gene sets.
4. Assign each nucleus to the panel with the highest score.
5. **Margin requirement.** If the top score minus the second score is below 0.05, the
   nucleus is labeled `ambiguous` and excluded. Ambiguous count is reported per donor.

No cluster is inspected. No UMAP is drawn before assignment. A UMAP may be drawn
afterward, as a diagnostic only, and it cannot change any label.

---

## 4. Prespecified failure conditions

Each of these stops the run and gets reported as a labeling failure, not as a result:

| Condition | Threshold |
|---|---|
| Fibroblast nuclei recovered for a donor | fewer than 200, that donor is excluded and reported |
| Donors surviving the floor | fewer than 6, the study stops |
| Detection sanity | fewer than 30 percent of called fibroblasts carry a nonzero COL1A1 or PDGFRA count, labeling failed, stop |
| Ambiguous fraction | above 40 percent of nuclei in a donor, labeling failed for that donor |
| Keratinocyte share | if fibroblasts outnumber keratinocytes by more than tenfold in whole skin, something is wrong, stop and investigate |

---

## 5. Paralog detectability, decided before looking

This is the condition most likely to end the study, so it is written down now rather
than discovered later.

For a given donor, a paralog is **detectable** if at least 5 percent of that donor's
fibroblast nuclei carry a nonzero count for it.

- Paralogs below the floor get no correlation computed for that donor.
- If RUNX2 is below the floor in more than half the surviving donors, the preregistered
  outcome is **"attribution impossible"**, which section 5 of the preregistration already
  names as the expected result. It is a limit of single nuclei sparsity, not a finding
  about biology, and it will not be written up as one.

RUNX2 is lowly expressed outside bone. Expect this.

---

## 6. A technical consequence that removes the fragment files entirely

Each donor's `filtered_feature_bc_matrix.h5` is Cell Ranger ARC output and should carry
both feature types, `Gene Expression` and `Peaks`, on the same barcodes. If it does,
**the ATAC fragment files are not needed.**

That matters. The fragments are the large objects here. The h5 files are not.

The reason this works is the design chosen in the preregistration for statistical
reasons, and it turns out to have a second benefit. The primary is a **within donor**
correlation, computed on that donor's own nuclei, then combined across donors by
Fisher z. Each donor can therefore be processed on its own per-sample peak set. No
common peak set across donors is required, so no re-counting from fragments is required.

**The cost, stated honestly.** Motif deviations computed on different peak sets are not
comparable across donors. The preregistration's secondary, the dcSSc versus HC group
comparison, therefore may not be executable at all without building a union peak set
and re-quantifying from fragments. Given that the same test already has a floor of
p = 0.00200 at 10 versus 4 and only one detectable outcome, that is a small loss. It is
recorded as a loss rather than quietly dropped.

**This is an assumption until the h5 is opened.** The first cell of the pipeline loads
one donor and raises if `Peaks` is absent from `feature_types`. It is not assumed true.

---

## 7. Per-sample download map

Per-GSM supplementary files, avoiding the full RAW tar and its spatial images:

```
https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM9338nnn/<GSM>/suppl/<GSM>_<SAMPLE>_filtered_feature_bc_matrix.h5
```

| Sample | GSM (multiome h5) | Sample | GSM (multiome h5) |
|---|---|---|---|
| SSC1 | GSM9338143 | SSC8 | GSM9338157 |
| SSC2 | GSM9338145 | SSC9 | GSM9338159 |
| SSC3 | GSM9338147 | SSC10 | GSM9338161 |
| SSC4 | GSM9338149 | HC1 | GSM9338163 |
| SSC5 | GSM9338151 | HC2 | GSM9338165 |
| SSC6 | GSM9338153 | HC3 | GSM9338167 |
| SSC7 | GSM9338155 | HC4 | GSM9338169 |

The even numbered GSMs in between are the `atac_fragments.tsv.gz` registrations and are
not downloaded.

---

## 8. Order of operations

1. Load **one** donor, SSC1. Confirm `Peaks` is present. This is cheap and decisive.
2. Run the labeling rule on that one donor. Report every diagnostic in section 4.
3. Only if section 4 passes on donor one, extend to all fourteen.
4. Only then score motifs.

Nothing is downloaded in bulk before step 1 answers.

---

# 9. Deviations and prespecified-check reporting

**Checks computed 2026-09-13. Appended to this document 2026-09-15.** After deposit of the
composition note (concept DOI 10.5281/zenodo.22693563, current version
10.5281/zenodo.22714302, v3). Appended so that sections 0 through 8 remain the document as
fixed on 2026-09-09.

The checks below were recomputed by re-running the labeling rule of sections 2 and 3
unchanged. That re-run reproduces every fibroblast count in the deposited note as an
identical integer and every donor's fibroblast share to within 0.04 percentage points,
so what follows is the record of what the deposited run computed, not a later
measurement reported in its place.

## 9.1 Deviation: the mitochondrial threshold moved from 5 percent to 20 percent

Section 3.1 fixes the cut at 5 percent, on the reasoning that this is single nuclei data
and high mitochondrial content indicates ambient contamination or damage.

Applied to the first donor, the 5 percent cut discarded 87.3 percent of barcodes. The
threshold was changed to 20 percent and all fourteen samples were processed under the
single changed threshold. No sample was processed under both.

The reasoning originally offered for the loss, that snRNA is shallow by construction, was
wrong and is retracted: the median here is 1,851 genes per nucleus. The loss was caused by
the mitochondrial cut itself.

The composition note discloses the change in its Method section. This entry exists so that
the protocol document, which is the document claiming things were fixed in advance, carries
the record of what moved.

## 9.2 Prespecified checks: what was run and what it returned

Section 3.5 requires the ambiguous count per donor. Section 4 prespecifies five stop
conditions. The composition note reports none of them by name. They are reported here.

| prespecified check | where fixed | result |
|---|---|---|
| Ambiguous count per donor | sec 3.5 | REPORTED, per-donor table below. 4 to 346 nuclei, 1.0 to 5.1 percent of QC nuclei. |
| Ambiguous fraction above 40 percent fails that donor | sec 4 | PASSES, 0 of 14 donors fail. Highest observed 5.1 percent. |
| Detection sanity: at least 30 percent of **called fibroblasts** carry nonzero COL1A1 or PDGFRA | sec 4 | PASSES, 0 of 14 donors fail. Range 83.3 to 100 percent. Two cells are low-n and are marked below. |
| Keratinocyte share: fibroblasts outnumbering keratinocytes more than tenfold | sec 4 | PASSES, 0 of 14 donors fail — vacuously. See 9.4. |
| Fibroblast nuclei per donor, floor of 200 | sec 4 | **FAILS on five of ten dcSSc donors.** SSC3, SSC4, SSC6, SSC7, SSC8. |
| Donors surviving the floor, minimum of six | sec 4 | **FAILS.** Five of ten dcSSc donors survive. This is the failure that stopped the study. |
| Paralog detectability floor, 5 percent of a donor's fibroblast nuclei | sec 5 | NOT RUN; the study stopped at the donor floor first. |

| donor | QC nuclei | ambiguous | ambiguous % | called fibroblasts | COL1A1 or PDGFRA in called fibroblasts | keratinocyte nuclei | fibroblast : keratinocyte |
|---|---|---|---|---|---|---|---|
| SSC1 | 5,761 | 251 | 4.4% | 1,296 | 98.7% | 2,139 | 0.61 |
| SSC2 | 3,444 | 90 | 2.6% | 274 | 100.0% | 1,035 | 0.27 |
| SSC3 | 2,269 | 31 | 1.4% | 139 | 97.1% | 1,958 | 0.07 |
| SSC4 | 2,760 | 33 | 1.2% | 5 | **100.0% (n = 5, not interpretable)** | 2,359 | 0.00 |
| SSC5 | 3,820 | 106 | 2.8% | 572 | 95.5% | 2,000 | 0.29 |
| SSC6 | 1,880 | 19 | 1.0% | 12 | **83.3% (n = 12, not interpretable)** | 1,551 | 0.01 |
| SSC7 | 267 | 4 | 1.5% | 89 | 98.9% | 48 | 1.85 |
| SSC8 | 2,801 | 45 | 1.6% | 170 | 95.3% | 1,830 | 0.09 |
| SSC9 | 6,835 | 346 | 5.1% | 419 | 93.6% | 3,775 | 0.11 |
| SSC10 | 5,000 | 204 | 4.1% | 239 | 99.6% | 3,752 | 0.06 |
| HC1 | 4,662 | 135 | 2.9% | 3,435 | 99.5% | 870 | 3.95 |
| HC2 | 3,993 | 190 | 4.8% | 2,061 | 97.2% | 1,171 | 1.76 |
| HC3 | 5,684 | 249 | 4.4% | 1,423 | 95.8% | 3,207 | 0.44 |
| HC4 | 5,923 | 268 | 4.5% | 1,027 | 96.1% | 2,472 | 0.41 |

**A detection percentage computed over 5 or 12 called fibroblasts is arithmetic, not a
check.** Those two cells carry their n and are not offered as evidence of anything.

## 9.3 A referent correction, and it changes what a reader concludes

The composition note's COL1A1 column is detection across **all QC nuclei**. The detection
sanity check at section 4 is detection across **called fibroblasts**. These are different
quantities and in this dataset they differ by more than twelvefold in the same donor:
SSC4 is 7.9 percent by the first and 100 percent by the second.

The consequence is not only that a prespecified check went unreported. A reader ticking
section 4 off against the note's column would read 7.9 percent against a 30 percent
threshold and conclude the check FAILED, when the check as written returns 100 percent.
The note's column should be read, and in any future version labeled, as "COL1A1 detected,
all QC nuclei."

## 9.4 The keratinocyte-share check passed and could not have failed

Section 4's keratinocyte condition is one-sided: it fires only if fibroblasts outnumber
keratinocytes by more than tenfold. Observed ratios run 0.00 to 3.95, the maximum being a
healthy control. No donor is within a factor of two and a half of the threshold, and the
threshold sits on the side of 1.0 that this tissue never approaches.

The check therefore passed on SSC4, a library that is 85 percent keratinocyte and yielded
five fibroblast nuclei. It would have passed on a donor with no fibroblasts at all.

What did catch that condition was the 200-nucleus per-donor floor and the six-donor study
floor, both of which fired as written. Nothing was missed. But this check did no work, and
recording it as "passed" without this paragraph would imply that it did.

A check is stated with the direction its failure would come from. A one-sided check
records the direction it cannot see. This one could not see the only direction the data
went.

