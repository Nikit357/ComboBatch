# Best approaches

What the source benchmark recommends, what qualifies each recommendation, and how to
express its data subsets with `--subset-query`.

Runnable versions of everything here live in `examples/configs/best_approaches/`.

> **Read the caveats before the table.** Three separate rankings of the same benchmark
> produce top-ten sets with **zero overlap**, and at least one recommended approach
> destroys biological signal while scoring well on batch mixing. A recommendation here
> means "this worked well on one dataset under one ranking", never "this is optimal".

## The decision tree

| Your data | Method | Imputation | Subset |
|---|---|---|---|
| Fresh-frozen only | `10_mnn` (alt `16_fsqn_r`) | `strict` | `RNA_BATCH.str.contains('_FF_')` |
| FFPE only ★ | `04_sva` | `softimpute` (alt `knn`) | `RNA_BATCH.str.contains('_FFPE_')` |
| RNA-seq only | `04_sva` (2nd `16_fsqn_r`) | `softimpute` (alt `knn`) | `RNA_BATCH.str.startswith('RNASeq')` |
| RNA-seq + Illumina arrays | `16_fsqn_r` + post-removal (alt `33_amdbnorm`, `13_fsmvn`) | `strict` | none |
| RNA-seq + mixed arrays | `10_mnn` + post-removal | `strict` | none |

★ best local batch mixing anywhere in the benchmark — which matters because local mixing is
where nearly everything else fails.

Column names in those queries are the source dataset's. Substitute your own: nothing in
ComboBatch knows what `RNA_BATCH` is.

## Eight caveats

**1. Global metrics are necessary but not sufficient.** Good PCReg/R² coexists with locally
separated batches. The global and local metric families are near-orthogonal — median
Spearman r ≈ 0.009 across the benchmark. Read at least one of each.

**2. Local metrics universally fail cross-platform.** kBET and iLISI rate roughly 60% of all
tested approaches as "bad". This is a property of the problem, not a defect in any method,
and nothing in the benchmark solved it. Do not go looking for the method that will.

**3. Three rankings disagree and the article has not settled which it reports.** The
equal-weight composite (the figure default, best group at 0.602), the literature-weighted
composite (top ten entirely `C_rnaseq_only__post1`, led by a *Shambhala* run at 0.732), and
a Borda count (top ten entirely `16_fsqn_r` on strategies A and E1–E3). The three top-ten
sets have **zero overlap** with each other or with the curated Best-15. Treat the Best-15 as
the manuscript's curated selection, not as the optimum.

**4. At least one Best-15 approach corrupts biology.** Marker-gene correlation before and
after harmonization was ~1 for the linear methods but fell as low as **0.05** for one of the
clustermap-best approaches. This is the caveat that most changes practice: a high
batch-mixing score is not evidence that the result is usable. Compute metric groups **L and
M** on your own outputs — every config in `examples/configs/best_approaches/` enables them
for exactly this reason.

**5. No method dominates.** 26 of 31 methods are Pareto-optimal somewhere in the
batch-mixing × biology-preservation plane. That is the argument for running a cross-product
rather than adopting a recommendation.

**6. `softimpute` beats `knn` by roughly 10 percentage points** on PCR. Where both are
listed, prefer softimpute.

**7. The factor-importance hierarchy is version-dependent.** The manuscript's η² analysis:
**method 36.3% > subset 25.6% > imputation 1.6% > post-removal 0.2%**. An earlier internal
note ordered them differently (subset > method > post-removal > imputation); quote the
manuscript numbers. Either way the practical conclusion is the same — which method you pick
matters far more than how you tune it.

**8. The Best-15 is a manual curation**, combining visual inspection with composite scoring.
It is not the output of an automated rule.

## Reproducing the benchmark's data subsets

The 14 filter strategies of the source pipeline are not shipped as code. Most are one
`--subset-query` expression. Counts are from the source dataset (7,238 samples).

| Strategy | `--subset-query` | n |
|---|---|---|
| `C_rnaseq_only` | `RNA_BATCH.str.startswith('RNASeq')` | 2,243 |
| `J_ff_only` | `RNA_BATCH.str.contains('_FF_')` | 3,167 |
| `K_ffpe_only` | `RNA_BATCH.str.contains('_FFPE_')` | 3,000 |
| `F_microarray_only` | `not RNA_BATCH.str.startswith('RNASeq')` | 4,931 |
| `G_affymetrix_only` | `RNA_BATCH.str.startswith('GPL570')` | 2,801 |
| `D_malignant_only` | `Diagnosis_cell_type_unified not in @NORMAL_GROUPS` | 6,285 |
| `S0_no_removal` | `Diagnosis_cell_type_unified not in @RARE_GROUPS` | 7,174 |

**`J_ff_only` and `K_ffpe_only` are not complementary:** 3,167 + 3,000 ≠ 7,174. Batches whose
middle token is `Unknown` — notably `GPL570_Unknown_Unknown`, 1,030 samples — fall into
neither. A query that looks like a partition usually is not one; check the counts.

Three strategies cannot be expressed as a static query, and are honest about it:

- **`H_affymetrix_extended`** is a hand-curated 16-label allow-list. It deliberately excludes
  `GPL14951_FFPE_Unknown`, which is Illumina rather than Affymetrix — a pattern match on the
  platform prefix would wrongly include it. Write it as an explicit `RNA_BATCH in [...]`.
- **`E1`/`E2`/`E3`** are iterative, data-dependent PCA-outlier removal. Approximate them with
  `--post-removal` applied repeatedly.
- **`I_rare_batches_removed`** is a data-dependent size threshold, not a fixed list.

## The practical advice

Run the cross-product on your own data and read the metrics:

```bash
combobatch dispatch --config examples/configs/full_crossproduct.yaml --n-workers 8
combobatch concat --out-dir ./tables
```

Then judge with at least one global metric (group A), one local metric (group B), and one
biology-preservation metric (group M) — because caveats 1, 2 and 4 each describe a way that
reading fewer than three would mislead you.
