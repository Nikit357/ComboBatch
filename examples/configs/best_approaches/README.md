# Best approaches

Five runnable configurations, one per branch of the benchmark's decision tree. Each
carries the subset expression, method and imputer that the dissertation benchmark
recommends for that data shape.

| Config | Scenario | Method | Imputation | On the example dataset |
|---|---|---|---|---|
| `ff_only_mnn.yaml` | Fresh-frozen only | `10_mnn` (alt `16_fsqn_r`) | `strict` | 104 samples, 3 batches |
| `ffpe_only_sva.yaml` | FFPE only ★ | `04_sva` | `softimpute` (alt `knn`) | 40 samples, **1 batch** |
| `rnaseq_only_sva.yaml` | RNA-seq only | `04_sva` (2nd `16_fsqn_r`) | `softimpute` (alt `knn`) | 34 samples, **1 batch** |
| `rnaseq_plus_illumina_fsqn.yaml` | RNA-seq + Illumina arrays | `16_fsqn_r` + post-removal | `strict` | 144 samples, 4 batches, **no Illumina** |
| `multiplatform_mnn.yaml` | RNA-seq + mixed arrays | `10_mnn` + post-removal | `strict` | 144 samples, 4 batches |

★ best local batch mixing anywhere in the benchmark.

**Three of the five are degenerate on the example dataset**, as the last column says: two
subsets leave a single batch, so there is no between-batch structure to correct, and the
example data contains no Illumina arrays. They demonstrate the invocation. Point them at
your own data for the scenario each is named for.

```bash
combobatch run --config examples/configs/best_approaches/multiplatform_mnn.yaml
```

Every config enables metric groups **L and M** — the marker-gene correlation and
cross-batch discrimination groups. That is deliberate: they are what would catch a
configuration destroying the biology it is supposed to preserve. See caveat 4.

---

## Read this before calling any of them "the best"

These are recommendations from one benchmark on one dataset. Eight findings qualify them,
and leaving them out would make this directory misleading.

**1. Global metrics are necessary but not sufficient.** A good PCReg/R² coexists with
locally separated batches. The global and local metric families are near-orthogonal —
median Spearman r ≈ 0.009 across the benchmark.

**2. Local metrics universally fail cross-platform.** kBET and iLISI score roughly 60% of
all tested approaches as "bad". This is a limitation of the problem, not of any method,
and no approach in the benchmark solved it.

**3. Three rankings disagree, and the article has not settled which it reports.** The
equal-weight composite, the literature-weighted composite and a Borda count produce top-10
sets with **zero overlap** — with each other or with the curated Best-15. Treat what is
here as "the manuscript's curated selection", never as "the optimum".

**4. At least one Best-15 approach corrupts biology.** Marker-gene correlation before and
after harmonization was ~1 for the linear methods but fell as low as **0.05** for one of
the clustermap-best approaches. Which one is recorded in `docs/BEST_APPROACHES.md`; until
you have checked, treat a high batch-mixing score with suspicion and read groups L and M
in your own output — which is why every config here computes them.

**5. No method dominates.** 26 of 31 methods are Pareto-optimal somewhere in the
batch-mixing × biology-preservation plane. There is no single right answer, which is the
reason this tool runs cross-products rather than one recommended pipeline.

**6. `softimpute` beats `knn` by roughly 10 percentage points** on PCR. Where both are
listed, prefer softimpute.

**7. The factor-importance hierarchy is version-dependent.** The manuscript's η² analysis
gives **method 36.3% > strategy 25.6% > imputation 1.6% > post-removal 0.2%**. An earlier
note ordered them differently; quote the manuscript numbers.

**8. The Best-15 is a manual curation** — visual inspection combined with composite
scoring, not the output of an automated rule.

The practical consequence: run the cross-product on *your* data
(`examples/configs/full_crossproduct.yaml`) and read the metrics, rather than adopting a
branch here because it won somewhere else.
