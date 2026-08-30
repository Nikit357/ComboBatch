# Plan: ComboBatch — a standalone imputation × harmonization × metrics tool

**Date:** 2026-08-28 (revised 2026-08-29, second review pass)
**Target repository:** `/Users/user890/Desktop/ComboBatch` → public GitHub repo `Nikit357/ComboBatch`
**Licence:** MIT
**Container registry:** GHCR only — `ghcr.io/nikit357/combobatch`
**Source repositories (read-only donors):**
- `../Follicular_lymphoma_disser/harmonization-scripts/` — 39 harmonization methods, 4 imputers, dispatcher pattern
- `../Follicular_lymphoma_disser/shambhala_adoption/Shambhala_containerized/` — Shambhala2 (pure-Python + Octave)
- `../Follicular_lymphoma_disser/harmonization-metrics-calculation/` — 14 metric groups (A–N)

---

## Overview

ComboBatch packages the FL dissertation's harmonization benchmark as a **general-purpose,
dataset-agnostic CLI tool**. A user points it at an expression matrix and an annotation table
(local or S3), names a biology column and a batch column, picks one or more imputation
strategies and one or more harmonization methods — with full hyperparameter control over both —
and gets back harmonized matrices plus, optionally, harmonization quality metrics. A dispatcher
runs the full cross-product in parallel. Everything ships in one Docker image named
`combobatch`, published publicly to GHCR, which the Karpenter pod then consumes instead of
installing R and Octave at startup.

**The three defining design decisions**, each of which the source code forces:

1. **Configuration replaces constants.** `bench_shared.py` hardcodes `BIO_COL`,
   `BATCH_COL`, `S3_BUCKET`, `S3_PREFIX`, the reference batch `"RNASeq_FF_PolyA"`, the 14 FL
   filter strategies, and `compute_batch_metrics.py` additionally hardcodes `BATCH_COLS`,
   `BIO_COLS`, the diagnosis literals `"Follicular_Lymphoma"` / `"Diffuse_Large_B_Cell_Lymphoma"`,
   and an S3-only storage path. All of this becomes explicit configuration.

2. **Single-source registries.** Today `METHODS` lives in `bench_shared.py` but `ALL_METHODS`
   is *retyped by hand* in two dispatchers; `ALL_STRATEGIES` and `ALL_IMPUTATION` likewise; and
   the metric groups have **no registry at all** — dispatch is a hardcoded sequence of
   `_run_group("A", …)` calls with a parallel `GROUP_SENTINEL_KEYS` dict in a different file.
   ComboBatch defines each registry exactly once and derives every CLI choice list from it.

3. **Every hyperparameter becomes reachable — which forces a new output key grammar.** 19 of the
   39 `normalize_*` functions declare Python-level hyperparameters and **all 39 already accept
   `**kw`**, but the only call site in production is frozen at
   `fn(exp, ann, batch_col=…, bio_col=…)`, so every one has always run at its default. Worse,
   the R-backed methods and imputers hardcode their R-side arguments in the interpolated R
   string, so parameters like `missForest(ntree=)` or `softImpute(rank.max=)` are unreachable
   even in principle. ComboBatch widens the Python call *and* threads parameters into the R
   calls. The consequence is that `{imp}__{method}__post{0|1}` is **no longer a unique output
   key** — two runs of `10_mnn` with `k=20` and `k=50` would collide — so hyperparameters are
   encoded into the output filename and recorded in `run_manifest.json`.

**Explicitly out of scope** (see §2.3): the 14 FL filter strategies, the FL marker-gene panel as
a *required* dependency, and the `01_raw` reference-download machinery.

---

## 1. Background — verified source inventory

Everything in this section was read from the source files, not inferred.

### 1.1 Harmonization methods (`bench_shared.py:2670`, constant `METHODS`)

`dict[str, tuple[callable, str]]` — key → `(function, harshness_tier)`. **39 entries.**
Tiers: `low` = 6, `medium` = 21, `high` = 12. Universal signature:

```python
def normalize_<name>(exp_df, ann_df, batch_col=BATCH_COL, [bio_col=BIO_COL,]
                     <hyperparameters with defaults>, **kw) -> pd.DataFrame
```

`exp_df` is **samples × genes** everywhere; R-backed methods transpose internally.

| Category | Method keys | ComboBatch disposition |
|---|---|---|
| Pure Python (14) | `01_raw`, `02_median_scaling`, `07_pycombat`, `08_inmoose_combatseq`, `11_harmony`, `12_scanorama`, `13_fsmvn`, `15_fsqn_py`, `17_quantile`, `18_rank`, `25_angel`, `26_xpn`, `30_recombat`, `39_procrustes` | ported, except `39_procrustes` (**removed**, §2.5) |
| R via `importr` (9) | `03_limma`, `04_sva`, `05_combat`, `06_combat_seq`, `09_ruv`, `10_mnn`, `14_qsmooth`, `22_tmm`, `23_vst` | ported |
| R via tempfile round-trip (13) | `16_fsqn_r`, `19_tdm`, `20_shambhala`, `21_harmonizr`, `27_dwd`, `28_npn`, `29_combat_ref`, `31_ruv3prps`, `33_amdbnorm`, `34_arsyn`, `36_explobatch`, `37_fabatch`, `38_harman` | ported, except `36_explobatch` (**removed**); `20_shambhala` **re-backed** by the vendored pure-Python + Octave implementation; `37_fabatch` **repaired** (§5.5a) |
| Always raise `NotImplementedError` (3) | `24_peer_k10`, `32_deepmnn`, `35_dasc` | **removed** (§2.5) |
| RNA-seq-only guard `_assert_rnaseq_only()` | `22_tmm`, `23_vst`, `39_procrustes` | kept as a capability flag on the two surviving methods |

**The 19 methods carrying Python-level hyperparameters** (name → params/defaults), all currently
unreachable:

| Key | Hyperparameters |
|---|---|
| `09_ruv` | `k=2` |
| `10_mnn` | `k=20` |
| `11_harmony` | `n_pcs=50` |
| `13_fsmvn`, `15_fsqn_py`, `16_fsqn_r`, `19_tdm`, `20_shambhala`, `29_combat_ref` | `target_group="RNASeq_FF_PolyA"` |
| `21_harmonizr` | `algorithm="ComBat"` \| `"limma"` |
| `24_peer_k10` | `n_factors=10` (method removed) |
| `25_angel` | `threshold=0.20` |
| `26_xpn` | `target_group`, `n_quantiles=50` |
| `27_dwd` | `target_group`, `min_batch_size=5` |
| `31_ruv3prps` | `k_factors=5`, `min_cell_size=2` |
| `33_amdbnorm` | `target_group` |
| `36_explobatch` | `maxdim=9` (method removed) |
| `38_harman` | `limit=0.1` |
| `39_procrustes` | `ec_batch_substr`, `coeffs_kit` (method removed) |

Note `target_group="RNASeq_FF_PolyA"` — an FL batch name — is the default of **nine** methods.
In ComboBatch this defaults to `None` and resolves to "largest batch" unless the user names one.

**A second, deeper layer of unreachable parameters.** Beyond the Python signatures above, every
R-backed method hardcodes its R-side arguments inside the interpolated R string — for example
`33_amdbnorm` fixes `polyFit(dis, 9)` and `genDistData(..., 500)`, `19_tdm` fixes
`log_target=FALSE`, `31_ruv3prps` fixes `k=5`. These are invisible to the Python signature and
cannot be reached at all today. §5.5 threads them through as well.

### 1.2 Imputation (`bench_shared.py:359` and `:384`)

Four strategies, but split across **two** functions:

- `prepare_dataset(ann_filter, exp_full)` → `strict`: `exp_full.loc[common].dropna(axis=1)`.
  Drops every gene with *any* NA. No threshold parameter.
- `prepare_dataset_imputed(ann_filter, exp_full, method="knn", knn_k=5, max_na_frac=0.20)` →
  `knn` | `missforest` | `softimpute`. Drops genes with `na_frac > max_na_frac` **before**
  imputing, then imputes. Passing `"strict"` here raises `ValueError`.

**Why the hyperparameters are unreachable.** Three independent causes, each needing a different fix:

| Imputer | Exposed in the signature | Actually reachable | Why not | ComboBatch fix |
|---|---|---|---|---|
| `strict` | — | — | Lives in a **different function** with no parameters at all | Becomes a registry entry with `max_na_frac=0.0` semantics |
| `knn` | `knn_k=5` | **none** | The parameter exists but **no caller ever passes it** — `run_prep_job.py:145` calls `prepare_dataset_imputed(ann_strat, comb_exp, method=imp)` and nothing else | Widen the call site; forward `**params` to `KNNImputer` |
| `softimpute` | none | none | The R call is the literal `softImpute(r_mat, type="svd")` — `rank.max`, `lambda`, `thresh`, `maxit` are **not in the Python signature at all**, so there is nothing to pass | Build the R argument list from `**params` |
| `missforest` | none | none | The R call is the literal `missForest(r_mat)` with **no arguments** — `ntree`, `maxiter`, `mtry` are unreachable in principle | Build the R argument list from `**params` |

`max_na_frac=0.20` is a parameter with a default that **no caller ever overrides**. The
requirement — "specify what maximum percentage of NA values should a gene have to be imputed" —
is exactly this parameter, promoted to the CLI.

### 1.3 Shambhala (`Shambhala_containerized/`)

- The reusable compute function is `shambhala.parallel.harmonize_parallel(input_df, p_df, rm, rs,
  k=5, n_workers=5, octave_scripts_dir=None, octave_bin="octave", timeout_s=6000,
  random_seed=None, disable_progress=False, fixed_clusters=None, python_cublock=False,
  skip_qn=False)`. **It takes and returns genes × samples**, unlike every other module in the
  package (samples × genes). The single transpose lives in `run_shambhala.py::main()`.
- There is **no single `harmonize()` entry point**. The ~40 lines of glue (gene intersection →
  NA strategy → P alignment → Q statistics → `harmonize_parallel` → NA restore) are **duplicated
  in four places** that have already drifted apart.
- P = quantile-normalization anchor (concatenated with each batch, then discarded).
  Q = target shape; only `rm = mean(log(Q+ε))` and `rs = std(log(Q+ε))` per gene are used.
- The 18 `shambhala_<P>_<Q>` variants are **benchmark-harness constants**
  (9 P × 2 Q in `shambhala_bench_shared.py`), not part of the library. The library takes P and Q
  as arbitrary file paths.
- **Scale contract (load-bearing):** `readExpressionData.m` applies `log2(x+1)` internally, so
  data piped to Octave must be **raw, non-log**. Output is `exp(rm + rs*log(x+1))` — also raw
  scale, and can reach ~100 000. There is no `--log-input` flag and no scale detection.
- Octave is required for everything except `--python-cublock`, which is **permanently xfailed**
  (ProcessPoolExecutor deadlock; 226 % mean relative difference vs Octave).
- **Packaging trap:** only `shambhala/` is packaged (`top_level.txt`). `octave/*.m` is **not**,
  so `harmonize_parallel`'s default `octave_scripts_dir` only resolves under `pip install -e .`.
- **Nested-parallelism trap:** total Octave processes = outer workers × `n_workers`. Documented
  guidance is `--n-workers 1 --n-shambhala-workers 16` on a 16-vCPU pod.

**P0 / Q0 calibration files — verified.** One canonical pair exists in three synced copies
(S3 `calibration_datasets/`, a local gzipped mirror at
`shambhala_adoption/Calibration_datasets/`, and the test fixtures). The `_small` fixture suffix
is a misnomer: `md5(P0_small.csv) == md5(P0_standard.csv)` and likewise for Q0 — they are the
same files.

| File | Shape | Orientation | Scale | Raw | Gzipped |
|---|---|---|---|---|---|
| `P0_standard.csv` | 39 samples × 11,768 genes | samples-as-rows | **raw linear** (min 4.44, max 15,596) | 7.47 MiB | **2.99 MiB** |
| `Q0_standard.csv` | 100 samples × 11,887 genes | samples-as-rows | **raw linear** (min 14.3, max 136,227) | 19.22 MiB | **7.16 MiB** |

P0 = 39 Affymetrix GPL570 healthy-tissue samples; Q0 = GTEx mixed-tissue RNA-seq. Their gene
sets differ (11,768 vs 11,887) and are intersected at load. **Total 10.1 MiB gzipped — trivially
shippable in the repo.**

### 1.4 Metric groups (`compute_batch_metrics.py`)

**There is no group registry to import.** The group universe is the literal string
`set("ABCDEFGHIJKLMN")` and dispatch is a hardcoded sequence of `_run_group(...)` calls at
`compute_batch_metrics.py:2356-2388`. The only group-keyed constants are
`GROUP_SENTINEL_KEYS` (14 entries) and `GROUPS_NEEDING_REFERENCE = {"L"}`, both in
`run_metrics_job.py`. **ComboBatch has to build this registry.**

Execution order is `E, K, A, J, B, C, D, G, H, I, L, M, N, F`.

| Grp | Function | Measures | Needs | Cost |
|---|---|---|---|---|
| E | `compute_group_e` | data quality, shapes, zero fractions, bimodality | ann | fast |
| K | `compute_group_k` | NA retention | — | fast |
| A | `compute_group_a` | PCA ANOVA R² per PC, PCR, DSC | ann, PCA, 999 perms | moderate |
| J | `compute_group_j` | % variance per PC 1–10 | PCA only (**no ann**) | fast |
| B | `compute_group_b` | kBET, iLISI, cLISI, ASW, CMS | ann, PCA | moderate |
| C | `compute_group_c` | UMAP/tSNE centroid dispersion, entropy | ann, UMAP+tSNE | **slow** |
| D | `compute_group_d` | pairwise KS, within-batch cohort effects, per-gene CV | ann | **slow** (own 3600 s timeout) |
| G | `compute_group_g` | kNN graph connectivity per biology group | ann, PCA | fast |
| H | `compute_group_h` | intra/inter Euclidean distance ratios | ann | moderate |
| I | `compute_group_i` | WaterMelon entropy score | ann, PCA, 200 perms | slow (`--skip-wm`) |
| L | `compute_group_l` | marker-gene correlation preservation | ann + **reference matrix** + panel | fast |
| M | `compute_group_m` | cross-batch rank agreement | ann + panel | moderate |
| N | `compute_group_n` | leave-one-batch-out prediction + perm control | ann, `n_perm` | slow |
| F | `compute_group_f` | variancePartition | ann, **R subprocess** | **very slow**, `--skip-slow` |

**Which columns each group consumes** — this is what makes the columns user-specifiable
(§5.8, and documented per-group in `docs/METRICS.md`):

| Constant | Value in the source | Consumed by |
|---|---|---|
| `BATCH_COLS` | `["RNA_BATCH", "PLATFORM_RNA", "RNASEQ_SOURCE", "COHORT_LABEL"]` | A, B (kBET/iLISI/ASW_batch), C, H |
| `BIO_COLS` | `["Major_group", "PLATFORM_RNA", "Diagnosis_cell_type_unified", "TUMOR_NORMAL"]` | B (cLISI/ASW_bio), G, H |
| `ALL_COLS` | union of the two (7 names) | PCR, F, H |
| `WM_BATCH_COLS` / `WM_BIO_COLS` | duplicated copies of the above | I |
| bare `"RNA_BATCH"` literal | — | E, D2/D3, I subsampling, M, N, and every sentinel key |
| bare `"COHORT_LABEL"` literal | — | E, D2, F (forced random effect), L |
| `PRED_CLASS_COL` / `PRED_BIO_COL` | `"Major_group"` / `"Diagnosis_cell_type_unified"` | N, M |

`compute_all_metrics()` runs each group in an isolated try/except, writes `error_<LETTER>` on
failure, flushes partial JSON after every group, and calls an `on_group_done` callback.
`compute_group_l` has the argument-order trap `(exp, ref, ann)` while all others are `(exp, ann, …)`.

**kBET and LISI are pure-Python reimplementations** — no R dependency. Only Group F touches R,
via `Rscript` subprocess (not rpy2). `rpy2>=3.5.5` in the metrics requirements is vestigial.

### 1.5 Dependency conflict analysis (verified)

| Package | harmonization | metrics | shambhala | Resolution |
|---|---|---|---|---|
| **pandas** | `>=1.3.4,<2.0.0` | `>=2.0` | `>=1.5` | **`>=1.5,<2.0`** — see below |
| numpy | `>=1.24,<2.0` | `>=1.24` | `>=1.23,<2.0` | `>=1.24,<2.0` |
| scikit-learn | `>=1.3,<1.6` | `>=1.3` | `>=1.3` | `>=1.3,<1.6` |
| awscli | `>=1.29` | `>=1.29` | (pod: `>=1.44`) | `>=1.44` |

The pandas conflict is the only hard blocker, and it is **resolvable**. `reComBat 0.1.4`
genuinely requires `pandas<2.0`, and pandas 1.5.3 ships a **py3.11 wheel only**. The metrics
`>=2.0` floor, however, is unjustified: I audited the API surface of all six metrics modules and
found only `pd.DataFrame` (85), `pd.Series` (6), `pd.read_csv` (5), `pd.notna` (2),
`pd.read_excel` (1), plus conservative methods (`astype`, `join`, `nunique`, `where`, `items`,
`dropna`, `value_counts`, `groupby`, …). **Zero pandas-2.x-only APIs** — no `applymap`,
`dtype_backend`, `ArrowDtype`, `future_stack`, `copy_on_write`, or `convert_dtypes`.

**Therefore: one image, Python 3.11 exactly, `pandas>=1.5,<2.0`.** Python 3.12 is ruled out by
pandas 1.5.3; Python 3.10 by Shambhala's `requires-python = ">=3.11"`.

### 1.6 Infrastructure facts

- Existing Dockerfile: `python:3.11-slim`, R from CRAN `bookworm-cran40` **unpinned** (drift
  hazard — see §7), no Octave, `CMD ["sleep", "infinity"]`, build context = repo root.
- `install_r_packages.R` ends with a `stop()` on any missing package, and currently installs and
  verifies **three packages that cannot succeed**: `FAbatch` (does not exist on Bioconductor or
  CRAN — the real package is `bapred`, on **CRAN**), R `reComBat` (it is a *Python* package), and
  `exploBATCH` (hard dep `fMM`, GitHub repo deleted). Any `docker build` reusing it as-is fails.
- Karpenter usage in live manifests is thin: `nodeSelector: karpenter.k8s.aws/instance-family: c6a`
  in three files, plus a `node-group=rnd-sandbox` toleration + required nodeAffinity pair in every
  pod. `capacity-type: on-demand` is documented as needed in an old plan but **never applied**.
- The donors publish to a private ECR
  (`028257207274.dkr.ecr.us-east-1.amazonaws.com/fl-batch-correction`), which requires node-IAM
  pull auth and an admin `ecr:CreateRepository` request. **ComboBatch does not use ECR** — see
  §5.10.
- Local machine verified 2026-08-29: `docker` 29.7.2, `kubectl`, `aws` CLI 2.34.33 with valid
  credentials, and **Python 3.11.16** installed via Homebrew at `/opt/homebrew/opt/python@3.11`.
  Note `python3` still resolves to the system 3.9.6, so scripts and docs invoke `python3.11`
  explicitly. Unit tests therefore run natively on this Mac; only the R- and Octave-backed
  integration tests need the container.
- **The Mac is arm64; the c6a nodes are amd64.** Every image build must be
  `docker buildx build --platform linux/amd64`, or the pod fails with `exec format error`.

---

## 2. Scope decisions

### 2.1 Ported nearly verbatim

- The 34 retained `normalize_*` implementations and their R interop helpers (`_py2rpy`,
  `_rpy2py`, `_r_gc`, `_assert_rnaseq_only`).
- All four imputation code paths.
- All 14 metric group functions and the embedding helpers.
- The dispatcher pattern: `ThreadPoolExecutor` supervising `subprocess.Popen` workers (this is
  what gives each rpy2 job a fresh R session — it is not incidental and must be preserved),
  the cgroup-aware `_free_gb()` memory guard, the timeout/exit-code contract, the failed-jobs log.
- Shambhala's `octave_bridge`, `parallel`, `q_rescale`, `na_handling`, `progress_display`,
  and the eight `octave/*.m` files.

### 2.2 Generalized

| Source hardcoding | ComboBatch |
|---|---|
| `BIO_COL`, `BATCH_COL` module constants | `--bio-col`, `--batch-col`, `--cohort-col` |
| `BATCH_COLS` / `BIO_COLS` / `ALL_COLS` lists | `--metric-batch-cols`, `--metric-bio-cols`, `--metric-cohort-col`, `--predict-class-col` (default: the single configured columns) |
| `S3_BUCKET` / `S3_PREFIX` constants, S3-only | URI-dispatched `StorageBackend` (local path or `s3://…`), `--endpoint-url` for MinIO |
| `target_group="RNASeq_FF_PolyA"` (9 methods) | `--reference-batch`; default `None` → largest batch |
| Group N literals `"Follicular_Lymphoma"` / `"Diffuse_Large_B_Cell_Lymphoma"` | `--predict-classes`, default = top-N most frequent levels of the configured class column |
| `marker_gene_annotation.csv` required at **import** time | `--marker-panel PATH`, optional; absent → L and M skipped with a clear reason. FL panel shipped under `examples/` |
| Group L needs `01_raw__post0` downloaded from S3 | The unharmonized input matrix is already in hand — pass it directly, no second download |
| `GROUP_SENTINEL_KEYS` naming tied to `RNA_BATCH` | Sentinels derived from the configured column names |
| `--ref-cache-dir /workspace/…`, `--tmp-dir /tmp/bench_metrics_jobs` | `--work-dir`, defaulting to a `tempfile` location |
| `POD_NAME` env for the failed-log key | `--run-id`, defaulting to `$POD_NAME` then `local` |
| Hardcoded R-side arguments inside interpolated R strings | Built from `**params` (§5.5) |

### 2.3 Deliberately NOT ported

- **The 14 FL filter strategies** (`build_filter_strategies`, `RARE_GROUPS`, `NORMAL_GROUPS`,
  `BAD_BATCHES_A/B`, `_AFFY_EXT_BATCHES`, `MIN_BATCH_SIZE`). These encode FL dataset knowledge and
  are meaningless for an arbitrary dataset. Subsetting is done with a generic `--subset-query`
  (§5.15), a pandas `DataFrame.query` string over the user's own annotation columns.
- **The `load_data()` FL fixups**: the derived `Diagnosis_with_coo` columns and the hardcoded
  `PUB_Suntsova_GSE120795 → RNASeq_FF_Unknown` batch relabel.
- **The 18 `shambhala_<P>_<Q>` variant keys.** ComboBatch exposes Shambhala as one method taking
  `P` and `Q` as parameters, defaulting to the shipped P0/Q0 (§5.6), which is strictly more general.
- **`make_slides.py`**, the notebook helpers, and `harmonization-tools-research/`.
- **`--skip-shambhala`** and the `01_raw` magic-baseline logic.

### 2.4 Kept as an option, defaulted off

**Post-removal** (`post0`/`post1`). A generic technique — drop the *n* most PCA-deviant batches
*after* harmonization via `identify_outlier_batches` — and the metrics layer already threads a
`post_rm` axis. It becomes:

| Flag | Default | Meaning |
|---|---|---|
| `--post-removal` | off | enable the `post1` variant at all |
| `--post-removal-n` | 1 | **how many batches to exclude** (source `n_outliers=1`) |
| `--post-removal-min-batch` | 20 | batches smaller than this are never candidates |
| `--post-removal-n-pcs` | 2 | PCs used for the centroid-distance ranking |

When post-removal is off, only the `post0` output is written and the `postN` segment is omitted
from output keys entirely.

### 2.5 Methods removed from the registry

The registry ships **34 methods**, down from 39. Original keys are preserved (so the numbering
has gaps at 24, 32, 35, 36, 39) because they are the identifiers used throughout the article and
the published `metrics_comprehensive.csv`; renumbering would break traceability.

| Removed | Why |
|---|---|
| `24_peer_k10` | `normalize_peer` unconditionally raises — PEER is unavailable for R 4.5 |
| `32_deepmnn` | unconditionally raises — scRNA-seq only, not applicable to bulk |
| `35_dasc` | unconditionally raises — DASC returns cluster assignments, not an expression matrix |
| `36_explobatch` | unconditionally raises — its hard dependency `fMM` had its GitHub repo deleted (HTTP 404), so it cannot be revived |
| `39_procrustes` | **Licensing.** The `BostonGene/Procrustes` repo carries a proprietary *"Limited License and Terms of Use for Research Use by Non-Profit and Government Institutions"* — no commercial use, eligibility restricted to students/post-docs/faculty at degree-granting or US-government institutions, and an explicit prohibition on altering or removing copyright notices. Vendoring it into an MIT-licensed public tool is not a redistribution ComboBatch is licensed to make, and a runtime-fetch workaround would still make an advertised capability depend on a licence most users cannot accept. Excluded outright. |

A tool should not advertise a method that can never run, nor one that most of its users are not
permitted to use. All five are documented in `docs/METHODS.md` under "not included, and why", so
the article's 39-method count stays explicable. The Procrustes entry additionally records what
was learned while evaluating it, so the decision is auditable: the repo is public and clones
without auth; it is a flat 2-module package (`procrustes_bg/transform.py`, `plot.py`) needing
only `numpy` and `pandas`; only `V4_coefficients.json` and `V7_UTR_coefficients.json` are
runtime-required (31.7 MB raw / 10.5 MB gzipped, the rest of the 319 MB tree being demo data);
`coeffs_kit="V7"` is invalid despite the source docstring advertising it; and the model is not
the per-gene slope+intercept the docstring claims but a multivariate fit with 1–373 predictor
genes per target gene. If a redistribution grant is ever obtained, that note is the starting
point for adding it back.

`20_shambhala` is **re-backed**: the source's key `20_shambhala` called the `Shambhala2` **R**
package through a tempfile round-trip. ComboBatch keeps the key but points it at the vendored
pure-Python + Octave implementation (§5.6), which is faster, better tested (69 tests), and drops
an R dependency. One Shambhala, not two. See §7.9 for the reproducibility consequence.

---

## 3. Target repository layout

Every subdirectory carries both a `README.md` (for users) and a `CLAUDE.md` (for agentic work).

```
ComboBatch/
├── README.md                          # install, quickstart, CLI reference, worked examples
├── CLAUDE.md                          # top-level architecture + constraints
├── LICENSE                            # MIT
├── pyproject.toml
├── requirements.txt
├── .dockerignore
├── .github/workflows/ci.yml
│
├── docs/
│   ├── README.md                      # index of the docs set
│   ├── HYPERPARAMETERS.md             # every method + imputer parameter, with examples
│   ├── METRICS.md                     # every metric group A–N, its columns and outputs
│   ├── METHODS.md                     # the 34 methods, citations, "not included, and why"
│   ├── OUTPUT_LAYOUT.md               # key grammar, run_manifest, sidecar schema
│   └── BEST_APPROACHES.md             # the benchmark's recommended configurations
│
├── combobatch/
│   ├── README.md  CLAUDE.md
│   ├── __init__.py                    # __version__, top-level re-exports
│   ├── config.py                      # dataclasses + YAML/CLI resolution
│   ├── storage.py                     # StorageBackend ABC, LocalBackend, S3Backend
│   ├── dataio.py                      # read/write expression + annotation, alignment
│   ├── logging_utils.py               # timestamped logging, thread pinning
│   ├── memory.py                      # cgroup-aware _free_gb, wait_for_memory, check_memory
│   ├── params.py                      # hyperparameter parsing, coercion, filename encoding
│   ├── subset.py                      # --subset-query evaluation and validation
│   │
│   ├── imputation.py                  # IMPUTER_REGISTRY: strict, knn, softimpute, missforest
│   ├── postremoval.py                 # identify_outlier_batches + apply_post_removal
│   │
│   ├── methods/
│   │   ├── README.md  CLAUDE.md
│   │   ├── __init__.py                # METHOD_REGISTRY (single source of truth, 34 entries)
│   │   ├── base.py                    # MethodSpec, HyperParam
│   │   ├── rinterop.py                # _py2rpy, _rpy2py, _r_gc, r_arglist, r_available
│   │   ├── python_methods.py          # the 13 pure-Python normalize_*
│   │   ├── r_methods.py               # the 20 R-backed normalize_*
│   │   └── shambhala_method.py        # adapter → combobatch.vendor.shambhala
│   │
│   ├── vendor/
│   │   └── shambhala/                 # README.md CLAUDE.md + library + octave/*.m (packaged!)
│   │
│   ├── data/
│   │   └── calibration/               # P0_standard.csv.gz, Q0_standard.csv.gz (10.1 MiB)
│   │
│   ├── harmonize.py                   # run_one_combination(): impute → harmonize → post-rm
│   │
│   ├── metrics/
│   │   ├── README.md  CLAUDE.md
│   │   ├── __init__.py                # METRIC_GROUP_REGISTRY (built here; none upstream)
│   │   ├── groups.py                  # the 14 compute_group_* functions
│   │   ├── embeddings.py              # compute_pca / compute_umap / compute_tsne
│   │   ├── panels.py                  # optional marker-panel loader
│   │   ├── runner.py                  # compute_metrics() orchestrator + on_group_done
│   │   └── concat.py                  # sidecars → summary CSV + long tables
│   │
│   └── cli/
│       ├── README.md  CLAUDE.md
│       ├── main.py                    # `combobatch` root parser
│       ├── run_cmd.py                 # combobatch run
│       ├── dispatch_cmd.py            # combobatch dispatch      ← the cross-product dispatcher
│       ├── metrics_cmd.py             # combobatch metrics
│       ├── concat_cmd.py              # combobatch concat
│       └── selftest_cmd.py            # combobatch selftest
│
├── examples/
│   ├── README.md                      # how to run every config below
│   ├── configs/
│   │   ├── minimal.yaml
│   │   ├── full_crossproduct.yaml
│   │   ├── hyperparameters.yaml       # worked hyperparameter-sweep example
│   │   ├── shambhala.yaml
│   │   └── best_approaches/           # §5.15 — the benchmark's recommendations
│   │       ├── README.md
│   │       ├── ff_only_mnn.yaml
│   │       ├── ffpe_only_sva.yaml
│   │       ├── rnaseq_only_sva.yaml
│   │       ├── rnaseq_plus_illumina_fsqn.yaml
│   │       └── multiplatform_mnn.yaml
│   ├── example_dataset/               # 5 GEO cohorts, ~13.9 MB gzipped (§5.12)
│   │   ├── README.md
│   │   ├── PROVENANCE.md
│   │   ├── example_exp.tsv.gz
│   │   └── example_ann.csv
│   └── marker_panel_fl.csv            # the FL panel, as an *example* not a dependency
│
├── tests/
│   ├── README.md  CLAUDE.md
│   ├── conftest.py
│   ├── unit/  integration/  e2e/  fixtures/
│
├── docker/
│   ├── README.md
│   ├── Dockerfile
│   ├── install_r_packages.R
│   └── entrypoint.sh
│
├── k8s/
│   ├── README.md  CLAUDE.md
│   ├── pod-combobatch.yaml
│   ├── job-combobatch.yaml
│   └── create_pvc.yaml
│
└── scripts/
    ├── README.md
    ├── build_example_dataset.py
    ├── build_and_push_image.sh
    └── smoke_test.sh
```

---

## 4. Core data contracts

### 4.1 Input

- **Expression**: rows = samples, columns = genes (FL convention). `.tsv`, `.tsv.gz`, `.csv`,
  `.csv.gz`, `.parquet`. Separator inferred from extension; `index_col=0`.
- **Annotation**: rows = samples, index must intersect the expression index.
  `.csv`, `.tsv`, `.csv.gz`.
- Both accept a local path or an `s3://bucket/key` URI, resolved by `storage.resolve()`.
- Alignment is by **index intersection**, never positional; duplicates dropped `keep="first"`;
  the tool asserts `len(exp) == len(ann)` afterwards and reports how many samples were dropped.

### 4.2 Output layout (identical shape local and S3)

```
{out_root}/
├── run_manifest.json                             # ← resolved params for every output
├── prepared/{imp}[__{ptag}]__exp.tsv.gz          # post-imputation, pre-harmonization
├── prepared/{imp}[__{ptag}]__ann.tsv.gz
├── exp/{imp}__{method}[__{ptag}][__post{0|1}].tsv.gz
├── metrics/{same-stem}_metrics.json
├── genes/{same-stem}_genes.json
├── metrics_summary.csv
├── marker_gene_correlations_long.csv             # only if Group L ran
├── prediction_folds_long.csv                     # only if Group N ran
└── failed_jobs_{run_id}.txt
```

**Hyperparameters are written into the filename.** `ptag` is:

1. **omitted entirely** when both the method and the imputer run at pure defaults, so default
   runs keep the familiar `{imp}__{method}__post0` names;
2. otherwise a **human-readable encoding** of the non-default parameters, in sorted key order,
   filesystem-safe: `k=50, limit=0.05` → `k50-limit0.05`. Floats keep their decimal point (safe
   on every filesystem ComboBatch targets); `/`, spaces and `__` are percent-escaped; `True` /
   `False` become `T` / `F`; a list becomes `+`-joined;
3. **truncated with a hash suffix** when the readable form would exceed 64 characters:
   `k50-limit0.05-…-a1b2c3d4`. This keeps names readable in the common case (one or two tuned
   parameters) without risking a filesystem name-length limit in the pathological case.
4. overridable with `--param-tag mytag` for full control.

`run_manifest.json` maps every `ptag` back to its complete resolved parameter dict — including
the defaults that were *not* encoded — so a result is always self-describing even when its name
is abbreviated.

> Note the source's filename grammar splits on `"__"` and requires **exactly 4 parts**; anything
> else is silently dropped from job enumeration. ComboBatch's parser must handle 2–5 segments and
> **fail loudly** on an unparseable key rather than skipping it silently.

### 4.3 Sidecar JSON

Flat dict, superset of the source schema: `imp`, `method`, `param_tag`, `params` (the full
resolved dict), `post_rm`, `harshness`, `status`, `compute_time_s`, `n_samples`, `n_genes`, plus
every metric key. `status ∈ {ok, skipped, failed, cached, upload_failed, post_removal_failed}`.
**Serialization fix:** the source uses bare `json.dump`, which emits `NaN` — invalid strict JSON
that non-Python parsers reject. ComboBatch uses a NaN-safe encoder (NaN/Inf → `null`) everywhere.

---

## 5. Files to create

### 5.1 `pyproject.toml`, `requirements.txt`

Package `combobatch`, version `0.1.0`, `requires-python = ">=3.11,<3.12"`, **MIT**.
Entry point: `[project.scripts] combobatch = "combobatch.cli.main:main"`.

**Core dependencies are R-free; system-dependent packages live in extras.** This was forced
during Phase 0: `rpy2` cannot be built without R on the `PATH`, and R is not installed on the
dev Mac, so a flat dependency list would make `pip install -e ".[dev]"` fail — and with it the
plan's own requirement (§5.13) that the unit suite run natively. The split is therefore:

| Target | Contents |
|---|---|
| `dependencies` | numpy, pandas, scikit-learn, scipy, statsmodels, umap-learn, psutil, pyyaml, boto3, botocore — all pure-Python-installable |
| `.[r]` | `rpy2` — the 20 R-backed methods, the R imputers, metric group F |
| `.[methods]` | harmonypy, scanorama, combat, inmoose, reComBat |
| `.[shambhala]` | qnorm |
| `.[panels]` | openpyxl, xlrd |
| `.[dev]` | pytest, pytest-timeout, moto[s3], black |
| `.[all]` | r + methods + shambhala + panels |

`requirements.txt` remains the **image's** complete environment (it assumes R and Octave are
present) and carries `awscli`, which the package itself never imports. The version constraints
are identical in both files; `tests/unit/test_docs_in_sync.py` will assert they do not drift.

**Critical packaging fix**: include the Octave `.m` files and the calibration data, which the
source package omits —

```toml
[tool.setuptools.package-data]
"combobatch.vendor.shambhala" = ["octave/*.m"]
"combobatch.data" = ["calibration/*.csv.gz"]
"combobatch.metrics" = ["data/*.csv"]
```

Unified `requirements.txt` (from §1.5), with a comment on every pin explaining *why*:

```
# --- core ---
numpy>=1.24,<2.0          # <2.0: rpy2/scanorama/numba ABI
pandas>=1.5,<2.0          # <2.0: reComBat 0.1.4 hard requirement; py3.11 wheel only
scikit-learn>=1.3,<1.6
scipy>=1.9
statsmodels>=0.14
umap-learn>=0.5
psutil>=5.9
pyyaml>=6.0
# --- storage ---
boto3>=1.28.17
botocore>=1.31.17
awscli>=1.44
# --- R bridge ---
rpy2>=3.5.5,<3.7          # ABI-bound to R 4.5.x
# --- harmonization ---
harmonypy>=0.0.9,<2.0
scanorama>=1.7
combat>=0.3.3
inmoose>=0.9.1
git+https://github.com/BorgwardtLab/reComBat
# --- shambhala ---
qnorm>=0.4.1              # pulls numba + llvmlite (~56 MB)
# --- panels ---
openpyxl>=3.1
xlrd>=2.0.1
```

### 5.2 `combobatch/config.py`

Frozen dataclasses, constructed from a YAML file, then overlaid by CLI flags (CLI wins):

```python
@dataclass(frozen=True)
class ColumnSpec:
    batch: str
    bio: str
    cohort: str | None = None
    # Metric-side column selection — the generalization of BATCH_COLS / BIO_COLS.
    metric_batch_cols: tuple[str, ...] = ()   # default: (batch,)
    metric_bio_cols: tuple[str, ...] = ()     # default: (bio,)
    predict_class_col: str | None = None      # default: bio

@dataclass(frozen=True)
class ImputationSpec:
    name: str                      # strict | knn | softimpute | missforest
    params: dict[str, object] = field(default_factory=dict)
    max_na_frac: float = 0.20

@dataclass(frozen=True)
class MethodSpec:
    name: str                      # a METHOD_REGISTRY key
    params: dict[str, object] = field(default_factory=dict)
    param_tag: str | None = None

@dataclass(frozen=True)
class MetricsSpec:
    enabled: bool = False
    groups: frozenset[str] = frozenset("ABEJK")   # fast default set
    skip_slow: bool = True
    skip_wm: bool = True
    n_perm: int = 20
    marker_panel: str | None = None
    predict_classes: tuple[str, ...] | None = None

@dataclass(frozen=True)
class PostRemovalSpec:
    enabled: bool = False
    n_batches: int = 1             # --post-removal-n
    min_batch_size: int = 20
    n_pcs: int = 2

@dataclass(frozen=True)
class RunConfig:
    exp_uri: str; ann_uri: str; out_uri: str
    columns: ColumnSpec
    imputations: tuple[ImputationSpec, ...]
    methods: tuple[MethodSpec, ...]
    metrics: MetricsSpec
    post_removal: PostRemovalSpec
    reference_batch: str | None = None
    subset_query: str | None = None
    random_seed: int = 42
    ...
```

`validate()` raises on: unknown method/imputer/group key, a column absent from the annotation,
a `--reference-batch` absent from the batch column, a `--subset-query` referencing a nonexistent
column or matching zero rows, an unknown hyperparameter name, and an unwritable `out_uri`.
**Validation happens before any expensive work** — the source discovers a bad strategy name only
after downloading 1.9 GB.

Example config (`examples/configs/full_crossproduct.yaml`):

```yaml
input:
  expression: ./examples/example_dataset/example_exp.tsv.gz
  annotation: ./examples/example_dataset/example_ann.csv
output: ./out/
columns:
  batch: RNA_BATCH
  bio: Diagnosis_cell_type_unified
  cohort: COHORT_LABEL
reference_batch: RNASeq_FF_rRNADepletion
subset_query: null
imputations:
  - name: strict
  - name: knn
    max_na_frac: 0.30
    params: {knn_k: 10, weights: distance}
  - name: softimpute
    params: {rank_max: 30, lambda: 1.0}
  - name: missforest
    params: {ntree: 50, maxiter: 5}
methods:
  - name: 10_mnn
    params: {k: 50}
  - name: 04_sva
  - name: 38_harman
    params: {limit: 0.05}
  - name: 20_shambhala          # P/Q default to the shipped P0/Q0
    params: {k: 5, n_workers: 1}
metrics:
  enabled: true
  groups: [A, B, E, J, K, L]
  marker_panel: ./examples/marker_panel_fl.csv
post_removal:
  enabled: true
  n_batches: 2
```

### 5.3 `combobatch/storage.py`

```python
class StorageBackend(ABC):
    def exists(self, key: str) -> bool: ...
    def read_bytes(self, key: str) -> bytes: ...
    def write_bytes(self, key: str, data: bytes) -> None: ...
    def list_keys(self, prefix: str) -> Iterator[str]: ...
    def delete(self, key: str) -> None: ...

class LocalBackend(StorageBackend): ...      # pathlib; mkdir -p on write
class S3Backend(StorageBackend): ...         # boto3; one client per backend instance

def resolve(uri: str, *, endpoint_url: str | None = None) -> tuple[StorageBackend, str]
```

Three fixes over the source:
1. **The source creates a fresh `boto3.client("s3")` at 15 separate call sites.** One client per
   backend, created once, reused.
2. **`s3_exists` swallows all `ClientError`s**, making a 403 indistinguishable from a 404.
   ComboBatch inspects the error code and re-raises anything that is not `404`/`NoSuchKey` — a
   permissions problem must not silently look like "output not yet computed".
3. `--endpoint-url` support, enabling MinIO and `moto`-based tests.

### 5.4 `combobatch/imputation.py`

`IMPUTER_REGISTRY: dict[str, ImputerSpec]` with **all four strategies as first-class registry
entries** behind one signature — removing the source's `strict`-lives-in-a-different-function
split:

```python
def impute(exp_df, *, method="strict", max_na_frac=0.20, params=None) -> tuple[pd.DataFrame, dict]
```

Returns the matrix **and a report** (`n_genes_in`, `n_genes_dropped_na_frac`, `n_genes_out`,
`n_cells_imputed`, `residual_na`, `elapsed_s`) which goes into the sidecar.

**Complete hyperparameter surface after the fix** (see §1.2 for why none of this is reachable
today; all are documented with worked examples in `docs/HYPERPARAMETERS.md`):

| Imputer | Parameters | Mechanism |
|---|---|---|
| `strict` | `max_na_frac` (fixed at `0.0`) | pure pandas; a gene with any NA is dropped |
| `knn` | `max_na_frac`, `knn_k`, `weights`, `metric` | forwarded to `sklearn.impute.KNNImputer` |
| `softimpute` | `max_na_frac`, `rank_max`, `lambda`, `thresh`, `maxit`, `type` | built into the R `softImpute()` call |
| `missforest` | `max_na_frac`, `ntree`, `maxiter`, `mtry`, `parallelize` | built into the R `missForest()` call |

Note `strict` **is** an imputer in the registry — it is the `max_na_frac = 0.0` corner of one
uniform rule, not a separate code path. This makes `--imputations strict,knn,softimpute,missforest`
a homogeneous list and lets the dispatcher treat all four identically.

Three behaviour fixes:
- **No silent fallback.** The source wraps imputation in `try/except` and falls back to strict
  **while still labelling the output `__knn__`**. ComboBatch fails the job loudly; a run labelled
  `knn` must have had KNN applied.
- **Warn on the no-op path.** If the `max_na_frac` filter already removed every NA, the source
  returns before imputing, so a `knn` run can be byte-identical to a `strict` run with no
  indication. ComboBatch records `imputation_was_noop: true` in the sidecar.
- **Residual NA is reported, not hidden.** KNN and softimpute can leave NaN where a whole batch
  is NaN for a gene; the count goes into the report.

### 5.5 `combobatch/methods/`

`METHOD_REGISTRY: dict[str, MethodSpec]` — the **single** source of truth, 34 entries.
`MethodSpec` carries what the source's bare `(fn, tier)` tuple cannot:

```python
@dataclass(frozen=True)
class MethodSpec:
    key: str; fn: Callable; harshness: str
    requires_r: bool
    r_packages: tuple[str, ...] = ()
    requires_octave: bool = False
    rnaseq_only: bool = False
    uses_reference_batch: bool = False
    hyperparams: dict[str, HyperParam] = ...   # name → (type, default, help, choices/bounds)
    citation: str = ""
    status: str = "available"
```

This drives, for free: `combobatch list-methods` (with each method's tunable parameters),
CLI validation of `--method-params` **before** a job launches, `--require-r/--no-r` filtering,
the generated `docs/HYPERPARAMETERS.md` tables, and a preflight that checks R/Octave availability
once at dispatcher startup instead of discovering it 400 jobs in.

**Reaching the R-side parameters.** `rinterop.r_arglist(params)` renders a validated R argument
list from a Python dict, which each R-backed method interpolates into its call:

```python
# before:  res <- polyFit(dis, 9)
# after:   res <- polyFit(dis, {r_arglist({"degree": degree})})
```

Every literal currently frozen inside an R string (`polyFit(dis, 9)`, `genDistData(..., 500)`,
`log_target=FALSE`, `RUVIII(k=5)`, `harman(limit=0.1)`, …) becomes a declared `HyperParam`.
`r_arglist` quotes strings, maps `True/False` → `TRUE/FALSE` and `None` → `NULL`, and **rejects
anything that is not a scalar or a list of scalars** — the R snippets are f-string-interpolated,
so this is the injection boundary and it is enforced in one place.

`normalize_*` bodies are ported verbatim except:
- `target_group` defaults change from `"RNASeq_FF_PolyA"` to `None`, with a shared
  `resolve_reference_batch(ann, batch_col, requested)` helper → requested, else largest batch,
  logging the choice.
- `normalize_dwd` keeps the removed-`penalty` fix (`bench_shared.py:1843`).
- The `exp_in = exp_df.dropna(axis=1) if …` NaN guard is factored into one decorator, and the
  number of genes it drops is **recorded** rather than silently discarded.

#### 5.5a Repairing `37_fabatch`

**The method is currently a silent no-op — verified against the `bapred` 1.1 manual.** The real
signature is `fabatch(x, y, batch, nbf=NULL, minerr=1e-06, probcrossbatch=TRUE, maxiter=100,
maxnbf=12)` with the adjusted matrix in **`$xadj`**. The source calls
`fabatch(xtr=, ytr=, batch=, type="among")` and reads `$adj.data`. There is no `xtr`, no `ytr`,
no `type` argument, and no `batchadjust()` function in the package at all. R raises
"unused arguments", the `tryCatch` handler fires, and `write.csv(exp_mat)` writes the
**uncorrected input**, which Python then returns as if it were corrected. Two independent
fail-silent bugs in the same direction. (No published number is contaminated: the benchmark's
`metrics.csv` has zero rows for `37_fabatch`.)

Required changes:

| Item | Now | Fixed |
|---|---|---|
| R call | `fabatch(xtr=…, ytr=…, batch=…, type="among")` | `fabatch(x = t(exp_mat), y = <2-level factor>, batch = batches)` |
| Result field | `res$adj.data` | `res$xadj` |
| Install | `BiocManager::install("FAbatch")` — package does not exist | `BiocManager::install("bapred")` — **CRAN** v1.1, GPL-2, needs Bioc deps `sva`, `affyPLM`, `affy`, `Biobase` |
| Verify list | `"FAbatch"` | `"bapred"` (`test_mock.py` already had this right) |
| Error handling | `tryCatch` writes the uncorrected matrix | re-raise; the worker records a real failure |
| Docs | `project_overview.md:164` says `batchadjust(y, batch, type="among")` | no such function — corrected in `docs/METHODS.md` |

**The `y` constraint — resolved.** `fabatch` requires `y` to be a **two-level factor** (a binary
target), so passing the batch vector only works when the dataset has exactly two batches.
ComboBatch adds `--fabatch-target-col` (default: the configured bio column) and:
- if the resolved target has exactly 2 levels → run;
- otherwise → raise `NotImplementedError` naming the column and its level count, so the job is
  recorded as an **honest SKIP** rather than a silent no-op.

Auto-binarizing (largest class vs rest) was considered and rejected: it would always run but
would silently change what the method is being asked to do. Constant-within-batch genes are
dropped up front (`fabatch` reports them in `res$badvariables`).

### 5.6 `combobatch/methods/shambhala_method.py` + `combobatch/vendor/shambhala/`

Vendor `shambhala/*.py` and `octave/*.m` unchanged, then write **the one glue function that the
source never had** (its logic currently exists in four drifted copies):

```python
def shambhala_harmonize(
    exp_df,                       # samples × genes, RAW (non-log) scale
    p_df=None, q_df=None,         # samples × genes; None → the shipped P0 / Q0
    *, k=5, n_workers=1, na_strategy="drop", knn_k=5, max_na_frac=0.20,
    q_pseudocount=1e-6, random_seed=None, octave_bin="octave",
    octave_scripts_dir=_VENDORED_OCTAVE_DIR, timeout_s=6000,
    disable_progress=True, precompute_qn_reference=False,
    synthetic_cublock_p=False, max_p_samples=None,
    precompute_cublock_clusters=False,
) -> pd.DataFrame
```

Registered as method key `20_shambhala`, **replacing** the source's R-package-backed
implementation. This is the containerized pure-Python + Octave pipeline: faster (the approved
A_E speed-up gives ~8.6×), better tested (69 tests, with a determinism baseline), and it removes
the `Shambhala2` R package from the dependency set entirely.

**Default P and Q ship with the package.** `combobatch/data/calibration/` holds
`P0_standard.csv.gz` (39 × 11,768) and `Q0_standard.csv.gz` (100 × 11,887) — 10.1 MiB gzipped
together, both samples-as-rows in raw linear scale, loaded through `importlib.resources`. Users
override with `params: {P: <path-or-s3-uri>, Q: <path-or-s3-uri>}`, so the 18 benchmark P/Q
variants remain reproducible by pointing at the other calibration files without any of them
being baked into the tool.

**Four source bugs fixed rather than copied:**

1. **`skip_qn` is not forwarded** in `run_shambhala_job.py:516-529`, so
   `--precompute-qn-reference` applies Python QN *and* Octave `quantilenorm` — double
   quantile normalization. `run_shambhala.py` does forward it. The glue function forwards
   `skip_qn=precompute_qn_reference`.
2. **`fixed_clusters` silently overrides `skip_qn`** in `octave_bridge.py`'s script selection,
   so the *approved* A_E speed-up combination also double-QNs. Fix: add a
   `Shambhala2_piped_preqn_fixed.m` variant, or refuse the combination with a clear error.
3. **`ProgressEvent` field mismatch** in `_run_python_cublock_batch` (`parallel.py:243-249`
   passes `done/total/speed_s_per_sample`; the NamedTuple declares
   `sample_done/sample_total/is_final`) → `TypeError` on any progress-enabled
   `python_cublock` run.
4. **`write_expression` ignores the extension** — `--output foo.tsv.gz` silently writes
   uncompressed CSV. ComboBatch routes all writes through `dataio`.

**Defaults chosen for the dispatcher context:** `n_workers=1` and `disable_progress=True`, because
total Octave processes multiply as outer × inner. `combobatch dispatch` warns if
`n_workers_outer × shambhala_n_workers > cpu_count`.

`python_cublock` is exposed but **documented as experimental and xfailed upstream**
(226 % mean relative difference, ProcessPoolExecutor deadlock).

### 5.7 `combobatch/harmonize.py`

```python
def run_one_combination(cfg, imp_spec, method_spec, *, exp_df=None, ann_df=None,
                        skip_if_exists=False) -> list[dict]
```

Pipeline: load/align → optional `subset_query` → impute (cached per `imp_spec`) → resolve
reference batch → resolve + validate hyperparameters → normalize → for each `post_rm` variant:
optional post-removal → compute inline metrics if requested → write matrix + sidecar.
Memory hygiene (`del` + `gc.collect()` + `_r_gc()`) is preserved from the source.

**Prepared-matrix caching**: as in the source's two-stage design, imputation output is written to
`prepared/` and reused across every method sharing that `(imp, ptag)` — imputation is far more
expensive than most harmonizations, and `missforest` in particular.

### 5.8 `combobatch/metrics/`

`METRIC_GROUP_REGISTRY: dict[str, MetricGroupSpec]` — **newly written**, since none exists upstream:

```python
@dataclass(frozen=True)
class MetricGroupSpec:
    letter: str; name: str; fn: Callable
    sentinel_key_template: str          # e.g. "r2_{batch_col}"
    needs_annotation: bool
    column_roles: tuple[str, ...]       # ("batch",) | ("bio",) | ("batch","bio") | ("cohort",)
    needs_reference: bool               # only L
    needs_panel: bool                   # L, M
    needs_embedding: tuple[str, ...]    # (), ("pca",), ("umap","tsne")
    needs_r: bool                       # only F
    speed: str                          # fast | moderate | slow | very_slow
    default_active: bool
    description: str                    # feeds docs/METRICS.md
```

This single table replaces the hardcoded `_run_group` sequence, the separate
`GROUP_SENTINEL_KEYS` dict in another file, and `GROUPS_NEEDING_REFERENCE` — and makes
`--groups` validation, `--skip-slow`, incremental sentinels, `combobatch list-metrics`, and the
generated `docs/METRICS.md` all derive from one place.

**User-specified metric columns.** `column_roles` declares which *roles* a group needs; the
concrete column names come from `ColumnSpec`. So `--metric-batch-cols RNA_BATCH,PLATFORM_RNA`
makes group A emit `r2_RNA_BATCH` and `r2_PLATFORM_RNA`, and `--metric-bio-cols` does the same
for the biology-preservation groups. Defaults are the single `--batch-col` / `--bio-col`, which
keeps a minimal invocation minimal; the FL four-column behaviour is reproduced by naming all four.
Every metric key is suffixed with the column it was computed on, exactly as upstream, so sentinel
keys stay consistent with whatever columns were configured.

`compute_metrics(exp_df, ann_df, cols, *, groups=None, ref_df=None, panel=None, …) -> dict`
keeps the source's per-group try/except isolation, `error_<LETTER>` reporting, partial flushing,
and `on_group_done` callback. Fixes applied:

- **`compute_group_l`'s argument order** `(exp, ref, ann)` is normalized to `(exp, ann, …, ref=…)`.
- **`--skip-slow` vs `--force-groups F`**: the source's comment claims force wins; it does not,
  because `skip_slow=True` still gates `_run_group("F", …)`. ComboBatch makes force actually win.
- **Group L's reference** is the in-memory pre-harmonization matrix — no `01_raw` download, no
  12.8 GB `ref_cache`, and no risk of a stale cache silently comparing against the wrong baseline.
- **Group N** generalizes from the FL diagnosis literals to `--predict-classes`, defaulting to the
  top-N most frequent levels; the reported metric names carry the actual class count.
- **Thread pinning** (`OMP_NUM_THREADS` etc. `setdefault("1")` before any numeric import) is kept
  — the source proved that pod-level env vars do not reach SSH-launched processes, so this must
  live in the Python entry point.
- Magic numbers currently unreachable (`n_permutations=200`, `MIN_COHORT_N=20`, `N_PCS=10`,
  `XB_MAX_SAMPLES=8000`, DSC's 999 perms, UMAP `n_neighbors=30`) become optional `MetricsSpec`
  fields with the same defaults.

`concat.py` ports `run_metrics_concat.py`: sidecars → `metrics_summary.csv` plus the three
long-format tables, with `_split_nested` / `_drop_remaining_containers` preserved.

### 5.9 CLI (`combobatch/cli/`)

```
combobatch run       --config c.yaml | --exp … --ann … --method … --imputation …
combobatch dispatch  --config c.yaml --imputations … --methods … --n-workers 8
combobatch metrics   --config c.yaml [--from-outputs s3://…/exp/]
combobatch concat    --out-dir ./tables --date-tag 260829
combobatch selftest  [--with-r] [--with-octave]
combobatch list-methods | list-metrics | list-imputers
```

**`dispatch` is the cross-product dispatcher.** Its two list flags name the axes explicitly:

```bash
combobatch dispatch --config base.yaml \
    --imputations strict,knn,softimpute,missforest \
    --methods 01_raw,04_sva,10_mnn,16_fsqn_r,20_shambhala \
    --n-workers 8 --skip-if-exists
# → 4 × 5 = 20 combinations; with --post-removal, 40 output matrices
```

Both accept `all` (expanding from the registries) and both are validated against the registries
before anything launches. `--dry-run` prints the resolved job matrix — every
`(imputation, method, params)` triple with its output key — plus the total output count and an
estimated wall time, without running anything. Other flags: `--memory-limit-gb`, `--timeout-s`,
`--retry-failed`, `--skip-if-exists`, `--run-id`.

`run` — a single `(imputation, method)` combination; also the worker subprocess invoked by
`dispatch`. `metrics` — standalone, satisfying "in a separate script if a user wants to calculate
the metrics later"; it enumerates existing outputs under `{out_root}/exp/` and computes the
requested groups, reusing the incremental sentinel logic so new groups can be added to finished
runs without recomputing the rest.

**Hyperparameters on the command line** (fully worked examples in `docs/HYPERPARAMETERS.md`
and `examples/configs/hyperparameters.yaml`):

```bash
# one method, one parameter
combobatch run … --method 10_mnn --method-params '10_mnn:k=50'

# several methods, several parameters each
combobatch dispatch … --methods 10_mnn,38_harman,31_ruv3prps \
    --method-params '10_mnn:k=50;38_harman:limit=0.05;31_ruv3prps:k_factors=8,min_cell_size=3'

# imputer parameters, including the NA-fraction ceiling
combobatch dispatch … --imputations knn,softimpute,missforest \
    --impute-params 'knn:knn_k=10,max_na_frac=0.3;softimpute:rank_max=30;missforest:ntree=50'

# a sweep — same method, three values, three distinct outputs
for k in 5 20 50; do
  combobatch run … --method 10_mnn --method-params "10_mnn:k=$k"
done
# → exp/strict__10_mnn__k5__post0.tsv.gz
#   exp/strict__10_mnn__k20__post0.tsv.gz
#   exp/strict__10_mnn__k50__post0.tsv.gz
```

Types are coerced from each `MethodSpec.hyperparams` declaration; unknown names are rejected with
a did-you-mean suggestion. YAML remains the recommended path for anything non-trivial.

### 5.10 `docker/Dockerfile` and the registry

Single image named **`combobatch`**, `python:3.11-slim`, built for `linux/amd64`.
Layer order: apt → R → R packages → Octave → pip → package install.

Changes from the source Dockerfile, each traceable to a documented failure:

| Change | Reason |
|---|---|
| **Pin R to 4.5.3 explicitly** | Unpinned `r-base` from `bookworm-cran40` drifts to 4.6.x; rpy2 3.6.x has a C-level ABI incompatibility → segfault in `33_amdbnorm`. **Verified 2026-08-29:** `https://cdn.posit.co/r/debian-12/pkgs/r-4.5.3_1_amd64.deb` returns HTTP 200 (67 MB), so the Posit standalone `.deb` is the primary pin. Fallback if it ever disappears: `apt-get install r-base=4.5.3-1~bookwormcran.0` + `apt-mark hold` — that exact version string is what `bookworm-cran40` serves today. |
| Add `octave octave-statistics` + the `/usr/local/bin/matlab` wrapper | Shambhala2 internally calls `system("matlab …")`; the wrapper strips MATLAB-only flags |
| Add `bapred` (**CRAN**, via `BiocManager::install` so its Bioc deps resolve) + `affy`, `affyPLM`, `Biobase` | Makes `37_fabatch` actually run (§5.5a) |
| **Remove `FAbatch` and R `reComBat`** from install *and* verify vectors | `FAbatch` exists on neither CRAN nor Bioconductor; `reComBat` is a Python package. Both currently trip the script's terminal `stop()` |
| **Remove `exploBATCH`** entirely | `36_explobatch` is dropped from the registry (§2.5); its `fMM` dependency is permanently gone |
| Add `variancePartition`, `BiocParallel`, `matrixStats` | Needed by metric group F; present only in the metrics pod today |
| Keep `libpng-dev zlib1g-dev`, pre-install `quantreg`/`Hmisc`/`lme4` | Without them `Hmisc` fails to compile and BiocManager **silently skips qsmooth** |
| Add `cmake` | Required by the metrics pod's build chain |
| **No Procrustes** | `39_procrustes` is removed from the registry (§2.5) — no clone, no `sys.path` hack, no proprietary dependency anywhere in the image |
| `ENV OMP_NUM_THREADS=1` … + write to `/etc/environment` | Pod env vars do not reach SSH-launched processes |
| `ENTRYPOINT ["combobatch"]`, `CMD ["--help"]` | Real entry point; `k8s` overrides with `sleep infinity` for interactive pods |

Per-package `|| echo "WARN"` for the optional/fragile GitHub installs, with **one** consolidated
verification step at the end that fails only on genuinely required packages.

Expected image size **6–9 GB**. *(Measured 2026-08-30: **5.07 GB** — see
`docker/remaining_selftest_failures_plan_260830.md`.)*

**Registry: GHCR only.** `ghcr.io/nikit357/combobatch`, marked public. The repo already lives
under `github.com/Nikit357`, so pushes reuse GitHub auth; a public package is anonymously
pullable with no Docker Hub rate limit; and the `rnd-sandbox` cluster can pull it directly, so no
mirror and no `imagePullSecrets` are needed. This removes the donors' node-IAM pull dependency
and the admin `ecr:CreateRepository` request entirely — **ComboBatch touches no AWS container
registry at all**, and there is exactly one place an image can be stale.

`scripts/build_and_push_image.sh`:

```bash
IMAGE="ghcr.io/nikit357/combobatch"
VERSION="$(date +%Y%m%d)"

docker buildx build --platform linux/amd64 \
    --tag "${IMAGE}:latest" --tag "${IMAGE}:${VERSION}" \
    --file docker/Dockerfile --push .
```

`--platform linux/amd64` is mandatory: the build host is arm64, the c6a nodes are amd64.
Pushing requires `echo $GITHUB_TOKEN | docker login ghcr.io -u nikit357 --password-stdin` with a
token carrying `write:packages`; the package must be flipped to public once, in the GitHub UI,
after the first push (see §7.13).

### 5.11 `k8s/pod-combobatch.yaml`

Same shape as the existing pods — namespace `rnd-sandbox`, `node-group` toleration + required
nodeAffinity, PVC `danya-nikitin-fl` at `/workspace`, `aws-credentials` Secret at `/root/.aws`,
`dshm` 4 GiB memory emptyDir, SSH for interactive work — with these differences:

- `image: ghcr.io/nikit357/combobatch:latest`, `imagePullPolicy: Always`. **The entire startup
  install script disappears** — startup drops from 30–40 min to an image pull. Because the GHCR
  package is public, **no `imagePullSecrets` is required**.
- SSH pubkey from the **Secret** `danya-nikitin-ssh-pubkey` (as the metrics pod does), not an
  inline ConfigMap with the key committed in the file (as `pod-ssh.yaml` does).
- `nodeSelector: karpenter.k8s.aws/instance-family: c6a` **plus
  `karpenter.sh/capacity-type: on-demand`** (corrected during Phase 7: the `karpenter.k8s.aws/`
  spelling this plan originally carried, inherited from the donor note, matches no node in
  the cluster — only *instance* attributes use the AWS-provider prefix, and a `nodeSelector`
  on a non-existent label leaves the pod `Pending` forever with no error) — documented as needed in
  `implementation-plans-old/nohup_job_failing.md` but never applied. Spot reclamation gives a
  2-minute SIGTERM, which cannot save a 6-hour missForest job.
- Annotation `karpenter.sh/do-not-disrupt: "true"`.
- Env: `OMP_NUM_THREADS=1` and friends, plus `COMBOBATCH_WORK_DIR=/workspace`.
- No ConfigMaps for requirements — the image *is* the environment. This removes the
  documented "keep `requirements.txt` in sync with the ConfigMap" hazard entirely.
- The `aws-credentials` Secret is still mounted, but now only for **data** access (reading an
  `s3://` input, writing an `s3://` output) — never for pulling the image.

`k8s/job-combobatch.yaml` — a `batch/v1` Job for unattended runs, `restartPolicy: Never`,
`backoffLimit: 0`, with the dispatcher as `command`.

**Note:** PVC `danya-nikitin-fl` is `ReadWriteOnce` — only one pod can mount it at a time.

### 5.12 Example dataset — verified and sized

The five cohorts were checked against `comb_ann_unified.csv` (7,238 × 484) and measured against
the real expression matrix on S3. **All five are present, and the subset fits GitHub
comfortably.**

| Accession | Exact `COHORT_LABEL` | n | `RNA_BATCH` | Platform | Source | Diagnoses |
|---|---|---|---|---|---|---|
| GSE64555 | `PUB_DLBCL_GSE64555` | 40 | `GPL570_FFPE_Unknown` | GPL570 microarray | FFPE | DLBCL 40 |
| GSE148070 | `PUB_FL_GSE148070` | 40 | `GPL17586_FF_Unknown` | GPL17586 microarray | FF | FL 40 |
| GSE119234 | `PUB_Health_B_cells_GSE119234` | 20 | `RNASeq_FF_rRNADepletion` | RNA-seq | FF | Memory 10, GC 5, Naive 5 |
| GSE69033 | `GSE69033` (bare, no prefix) | 30 | `GPL20188_FF_Unknown` | GPL20188 microarray | FF | Naive/Centroblast/Centrocyte/Memory/Plasma, 6 each |
| GSE62241 | `PUB_FL_GSE62241` | 28 → **14** | `RNASeq_FF_rRNADepletion` | RNA-seq | FF | FL 24, Bone_marrow_CD19+ 2, Centrocyte 2 |

**Totals: 158 annotation rows → 144 unique biological samples.** Combined diagnoses (deduped):
FL 50, DLBCL 40, Memory 16, Naive 11, Centrocyte 6, Centroblast 6, Plasma 6, GC 5.

All four required properties hold: **FL and DLBCL both present**; **RNA-seq present**
(`RNASeq_FF_rRNADepletion`, 48 rows); **microarray present** (3 distinct platforms); **4 distinct
platforms and batches**, spanning both FF and FFPE. This is a genuinely good demonstration set —
it exercises cross-platform harmonization, preservation-type mixing, and malignant-vs-normal
biology at once.

**Measured sizes** (exact, from range-fetching the 158 rows out of the 1.8 GiB `comb_exp.tsv`,
not extrapolated). Total gene count 21,890; mean 216,024 uncompressed bytes per selected row;
measured gzip ratio 0.426:

| Variant | Shape | Uncompressed | **Gzipped** |
|---|---|---|---|
| All 158 rows, all genes | 158 × 21,890 | 34.28 MB | 14.61 MB |
| **Deduped 144 rows, all genes** | **144 × 21,890** | **32.63 MB** | **13.89 MB** |
| 158 rows, all-NA genes dropped | 158 × 16,013 | 33.30 MB | 14.54 MB |
| 158 rows, 5-cohort gene intersection | 158 × 6,797 | 17.22 MB | 7.70 MB |

**Decision: ship the deduped 144 × 21,890 matrix, ~13.9 MB gzipped, committed to the repo.**
That is 7× under the 100 MB budget and under GitHub's 50 MB soft warning. No Git LFS, no download
script, no `.gitattributes`.

Two data-quality findings that `scripts/build_example_dataset.py` must act on:

1. **GSE62241 is duplicated.** Its 28 rows are the same 14 SRA experiments entered twice — once
   as `SRX730599-PUB_FL_GSE62241` (with full metadata) and once as bare `SRX730599` (metadata all
   NaN). Verified at the expression level: r = 1.000000 across 6,797 shared genes for three
   spot-checked samples. The two copies also **disagree on `Diagnosis_cell_type_unified`** for 4
   samples. The build script drops the bare-ID block and keeps the annotated one.
2. **Gene coverage is bimodal, and this is a feature.** GSE69033 and GSE119234 cover an identical
   6,925-gene subspace; the other three cover ~15.9–16.0k. Only 6,797 genes are present in all
   five cohorts, and 5,877 of the 21,890 columns are entirely NA across these samples. Keeping
   the ragged matrix is deliberate: **the missingness is platform-driven, which is exactly what
   makes this dataset a real test of the imputation axis.** A dense intersection matrix would
   make `strict`, `knn`, `softimpute` and `missforest` produce near-identical results and defeat
   the demo.

`scripts/build_example_dataset.py` (run once by a maintainer, output committed) reads
`comb_ann_unified.csv` → filters to the five `COHORT_LABEL` values → drops the 14 duplicate
bare-ID rows → range-fetches those samples from
`s3://FL_batch_correction/exp/comb_exp.tsv` → writes
`example_exp.tsv.gz` + `example_ann.csv` (a reduced column set, not all 484) → emits
`PROVENANCE.md` with accessions, platforms, counts, the exact filter expression, generation date,
and per-series GEO attribution.

### 5.13 Tests

**`tests/unit/`** — pytest, no network, no R, no Octave. Runs natively on the dev Mac (Python
3.11.16) in < 60 s.

| Module | Covers |
|---|---|
| `test_config.py` | YAML→dataclass, CLI overlay precedence, validation failures (unknown method, missing column, bad group letter) |
| `test_storage.py` | Local backend round-trip; S3 via `moto`; URI dispatch; **403 vs 404 distinction** |
| `test_params.py` | `--method-params` parsing, coercion, unknown-name rejection; **filename encoding** — readability, escaping, the 64-char truncation-with-hash path, and round-tripping |
| `test_subset.py` | `--subset-query` evaluation, unknown-column and zero-row rejection |
| `test_imputation.py` | All 4 imputers on synthetic NA patterns; `max_na_frac` boundaries; residual-NA reporting; **no-op detection**; **no silent strict fallback**; every declared hyperparameter changes the result |
| `test_registry.py` | Every `METHOD_REGISTRY` value is callable with the documented signature; no duplicate keys; every declared hyperparameter exists in the function's signature (`inspect.signature`); **the 5 removed keys are absent** |
| `test_rinterop.py` | `r_arglist` rendering, quoting, `TRUE/FALSE/NULL` mapping, and **rejection of non-scalar input** (the R-injection boundary) |
| `test_methods_all.py` | **All 34 methods** on synthetic data — see below |
| `test_hyperparams.py` | Passing a non-default hyperparameter **changes the output**, per method — the regression guard for the bug that all hyperparameters were unreachable |
| `test_postremoval.py` | Outlier identification with `n_batches` 1 and 2; post-removal failure does not silently emit a post0-identical post1 |
| `test_metrics_groups.py` | Port the source's 21 assertions: R² low/high, group isolation, J ordering, K NA counts, L identity-ρ=1 and per-batch-shift ceiling, M margin, N separability |
| `test_metrics_columns.py` | User-specified metric columns produce the expected metric key names and sentinels |
| `test_metrics_registry.py` | Registry ↔ function agreement; sentinel templates resolve; `--force-groups F` **overrides** `--skip-slow` |
| `test_keys.py` | Output-key grammar round-trips; unparseable keys **raise** rather than being silently skipped |
| `test_docs_in_sync.py` | Generated docs match the registries |

**`test_methods_all.py` covers every method, not just the pure-Python ones.** It is parametrized
over the whole registry and asserts the shared contract for each: returns a `DataFrame`, index
preserved and in order, no all-NaN output, no unexpected gene loss beyond what the NaN guard
reports, and determinism under a fixed seed. Methods whose backend is unavailable in the current
environment are marked SKIP via the `MethodSpec` capability flags rather than silently passing —
so running the same file inside the Docker image exercises all 34 for real, and running it on the
Mac exercises the 13 pure-Python ones and honestly reports the rest as skipped. A companion
`test_methods_all_backends.py` in `tests/integration/` asserts that **inside the image, zero
methods skip** — that is what turns "we have 34 methods" into a tested claim.

**`tests/integration/`** — gated by `@pytest.mark.skipif` on availability, mirroring the source's
`octave_required` guard: the 20 R-backed methods, the R imputers, group F, `37_fabatch` against
real `bapred` (including the ≠2-level SKIP path), and Shambhala on the 10-sample fixture with the
determinism check against a committed baseline.

**`tests/e2e/`** — full CLI against `examples/example_dataset/`:
1. `combobatch run` with one method → output exists, correct shape, sidecar valid strict JSON.
2. `combobatch dispatch` over a 2 × 3 grid → 6 outputs, `metrics_summary.csv` has 6 rows.
3. `--skip-if-exists` on a second run → all jobs cached in < 5 s.
4. Metrics inline vs. standalone `combobatch metrics` → **identical values**.
5. Local backend and S3 backend (`moto`) produce identical layouts.
6. Two runs of the same method with different hyperparameters → **distinct, readable** output
   keys, both recorded in `run_manifest.json`.
7. Each `examples/configs/best_approaches/*.yaml` runs end to end on the example dataset.
8. The imputation axis is non-degenerate on the example dataset — `strict` and `softimpute`
   produce different gene counts and different metric values, confirming the bimodal-coverage
   property §5.12 relies on.

**`combobatch selftest`** replaces `test_mock.py`: runs all 34 methods + 4 imputers on an
80 × 200 synthetic matrix, classifies PASS / SKIP / FAIL with the source's
`_is_rpy2_conversion_error` heuristic (an rpy2 conversion failure disguised as
`NotImplementedError` is a FAIL, not a SKIP), and exits 0 iff no FAIL. This is the image's
acceptance gate.

### 5.14 Documentation

Every subdirectory gets a `README.md` (for users) and a `CLAUDE.md` (for agentic work), per §3.
The `CLAUDE.md` files are deliberately short and non-duplicating: each states what its directory
owns, the invariants that hold inside it, and which registry or contract governs edits there.

Top-level `docs/`:

| File | Contents |
|---|---|
| `HYPERPARAMETERS.md` | **Generated from the registries.** Every method and imputer, each parameter's type, default, bounds/choices, what it does, and the R-side argument it maps to where applicable. Worked CLI and YAML examples for single runs, multi-method dispatch, and sweeps. A "which parameters actually matter" section drawn from the benchmark. |
| `METRICS.md` | **Generated from `METRIC_GROUP_REGISTRY`.** Each group A–N: what it measures, the exact metric keys emitted, which annotation columns it consumes and how to change them, cost, whether it needs a reference matrix or a marker panel, its polarity (higher-is-better vs lower-is-better), and its citation. Includes the Group L ceiling caveat (§7.7). |
| `METHODS.md` | The 34 methods with citations and harshness tiers; **the five removed methods and why**, including the Procrustes licensing rationale and the technical notes gathered while evaluating it; the RNA-seq-only guard; the reference-batch resolution rule; the `20_shambhala` implementation substitution. |
| `OUTPUT_LAYOUT.md` | Key grammar, `ptag` encoding rules, `run_manifest.json` schema, sidecar schema, status values. |
| `BEST_APPROACHES.md` | §5.15 — the benchmark's recommendations, with the caveats that qualify them. |

A CI check asserts that `HYPERPARAMETERS.md`, `METRICS.md` and `METHODS.md` match what the
registries currently declare, so the docs cannot silently go stale — the exact failure mode that
produced the donors' four contradictory "39 methods / 31 methods / 12 strategies / 14 strategies"
counts.

### 5.15 `--subset-query` and the benchmark's best approaches

**`--subset-query`** takes a pandas `DataFrame.query` string evaluated against the user's own
annotation columns. No FL-specific strategies ship as code. It is validated before any download:
unknown column → error naming the available columns; zero matching rows → error.

The FL strategies are reproducible as one-liners, which `docs/BEST_APPROACHES.md` documents as
worked examples of the flag:

| FL strategy | `--subset-query` equivalent | n |
|---|---|---|
| `C_rnaseq_only` | `RNA_BATCH.str.startswith('RNASeq')` | 2,243 |
| `J_ff_only` | `RNA_BATCH.str.contains('_FF_')` | 3,167 |
| `K_ffpe_only` | `RNA_BATCH.str.contains('_FFPE_')` | 3,000 |
| `F_microarray_only` | `not RNA_BATCH.str.startswith('RNASeq')` | 4,931 |
| `G_affymetrix_only` | `RNA_BATCH.str.startswith('GPL570')` | 2,801 |
| `D_malignant_only` | `Diagnosis_cell_type_unified not in @NORMAL_GROUPS` | 6,285 |
| `S0_no_removal` | `Diagnosis_cell_type_unified not in @RARE_GROUPS` | 7,174 |

Three of the 14 are **not** expressible as a static query and are documented as such:
`H_affymetrix_extended` is a hand-curated 16-label allow-list (it deliberately excludes
`GPL14951_FFPE_Unknown`, which is Illumina, not Affymetrix); `E1/E2/E3` are data-dependent
iterative PCA-outlier removal; `I_rare_batches_removed` is a data-dependent size threshold.
A user can express the first as an explicit `in [...]` list; the others correspond to
`--post-removal` applied iteratively.

**A caveat worth stating in the docs:** `J_ff_only` and `K_ffpe_only` are *not* complementary —
3,167 + 3,000 ≠ 7,174, because batches whose middle token is `Unknown` (notably
`GPL570_Unknown_Unknown`, 1,030 samples) fall into neither.

**`examples/configs/best_approaches/`** ships the benchmark's recommendations as runnable configs,
one per decision-tree branch, each carrying the `--subset-query` and hyperparameters that
reproduce it. From the manuscript's Figure 11 / Table 1 and the curated Best-15:

| Config | Scenario | Method | Imputation | Subset |
|---|---|---|---|---|
| `ff_only_mnn.yaml` | FF only | `10_mnn` (alt `16_fsqn_r`) | strict | `RNA_BATCH.str.contains('_FF_')` |
| `ffpe_only_sva.yaml` | FFPE only ★ best local mixing in the benchmark | `04_sva` | **softimpute or knn** | `RNA_BATCH.str.contains('_FFPE_')` |
| `rnaseq_only_sva.yaml` | RNA-seq only | `04_sva` (2nd `16_fsqn_r`) | **softimpute or knn** | `RNA_BATCH.str.startswith('RNASeq')` |
| `rnaseq_plus_illumina_fsqn.yaml` | RNA-seq + Illumina arrays | `16_fsqn_r`, alt `33_amdbnorm` / `13_fsmvn` + post-removal | strict | none |
| `multiplatform_mnn.yaml` | RNA-seq + mixed arrays | `10_mnn` + post-removal | strict | none |

The curated **Best-15** (5 methods × 6 strategies, all `post_rm=False`) is tabulated in
`docs/BEST_APPROACHES.md` with each run's original benchmark `run_id`, so a reader can trace any
recommendation back to the published `metrics_comprehensive.csv`.

`docs/BEST_APPROACHES.md` must carry the caveats alongside the recommendations, or it becomes
misleading:

1. **Global metrics are necessary but not sufficient.** Good PCReg/R² coexists with locally
   separated batches; local and global metric families are near-orthogonal (median Spearman
   r ≈ 0.009).
2. **Local metrics (kBET, iLISI) universally fail cross-platform** — a limitation of the problem,
   not of any method. ~60 % of all tested approaches score "bad".
3. **Three rankings disagree and the article has not settled which it reports**: equal-weight
   composite (the figure default, favours the clustermap-best group at 0.602), literature-weighted
   composite (top-10 is entirely `C_rnaseq_only__post1`, led by a *Shambhala* run at 0.732), and
   Borda count (top-10 is entirely `16_fsqn_r` on A/E1–E3). The three top-10 sets have **zero
   overlap** with each other and with the Best-15. Ship the Best-15 as "the manuscript's curated
   selection", not as "the optimum".
4. **At least one Best-15 approach corrupts biology**: marker-gene correlation before/after
   harmonization was ~1 for linear methods but **as low as 0.05** for one of the
   clustermap-best approaches. The affected config carries an explicit warning.
5. **No method dominates** — 26 of 31 methods are Pareto-optimal somewhere in the
   batch-mixing × biology-preservation plane.
6. **Softimpute ≈ 10 pp better than KNN** by PCR; prefer it where both are listed.
7. **The factor-importance hierarchy is version-dependent.** The June-2026 note says
   *strategy > method > post-removal > imputation*; the current manuscript's η² analysis says
   **method 36.3 % > strategy 25.6 % > imputation 1.6 % > post-removal 0.2 %**. Quote the
   manuscript numbers.
8. **The Best-15 is a manual curation**, combining visual inspection with composite scoring — not
   the output of an automated rule, and Figure 3 itself is still unimplemented as a publication
   figure.

---

## 6. Source bugs to fix, not copy

Consolidated for the implementation phase. Each was verified in the source.

| # | Location | Bug | Fix |
|---|---|---|---|
| 1 | `run_prep_job.py:152-159` | Imputation failure falls back to strict but keeps the `__knn__` label | Fail loudly |
| 2 | `run_norm_job.py:161-198` | Post-removal exception still writes the `post1` key (identical to post0) | Mark `status: "post_removal_failed"` |
| 3 | `bench_shared.py:s3_exists` | Swallows all `ClientError`s — 403 looks like 404 | Re-raise non-404 |
| 4 | `run_shambhala_job.py:516` | `skip_qn` not forwarded → double quantile normalization | Forward it |
| 5 | `octave_bridge.py` | `fixed_clusters` overrides `skip_qn` → the approved A_E combo double-QNs | Add a combined `.m` variant or refuse |
| 6 | `parallel.py:243-249` | `ProgressEvent` constructed with wrong field names → `TypeError` | Match the NamedTuple |
| 7 | `io_utils.write_expression` | Ignores the file extension; `.tsv.gz` silently written as plain CSV | Route through `dataio` |
| 8 | `run_metrics_job.py` | Comment claims `--force-groups F` beats `--skip-slow`; it does not | Make force win |
| 9 | `run_norm_job.py` | `json.dump` emits bare `NaN` — invalid strict JSON | NaN-safe encoder everywhere |
| 10 | `install_r_packages.R` | Installs and verifies `FAbatch` (exists nowhere) and R `reComBat` (a Python package) → terminal `stop()` fails the build | Remove both; add `bapred` from CRAN |
| 11 | `run_metrics_parallel.py:_list_exp_files` | Keys not splitting into exactly 4 parts are **silently dropped** | Fail loudly on unparseable keys |
| 12 | `bench_shared.py` / dispatchers | `ALL_METHODS`/`ALL_STRATEGIES`/`ALL_IMPUTATION` retyped by hand in 2–3 files | Derive from the registries |
| 13 | `compute_batch_metrics.py` | Importing the module requires `marker_gene_annotation.csv` to exist, even when L/M are not requested | Lazy load |
| 14 | `bench_shared.py:normalize_fabatch` | `fabatch(xtr=, ytr=, type=)` + `$adj.data` match no `bapred` API; the `tryCatch` then writes the *uncorrected* matrix, which is scored as if corrected — a silent no-op | `fabatch(x, y, batch)` + `$xadj`; re-raise instead of falling back (§5.5a) |

Two further defects were found in `normalize_procrustes` — the docstring advertises a
`coeffs_kit="V7"` for which no coefficient file exists, and it describes the model as per-gene
slope+intercept when it is actually multivariate. Both are **moot for ComboBatch** since
`39_procrustes` is not included (§2.5), but they are recorded in `docs/METHODS.md` so the donor
repo can fix them independently.

---

## 7. Side effects and caveats

1. **Output keys change shape.** Encoding hyperparameters into filenames means ComboBatch keys
   are not byte-compatible with the existing `FL_batch_correction/` S3 layout. Intentional — the
   tool writes to a user-specified `out_root` and never to the dissertation's prefix. **No
   existing dissertation artefact is read, moved, or overwritten by anything in this plan.**

2. **R 4.5.3 pinning — resolved.** Verified 2026-08-29: the Posit Debian-12 build
   (`cdn.posit.co/r/debian-12/pkgs/r-4.5.3_1_amd64.deb`) exists and is 67 MB, so the primary
   pin is available and this is no longer the plan's highest-risk step. `bookworm-cran40`
   currently serves `r-base 4.5.3-1~bookwormcran.0`, which is the exact string for the
   `apt-mark hold` fallback — but it *will* drift to 4.6.x, so the pin is still mandatory.
   Whichever path is used, the built image must be checked with
   `Rscript -e "cat(R.version$version.string)"` before anything else.

3. **Image size 6–9 GB.** Pulls are slow on a cold node. A `combobatch-slim` variant (no Octave,
   no R) is **deferred past v0.1**.

4. **`pandas<2.0` is a long-term liability.** It is EOL and will eventually block a dependency
   upgrade. The escape hatch is to drop `30_recombat` (the sole `<2.0` constraint) and move to
   pandas 2.x + Python 3.12.

5. **Nested parallelism with Shambhala.** Outer workers × inner `n_workers` multiply. Default
   inner to 1 and warn when the product exceeds the CPU count.

6. **Group L semantics differ from the source** — better, but different. ComboBatch compares
   against the in-memory pre-harmonization matrix; the source compares against `01_raw__post0`
   from S3. Numbers will not be identical when the source's `01_raw` differs from the user's
   actual input.

7. **Group L has a known ceiling.** Any harmonizer applying a per-batch *monotone* transform
   (centering, scaling, z-score, median scaling, and to first order ComBat / limma) leaves
   within-cohort Spearman correlations unchanged, so Group L sits at ~1.0 for that whole method
   class by construction. It is a "nothing broke" guard rail; Group M's margin carries the
   discriminative signal. Stated in `docs/METRICS.md`.

8. **`--platform linux/amd64` is mandatory** on every build from this arm64 Mac. Python 3.11.16
   is now installed locally, so unit tests run natively; only R- and Octave-backed tests need the
   container.

9. **The tool ships 34 of the article's 39 methods.** Four cannot run at all and one is
   licence-restricted (§2.5). Anyone reproducing the article from ComboBatch will therefore not
   be able to reproduce the `39_procrustes` rows — which is moot in practice, since the published
   `metrics.csv` has **zero** `status=ok` rows for it (the pod's clone step was never actually
   present). `docs/METHODS.md` states this explicitly so the gap is not mistaken for an omission.

10. **`20_shambhala` results will not be bit-identical to the published benchmark.** The key now
    points at the containerized pure-Python + Octave implementation instead of the `Shambhala2`
    R package. This is the intended replacement — faster and better tested — but any claim of
    exact reproducibility against the article's Shambhala rows must be qualified.
    `docs/METHODS.md` records both implementations and the substitution.

11. **`37_fabatch` may still SKIP on many datasets** — `bapred::fabatch` structurally requires a
    binary target, so a dataset whose bio column has 3+ levels cannot run it. The SKIP is now
    honest rather than a silent no-op, but the method will not be universally available.

12. **The example dataset ships real GEO-derived expression values.** Attribution per series is
    in `PROVENANCE.md`. If any series turns out to carry a restrictive licence, substitute a
    backup cohort — the build script is parametrized by the cohort list.

13. **GHCR publication is a manual one-time step.** A newly pushed GHCR package is **private by
    default**; it must be flipped to public in the GitHub UI after the first push, or the pod
    will fail to pull with a 403 that looks like a missing image.

14. **The dissertation repo's `k8s/README.md` and `k8s/CLAUDE.md` describe a pod that does not
    exist on disk** (a `/app/venv`, a Posit R install, a Procrustes clone — none present in
    `pod-ssh.yaml`). Do not use them as a specification; the manifests were read directly.

---

## 8. Verification commands

```bash
# ── Local, no R/Octave (Python 3.11.16 installed via Homebrew) ───────────────
# Note: `python3` is still the system 3.9.6 — always invoke python3.11 explicitly.
cd /Users/user890/Desktop/ComboBatch
python3.11 -m venv .venv && source .venv/bin/activate
python --version                          # → Python 3.11.16 inside the venv
pip install -e ".[dev]"
pytest tests/unit -v                      # must be green with no network
combobatch list-methods                   # 34 rows, with hyperparameters
combobatch list-metrics                   # 14 rows A–N
combobatch list-imputers                  # 4 rows incl. strict

# docs match the registries (CI enforces this too)
pytest tests/unit/test_docs_in_sync.py -v

# the removed methods must be absent, not merely skipped
combobatch list-methods | grep -E '24_peer|32_deepmnn|35_dasc|36_explobatch|39_procrustes' \
    && echo "FAIL: a removed method is still registered" || echo "OK: 5 methods removed"

# ── Image build (amd64 from an arm64 Mac) ────────────────────────────────────
docker buildx build --platform linux/amd64 -t combobatch:test -f docker/Dockerfile .
docker run --rm --platform linux/amd64 combobatch:test \
    python -c "import rpy2.robjects as ro; print(ro.r('R.version.string')[0])"
# MUST print R version 4.5.x — 4.6.x means the pin failed and rpy2 will segfault

docker run --rm --platform linux/amd64 combobatch:test combobatch selftest
# exit 0 = all 34 methods PASS or SKIP

docker run --rm --platform linux/amd64 combobatch:test \
    Rscript -e "library(variancePartition); library(bapred); cat('R deps OK\n')"
docker run --rm --platform linux/amd64 combobatch:test octave --version

# every method must actually run inside the image — zero skips
docker run --rm --platform linux/amd64 combobatch:test \
    pytest tests/integration/test_methods_all_backends.py -v

# ── Publish to GHCR (public, the only registry) ──────────────────────────────
echo "$GITHUB_TOKEN" | docker login ghcr.io -u nikit357 --password-stdin
bash scripts/build_and_push_image.sh
# then flip the package to public once in the GitHub UI, and verify anonymously:
docker logout ghcr.io && docker pull ghcr.io/nikit357/combobatch:latest

# ── End-to-end on the example dataset ────────────────────────────────────────
docker run --rm --platform linux/amd64 -v "$PWD:/data" combobatch:test \
    combobatch dispatch --config /data/examples/configs/full_crossproduct.yaml --dry-run
# prints the 4 × 4 job matrix and the resulting output keys, runs nothing

docker run --rm --platform linux/amd64 -v "$PWD:/data" combobatch:test \
    combobatch dispatch --config /data/examples/configs/full_crossproduct.yaml --n-workers 2
ls out/exp/ out/metrics/ && head -3 out/metrics_summary.csv
python -c "import json; json.load(open('out/run_manifest.json'))"   # valid strict JSON

# ── Hyperparameters reach the method AND appear in the filename ──────────────
combobatch run … --method 10_mnn --method-params '10_mnn:k=5'
combobatch run … --method 10_mnn --method-params '10_mnn:k=50'
ls out/exp/ | grep 10_mnn
#   strict__10_mnn__k5__post0.tsv.gz
#   strict__10_mnn__k50__post0.tsv.gz
python - <<'PY'
import pandas as pd
a = pd.read_csv("out/exp/strict__10_mnn__k5__post0.tsv.gz",  sep="\t", index_col=0)
b = pd.read_csv("out/exp/strict__10_mnn__k50__post0.tsv.gz", sep="\t", index_col=0)
assert not a.equals(b), "hyperparameters had no effect — the forwarding is broken"
print("OK: hyperparameters reach the method and are visible in the key")
PY

# imputer hyperparameters, incl. the NA ceiling
combobatch run … --imputation knn --impute-params 'knn:knn_k=10,max_na_frac=0.3'

# ── Metrics: user-specified columns, computed later, separately ──────────────
combobatch metrics --config examples/configs/minimal.yaml --groups A,B,E,J,K \
    --metric-batch-cols RNA_BATCH,PLATFORM_RNA --metric-bio-cols Diagnosis_cell_type_unified
python -c "import json;d=json.load(open('out/metrics/strict__04_sva__post0_metrics.json'));\
print([k for k in d if k.startswith('r2_')])"   # → ['r2_RNA_BATCH', 'r2_PLATFORM_RNA']
combobatch concat --out-dir ./tables --date-tag 260829

# ── subset-query and a best-approaches config ────────────────────────────────
combobatch run --config examples/configs/best_approaches/ffpe_only_sva.yaml --dry-run
combobatch dispatch --config examples/configs/minimal.yaml \
    --subset-query "RNA_BATCH.str.startswith('RNASeq')" --dry-run

# ── K8s (public image — no imagePullSecrets) ─────────────────────────────────
kubectl apply -f k8s/pod-combobatch.yaml -n rnd-sandbox
kubectl get pod danya-nikitin-combobatch -n rnd-sandbox -w   # Running in ~2 min, not ~40
kubectl exec -it danya-nikitin-combobatch -n rnd-sandbox -- combobatch selftest
```

---

## 9. Decisions on record

All previously open questions are settled:

| Question | Decision |
|---|---|
| Which 5 GEO cohorts | GSE64555, GSE148070, GSE119234, GSE69033, GSE62241 — verified, 144 unique samples, ~13.9 MB gzipped, committed |
| Repository | New public GitHub repo `Nikit357/ComboBatch` |
| Licence | MIT |
| Container registry | **GHCR only**, `ghcr.io/nikit357/combobatch`, public. No ECR, no mirror. |
| FL filter strategies | Not shipped as code — `--subset-query` over the user's own annotation columns |
| `24_peer_k10`, `32_deepmnn`, `35_dasc` | Removed |
| `36_explobatch` | Removed |
| `37_fabatch` | Repaired; honest SKIP when the target column has ≠2 levels (no auto-binarizing) |
| `39_procrustes` | **Removed** — proprietary non-commercial licence, incompatible with an MIT public tool |
| `20_shambhala` | Replaced with the containerized pure-Python + Octave implementation |
| `combobatch-slim` | Deferred past v0.1 |
| Local Python | 3.11.16 installed via Homebrew; scripts invoke `python3.11` explicitly |

---

## 10. TODO

### Phase 0 — scaffold ✅ complete (2026-08-29)
- [x] `git init`; create the §3 directory tree with `__init__.py`, `README.md` and `CLAUDE.md` in every subdirectory
- [x] `pyproject.toml` — metadata, `requires-python = ">=3.11,<3.12"`, MIT, entry point, **`package-data` for `octave/*.m` and `calibration/*.csv.gz`**
- [x] `requirements.txt` — unified pins with a justifying comment per pin
- [x] `LICENSE` (MIT), `.dockerignore`, `.gitignore`
- [x] `python3.11 -m venv .venv` and confirm `pip install -e ".[dev]"` resolves — installed cleanly; pandas 1.5.3, numpy 1.26.4, scikit-learn 1.5.2 on Python 3.11.16, all pins holding
- [x] **Verify the Posit R 4.5.3 Debian-12 `.deb` URL exists** — **confirmed present** (HTTP 200, 67,072,216 bytes). Fallback version string also captured; see §5.10.
- [x] Minimal `combobatch/cli/main.py` so the console script resolves — `combobatch --version` and `--help` work; every subcommand is a stub exiting 2 and naming its phase

### Phase 1 — foundation ✅ complete (2026-08-29)
- [x] `combobatch/logging_utils.py` — thread pinning before numeric imports; timestamped logs. Verified the module imports no numeric library, so `pin_threads()` can still take effect.
- [x] `combobatch/memory.py` — port cgroup-aware `free_gb`, `wait_for_memory`, `check_memory`. Also guards against the cgroup-v1 "unlimited" sentinel near 2^63, which would otherwise read as vast free memory.
- [x] `combobatch/storage.py` — ABC + Local + S3; **403/404 conflation fixed**; one client per backend
- [x] `combobatch/dataio.py` — readers/writers, extension-driven compression, index alignment
- [x] `combobatch/params.py` — parsing, coercion, **readable filename encoding + 64-char truncation-with-hash**, key grammar (**fails loudly**)
- [x] `combobatch/subset.py` — `--subset-query` evaluation and pre-flight validation
- [x] `combobatch/config.py` — dataclasses (incl. `PostRemovalSpec`, metric column roles), YAML loader, CLI overlay, `validate()` before expensive work
- [x] `tests/unit/test_storage.py`, `test_config.py`, `test_params.py`, `test_keys.py`, `test_subset.py` — **163 tests, all passing**, no network and no R
- [x] `tests/unit/test_dataio.py` *(added beyond the plan's list: `dataio` implements bug #7, and an untested bug fix is not a fixed bug)*

**Naming deviation, deliberate.** §5.2 named the config-side selection dataclasses
`ImputationSpec` / `MethodSpec`, but §5.5 gives the *registry* entry the same name
`MethodSpec`. Two different classes sharing a name across modules is exactly the kind of
confusion this plan exists to remove, so the config-side ones are **`ImputationSelection`**
and **`MethodSelection`** — the user's *choice* of a method — leaving `MethodSpec` free for
the registry entry that *describes* a method in Phase 2.

### Phase 2 — imputation and methods ✅ complete (2026-08-29)
- [x] `combobatch/imputation.py` — **all 4 imputers as registry entries incl. `strict`**; full hyperparameter surface (softImpute's `rank.max`/`lambda` and missForest's `ntree`/`maxiter`/`mtry` were not in the donor's Python signature at all); **no silent strict fallback**; no-op detection
- [x] `combobatch/methods/rinterop.py` — converters, `r_gc`, **`r_arglist` (the R-injection boundary)**, availability probes. Importable with no rpy2 present, which is what lets the unit suite run on the Mac.
- [x] `combobatch/methods/base.py` — `MethodSpec` + `HyperParam`, plus `missing_backends()` so an absent backend is an honest skip
- [x] `combobatch/methods/python_methods.py` — the **13** pure-Python methods (no Procrustes)
- [x] `combobatch/methods/r_methods.py` — the 20 R-backed methods; **every hardcoded R-side argument now a declared `HyperParam`** (`polyFit(dis, 9)`, `genDistData(…, 500)`, `log_target=FALSE`, `RUVIII(k=5)`, `cpm(prior.count=)`, `huge.npn(npn.func=)`, `ARSyNseq(norm=)`, and RUVg's housekeeping-gene list); `normalize_dwd` keeps the removed-`penalty` fix
- [x] **`37_fabatch` repair** — `fabatch(x, y, batch)` + `$xadj`; `target_col` hyperparameter; honest SKIP on ≠2 levels; **the `tryCatch` that wrote uncorrected data is gone**
- [x] `resolve_reference_batch()` — replaces the dataset-specific default in 9 methods; falls back to the largest batch and logs the choice
- [x] `combobatch/methods/__init__.py` — `METHOD_REGISTRY`, **33 entries** (`20_shambhala` joined in Phase 3 → 34), full hyperparameter declarations, plus `EXCLUDED_METHODS` recording why each of the 5 is absent
- [x] `combobatch/postremoval.py` — `identify_outlier_batches` with **`n_batches`**; **the silent post1==post0 failure now raises `PostRemovalError`**
- [x] `tests/unit/test_imputation.py`, `test_registry.py` (**asserts the 5 removed keys are absent**), `test_rinterop.py`, `test_methods_all.py`, `test_hyperparams.py`, `test_postremoval.py` — **470 passing, 150 skipped**
- [x] `tests/conftest.py` — synthetic fixtures *(pulled forward from Phase 9; the Phase 2 tests need them)*
- [x] `tests/integration/test_methods_all_backends.py` — *(pulled forward from §5.13: asserts zero methods skip inside the image, which is what makes the method count a tested claim rather than an assertion)*

**Skip accounting.** On the dev Mac 8 methods run natively (`01_raw`, `02_median_scaling`,
`13_fsmvn`, `15_fsqn_py`, `17_quantile`, `18_rank`, `25_angel`, `26_xpn`); the other 25 skip
with a named missing backend — 150 skips = 25 methods × 6 contract tests, exactly as intended.

### Phase 3 — Shambhala ✅ complete (2026-08-29)
- [x] Vendor `shambhala/*.py` + `octave/*.m` into `combobatch/vendor/shambhala/` — **`io_utils.py` deliberately not vendored** (its writer ignored the file extension, §6 bug 7); everything routes through `dataio`. `octave/` now lives *inside* the package, fixing the upstream default that only resolved under `pip install -e .`
- [x] Copy `P0_standard.csv.gz` + `Q0_standard.csv.gz` (10.1 MiB) into `combobatch/data/calibration/` — shapes verified as 39 × 11,768 and 100 × 11,887, both raw linear
- [x] Write the single `shambhala_harmonize()` glue function (replacing 4 drifted copies), defaulting P/Q to the shipped files
- [x] **Fix:** forward `skip_qn`; resolve the `fixed_clusters`/`skip_qn` override (one `_PIPELINE_SCRIPTS` table + the previously missing `Shambhala2_piped_preqn_fixed.m`); fix `ProgressEvent`
- [x] Register as `20_shambhala`, **replacing the R-package implementation**; default `n_workers=1`, `disable_progress=True` — registry now **34 methods**, 16 declared hyperparameters
- [x] Nested-parallelism warning — `warn_nested_parallelism()` written and unit-tested here; *the dispatcher call site is wired in Phase 4*
- [x] `tests/integration/test_shambhala.py` — Octave-gated. **Deviation:** determinism is asserted seed-to-seed (same seed twice) rather than against a committed numeric baseline, because a baseline can only be generated where Octave exists and one pinned from a different Octave build would fail for the wrong reason. Generate it in the image in Phase 6 if a stronger guarantee is wanted.
- [x] `tests/unit/test_shambhala_glue.py` + `test_shambhala_vendor.py` — regression tests for all four bugs above, none of which need Octave. **510 passing, 156 skipped**

**Two defects found and fixed beyond the plan.** Neither was in §6 because both are only
reachable through options ComboBatch newly exposes:

1. **`q_pseudocount=0` crashed with a `KeyError`.** That setting excludes zero-count genes
   from the Q statistics, the rescale step then drops them, and upstream's gene restoration
   indexes by the *original* gene list. The glue now records them as unusable so they come
   back as NaN, with the count logged.
2. **No scale detection anywhere.** Log-scale input is silently double-logged into
   plausible-but-wrong numbers. `_warn_if_log_scale()` warns below a maximum of 30. A
   warning rather than an error: it is a heuristic over user data.

`qnorm` was added to the `dev` extra, so the `skip_qn` regression test runs unconditionally
instead of self-skipping.

### Phase 4 — harmonize + dispatcher ✅ complete (2026-08-29)
- [x] `combobatch/harmonize.py` — `run_one_combination()`, prepared-matrix caching, memory hygiene. **The post-removal failure now yields `status: post_removal_failed` and writes nothing**, rather than a `post1` file identical to `post0`
- [x] `combobatch/cli/run_cmd.py` — also the worker `dispatch` launches, and the owner of the argument surface both commands share, so the two cannot disagree about what a configuration means
- [x] `combobatch/cli/dispatch_cmd.py` — **`--imputations` / `--methods` cross-product lists (accepting `all`)**, ThreadPoolExecutor + subprocess, memory guard, timeouts, failed-jobs log, `--dry-run` job matrix, `--skip-if-exists`, `--retry-failed`, SIGTERM handler
- [x] `combobatch/cli/main.py` + `list-methods` / `list-imputers` (`list-metrics` still names Phase 5), plus `combobatch/__main__.py` so workers launch as `sys.executable -m combobatch`, which resolves where a console script may not
- [x] `run_manifest.json` writer — full resolved params per output, merged rather than replaced across invocations
- [x] NaN-safe JSON in `dataio` (`json_ready` / `dumps_json` / `write_json` / `read_json`) — plan §6 bug 9 called for it "everywhere" without naming a home; `dataio` already owns serialization
- [x] `tests/unit/test_harmonize.py`, `test_cli.py`, `test_dispatch.py` — **586 passing, 165 skipped**

**Two design points the plan left open, resolved here.**

1. **One tag slot, two parameter sources.** §4.2 says the tag covers "both the method and
   the imputer", but the grammar has a single slot. `combined_param_tag()` merges the two
   non-default dicts and **raises** if a name is tuned on both sides, rather than letting
   one silently overwrite the other. No real pair collides today; the guard is tested with
   a synthetic one.
2. **`max_na_frac` had a wrong shared default.** §4 gave `ImputationSpec.max_na_frac` a
   literal `0.20` while §5.4 fixes `strict` at `0.0` — so a plain `--imputation strict`
   would have run with a 0.20 ceiling and not been strict at all. The field is now
   `None` = "whatever this imputer declares". This is a Phase 1 correction, made here
   because Phase 4 is what first resolves the value.

**Three defects found by the new tests.** A worker validated its configuration *without*
the registries, so an unknown method key surfaced as a bare `KeyError` after the input had
already been downloaded — the exact "validate before expensive work" failure §5.2 exists to
prevent. `--method-params` names were checked only inside the job, so a typo could reach an
output filename; they are now rejected at parse time with a did-you-mean. And the memory
guard could wait forever in silence on a machine that never has the threshold free, so it
is now bounded by `--memory-wait-s` (default 600) and launches anyway on expiry, leaving
the worker's own guard to make the honest call.

**Known and bounded:** the `prepared/` cache is advisory, not a lock. Workers that start
together all miss it and all impute, writing identical content — at most one redundant
imputation per concurrent worker. Correctness does not depend on who wins.

### Phase 5 — metrics ✅ complete (2026-08-29)
- [x] `combobatch/metrics/embeddings.py` — PCA/UMAP/tSNE helpers, every buried constant now a parameter
- [x] `combobatch/metrics/groups.py` — all 14 groups ported; **`compute_group_l`'s reference is now a keyword** so it can never be transposed with the annotation silently; **every hardcoded column literal replaced by `MetricColumns`**, resolved from the user's `ColumnSpec`
- [x] Generalize Group N's class labels → the classes the user names, else the N most frequent levels of the configured column; the result reports `pv_lobo{n}_n_classes_used`, the count actually used
- [x] `combobatch/metrics/panels.py` — **lazy** optional loader; a user CSV, not a shipped dataset-specific table. No panel means every shared gene, which is the honest generic default
- [x] `combobatch/metrics/__init__.py` — `METRIC_GROUP_REGISTRY`, **14 entries**, with `column_roles`, sentinel templates, speeds and backend declarations
- [x] `combobatch/metrics/runner.py` — isolation, partial flush, `on_group_done`, sentinels resolved against the configured columns; **`--force-groups` beats `--skip-slow`**, verified end to end
- [x] Group L takes the in-memory reference — no baseline download, no `ref_cache`
- [x] `combobatch/metrics/concat.py` + `cli/concat_cmd.py` — wide summary plus the long tables; **an unparseable sidecar name is an error, not a silent drop**
- [x] Wire `--metrics` into both `run` and `dispatch`; standalone `cli/metrics_cmd.py` with `--dry-run`, `--force-groups` and incremental sentinels
- [x] `tests/unit/test_metrics_groups.py`, `test_metrics_columns.py`, `test_metrics_registry.py`, `test_metrics_runner.py`, `test_metrics_concat.py` — **698 passing, 165 skipped**

**Three registries, all now real.** `METRIC_GROUP_REGISTRY` replaces what upstream spread
over three files — a hardcoded `_run_group` sequence, a `GROUP_SENTINEL_KEYS` dict in a
*different module*, and a third constant for "needs the reference". `test_metrics_registry.py`
cross-checks the declarations against the functions by `inspect.signature`, and asserts
`config.DEFAULT_METRIC_GROUPS` still equals the registry's own default set rather than
drifting into a second source of truth.

**The magic numbers are now `MetricsSpec` fields**, as §5.8 required: `n_pcs`, `min_test_n`,
`min_cohort_n`, `wm_permutations`, `wm_max_samples`, `xb_max_samples`, `dsc_permutations`,
`asw_max_samples`, `ks_genes`, `collect_detail`. Defaults reproduce the donor exactly.

**Two behaviours worth knowing.** `--skip-slow` now excludes every group the registry marks
`very_slow` (F, I and N), not just F as upstream did — and it **logs each exclusion by name
with the flag that undoes it**, because a group the user explicitly asked for must never
vanish without a word. A group whose backend is absent is an honest named skip
(`error_F: "skipped: needs R (Rscript)"`), matching the method registry's convention rather
than surfacing as a failed computation.

**Group L's ceiling is now tested, not just documented:** `test_metrics_groups.py` asserts
that a monotone per-gene transform scores exactly 1.0, which is precisely why §7.7 says the
group is a guard rail and group M's margin carries the discriminative signal.

### Phase 6 — Docker and GHCR ✅ complete (2026-08-29, except the two live steps)
- [x] `docker/install_r_packages.R` — **dropped `FAbatch`, R `reComBat`, `exploBATCH`; added `bapred` (CRAN, via `BiocManager::install` so its Bioc deps resolve) + `affy`/`affyPLM`/`Biobase`, `variancePartition`, `BiocParallel`, `matrixStats`**; pre-installs `quantreg`/`Hmisc`/`lme4`. **Also dropped `DASC`** — the donor installed it, but `35_dasc` is not a registered method, so it would have been a dependency nothing could reach.
- [x] `docker/Dockerfile` — `python:3.11-slim-bookworm`, **R 4.5.3 pinned**, Octave + `matlab` wrapper, `cmake`, thread env vars, **no Procrustes**
- [x] `docker/entrypoint.sh`, `docker/README.md`
- [x] `scripts/build_and_push_image.sh` — `buildx --platform linux/amd64`, **GHCR tags only**
- [x] `tests/unit/test_docker_assets.py` — **new, not in the original checklist.** The R script is the only place a registry is restated in another language, and a drift there does not fail the build: the method simply SKIPs at run time and disappears from the benchmark. This asserts the script's `required` vector equals the registry union exactly, in both directions, and pins the Dockerfile invariants.
- [x] `scripts/smoke_test.sh` — the "build; verify R; zero skips" gate below, as a runnable script
- [ ] Build; verify R is 4.5.x; `selftest`; `test_methods_all_backends.py` with zero skips *(requires an hour-long amd64 build with network access — run `bash scripts/smoke_test.sh`, **with the FortiClient VPN disconnected**: see `docker/tls_inspection_ca_plan_260829.md`, which diagnoses the `curl: (60)` failure of the first attempt on 2026-08-29 and adds the preflight probe and the exit-code classification that replaced its misleading message; then `docker/r_package_deps_plan_260830.md`, which diagnoses the second attempt's failure 36 minutes in — `libuv1-dev` missing from the apt layer, surfacing five dependency levels away as a missing `qsmooth` — and drops the `Suggests` closure that `dependencies = TRUE` was pulling; then `docker/mid_build_network_loss_plan_260830.md`, which diagnoses the third attempt — the VPN reconnected 14 minutes in, discarding 857 s of finished compiles — and splits the R install into six self-verifying, individually cached layers so a network loss costs one stage; then `docker/recombat_install_plan_260830.md`, which diagnoses the fourth — reComBat declares the deprecated `sklearn` stub and `python <3.11`, so pip refuses it twice over — and installs it alone at a pinned commit with `--no-deps --ignore-requires-python`; then `docker/image_verification_failures_plan_260830.md`, for the first image that **built** — its gate failed 3 of 8 checks on one rpy2 invisible-result bug, and fixing that exposed `27_dwd` silently returning its input unchanged, `19_tdm`'s missing `binr`, `37_fabatch`'s mangled annotation column, and 4.5 GB of CUDA pulled in by `harmonypy` 0.2.0; finally `docker/remaining_selftest_failures_plan_260830.md`, which took the gate from 7/8 to a clean pass — HarmonizR's default `ComBat_mode` writes an empty file, and Shambhala needs several hundred genes where the fixture gave it 200)*
- [ ] `docker login ghcr.io`; push; **flip the GHCR package to public**; verify an anonymous `docker pull` *(requires a `write:packages` token and the GitHub UI)*

**Three deviations from the plan text, each forced by something the plan itself asks for.**

*`ENTRYPOINT ["combobatch"]` was not taken literally.* Five of §8's own verification commands
run something else inside the image (`pytest`, `Rscript -e`, `octave --version`,
`python -c …`), and the interactive pod overrides the command with `sleep infinity`; a literal
`combobatch` entry point makes all six impossible. `entrypoint.sh` execs a real executable as
given and hands everything else to the CLI, so `docker run IMAGE run --config c.yaml` works
and a mistyped subcommand gets argparse's did-you-mean instead of `exec: nonsense: not found`.
Six dispatch paths are tested by running the script.

*The base image is `python:3.11-slim-bookworm`, not `python:3.11-slim`.* The R pin is a
Debian 12 package; when the `3.11-slim` tag moves to the next Debian release, the pin breaks.

*`requirements.txt` gained `pytest`, `pytest-timeout` and `moto[s3]`.* The image runs
`tests/integration/test_methods_all_backends.py` as its acceptance gate — it is the only
environment where every backend is actually present — so the runner has to be in it. Floors
match the `dev` extra. The Dockerfile also sets `COMBOBATCH_FULL_ENV=1`, without which that
file self-skips and an image with a broken backend would pass its own gate.

**Verified without building:** `docker buildx build --check` reports no warnings; the Posit
`r-4.5.3_1_amd64.deb` still returns HTTP 200 (67 MB, re-checked 2026-08-29); all three shell
scripts pass `sh -n`/`bash -n`; the R script's `required` vector and the registry union are
identical 21-element sets.

### Phase 7 — K8s ✅ complete (2026-08-29, except the live launch)
- [x] `k8s/pod-combobatch.yaml` — **public GHCR image, no `imagePullSecrets`**, SSH key from **Secret**, `capacity-type: on-demand`, `karpenter.sh/do-not-disrupt`, thread env vars, no ConfigMaps
- [x] `k8s/job-combobatch.yaml`, `k8s/create_pvc.yaml`
- [x] `k8s/README.md`, `k8s/CLAUDE.md`
- [x] `tests/unit/test_k8s_manifests.py` — **new, not in the original checklist.** Every constraint in these manifests fails *silently*: a `nodeSelector` on a non-existent label leaves the pod `Pending` with no error, a missing `capacity-type` puts a six-hour job on a reclaimable Spot node, and a reintroduced requirements ConfigMap drifts from `requirements.txt` unnoticed. `kubectl apply` validates none of it. 32 assertions do.
- [ ] Launch; confirm startup is an image pull, not a 40-minute install; run `selftest` in-pod *(blocked on the Phase 6 live steps — the image must be pushed and made public first)*

**One plan-level correction, verified against the live cluster.** §5.11 specified
`karpenter.k8s.aws/capacity-type: on-demand`, inherited from the donor's
`nohup_job_failing.md`. That label does not exist: only *instance* attributes carry the
AWS-provider prefix, and capacity type is a core Karpenter label. `kubectl get nodes -L
karpenter.sh/capacity-type,karpenter.k8s.aws/instance-family` shows the former populated
(`on-demand`/`spot`) and the latter empty. The manifests use `karpenter.sh/capacity-type`,
§5.11 above is corrected, and a test pins both the right spelling and the absence of the
wrong one. Had it shipped as planned, the pod would have sat `Pending` indefinitely with
nothing in the events explaining why.

**One Phase 6 amendment, forced by this phase.** The pod offers SSH, but the image had no
`openssh-server` — so the pod would have had to `apt-get` it at startup, which is exactly
the boot-time-install pattern `k8s/CLAUDE.md` forbids and which needs network egress the
node may not have. `openssh-server`, `rsync`, `tmux` and `procps` are now in the image's
apt layer (~10 MB in a 6–9 GB image). Host keys are still generated per pod, so two pods
never share an SSH identity.

**Verified without launching:** `kubectl apply --dry-run=client` accepts all three
manifests against the live cluster's schema; both referenced Secrets (`aws-credentials`,
`danya-nikitin-ssh-pubkey`) exist in `rnd-sandbox`; the PVC `danya-nikitin-fl` is Bound at
100 GiB on `ebs-gp3-delete`, matching `create_pvc.yaml` exactly; the cluster runs Karpenter
v1 (`nodepools.karpenter.sh`), for which `do-not-disrupt` is the correct annotation name.

### Phase 8 — example dataset and example configs ✅ complete (2026-08-29)
- [x] `scripts/build_example_dataset.py` — the 5 verified `COHORT_LABEL` values; **drops the 14 duplicate bare-ID GSE62241 rows**; keeps the ragged 21,890-gene matrix
- [x] Generate + commit `example_exp.tsv.gz` (**144 × 21,890, 13.82 MB** — the estimate was 13.89 MB) and `example_ann.csv` (144 × 12)
- [x] `examples/example_dataset/PROVENANCE.md` — accessions, platforms, counts, diagnosis breakdown, both judgement calls, per-series GEO attribution
- [x] `examples/configs/{minimal,full_crossproduct,hyperparameters,shambhala}.yaml`
- [x] `examples/configs/best_approaches/` — 5 configs + README carrying all 8 caveats
- [x] Copy the FL marker panel to `examples/marker_panel_fl.csv` (633 genes, 617 present in the example matrix)
- [x] `tests/unit/test_example_dataset.py` — **new.** Every shipped config parses and validates *against the shipped annotation*, so a renamed method or column breaks the suite rather than the user's first command; and the numbers the READMEs quote are asserted rather than remembered.

**The plan's predictions held, with one arithmetic slip.** Measured: 158 rows → 144 samples,
13.82 MB gzipped, 5,877 all-NA genes, 6,797 genes present in all five cohorts, 41.7 % NA
overall, and GSE69033 + GSE119234 sharing exactly the predicted 6,925-gene subspace. The
14 duplicate accessions and the 4 diagnosis disagreements are all real and all confirmed.

§5.12's *deduped* diagnosis breakdown is wrong, though: it lists FL 50 and sums to 140, not
144. The truth is FL 54, DLBCL 40, Memory 16, Naive 11, Centroblast 6, Centrocyte 6,
Plasma 6, GC 5. `PROVENANCE.md` carries the measured table.

**A design conflict this phase surfaced, and the fix.** `RunConfig.validate()` rejected any
method listed twice — which forbids the parameter sweep the tool exists for. `10_mnn` at
default `k` and `10_mnn` at `k=50` write to *different* files (`strict__10_mnn` and
`strict__10_mnn__k50`), so there was never a collision to prevent. Validation now keys on
`(name, tag)`. And because the effective tag is only knowable with the registry in hand,
`build_job_matrix` additionally rejects two jobs that resolve to the same output key — which
catches the case validation cannot see, a parameter restated at its default value producing
no tag at all. Both directions are tested.

**Three of the five best-approach configs are degenerate on the example dataset**, and their
README says so in the table rather than leaving it to be discovered: `_FFPE_` leaves 40
samples in one batch, `RNASeq` leaves 34 in one batch, and the dataset contains no Illumina
arrays at all. They are faithful to the benchmark's recommendations and demonstrate the
invocation; the scenarios they are named for need the user's own data.

**Two deviations.** `minimal.yaml` uses `15_fsqn_py`, not `05_combat`: it is the first
command a new user runs and must work on a bare `pip install -e .`, which `05_combat`'s
R dependency would break (a test now enforces that). And `PROVENANCE.md` does not name the
internal source bucket — the file is destined for a public repository, and the URI would
have contradicted its own first line.

**Verified end to end on the real data:** `combobatch run --config examples/configs/minimal.yaml`
completes in 12 s (`strict` keeps 6,797 of 21,890 genes, exactly as predicted), and a
4-job dispatch produced a result worth putting in the README — `02_median_scaling` cuts
`r2_batch` from 0.297 to 0.038 while cutting `r2_bio` from 0.313 to 0.133. A single
batch-effect number would have ranked it first; that is caveat 1 reproducing itself on the
shipped dataset.

### Phase 9 — tests and docs ✅ complete (2026-08-29, except the full-scale validation)
- [x] `tests/conftest.py` + synthetic fixtures *(landed in Phase 1; no edit needed)*
- [x] `tests/e2e/` — all 8 scenarios from §5.13, 18 tests, 127 s against the real example dataset
- [x] `combobatch selftest` — 34 methods + 4 imputers, PASS/SKIP/FAIL, `--with-r` / `--with-octave` to turn a missing backend from a skip into a failure
- [x] `docs/HYPERPARAMETERS.md`, `docs/METRICS.md`, `docs/METHODS.md` (with the 5 removed methods, the Procrustes licensing rationale and the `20_shambhala` substitution), `docs/OUTPUT_LAYOUT.md`, `docs/BEST_APPROACHES.md` (all 8 caveats)
- [x] `tests/unit/test_docs_in_sync.py` — generated docs match the registries, byte for byte
- [x] `README.md` + `README.md`/`CLAUDE.md` pairs — added `combobatch/vendor/{README,CLAUDE}.md`, `docker/CLAUDE.md`, `examples/configs/README.md`, and removed six stale "*Scaffold*" notices
- [x] `.github/workflows/ci.yml` — Black, unit tests on 3.11, docs-in-sync, `selftest`, e2e; image build + push to GHCR on a `v*` tag, gated on `selftest --with-r --with-octave` and zero method skips
- [ ] Full-scale validation: reproduce one known dissertation result *(needs the pod and the 1.9 GB matrix; blocked on the Phase 6 and 7 live steps)*

**The three generated docs are generated, not written.** `combobatch/docsgen.py` renders
`METHODS.md`, `HYPERPARAMETERS.md` and `METRICS.md` from the registries, and the sync test
compares byte for byte. `METRICS.md` goes further: it **runs** each group on a fixture and
lists the keys that actually come back, so the key list is an observation rather than a
promise. Group F is excluded from probing on purpose — it needs R, and probing it would
make the file differ between a laptop and the image.

**Two fields were added to `MetricGroupSpec`.** §5.14 requires METRICS.md to state each
group's polarity and citation, and neither existed in the registry. Putting them anywhere
else would have been the parallel-table pattern this project exists to avoid, so all 14
entries now carry `polarity` and `citation` inline.

**A `--skip-slow`-shaped bug in my own selftest fixture, found by the fixture failing.**
The donor punches NAs uniformly at 10%; at 200 genes that leaves *every* gene with at least
one NA, so `strict` correctly keeps zero genes and the imputation axis is untestable.
`punch_holes` now confines NAs to a quarter of the genes, matching how real
platform-driven missingness actually behaves, and `_run_one` treats an empty gene axis as a
failure distinct from an all-NaN one.

**A correction to Phase 8, found by e2e scenario 8.** §5.12 claims the ragged coverage
makes the dataset "a real test of the imputation axis". True — but only above a threshold,
and the default is below it. Missingness here is *banded*, because a gene is measured by a
whole cohort or by none of it, so the per-gene NA fraction takes five values: 0.0 (6,797
genes), 0.3472 (9,080), 0.625 (8), 0.6528 (128) and 1.0 (5,877). `max_na_frac` only matters
where it crosses a band, and **the default 0.20 falls below every one of them** — so at
defaults `knn`, `softimpute` and `missforest` keep exactly the genes `strict` keeps and
return an identical matrix. `full_crossproduct.yaml` now sets 0.4, `PROVENANCE.md`
tabulates the bands, and a unit test pins them.

**Verified:** 840 unit tests, 18 e2e, 9 integration auto-skipped without R/Octave, Black
clean, `scripts/generate_docs.py --check` in sync, `combobatch selftest` exit 0.

The full-scale validation named above is `C_rnaseq_only` + `softimpute` + `04_sva`,
expressed as `--subset-query "RNA_BATCH.str.startswith('RNASeq')"`, checked against the
published `metrics_comprehensive.csv` row. It needs the pod, R, and the 1.9 GB matrix.
```

---

**Please review and approve before I implement.**
