# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Status: pre-implementation

**This repository contains no code yet.** The only file is
`combobatch_tool_plan_260828.md` — an approved-pending-review implementation plan. It is not a
git repository yet either.

Do not treat the plan as a description of existing behaviour. Every command in the "Commands"
section below is **planned**, not working. Before writing any code, read the plan in full: it
carries per-file specifications, a 10-phase TODO checklist, and — most importantly — the list of
source bugs that must be *fixed rather than copied*.

---

## What ComboBatch is

A standalone, dataset-agnostic CLI tool for benchmarking bulk transcriptomic batch-effect
correction. The user supplies an expression matrix and an annotation table (local or S3), names a
biology column and a batch column, and selects any combination of:

- **4 imputation strategies** — `strict`, `knn`, `softimpute`, `missforest`, with a configurable
  maximum per-gene NA fraction and per-imputer hyperparameters
- **34 harmonization methods** — the source's 39-method registry minus four that can never run
  (`24_peer_k10`, `32_deepmnn`, `35_dasc`, `36_explobatch`) and one excluded on licensing grounds
  (`39_procrustes` — BostonGene proprietary, non-commercial, incompatible with an MIT public
  tool), with `20_shambhala` re-backed by the vendored pure-Python + Octave implementation.
  Original keys are kept, so the numbering has gaps; each method accepts hyperparameters as
  keywords
- **14 metric groups (A–N)** — computed inline with the harmonization, or later in a separate run

A dispatcher runs the full (imputation × method) cross-product in parallel. Outputs and metrics
are written to a local directory or an S3 prefix, in the same layout either way. Everything ships
in one Docker image on a public registry, which a Karpenter-provisioned K8s pod consumes instead
of installing R and Octave at startup.

**Why it exists.** The FL dissertation's harmonization benchmark is genuinely reusable science,
but it is welded to one dataset: hardcoded column names, an FL-specific S3 bucket and prefix, 14
FL filter strategies, FL diagnosis string literals in the metrics, and a marker-gene panel
required at *import* time. ComboBatch is that pipeline with the FL assumptions lifted into
configuration. It is the tool named in the Article 1 title
("*…a novel computational pipeline ComboBatch…*"), so it is intended to become a public,
citable repository.

---

## Source repositories (read these before porting)

ComboBatch ports code from three sibling directories under
`../Follicular_lymphoma_disser/`. They are **read-only donors** — never write to them, and never
write to the dissertation's S3 prefix `FL_batch_correction/`.

| Donor | Supplies | Read first |
|---|---|---|
| `harmonization-scripts/` | `bench_shared.py` (39 methods, 4 imputers, S3 I/O, PCA R²), the two-stage dispatcher/worker pattern, `Dockerfile`, `install_r_packages.R`, `k8s/` | `CLAUDE.md`, then `project_overview.md` |
| `shambhala_adoption/Shambhala_containerized/` | `shambhala/*.py` + `octave/*.m`, the Octave bridge, Q-rescaling, NA handling | `CLAUDE.md` |
| `harmonization-metrics-calculation/` | `compute_batch_metrics.py` (groups A–N), the metrics worker/dispatcher/concat trio, `marker_panels.py` | `CLAUDE.md` |

**Do not trust the donors' `k8s/README.md` and `k8s/CLAUDE.md`.** They describe a pod with a
`/app/venv`, a Posit R install, and a Procrustes git clone — none of which exist in the manifests
on disk. Read the YAML directly.

---

## Architecture

Three registries are the spine of the design. Each is defined **exactly once**, and every CLI
choice list, validation check, and preflight probe derives from it.

| Registry | Location | Replaces |
|---|---|---|
| `METHOD_REGISTRY` | `combobatch/methods/__init__.py` | `METHODS` in `bench_shared.py` **plus** `ALL_METHODS` retyped by hand in two dispatchers |
| `IMPUTER_REGISTRY` | `combobatch/imputation.py` | `ALL_IMPUTATION`, retyped in two dispatchers; and the `strict`-lives-in-a-different-function split |
| `METRIC_GROUP_REGISTRY` | `combobatch/metrics/__init__.py` | **Nothing — this does not exist upstream.** Group dispatch there is a hardcoded `_run_group(...)` sequence with a parallel `GROUP_SENTINEL_KEYS` dict in a *different file*. It must be written from scratch. |

`MethodSpec` carries what the donor's bare `(fn, tier)` tuple cannot: `requires_r`,
`r_packages`, `requires_octave`, `rnaseq_only`, `uses_reference_batch`, `status`, and a
declaration of each tunable hyperparameter. That metadata is what lets the dispatcher validate
`--method-params` and probe for R/Octave *once at startup*, instead of discovering a missing
dependency 400 jobs into a run.

Data flow for one job:

```
load + align (index intersection, never positional)
  → optional --subset-query
  → impute            [cached per (imputation, param-tag) under prepared/]
  → resolve reference batch
  → resolve + validate hyperparameters
  → normalize
  → optional post-removal (post0 / post1)
  → optional metrics
  → write matrix + JSON sidecar
```

**The dispatcher must stay `ThreadPoolExecutor` supervising `subprocess.Popen` workers.** This is
not incidental: it is what gives every rpy2 job a fresh R session. Threads-around-R and
`multiprocessing` fork both break rpy2. The threads only supervise pipes.

### Concepts inherited from the donors

- **Post-removal** (`post0`/`post1`) — drop the most PCA-deviant batch *after* harmonization.
  Generic, so it is kept, but defaulted **off** (`--post-removal`).
- **Harshness tier** (`low`/`medium`/`high`) — per-method metadata, carried through to sidecars.
- **Incremental metric sentinels** — a group is recomputed only if its sentinel key is absent, so
  new groups can be added to finished runs without recomputing the rest.

### Deliberately not ported

The 14 FL filter strategies (`build_filter_strategies`, `RARE_GROUPS`, `BAD_BATCHES_A/B`,
`_AFFY_EXT_BATCHES`), the `load_data()` FL fixups (the `Diagnosis_with_coo` derived columns, the
`PUB_Suntsova_GSE120795` batch relabel), the 18 `shambhala_<P>_<Q>` variant keys, and the
`01_raw` magic-baseline machinery. Subsetting is replaced by a generic `--subset-query`
(a pandas `DataFrame.query` string).

---

## Non-negotiable constraints

Each of these was established by a documented failure in the donor repos. Changing one without
reading the plan's §7 will break the build or corrupt results.

**Python 3.11 exactly.** `reComBat 0.1.4` requires `pandas<2.0`; pandas 1.5.3 ships a py3.11
wheel only (3.12 source builds die on Cython); Shambhala requires `>=3.11`. Both ends are pinned.

**`pandas>=1.5,<2.0`.** The donor metrics `requirements.txt` says `>=2.0`, but that floor is
unjustified — its modules use only `pd.DataFrame`, `pd.Series`, `read_csv`, `notna`, `read_excel`
and conservative methods, with zero pandas-2.x-only APIs. This is what makes a single unified
image possible. If it ever binds, the escape is to drop `30_recombat` and move to pandas 2.x.

**R must be 4.5.x, pinned.** rpy2 3.6.x has a C-level ABI incompatibility with R 4.6.0 — the
symptom is a segfault around `33_amdbnorm`. The donor Dockerfile installs `r-base` **unpinned**
from `bookworm-cran40` and will drift. Verify with
`Rscript -e "cat(R.version\$version.string)"` in the built image.

**Build for `linux/amd64`.** The dev Mac is arm64; the c6a nodes are amd64. Always
`docker buildx build --platform linux/amd64`, or the pod fails with `exec format error`.

**Expression orientation is samples × genes** throughout the public API. Two internal exceptions,
both easy to get wrong: R-backed methods transpose to genes × samples internally, and Shambhala's
`harmonize_parallel` takes *and returns* genes × samples while every other Shambhala module uses
samples × genes.

**Shambhala input must be raw (non-log) scale.** `readExpressionData.m` applies `log2(x+1)`
itself, and the output is `exp(rm + rs·log(x+1))` — also raw scale, values reaching ~100 000.
There is no scale detection and no `--log-input` flag anywhere in the donor code.

**Shambhala parallelism multiplies.** Total Octave processes = outer workers × inner `n_workers`.
Default the inner count to 1 and pass `disable_progress=True` from the dispatcher.

**Thread pinning belongs in the Python entry point**, not the pod manifest. The donors proved
that pod-level `OMP_NUM_THREADS` env vars do not reach SSH-launched processes; unpinned BLAS
threads then oversubscribe the node and cost more than the parallelism gains.

**Hyperparameters make output keys non-unique.** Since two runs of `10_mnn` with different `k`
would otherwise collide, keys carry a parameter tag and `run_manifest.json` maps every tag back
to its resolved parameters. Omit the tag when a method runs at pure defaults.

---

## Commands (planned — none work yet)

```bash
# Setup. Python 3.11.16 is installed via Homebrew; `python3` is still the system 3.9.6,
# so always invoke python3.11 explicitly.
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

# Tests
pytest tests/unit -v                       # no network, no R, no Octave; must stay green
pytest tests/unit/test_registry.py -v      # a single module
pytest tests/unit/test_hyperparams.py::test_mnn_k_changes_output -v   # a single test
pytest tests/integration -v                # auto-skips without R / Octave
pytest tests/e2e -v                        # runs the CLI against examples/example_dataset/

# Dependency smoke test — all 39 methods + 4 imputers on synthetic data.
# Exit 0 iff every method PASSes or intentionally SKIPs. This is the image's acceptance gate.
combobatch selftest

# Discovery — these read the three registries
combobatch list-methods    # 34 rows, with each method's tunable hyperparameters
combobatch list-metrics    # 14 rows, A-N
combobatch list-imputers

# One combination (also the worker subprocess that `dispatch` launches)
combobatch run --config examples/configs/minimal.yaml
combobatch run --exp ./exp.tsv.gz --ann ./ann.csv --out ./out/ \
    --batch-col RNA_BATCH --bio-col Diagnosis_cell_type_unified \
    --imputation knn --method 10_mnn --method-params '10_mnn:k=50'

# Cross-product; --dry-run prints the job matrix and output count without running
combobatch dispatch --config examples/configs/full_crossproduct.yaml --n-workers 8 --dry-run
combobatch dispatch --config examples/configs/full_crossproduct.yaml \
    --n-workers 8 --skip-if-exists --memory-limit-gb 8.0

# Metrics computed later, against outputs that already exist
combobatch metrics --config examples/configs/minimal.yaml --groups A,B,E,J,K
combobatch concat --out-dir ./tables --date-tag 260829

# Image (amd64 is mandatory from an arm64 Mac)
docker buildx build --platform linux/amd64 -t combobatch:test -f docker/Dockerfile .
docker run --rm --platform linux/amd64 combobatch:test \
    python -c "import rpy2.robjects as ro; print(ro.r('R.version.string')[0])"   # MUST be 4.5.x
docker run --rm --platform linux/amd64 combobatch:test combobatch selftest
bash scripts/build_and_push_image.sh       # buildx → ghcr.io/nikit357/combobatch

# K8s
kubectl apply -f k8s/pod-combobatch.yaml -n rnd-sandbox
kubectl exec -it danya-nikitin-combobatch -n rnd-sandbox -- combobatch selftest
```

---

## Conventions

Inherited from the dissertation project: Python 3.11, `snake_case` / `PascalCase` /
`UPPER_SNAKE_CASE`, type hints on all new functions, docstrings with `Parameters` / `Returns`
sections, Black at line length 88, `pip` (not `uv`).

Specific to this repo:

- **Fail loudly.** The donors silently swallow several failures, and each one produces output that
  looks successful but is not: imputation falling back to `strict` while keeping the `__knn__`
  label, post-removal failing and still writing a `post1` key identical to `post0`, `s3_exists`
  treating a 403 as a 404, and job enumeration silently dropping any key that does not split into
  exactly four parts. Plan §6 lists 13 such bugs. Fix them; do not port them.
- **Validate before expensive work.** `RunConfig.validate()` must reject unknown method/imputer/
  group keys, missing annotation columns and an unwritable output root *before* anything is
  downloaded — the donors discover a bad strategy name only after loading 1.9 GB.
- **Serialize with the NaN-safe encoder.** Plain `json.dump` emits bare `NaN`, which is invalid
  strict JSON and rejected by non-Python parsers.
- **Never hardcode a column name, bucket, prefix, or reference batch.** `target_group=
  "RNASeq_FF_PolyA"` is the default of nine donor methods; here it resolves to the largest batch
  unless `--reference-batch` is given.
- Optional heavy dependencies load **lazily** — importing the metrics module must not require the
  marker-panel CSV to exist, as it does upstream.

## Infrastructure

- **Registry:** GHCR only — `ghcr.io/nikit357/combobatch`, public. The pod pulls it directly, so
  no `imagePullSecrets` and no AWS container registry are involved anywhere. Note a new GHCR
  package is **private by default** and must be flipped to public once in the GitHub UI, or the
  pod fails to pull with a 403 that looks like a missing image.
- **K8s:** namespace `rnd-sandbox`, `node-group=rnd-sandbox` toleration + required nodeAffinity,
  `karpenter.k8s.aws/instance-family: c6a`, plus `capacity-type: on-demand` — documented as needed
  in the donors but never applied, and necessary because Spot's 2-minute SIGTERM cannot save a
  6-hour missForest job.
- **PVC** `danya-nikitin-fl` at `/workspace` is `ReadWriteOnce` — one pod at a time.
- Image size **5.07 GB**, measured 2026-08-30 (Python + R 4.5.3 + ~25 R packages + Octave).
  Keep `harmonypy<0.2`: 0.2.0 depends on `torch` and took the image to 12.8 GB.
- Verified present on the dev Mac (2026-08-29): `docker` 29.7.2, `kubectl`, `aws` CLI 2.34.33
  with valid credentials, and **Python 3.11.16** via Homebrew at `/opt/homebrew/opt/python@3.11`.
  `python3` still resolves to the system 3.9.6, so invoke `python3.11` explicitly. Unit tests
  therefore run natively; only the R- and Octave-backed tests need the container.
