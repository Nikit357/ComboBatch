# ComboBatch

Benchmark **imputation × harmonization** combinations for bulk transcriptomics, and score
the results with batch-effect quality metrics.

Point ComboBatch at an expression matrix and an annotation table — local or on S3 — name a
biology column and a batch column, and pick any combination of imputation strategies and
harmonization methods. A dispatcher runs the full cross-product in parallel and writes
harmonized matrices plus, optionally, quality metrics. Everything runs from one Docker image.

> **Status: v0.1.0, feature-complete and untagged.** Every subcommand works and the
> example dataset ships with the repository, so the quickstart below runs on a fresh clone
> with no R, no Octave and no S3 access. The published image has not been built yet — see
> [`docker/BUILD_AND_PUSH.md`](docker/BUILD_AND_PUSH.md).

## What it does

| Axis | Options |
|---|---|
| **Imputation** | `strict`, `knn`, `softimpute`, `missforest` — each with a configurable maximum per-gene NA fraction and full hyperparameter control |
| **Harmonization** | 34 methods, from `01_raw` and `03_limma` through `10_mnn`, `04_sva`, `16_fsqn_r` and `20_shambhala`, each accepting hyperparameters as keywords |
| **Metrics** | 14 groups (A–N): PCA variance, kBET/LISI/ASW, embedding dispersion, distribution tests, variancePartition, graph connectivity, NA retention, marker-correlation preservation, cross-batch rank agreement, predictive validation |

Metrics can be computed inline with the harmonization or later, in a separate run, over
outputs that already exist.

## Install

```bash
# From a clone. Python 3.11 is required (see "Why 3.11" below).
python3.11 -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"

combobatch --version
combobatch --help
```

The core install is deliberately **R-free and Octave-free**, so it works on any laptop and
runs the whole unit-test suite. Methods that need a system dependency live behind extras:

| Extra | Brings | Needed for |
|---|---|---|
| `.[r]` | `rpy2` | the 20 R-backed methods, the R imputers, metric group F. **Requires R 4.5.x already installed.** |
| `.[methods]` | harmonypy, scanorama, combat, inmoose, reComBat | `11_harmony`, `12_scanorama`, `07_pycombat`, `08_inmoose_combatseq`, `30_recombat` |
| `.[shambhala]` | qnorm | `20_shambhala` (also needs Octave) |
| `.[panels]` | openpyxl, xlrd | rebuilding a marker panel for metric groups L/M |
| `.[all]` | all of the above | — |

For everything at once, use the image instead:

```bash
docker run --rm ghcr.io/nikit357/combobatch:latest combobatch --help
```

## Quickstart

The repository ships a 144-sample, 5-cohort example dataset, so this runs on a fresh clone
with nothing else installed:

```bash
combobatch selftest                                     # what can run here
combobatch run --config examples/configs/minimal.yaml   # ~12 s, no R needed
```

```bash
# One combination
combobatch run --exp ./exp.tsv.gz --ann ./ann.csv --out ./out/ \
    --batch-col RNA_BATCH --bio-col Diagnosis \
    --imputation knn --method 10_mnn --method-params '10_mnn:k=50'

# The cross-product, in parallel
combobatch dispatch --config examples/configs/full_crossproduct.yaml \
    --imputations strict,knn,softimpute \
    --methods 01_raw,04_sva,10_mnn,16_fsqn_r \
    --n-workers 8 --dry-run

# Metrics later, over outputs that already exist
combobatch metrics --config examples/configs/minimal.yaml --groups A,B,E,J,K
```

Inputs are **samples × genes**. Outputs land under `{out}/exp/`, metrics under
`{out}/metrics/`, and every run writes a `run_manifest.json` recording the exact resolved
hyperparameters behind each output filename.

## Documentation

| Document | Contents |
|---|---|
| [`docs/METHODS.md`](docs/METHODS.md) | The 34 methods, citations, and which are excluded and why |
| [`docs/HYPERPARAMETERS.md`](docs/HYPERPARAMETERS.md) | Every method and imputer parameter, with worked examples |
| [`docs/METRICS.md`](docs/METRICS.md) | Metric groups A–N, their columns, outputs and cost |
| [`docs/OUTPUT_LAYOUT.md`](docs/OUTPUT_LAYOUT.md) | Output key grammar, manifest and sidecar schemas |
| [`docs/BEST_APPROACHES.md`](docs/BEST_APPROACHES.md) | Recommended configurations from the source benchmark, with caveats |

## Why Python 3.11 exactly

`reComBat` requires `pandas<2.0`; pandas 1.5.3 ships a Python 3.11 wheel only and fails to
build from source on 3.12; the vendored Shambhala requires `>=3.11`. Both ends are pinned,
and `requires-python` enforces it.

## Provenance

ComboBatch generalizes the harmonization benchmark built for a follicular lymphoma
transcriptomics study, lifting its dataset-specific assumptions into configuration. It is
the pipeline named in that work's methods.

## Licence

MIT — see [`LICENSE`](LICENSE), including the notes on invoked third-party code.
