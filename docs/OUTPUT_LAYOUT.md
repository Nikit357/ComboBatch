# Output layout

Every run writes the same tree, whether the destination is a local directory or an `s3://`
prefix. Nothing about the layout changes between the two backends, so a run can be moved
between them without rewriting a path.

```
out/
├── run_manifest.json                       # what was run, with fully resolved parameters
├── prepared/
│   ├── strict__exp.tsv.gz                  # imputed matrix, cached per imputation
│   └── strict__ann.tsv.gz                  # its aligned annotation
├── exp/
│   ├── strict__10_mnn.tsv.gz               # harmonized matrices
│   ├── strict__10_mnn__k50.tsv.gz
│   └── knn__04_sva__post1.tsv.gz
├── genes/
│   └── strict__10_mnn_genes.json           # gene list per output, for cheap comparison
├── metrics/
│   └── strict__10_mnn_metrics.json         # one sidecar per output
├── metrics_summary.csv                     # written by `combobatch concat`
└── failed_jobs_<run_id>.txt                # only when something failed
```

## Key grammar

```
{imputation}__{method}[__{ptag}][__post{0|1}]
```

Segments are separated by a double underscore, and neither an imputer nor a method key may
contain one — that is checked when the key is built, not when it is parsed.

| Segment | Present when | Example |
|---|---|---|
| `imputation` | always | `strict`, `knn`, `softimpute`, `missforest` |
| `method` | always | `10_mnn`, `04_sva` |
| `ptag` | any parameter is not at its default | `k50`, `max_iter50-n_pcs30` |
| `post0` / `post1` | `--post-removal` is enabled | `post1` = deviant batch dropped |

With post-removal disabled the segment is **omitted entirely** rather than written as
`post0`, so an ordinary run's filenames stay short.

The prepared cache uses a shorter grammar, since it does not depend on the method:

```
{imputation}[__{ptag}]__{exp|ann}
```

**Unparseable keys raise.** The pipeline this was ported from silently dropped any key that
did not split into exactly four parts, which quietly excluded results from the summary
tables with nothing in the log.

## Parameter tags

`{"k": 50, "limit": 0.05}` encodes to `k50-limit0.05`: names sorted, values appended, joined
by single hyphens. Only **non-default** parameters appear, which is what keeps a default run
tag-free.

Tags longer than 64 characters are truncated and a stable hash of the *full* parameter dict
is appended, so a long sweep stays unique without producing an unusable filename. The
readable prefix survives, so the file is still recognizable.

The imputer and the method share this one slot. Tuning the same parameter name on both
sides raises rather than letting one silently overwrite the other; `param_tag` sets an
explicit tag when that happens or when the generated one is unwieldy.

## `run_manifest.json`

Written once per run, before any job starts. Filenames carry only non-default parameters;
the manifest is what makes a result self-describing.

```json
{
  "combobatch_version": "0.1.0",
  "run_id": "example",
  "created_utc": "2026-08-29T18:42:11+00:00",
  "input":   {"expression": "…/example_exp.tsv.gz", "annotation": "…/example_ann.csv"},
  "columns": {"batch": "RNA_BATCH", "bio": "Diagnosis_cell_type_unified", "cohort": "COHORT_LABEL"},
  "reference_batch": null,
  "subset_query": null,
  "random_seed": 42,
  "post_removal": {"enabled": false, "n_batches": 1, "min_batch_size": 20, "n_pcs": 2},
  "outputs": {
    "strict__10_mnn__k50": {
      "imputation": "strict",
      "imputation_params": {"max_na_frac": 0.0},
      "max_na_frac": 0.0,
      "method": "10_mnn",
      "method_params": {"k": 50},
      "param_tag": "k50",
      "post_rm": null,
      "harshness": "medium",
      "output": "exp/strict__10_mnn__k50.tsv.gz"
    }
  }
}
```

`imputation_params` and `method_params` are **fully resolved** — defaults included — so a
tag never has to be decoded by eye, and a result stays interpretable after the defaults
change in a later version.

## Sidecars

One JSON per output key, in `metrics/`, written whatever happened:

| Field | Meaning |
|---|---|
| `imp`, `method`, `param_tag`, `post_rm` | the job's identity, matching the key |
| `key` | the output key |
| `harshness` | how aggressively the method reshapes the matrix |
| `status` | see below |
| `run_id` | which run produced it |
| `error` | present only on failure; the exception text |
| `imputation_*` | what imputation did — genes kept, cells filled, residual NAs |
| metric keys | every group that ran, flat, e.g. `r2_RNA_BATCH` |
| `error_<LETTER>` | a group that failed or was skipped, with the reason |
| `metrics_groups_run` | which groups actually ran |

### Status values

| Status | Meaning | Exit code |
|---|---|---|
| `ok` | ran and wrote an output | 0 |
| `cached` | output already present, `--skip-if-exists` | 0 |
| `skipped` | the method declined this input, e.g. RNA-seq-only guard | 0 |
| `failed` | the method raised | 1 |
| `post_removal_failed` | harmonization worked, post-removal did not | 1 |
| `upload_failed` | computed, but the result could not be written | 1 |

A skip is never written as an output matrix, and `post_removal_failed` is never written as
a `post1` file identical to `post0` — the donor did exactly that, producing a result that
looked like a successful second variant.

## Exit codes

| Code | Meaning |
|---|---|
| 0 | success, or cached |
| 1 | the job failed |
| 2 | not enough memory to start |
| 3 | a required input is missing |

The dispatcher distinguishes all four and records failures in
`failed_jobs_<run_id>.txt`, which `--retry-failed` reads back. The file is flushed on
SIGTERM, so a Karpenter eviction or a `kubectl delete` still leaves a resumable run.

## JSON is strict JSON

Everything is written through a NaN-safe encoder: `NaN` and `±Inf` become `null`. Plain
`json.dump` emits bare `NaN`, which is invalid JSON and rejected by every non-Python
parser — the kind of defect found only when someone tries to read the results in R.
