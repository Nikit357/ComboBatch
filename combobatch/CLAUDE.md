# CLAUDE.md — `combobatch/`

## What this package owns

The whole tool. Three registries are its spine, each defined exactly once:
`METHOD_REGISTRY` (`methods/__init__.py`), `IMPUTER_REGISTRY` (`imputation.py`), and
`METRIC_GROUP_REGISTRY` (`metrics/__init__.py`). Every CLI choice list, validation check,
preflight probe and generated doc derives from them. Never hardcode a second copy of a
method, imputer or group list anywhere — that duplication is the single most common defect
in the code this package was ported from.

## Invariants

- **Expression matrices are samples × genes** in every public signature. R-backed methods
  transpose internally; the vendored Shambhala's `harmonize_parallel` takes and returns
  genes × samples. Those are the only two exceptions and both are local.
- **No hardcoded column names, buckets, prefixes or reference batches.** Everything comes
  from `ColumnSpec` / `RunConfig`. If you find yourself typing `"RNA_BATCH"`, stop.
- **Fail loudly.** The donor code silently fell back to a different imputer, wrote a
  post-removal output identical to its input, treated S3 403 as 404, and dropped
  unparseable keys from job enumeration — each producing results that look successful.
  Prefer a raised exception with a named cause over a plausible-looking result.
- **Validate before expensive work.** `RunConfig.validate()` runs before any download.
- **Serialize with the NaN-safe encoder.** Plain `json.dump` emits bare `NaN`, which is
  invalid strict JSON.
- **Optional heavy dependencies load lazily.** Importing `metrics` must not require a
  marker-panel CSV; importing the package must not require rpy2.

## Style

Python 3.11, Black at 88 columns, type hints on new functions, docstrings with
`Parameters` / `Returns`. Comments explain *why*, not *what*.
