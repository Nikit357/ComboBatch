# CLAUDE.md — `combobatch/metrics/`

## What this subpackage owns

`METRIC_GROUP_REGISTRY` and the 14 group functions. The registry drives `--groups`
validation, `--skip-slow`, the incremental sentinel logic, `combobatch list-metrics` and the
generated `docs/METRICS.md`.

## Invariants

- **Group failures are isolated.** Every group runs in its own try/except; a failure records
  `error_<LETTER>` and the others still run. Results are flushed after each group, so a
  crash or timeout in a later group never discards earlier ones.
- **Columns come from `ColumnSpec`, never from literals.** `column_roles` declares which
  roles a group needs (`batch`, `bio`, `cohort`); the concrete names are the user's. Metric
  keys are suffixed with the column they were computed on, so sentinel key templates must be
  resolved against the configured columns, not a fixed string.
- **Group L takes the pre-harmonization matrix in memory.** No baseline is downloaded and no
  reference cache exists — a stale cache would silently compare against the wrong baseline.
- **Group L has a ceiling.** Any per-batch monotone transform leaves within-cohort Spearman
  correlations unchanged, so group L sits near 1.0 for that whole class of method by
  construction. It is a "nothing broke" guard rail; group M's margin carries the
  discriminative signal. Say so wherever group L is reported.
- **Group N's class labels are configurable.** Never reintroduce dataset-specific diagnosis
  strings; default to the most frequent levels of the configured class column and let the
  metric names carry the actual class count.
- **Thread pinning happens in the Python entry point**, before any numeric import. Pod-level
  environment variables do not reach processes launched over SSH.
- **`panels.py` loads lazily.** Importing this subpackage must not require a marker CSV.

## Adding a group

Add the function to `groups.py` and one `MetricGroupSpec` to the registry, including its
sentinel key template and `column_roles`. Docs and CLI follow automatically.
