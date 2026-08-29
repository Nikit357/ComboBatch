# CLAUDE.md — `combobatch/cli/`

## Invariants

- **Choice lists come from the registries.** `--methods`, `--imputations` and `--groups`
  validate against `METHOD_REGISTRY`, `IMPUTER_REGISTRY` and `METRIC_GROUP_REGISTRY`. Never
  hand-maintain a parallel list of valid values.
- **Validate before spending.** Config validation, hyperparameter-name checking and the
  R/Octave preflight all run before the first download or subprocess launch.
- **`dispatch` supervises subprocesses, not threads.** A `ThreadPoolExecutor` pumps pipes;
  the real work is `subprocess.Popen` of `run`. This is not incidental — threads around rpy2
  break, and `multiprocessing` fork breaks R. Do not "simplify" it.
- **`--dry-run` must resolve everything and run nothing**: full job matrix, output keys,
  totals.
- **Exit codes are a contract**: `0` success or cached, `1` failure, `2` insufficient memory,
  `3` a required input is missing. The dispatcher distinguishes them; keep them stable.
- **A SIGTERM handler flushes the failed-jobs log** before exiting — Karpenter gives a
  two-minute warning and an unflushed log loses the run's resume state.
