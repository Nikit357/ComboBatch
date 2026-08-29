# `combobatch/` — the package

| Module | Owns |
|---|---|
| `config.py` | Frozen dataclasses (`RunConfig`, `ColumnSpec`, `ImputationSpec`, `MethodSpec`, `MetricsSpec`, `PostRemovalSpec`), YAML loading, CLI overlay, and `validate()`. |
| `storage.py` | `StorageBackend` ABC with local and S3 implementations, selected by URI scheme. |
| `dataio.py` | Reading and writing expression/annotation matrices; extension-driven compression; index alignment. |
| `params.py` | Hyperparameter parsing and coercion, the output-key grammar, and the readable filename encoding. |
| `subset.py` | `--subset-query` evaluation and pre-flight validation. |
| `imputation.py` | `IMPUTER_REGISTRY` — `strict`, `knn`, `softimpute`, `missforest` behind one signature. |
| `postremoval.py` | Post-harmonization outlier-batch removal. |
| `harmonize.py` | `run_one_combination()` — the impute → harmonize → post-remove → score pipeline. |
| `memory.py` | cgroup-aware free-memory checks used by the dispatcher and workers. |
| `logging_utils.py` | Timestamped logging and BLAS/OpenMP thread pinning. |
| `methods/` | The harmonization method registry and implementations. |
| `metrics/` | The metric-group registry and implementations. |
| `cli/` | Argument parsing and subcommands. |
| `vendor/` | Third-party code vendored verbatim. |
| `data/` | Packaged data files shipped with the wheel. |
