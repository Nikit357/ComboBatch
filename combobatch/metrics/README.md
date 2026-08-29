# `combobatch/metrics/` — harmonization quality metrics

`METRIC_GROUP_REGISTRY` in `__init__.py` describes all 14 groups (A–N). It has no upstream
equivalent — the donor code dispatched groups through a hardcoded call sequence with a
separate sentinel dict in another file. See `docs/METRICS.md` for what each group measures.

| Module | Contents |
|---|---|
| `groups.py` | The 14 `compute_group_*` functions. |
| `embeddings.py` | Shared PCA / UMAP / t-SNE helpers, computed once and reused. |
| `panels.py` | Optional marker-panel loader (groups L and M only). |
| `runner.py` | `compute_metrics()` — group isolation, partial flushing, sentinels. |
| `concat.py` | Sidecars → `metrics_summary.csv` plus the long-format detail tables. |
