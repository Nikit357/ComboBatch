# ComboBatch documentation

| Document | Contents |
|---|---|
| `METHODS.md` | The 34 harmonization methods, their citations and harshness tiers; the five excluded methods and why; the RNA-seq-only guard; reference-batch resolution. |
| `HYPERPARAMETERS.md` | Every method and imputer parameter — type, default, bounds, the R-side argument it maps to — with worked CLI and YAML examples. |
| `METRICS.md` | Metric groups A–N: what each measures, the keys it emits, which annotation columns it consumes, cost, polarity, citation. |
| `OUTPUT_LAYOUT.md` | Output key grammar, the `ptag` hyperparameter encoding, `run_manifest.json` and sidecar schemas. |
| `BEST_APPROACHES.md` | Recommended configurations carried over from the source benchmark, with the caveats that qualify them. |

`METHODS.md`, `HYPERPARAMETERS.md` and `METRICS.md` are **generated from the registries** in
`combobatch/methods/`, `combobatch/imputation.py` and `combobatch/metrics/`. Do not hand-edit
their generated sections — change the registry and regenerate. `tests/unit/test_docs_in_sync.py`
fails the build if they drift.
