# `combobatch/methods/` — harmonization methods

`METHOD_REGISTRY` in `__init__.py` is the single source of truth: 34 methods, keyed by their
original benchmark identifiers (`01_raw` … `38_harman`). The numbering has gaps at 24, 32, 35,
36 and 39 — those methods are excluded, and keeping the original keys preserves traceability
to the published benchmark. See `docs/METHODS.md`.

| Module | Contents |
|---|---|
| `base.py` | `MethodSpec` and `HyperParam` — capability flags, hyperparameter declarations, citation. |
| `rinterop.py` | rpy2 helpers: converters, R garbage collection, availability probes, and `r_arglist`. |
| `python_methods.py` | The 13 pure-Python methods. |
| `r_methods.py` | The 20 R-backed methods. |
| `shambhala_method.py` | Adapter over the vendored Shambhala implementation (`20_shambhala`). |
