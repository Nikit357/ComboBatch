# `combobatch/cli/` — command-line interface

```
combobatch run              one (imputation, method) combination
combobatch dispatch         the imputation x method cross-product, in parallel
combobatch metrics          compute metrics over outputs that already exist
combobatch concat           aggregate sidecars into summary tables
combobatch selftest         smoke-test every method and imputer on synthetic data
combobatch list-methods | list-metrics | list-imputers
```

`run` is also the worker subprocess that `dispatch` launches — one process per job, which is
what gives every rpy2 job a fresh R session.
