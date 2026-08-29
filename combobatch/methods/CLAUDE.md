# CLAUDE.md — `combobatch/methods/`

## Contract every method obeys

```python
normalize_<name>(exp_df, ann_df, *, batch_col, bio_col=None, **hyperparams) -> pd.DataFrame
```

`exp_df` is samples × genes; the return value is the same shape and index order. Registering
a method means adding one `MethodSpec` to `METHOD_REGISTRY` — nothing else, anywhere.

## Invariants

- **Declare every hyperparameter.** `MethodSpec.hyperparams` is what the CLI validates
  against, what `docs/HYPERPARAMETERS.md` is generated from, and what
  `tests/unit/test_registry.py` cross-checks against the function signature. A parameter the
  function accepts but does not declare is invisible and will never be reachable — that is
  precisely the bug this package exists to fix.
- **`r_arglist()` is the R-injection boundary.** R snippets are f-string-interpolated, so
  every value crossing into R goes through `r_arglist`, which rejects anything that is not a
  scalar or a list of scalars. Never interpolate a raw user value into an R string.
- **No hardcoded R-side arguments.** The donor code froze values like `polyFit(dis, 9)` and
  `RUVIII(k=5)` inside the R string, making them unreachable. Expose them as `HyperParam`s.
- **`target_group` defaults to `None`**, resolved by `resolve_reference_batch()` to the
  largest batch. Never default it to a dataset-specific batch name.
- **Availability is declared, not discovered.** `requires_r`, `r_packages`,
  `requires_octave` and `rnaseq_only` let the dispatcher preflight once at startup instead
  of failing hundreds of jobs in.

## Adding a method

1. Implement `normalize_<name>` in `python_methods.py` or `r_methods.py`.
2. Add a `MethodSpec` to `METHOD_REGISTRY` with full `hyperparams` and a citation.
3. That is all — CLI lists, docs and tests pick it up automatically.
