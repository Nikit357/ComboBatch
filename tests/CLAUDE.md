# CLAUDE.md — `tests/`

## Invariants

- **`tests/unit/` must run with no network, no R and no Octave.** It is the suite that runs
  on the dev laptop and in CI. If a test needs a backend, it belongs in `integration/`.
- **Skips must be honest.** Methods whose backend is missing are skipped via the
  `MethodSpec` capability flags, never quietly passed. `tests/integration/test_methods_all_backends.py`
  asserts that **inside the image, zero methods skip** — that is what turns "34 methods" from
  a claim into a tested fact.
- **`test_hyperparams.py` is the anti-regression guard.** For each method with a declared
  hyperparameter, changing it must change the output. The bug this whole tool exists to fix
  was that every hyperparameter was silently unreachable, so a test asserting only "it runs"
  would have passed against the broken code.
- **`test_registry.py` cross-checks declarations against reality** via `inspect.signature`,
  and asserts the five excluded method keys are absent.
- **`test_docs_in_sync.py` fails the build when generated docs drift** from the registries.
- Use `moto` for S3; never touch a real bucket from a test.

## Fixtures

Synthetic data is built in `conftest.py` with a fixed seed. Committed fixtures stay small;
the example dataset (~13.9 MB) is shared with `examples/` and is not duplicated here.
