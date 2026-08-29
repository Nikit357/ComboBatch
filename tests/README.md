# Tests

```bash
pytest tests/unit -v                    # no network, no R, no Octave; runs on a laptop
pytest tests/unit/test_registry.py -v   # a single module
pytest tests/integration -v             # auto-skips what the environment cannot run
pytest tests/e2e -v                     # drives the CLI over examples/example_dataset/
```

| Directory | Scope |
|---|---|
| `unit/` | Pure logic. Must stay green with no network and no system dependencies, in under a minute. |
| `integration/` | R-backed methods, R imputers, metric group F, Shambhala. Skipped unless the backend is present; inside the Docker image nothing may skip. |
| `e2e/` | Full CLI runs against the committed example dataset, local and S3 (`moto`) backends. |
| `fixtures/` | Small committed inputs and determinism baselines. |

`combobatch selftest` is a separate, coarser gate: it runs every method and imputer on a
synthetic matrix and is the acceptance check for a built image.
