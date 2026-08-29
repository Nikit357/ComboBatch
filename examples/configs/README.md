# Example configurations

Every file here runs against the shipped example dataset, from the repository root.
`tests/unit/test_example_dataset.py` validates all of them against the real annotation, so
a renamed method or column breaks the test suite rather than your first command.

| Config | Shows | Command |
|---|---|---|
| `minimal.yaml` | One imputer, one method, default metrics. No R required. | `combobatch run --config …` |
| `full_crossproduct.yaml` | 4 imputers × 12 methods, for the dispatcher. | `combobatch dispatch --config … --dry-run` |
| `hyperparameters.yaml` | Tuning both axes, including one method swept at two values of `k`. | `combobatch dispatch --config … --dry-run` |
| `shambhala.yaml` | `20_shambhala` and its four constraints. Needs Octave. | `combobatch run --config …` |
| `best_approaches/` | The source benchmark's recommendations, one per scenario. | see its README |

Anything on the command line overrides the file, so a config is a starting point rather
than a commitment:

```bash
combobatch run --config examples/configs/minimal.yaml --method 17_quantile --out ./elsewhere/
```

Two things worth reading before adapting one to your own data: the `columns:` block is the
only part that is genuinely required, and `max_na_frac` needs setting deliberately — on
this dataset the default leaves the imputers with nothing to do, for reasons
`../example_dataset/README.md` explains.
