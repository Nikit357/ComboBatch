# `combobatch/data/` — packaged data

Data files shipped inside the wheel and read through `importlib.resources`.

| Path | Contents |
|---|---|
| `calibration/P0_standard.csv.gz` | Shambhala P reference: 39 samples × 11,768 genes, samples-as-rows, **raw linear scale**. 3.0 MiB. |
| `calibration/Q0_standard.csv.gz` | Shambhala Q reference: 100 samples × 11,887 genes, samples-as-rows, **raw linear scale**. 7.2 MiB. |

These are the defaults for `20_shambhala` when the user supplies no `P`/`Q`. P is the
quantile-normalization anchor; Q defines the target expression shape. Their gene sets differ
and are intersected at load. Override with `params: {P: <path-or-uri>, Q: <path-or-uri>}`.

Both are loaded by `combobatch.methods.shambhala_method.load_calibration`, which caches
them — treat the returned frames as read-only.
