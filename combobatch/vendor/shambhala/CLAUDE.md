# CLAUDE.md — vendored Shambhala

## Rules

1. **Vendored verbatim.** Adaptations belong in `combobatch/methods/shambhala_method.py`.
   The only edits permitted here are the specific upstream bug fixes listed below.
2. **`octave/*.m` must ship as package data.** Upstream packaged only the Python modules, so
   a non-editable install silently lacked the `.m` files it needs at runtime. The
   `package-data` entry in `pyproject.toml` is what prevents that; do not remove it.

## Load-bearing contracts

- **Orientation flips here.** `harmonize_parallel` takes and returns **genes × samples**,
  while every other module in this package uses samples × genes. The wrapper owns the
  transpose.
- **Input must be raw, non-log.** `readExpressionData.m` applies `log2(x+1)` itself, and the
  output is `exp(rm + rs*log(x+1))` — also raw scale, reaching ~100 000. There is no scale
  detection anywhere. Passing log-transformed data produces silently wrong numbers.
- **CuBlock requires NaN-free input**, so NA handling must run before the Octave call.
- **`ProcessPoolExecutor`, never threads.** Each worker owns an Octave subprocess; threads
  would share stdin/stdout file descriptors and corrupt the pipe.
- **Parallelism multiplies.** Total Octave processes = outer dispatcher workers × inner
  `n_workers`. The wrapper defaults inner to 1 and passes `disable_progress=True`.

## Upstream bugs — all four fixed, do not reintroduce

1. `skip_qn` not forwarded in the benchmark path → double quantile normalization. Fixed in
   `combobatch/methods/shambhala_method.py`, which is the *only* copy of that glue.
2. `fixed_clusters` overrode `skip_qn` in Octave script selection → the approved speed-up
   combination also double-normalized. Fixed by `_PIPELINE_SCRIPTS`, a total function over
   all four combinations, plus the previously missing
   `octave/Shambhala2_piped_preqn_fixed.m`. Add a combination only by adding a row *and* a
   script; never by re-introducing a second `if`.
3. `ProgressEvent` constructed with the wrong field names in the Python-CuBlock path →
   `TypeError` on any progress-enabled run. Fixed in `parallel.py`.
4. `write_expression` ignored the file extension. Fixed by not vendoring `io_utils.py` at
   all; `combobatch.dataio` handles every read and write.

Regression tests: `tests/unit/test_shambhala_vendor.py` (2, 3) and
`tests/unit/test_shambhala_glue.py` (1, 4).
