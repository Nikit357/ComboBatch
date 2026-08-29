# Vendored Shambhala2

Pure-Python + Octave implementation of the Shambhala2 harmonization algorithm
(Borisov et al.), vendored here so that `20_shambhala` works from a plain install.

**Vendored verbatim**, apart from the deviations listed below. Do not refactor these modules
to match ComboBatch conventions — that would make it impossible to diff against upstream.
The ComboBatch-facing wrapper lives in `combobatch/methods/shambhala_method.py`; put
adaptations there.

| Module | Role |
|---|---|
| `octave_bridge.py` | The only interface to Octave: formats TSV on stdin, parses stdout. |
| `parallel.py` | Splits samples into batches, dispatches via `ProcessPoolExecutor`, rescales. |
| `q_rescale.py` | Q reference statistics and the final rescaling. |
| `na_handling.py` | NA drop/impute before Octave, and restoration afterwards. |
| `cublock_python.py` | Experimental pure-Python CuBlock; not equivalent to the Octave one. |
| `progress_display.py` | Live per-worker progress. |
| `octave/*.m` | Octave sources. Shipped as package data. |

## Deviations from upstream

Each is marked `VENDOR FIX` at the site.

| Change | Why |
|---|---|
| Imports rewritten to `combobatch.vendor.shambhala.*` | The package moved. |
| `octave/` moved inside the package | Upstream left it a sibling of the package and never packaged it, so the default script directory only resolved under `pip install -e .`. |
| One script-selection table (`_PIPELINE_SCRIPTS`) | Upstream used two independent `if` blocks, so `fixed_clusters` silently overrode `skip_qn`. |
| New `octave/Shambhala2_piped_preqn_fixed.m` | The fourth combination had no script, which is what forced the override above. |
| `ProgressEvent` built with its real field names | Upstream passed `done`/`total`/`speed_s_per_sample`; the NamedTuple declares `sample_done`/`sample_total`/`is_final`, so any progress-enabled `python_cublock` run raised `TypeError`. |
| `io_utils.py` **not vendored** | Its writer ignored the file extension, so `.tsv.gz` was written as plain CSV. All I/O goes through `combobatch.dataio`. |
