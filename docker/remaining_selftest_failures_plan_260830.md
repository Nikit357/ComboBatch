# Plan — the last two selftest failures, and two things the report got wrong

**Status:** awaiting approval. No code has been changed.
**Date:** 2026-08-30
**Trigger:** `bash scripts/smoke_test.sh`, log at `docker/docker_build_logs_260830_3.txt`.
Build **FINISHED** in 2738 s; **7 of 8** gate checks pass.

---

## Overview

Everything from the previous four plans held. `19_tdm` passes (0.8 s) now that `binr` is
installed, `27_dwd` passes (20.6 s) and genuinely corrects, `37_fabatch` passes (3.8 s),
`list-methods` prints all 34 with their hyperparameters, the zero-skip integration suite
passes, and the Octave X11 noise is gone — `octave-cli` did its job, leaving only the
benign shutdown line.

**Selftest: 36 passed, 0 skipped, 2 failed.** The two are the ones deliberately left open,
and both are now diagnosed:

| | reported | actual cause |
|---|---|---|
| `21_harmonizr` | `returned 0 genes` | HarmonizR's **default** `ComBat_mode=1` writes an empty file. We never pass the parameter |
| `20_shambhala` | `returned an all-NaN matrix` | **200 genes is too few.** Not a defect in the image |

Neither is a build problem. One is a hardcoded R-side default we never exposed — the defect
class `combobatch/methods/CLAUDE.md` exists to prevent — and one is a limit of the selftest
fixture that the method should have refused loudly instead of returning NaNs.

Two further defects, both in what the run *reported* rather than what it did:

- **`rpy2  not installed`** appears in the ENVIRONMENT block of a run in which every R
  method passed. It is false.
- **The image is 5.07 GB**, and six documents still say 6–9 GB. The `harmonypy` pin removed
  **7.7 GB**; the claim is now wrong in the other direction.

---

## Background — what was measured

### `21_harmonizr`: the mode decides, not the data

`ComBat_mode` selects ComBat's two switches inside `HarmonizR:::splitting`:

| mode | `par.prior` | `mean.only` | 200×80, 3 batches, exponential | same, normal |
|---|---|---|---|---|
| **1 (default)** | TRUE | FALSE | **1 line — empty** | **1 line — empty** |
| 2 | TRUE | TRUE | 201 lines | 201 lines |
| 3 | FALSE | FALSE | **1 line — empty** | **1 line — empty** |
| 4 | FALSE | TRUE | 201 lines | 201 lines |

The two modes that fail are exactly the two with `mean.only = FALSE`, and they fail on both
distributions. So it is not our fixture: **ComBat's variance adjustment produces nothing
through HarmonizR's block machinery here, and HarmonizR reports that by writing a
three-byte file rather than by raising.**

`algorithm="limma"` is also empty. `ur=FALSE` makes no difference.

Our call passes `data`, `description`, `algorithm`, `output_file`, `plot` and `verbosity` —
never `ComBat_mode`. It is not in the `MethodSpec` either, so it is unreachable from a
config, which is the precise thing `methods/CLAUDE.md` forbids: *"A parameter the function
accepts but does not declare is invisible and will never be reachable."*

The description layout was **not** the problem — both our layout and the swapped one
produce full output at mode 2, so that hypothesis is discarded.

### `20_shambhala`: a gene-count floor, and a silent partial failure

Same fixture shape, only the gene count varying:

| genes | NaN fraction | wall clock |
|---|---|---|
| 200 (the current fixture) | **1.000** | 6.7 s |
| 400 | **0.150** | 18.6 s |
| 800 | 0.000 | 43.3 s |
| 1200 | 0.000 | 66.0 s |
| 2000 | 0.000 | — |

The gene-symbol fix from the previous plan worked: the log now reads *"Shambhala gene
intersection: 200 genes (input had 200, P 11768, Q 11887)"* and the Octave pipeline runs to
completion. It simply has too little to work with, and degrades gradually rather than
failing.

**The 400-gene row is the important one.** A 15% NaN matrix is the failure mode this
project cares about most: it is not all-NaN, so the selftest's guard would not catch it, and
a real run would write a partly-empty matrix labelled as harmonized.

### The two reporting defects

```
  Python        3.11.16
  rpy2          not installed      ← false; every R method in the same run passed
  R             R version 4.5.3
```

`describe_environment()` reads `rpy2.__version__`, which **rpy2 3.6 does not define**; the
`except Exception` then turns that `AttributeError` into "not installed".
`importlib.metadata.version("rpy2")` returns `3.6.7` — verified in the image.

Size, measured on the built image:

| | before the `harmonypy` pin | now |
|---|---|---|
| image | 12.8 GB | **5.07 GB** |
| `site-packages` | 5804 MB | **1124 MB** |
| `torch` / `nvidia` / `triton` | 4546 MB | absent |

---

## Files to change

### 1. `combobatch/methods/r_methods.py` — `normalize_harmonizr`

#### 1a. Accept and forward `ComBat_mode`

```python
def normalize_harmonizr(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    algorithm: str = "ComBat",
    combat_mode: int = 2,
    **kw,
) -> pd.DataFrame:
```

and in the R call:

```r
            harmonizR(data_as_input={r_literal(data_path)},
                      description_as_input={r_literal(desc_path)},
                      algorithm={r_literal(algorithm)},
                      ComBat_mode={r_literal(int(combat_mode))},
                      output_file={r_literal(out_base)},
                      plot=FALSE, verbosity=0)
```

The docstring records why the default is 2 rather than HarmonizR's own 1: modes 1 and 3
(`mean.only = FALSE`) write an empty file on data with no missing values, measured on two
distributions.

#### 1b. Refuse an empty result

HarmonizR signals failure by writing a three-byte file. Existing code reads it and returns
an empty frame, which reached the selftest as a generic "0 genes":

```python
        if result.empty or result.shape[0] == 0:
            raise RuntimeError(
                "21_harmonizr produced an empty result. HarmonizR writes an empty file "
                "rather than raising when its ComBat blocks yield nothing; "
                f"ComBat_mode={combat_mode} may not suit this matrix."
            )
```

### 2. `combobatch/methods/__init__.py` — declare it

```python
            "combat_mode": HyperParam(
                "combat_mode",
                int,
                2,
                "HarmonizR's ComBat parameter set: 1 par.prior/scale, 2 par.prior/"
                "mean-only, 3 non-parametric/scale, 4 non-parametric/mean-only. "
                "Modes 1 and 3 write an empty result on data without missing values.",
                choices=[1, 2, 3, 4],
            ),
```

`docs/HYPERPARAMETERS.md` and `docs/METHODS.md` regenerate.

### 3. `combobatch/methods/shambhala_method.py` — never return a NaN matrix quietly

The 400-gene case returns 15% NaN and nothing notices. After harmonization:

```python
    nan_fraction = float(np.isnan(out.to_numpy(dtype=float)).mean())
    if nan_fraction > _MAX_NAN_FRACTION:
        raise RuntimeError(
            f"20_shambhala returned {nan_fraction:.1%} NaN. Shambhala needs a few "
            "hundred genes at minimum — measured: 200 genes gives 100% NaN, 400 gives "
            "15%, 800 gives none. Widen the input or raise --max-na-frac."
        )
```

with `_MAX_NAN_FRACTION = 0.01` — anything above a rounding-level trace is a failure, not a
result.

### 4. `combobatch/cli/selftest_cmd.py`

#### 4a. Report rpy2's version correctly

```python
    try:
        import importlib.metadata as _md

        lines.append(f"rpy2          {_md.version('rpy2')}")
    except Exception:
        lines.append("rpy2          not installed")
```

`rpy2.__version__` was removed in 3.6; the distribution metadata is the supported route.

#### 4b. Give `20_shambhala` enough genes

`N_GENES = 200` is right for the other 33 methods and wrong for this one. Rather than
slowing every method down tenfold, the selftest builds a **wider matrix for Shambhala
only**:

```python
SHAMBHALA_N_GENES = 800   # 200 -> 100% NaN, 400 -> 15%, 800 -> clean (43 s)
```

`_run_one` uses the wide fixture when the spec is `20_shambhala`. The alternative — skipping
it — would stop the selftest proving that the Octave backend works, which is most of why
Octave is in the image at all.

### 5. `tests/unit/`

| Test | Guards against |
|---|---|
| `21_harmonizr` declares `combat_mode`, and its default is not 1 or 3 | silently reinstating a mode that writes an empty file |
| the `harmonizR(` call passes `ComBat_mode` | the parameter becoming unreachable again |
| `normalize_harmonizr` raises on an empty result | an empty matrix being returned as a success |
| `shambhala_method` has an all-NaN guard | a partly-NaN matrix passing as harmonized |
| `describe_environment` does not use `rpy2.__version__` | the false "not installed" returning |

### 6. Documentation — the size, with the measured number

`CLAUDE.md:251`, `docker/README.md:64`, `docker/BUILD_AND_PUSH.md:19,41,215`, and
`combobatch_tool_plan_260828.md:1050` say 6–9 GB. Measured: **5.07 GB**, `site-packages`
1124 MB. The disk requirement in step 0 (≥ 25 GB) stays as it is — the build cache, not the
image, sets that.

`BUILD_AND_PUSH.md` also gains a troubleshooting row for `21_harmonizr` returning nothing.

### 7. `combobatch_tool_plan_260828.md`

One clause under Phase 6 recording this as the fifth and — with these two fixed — final
diagnosis before a clean gate.

---

## Files that do NOT need to change

| File | Why not |
|---|---|
| `docker/Dockerfile`, `install_r_packages.R` | Nothing here needs a new package or layer; `binr`, `libuv1-dev` and the staged installs are all confirmed working |
| `requirements.txt`, `pyproject.toml` | The `harmonypy` pin did what it was meant to; nothing else changes |
| `combobatch/methods/rinterop.py` | `isinstalled` works — 34 methods enumerated, zero skips |
| `combobatch/vendor/shambhala/` | Vendored code is not at fault: the Octave pipeline ran to completion and returned a matrix. The guard belongs in the glue, per `vendor/CLAUDE.md` — *"Adapt at the boundary, not inside."* |
| `scripts/smoke_test.sh` | Reported the state accurately |

---

## Side effects and caveats

- **Only the source layer rebuilds.** Everything through pip stays cached; expect a rebuild
  of a couple of minutes, not an hour.
- **`21_harmonizr`'s output changes** — from an empty matrix to a real one, and from
  ComBat's full location-and-scale adjustment to mean-only. That is a genuine methodological
  choice and it is now visible and settable, where before it was neither. Nothing has been
  published from this method, since it has never produced output.
- **The selftest gets slower by roughly 40 seconds** for the wide Shambhala fixture. The
  suite already spends 63 s in `38_harman`.
- **The Shambhala NaN guard may fire on real user data** with few genes. That is the point;
  the message names the measured thresholds so the cause is actionable.
- **`20_shambhala` remains unproven below ~800 genes.** The guard turns silence into a clear
  error but does not make the method work on narrow matrices, and nothing here claims it
  does.
- **The remaining Octave line** (`error: ignoring const execution_exception&`) is unchanged
  and documented: exit code 0, and the bridge gates on return codes.

---

## Verification commands

```bash
# 1. The HarmonizR mode is the whole story — empty at 1 and 3, full at 2 and 4
docker run --rm --platform linux/amd64 -w /tmp combobatch:test Rscript -e '
 suppressMessages(library(HarmonizR)); set.seed(42); g <- 200; n <- 80
 m <- matrix(rexp(g*n, 1/5)+2, nrow=g, dimnames=list(paste0("G",1:g), sprintf("S%04d",0:(n-1))))
 ids <- colnames(m); bi <- (0:(n-1))%%3+1
 write.table(m,"/tmp/d.tsv",sep="\t",quote=FALSE,col.names=NA)
 write.csv(data.frame(ID=ids,batch=paste0("Batch",bi-1),sample=bi),"/tmp/desc.csv",row.names=FALSE)
 for (md in 1:4) { unlink("/tmp/o.tsv")
   harmonizR("/tmp/d.tsv","/tmp/desc.csv",ComBat_mode=md,output_file="/tmp/o",plot=FALSE,verbosity=0)
   cat("mode",md,"->",length(readLines("/tmp/o.tsv")),"lines\n") }'

# 2. rpy2 reports a version
docker run --rm --platform linux/amd64 -w /tmp combobatch:test \
    python -c "import importlib.metadata as m; print(m.version('rpy2'))"      # 3.6.7

# 3. Laptop suite, formatting, generated docs
python -m pytest tests/unit -q
black --check .
python scripts/generate_docs.py --check
docker buildx build --check -f docker/Dockerfile .

# 4. The size claim, measured not assumed
docker images combobatch:test --format '{{.Size}}'

# 5. Rebuild and the whole gate, VPN off
bash scripts/smoke_test.sh        # expect: passed all checks (0 skipped)

# 6. Both methods, against the built image
docker run --rm --platform linux/amd64 combobatch:test combobatch selftest \
    | grep -E "20_shambhala|21_harmonizr|RESULT"
```

---

## Verified in the image before the rebuild

All three fixes were exercised against the existing image with the corrected source
mounted over the installed package:

```
ENV: rpy2          3.6.7          (was "not installed")
RES 21_harmonizr:  (80, 200)  finite: True
RES 20_shambhala:  (80, 800)  nanfrac: 0.0
```

**One recurring self-inflicted pattern, noted for next time.** Three of the guards written
across these plans failed on their first run because the forbidden string —
`requireNamespace`, `read.csv(... ann_path`, `rpy2.__version__` — appeared in the *comment
explaining why it is forbidden*. Source-shape assertions need to strip comments before
matching; the ones here now do.

---

## TODO

- [x] `combobatch/methods/r_methods.py` — `normalize_harmonizr` accepts `combat_mode`,
      defaults to 2, forwards it as `ComBat_mode`
- [x] `combobatch/methods/r_methods.py` — raise on an empty HarmonizR result instead of
      returning an empty frame
- [x] `combobatch/methods/r_methods.py` — docstring: why the default is not HarmonizR's
- [x] `combobatch/methods/__init__.py` — declare `combat_mode` on `21_harmonizr`
- [x] `combobatch/methods/shambhala_method.py` — refuse a result that is more than 1% NaN,
      naming the measured gene-count thresholds
- [x] `combobatch/cli/selftest_cmd.py` — `importlib.metadata.version("rpy2")`
- [x] `combobatch/cli/selftest_cmd.py` — a wide (800-gene) fixture for `20_shambhala` only
- [x] `tests/unit/` — the five guards above
- [x] `docs/HYPERPARAMETERS.md`, `docs/METHODS.md` — regenerate
- [x] Documentation — the measured 5.07 GB in six places
- [x] `docker/BUILD_AND_PUSH.md` — troubleshooting row for an empty HarmonizR result
- [x] `combobatch_tool_plan_260828.md` — record this under Phase 6
- [x] Run verification steps 1–4 — 872 unit tests pass, Black clean, docs in sync,
      `buildx --check` clean; all three fixes confirmed in the image
- [ ] **VPN off**, rebuild and run `bash scripts/smoke_test.sh` to a clean pass
- [ ] Run verification step 6 against the rebuilt image
