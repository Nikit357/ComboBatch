# Plan — the image builds; the acceptance gate does not pass

**Status:** awaiting approval. No code has been changed.
**Date:** 2026-08-30
**Trigger:** `bash scripts/smoke_test.sh`, log at `docker/docker_build_logs_260830_2.txt`.
Build **FINISHED** in 3990 s. The gate then failed **3 of 8** checks.

---

## Overview

**The image is real.** R 4.5.3 answers through rpy2, `variancePartition` and `bapred`
load, Octave and the `matlab` wrapper work, all 14 metric groups import without a panel
file. Four builds of work landed.

Three checks fail — `list-methods`, `test_methods_all_backends.py` and `selftest` — and all
three are **one bug, one line**:

```python
result = ro.r(f"requireNamespace('{name}', quietly=TRUE)")
return bool(result[0])          # TypeError: 'NoneType' object is not subscriptable
```

`requireNamespace()` returns its logical **invisibly**, and rpy2 3.6 maps an invisible
result to `None`. Measured in the built image:

| expression | returns |
|---|---|
| `requireNamespace('limma', quietly=TRUE)` | **None** |
| `requireNamespace('nosuchpkg', quietly=TRUE)` | **None** |
| `isTRUE(requireNamespace('limma', quietly=TRUE))` | `BoolVector` → True |
| `rpy2.robjects.packages.isinstalled('limma')` | `True` |

`None` for the installed package **and** the missing one, so this check could never have
worked — it can only crash. Every R-backed method and both R imputers go through it.

**I verified the fix without rebuilding**, by patching the installed file inside the image
that is already on your machine:

| check | before | after the one-line fix |
|---|---|---|
| `combobatch list-methods` | TypeError | **34 rows** |
| `pytest tests/integration/test_methods_all_backends.py` | 2 failed | **3 passed** — zero skips |
| `combobatch selftest` | TypeError, exit 1 | runs: **34 passed, 0 skipped, 4 failed** |

That last line is the point of this plan. The crash was hiding a second layer, and one item
in it is worse than the crash.

---

## The findings, in order of severity

### A — `27_dwd` silently returns its input unchanged, and reports success

Not in the failure list. It **passed**. While passing, it printed:

```
DWD failed for batch 'Batch1': argument "expon" is missing, with no default
DWD failed for batch 'Batch2': argument "expon" is missing, with no default
```

Every batch fell through this handler in `r_methods.py`:

```r
tryCatch({
    sol <- genDWD(X = X_combined, y = y)
    ...
    result[, batch_mask] <- X_batch + shift * w
}, error = function(e) {
    message("DWD failed for batch '", b, "': ", conditionMessage(e))
})
```

When it fails, `result[, batch_mask]` is never assigned, so `result` is returned as the
**uncorrected input**. The selftest checks a DataFrame with rows, genes and not-all-NaN — an
unchanged matrix satisfies all four. `27_dwd` would have produced a full benchmark column of
"harmonized" results that are the raw input.

The cause is exact. In the image:

```
genDWD formals:  X, y, C, expon, tol, maxIter, method, printDetails, rmzeroFea, scaleFea
no default:      X, y, C, expon
penaltyParameter formals: X, y, expon, rmzeroFea, scaleFea
```

We pass `X` and `y`. `C` and `expon` are required and absent, so the call **always** raises.

**The docstring above it already records this happening once:** *"Note `genDWD` is called
without a `penalty` argument: DWDLargeR removed it, and passing it made every batch fall
through to the error handler uncorrected."* Someone met the same symptom, removed an
argument, saw the message stop being fatal, and did not notice the batches were still
uncorrected. The swallow is what let both rounds pass.

This is the donor defect class named in `combobatch/CLAUDE.md` — *"post-removal failing and
still writing a `post1` key identical to `post0`"* — reproduced in our own code.

### B — four methods fail, now that they can be seen

```
19_tdm         FAIL  (22.9s)  Error in if (answer %in% allowed) break : argument is of length zero
20_shambhala   FAIL  ( 1.3s)  Shambhala found no genes shared by the input, P and Q
21_harmonizr   FAIL  ( 3.4s)  returned 0 genes
37_fabatch     FAIL  ( 0.8s)  Error in fabatch(x = t(exp_mat), y = y, batch = batches) : …
```

Two have established causes; two need diagnosis (§4 below).

### C — the image is 12.8 GB, documented as 6–9

```
/usr/local/lib/python3.11/site-packages/nvidia    2724 MB
/usr/local/lib/python3.11/site-packages/torch     1131 MB
/usr/local/lib/python3.11/site-packages/triton     691 MB
                                          total   4546 MB
```

**4.5 GB of CUDA and GPU runtime, in a CPU-only image for `c6a` nodes.** `combobatch` never
imports torch. It arrives from one line:

```
169.4 Collecting torch (from harmonypy<2.0,>=0.0.9->-r /tmp/requirements.txt (line 29))
```

`harmonypy>=0.0.9,<2.0` resolves to 0.2.0, which is the **only** release that hard-depends on
torch:

| harmonypy | requires |
|---|---|
| 0.1.0 | numpy, pandas, scikit-learn, scipy (`jax` optional) |
| **0.2.0** | numpy, pandas, scikit-learn, scipy, **torch** |
| 2.0.0 | numpy (torch dropped again; outside the `<2.0` pin) |

Removing it takes the image to roughly **8.3 GB**, inside the documented range.

### D — Octave prints an error and exits 0

```
octave: X11 DISPLAY environment variable not set
octave: disabling GUI features
error: ignoring const execution_exception& while preparing to exit
```

Cosmetic, and confirmed harmless: **exit code 0**, stdout correct, and the same three lines
appear from plain `octave --no-gui` — not the wrapper's doing. The Shambhala bridge gates on
`returncode != 0`, never on stderr text (`octave_bridge.py` lines 286, 515), so nothing
misreads it. It is noise that reads like a failure in every Shambhala run.

---

## Files to change

### 1. `combobatch/methods/rinterop.py` — the blocking one line

**Before** (lines 39–46):

```python
    import rpy2.robjects as ro

    try:
        result = ro.r(f"requireNamespace('{name}', quietly=TRUE)")
    except Exception:
        return False
    return bool(result[0])
```

**After:**

```python
    # rpy2's own API, not `ro.r("requireNamespace(...)")`: requireNamespace returns its
    # logical *invisibly*, and rpy2 3.6 maps an invisible result to None — for the
    # installed package and the missing one alike, so that form can only ever crash.
    # It also removes the last interpolation of a name into R source text in this file.
    from rpy2.robjects.packages import isinstalled

    try:
        return bool(isinstalled(name))
    except Exception:
        return False
```

**A deliberate non-change:** the crash is not moved inside the `try`. Catching it would make
every R package report "unavailable", turning a loud failure into 21 silently skipped
methods — the opposite of the `CLAUDE.md` "fail loudly" rule, and the same trade that
produced finding A.

### 2. `combobatch/methods/r_methods.py` — `27_dwd` must pass the arguments it needs

#### 2a. Supply `C` and `expon`

```r
                tryCatch({{
                    C   <- penaltyParameter(X_combined, y, expon = {r_literal(expon)})
                    sol <- genDWD(X = X_combined, y = y, C = C,
                                  expon = {r_literal(expon)})
```

`penaltyParameter` is DWDLargeR's own helper for `C`; both take the same `expon`.

#### 2b. Stop the handler from returning an uncorrected matrix

A per-batch failure must not look like success. The handler records the batch, and the
script fails if any batch was left uncorrected:

```r
            failed <- character(0)
            ...
                }}, error = function(e) {{
                    failed <<- c(failed, b)
                    message("DWD failed for batch '", b, "': ", conditionMessage(e))
                }})
            }}
            if (length(failed)) {{
                stop("27_dwd left ", length(failed), " batch(es) uncorrected: ",
                     paste(failed, collapse=", "))
            }}
```

#### 2c. Declare `expon` as a hyperparameter

Per `combobatch/methods/CLAUDE.md`: *"No hardcoded R-side arguments … Expose them as
`HyperParam`s."* `expon = 1.0` is DWDLargeR's standard DWD exponent. It goes in the
`MethodSpec` for `27_dwd`, and `docs/HYPERPARAMETERS.md` regenerates.

### 3. `requirements.txt` and `pyproject.toml` — drop 4.5 GB of GPU stack

```
harmonypy>=0.0.9,<0.2     # 11_harmony; <0.2 — 0.2.0 hard-depends on torch, which drags
                          # in nvidia/* and triton: 4.5 GB of CUDA in a CPU-only image.
                          # 0.1.0 has the same run_harmony API and self-describes as
                          # "NumPy, R package compatible".
```

Both files, identically — `tests/unit/test_requirements_match_pyproject.py` now enforces
that, and this is its first real use.

**Verified equivalent, not assumed.** Same input, same seed, both versions:

| harmonypy | backend | batch gap after |
|---|---|---|
| 0.1.0 | NumPy | 0.1253 |
| 0.2.0 | torch | 0.1249 |

A float-accumulation difference between backends, not a different algorithm — unlike the
reComBat `linear`/`elastic_net` case, where the plan chose to preserve behaviour precisely
because it *was* different.

### 4. `combobatch/cli/selftest_cmd.py` and the four failing methods

- **`20_shambhala` — a fixture limit, not a defect.** Shambhala matches the input against
  real gene symbols in the P and Q calibration panels; `make_synthetic()` emits `g0…g199`,
  so the intersection is empty and it cannot run by construction. It should **SKIP with that
  reason**, not FAIL. This adds a skip category beyond "missing backend"; `tests/CLAUDE.md`
  requires skips be honest, and "the synthetic fixture has no real gene symbols" is honest
  where "failed" is misleading. The alternative — give `make_synthetic()` real symbols drawn
  from the shipped panel — is better still and is the recommendation, with the skip as
  fallback if the panels do not overlap enough.
- **`19_tdm` — an interactive prompt with no stdin.** `if (answer %in% allowed) break` with
  `answer` of length zero is TDM asking a question and reading EOF. Fix at the call site by
  supplying the answer non-interactively; diagnosis step in the TODO.
- **`21_harmonizr` — "returned 0 genes"** and **`37_fabatch` — an `fabatch()` error whose
  message the summary truncates.** Both need one diagnostic run each inside the image before
  a fix can be written honestly. Listed as diagnosis items, not as fixes.

### 5. `docker/Dockerfile` — the `matlab` wrapper calls `octave-cli`

`octave-cli` removes the two X11 lines (measured); the `execution_exception` line survives
every variant and is an Octave 7.3 shutdown artefact. `QT_QPA_PLATFORM=offscreen` changes
nothing (measured). So: switch the binary, and document the remaining line rather than
pretend it can be suppressed.

```bash
exec octave-cli "${args[@]}"
```

### 6. Documentation — the size claim, in six places

`CLAUDE.md:251`, `docker/README.md:64`, `docker/BUILD_AND_PUSH.md:19,41,215`, and
`combobatch_tool_plan_260828.md:1050`. After §3 the image should be ~8.3 GB, so 6–9 GB
becomes true again — but the number must be **re-measured after the rebuild**, not assumed.
A troubleshooting row for the Octave line, and one for `TypeError: 'NoneType' object is not
subscriptable` pointing at the invisible-result trap.

### 7. `tests/unit/` — make the blind spot smaller

The unit suite could not have caught finding A: `r_available()` returns `False` on a laptop
without rpy2, so `r_package_available` returns before reaching the broken line. Two cheap
guards that do run there:

| Test | Guards against |
|---|---|
| stub `isinstalled` to return `True`/`False`, assert both pass through | the wrapper regressing to a form that returns `None` |
| assert `rinterop.py` contains no `requireNamespace` | the exact expression that cannot work coming back |
| assert `harmonypy` is pinned `<0.2` in both files | 4.5 GB of CUDA returning silently |
| an R-snippet guard: no `tryCatch(..., error = message(...))` that leaves a result unassigned | finding A's shape, in the other R methods |

The last is the valuable one and the hardest to phrase; a first cut is to require that every
`error = function(e)` handler in `r_methods.py` either re-raises or records into a variable
that is checked afterwards.

---

## Files that do NOT need to change

| File | Why not |
|---|---|
| `docker/install_r_packages.R`, the staged layers | All cached and correct; every R package the registries need is present, which the 3-passed integration run now proves |
| the reComBat install | Works: `[15/21] … 67.5s`, and its import check passed |
| `combobatch/vendor/shambhala/octave_bridge.py` | Gates on return code, not stderr, so the Octave noise cannot mislead it |
| `combobatch/imputation.py`, `selftest_cmd.py:239` | The other `ro.r()` call sites either discard the result or evaluate a visible expression — audited, all seven |
| `scripts/smoke_test.sh` | It reported exactly what happened; nothing to fix in the gate itself |

---

## Side effects and caveats

- **§3 invalidates the pip layer**, which is 62 minutes. §1, §2 and §4 are in the source
  layer and rebuild in seconds; §5 invalidates Octave and everything after. Order the work so
  the pip layer is rebuilt once.
- **`27_dwd`'s results change from "identity" to "corrected".** No published run used it —
  this is the first successful image — so nothing is invalidated. Its `hyperparams` entry
  also adds `expon` to the output key grammar when set away from the default.
- **`27_dwd` may now fail loudly** where it silently passed. That is the intent; a method
  that cannot run should not occupy a benchmark column.
- **The harmonypy downgrade changes `11_harmony` at float precision** (0.1253 vs 0.1249 on
  the fixture). Nothing has been published from either.
- **Only the smoke test is affected by §5** — the `matlab` name and flag-stripping behaviour
  are unchanged, which is what Shambhala2's internal `system("matlab …")` needs.
- **The four method failures are not all fixable in one pass.** Two are diagnosis items and
  may turn out to be fixture limits rather than defects.

---

## Verification commands

```bash
# 1. The rpy2 trap, in the image as built (None for installed AND missing)
docker run --rm --platform linux/amd64 combobatch:test python -c \
 "import rpy2.robjects as ro
from rpy2.robjects.packages import isinstalled
print('bare :', ro.r(\"requireNamespace('limma', quietly=TRUE)\"))
print('isinstalled:', isinstalled('limma'), isinstalled('nosuchpkg'))"

# 2. genDWD really does require C and expon
docker run --rm --platform linux/amd64 combobatch:test Rscript -e \
 'library(DWDLargeR); f <- formals(genDWD);
  print(names(f)[sapply(f, function(a) is.name(a) && !nzchar(as.character(a)))])'

# 3. Laptop suite and lint after the edits
python -m pytest tests/unit -q
black --check .
python scripts/generate_docs.py --check      # 27_dwd gains a hyperparameter
docker buildx build --check -f docker/Dockerfile .

# 4. Rebuild, VPN off
bash scripts/smoke_test.sh

# 5. Re-measure the size claim rather than assuming it
docker images combobatch:test --format '{{.Size}}'
docker run --rm --platform linux/amd64 combobatch:test \
    du -sm /usr/local/lib/python3.11/site-packages | tail -1

# 6. 27_dwd must now change its input
docker run --rm --platform linux/amd64 combobatch:test python -c \
 "import numpy as np, pandas as pd
from combobatch.methods.r_methods import normalize_dwd
rng = np.random.default_rng(0)
exp = pd.DataFrame(rng.normal(size=(30, 40)) + 5, index=[f's{i}' for i in range(30)],
                   columns=[f'g{j}' for j in range(40)])
exp.iloc[15:] += 3
ann = pd.DataFrame({'b': ['a']*15 + ['b']*15}, index=exp.index)
out = normalize_dwd(exp, ann, batch_col='b')
print('changed:', not np.allclose(out.values, exp.values))"
```

---

## What the diagnosis items turned out to be

Each of the three was bigger or different than the plan assumed.

**`19_tdm` — a missing package, not a prompt to answer.** TDM's `load_it()` helper calls
`BiocManager::install()` for anything missing *at call time*; `binr` is absent from the
image, so it tried to install mid-run, prompted, and read EOF. Declaring `binr` in
`19_tdm`'s `r_packages` fixes it at build time, and the registry cross-check carries it
into `install_r_packages.R` automatically. Until the image is rebuilt it now reports
**`SKIP — needs R package binr`**, which is the honest state rather than a confusing R
error.

**`37_fabatch` — not a fabatch problem, and not confined to `37_fabatch`.** The annotation
was read without `check.names=FALSE`, so R renamed `__target__` to `X__target__`. The
lookup then returned an **empty vector instead of raising**, and `fabatch` complained
about `y` several steps later. **Eight** annotation reads had the same omission, and
twelve lookups index a column by name — so *any* user column that is not a syntactic R
name (`Cell type`, `Diagnosis (WHO)`, `1st_batch`) fails the same silent way. All eight
fixed; `37_fabatch` now **passes**.

**`27_dwd` needed three fixes, not one.** Supplying `C` and `expon` was necessary and not
sufficient:

1. `genDWD` reaches into its argument with `@`, so a base R matrix dies with *"no
   applicable method for `@`"*. It requires an S4 `Matrix`.
2. The result's `w` is the separating direction, one entry per gene; `beta` is the scalar
   intercept. The code preferred `beta`, which would have made a length-1 "direction" and
   a uniform shift.

With all three, `27_dwd` returns a finite 80×200 matrix that **differs from its input** —
the method working for the first time.

**Still open:** `20_shambhala` moved from *"no genes shared with P and Q"* (the fixture
had placeholder names) to *"returned an all-NaN matrix"*, and `21_harmonizr` still returns
0 genes. Both are real and neither is diagnosed; they belong in their own plan.

Selftest inside the image, with the fixed source mounted: **35 passed, 1 skipped, 2 failed**
— against 34 passed / 0 skipped / 4 failed before, and a hard crash before that.

---

## TODO

- [x] `combobatch/methods/rinterop.py` — `isinstalled` in place of the `requireNamespace`
      round trip
- [x] `tests/unit/test_rinterop.py` — stub-based test that both outcomes survive, and a
      source guard against `requireNamespace`
- [x] `combobatch/methods/r_methods.py` — pass `C` (via `penaltyParameter`) and `expon` to
      `genDWD`
- [x] `combobatch/methods/r_methods.py` — *(beyond the plan)* convert to a sparse `Matrix`,
      and take `sol$w` rather than `sol$beta`
- [x] `combobatch/methods/r_methods.py` — the DWD error handler must fail the run rather
      than return an uncorrected matrix
- [x] `combobatch/methods/__init__.py` — declare `expon` as a `27_dwd` hyperparameter
- [x] `combobatch/methods/r_methods.py` — update the stale `penalty` note in the docstring
- [x] `docs/HYPERPARAMETERS.md` and `docs/METHODS.md` — regenerated
- [x] `requirements.txt` + `pyproject.toml` — `harmonypy>=0.0.9,<0.2`, with the reason
- [x] `tests/unit/test_docker_assets.py` — the harmonypy ceiling *(covered by the existing
      requirements/pyproject drift test, which now compares the identical `<0.2` pin)*
- [x] `docker/Dockerfile` — `matlab` wrapper execs `octave-cli`
- [x] `combobatch/cli/selftest_cmd.py` — `make_synthetic()` draws real gene symbols from
      the P∩Q calibration panels, with a placeholder fallback
- [x] **`19_tdm`** — `binr` declared in the spec and installed in the `cran` stage
- [x] **`37_fabatch`** — `check.names=FALSE` on all eight annotation reads
- [ ] **`21_harmonizr`** — still returns 0 genes; undiagnosed, needs its own plan
- [ ] **`20_shambhala`** — now returns an all-NaN matrix; undiagnosed, needs its own plan
- [x] `tests/unit/test_r_snippets.py` — new: no error handler may only print; every
      annotation read disables name mangling; `genDWD` gets all four arguments
- [x] Documentation — two `CLAUDE.md` files, three troubleshooting rows, and a
      `site-packages` size check in step 4
- [x] `combobatch_tool_plan_260828.md` — recorded under Phase 6
- [x] Run verification steps 1–3 — 867 unit tests pass, 18 e2e, Black clean, docs in sync,
      `buildx --check` clean
- [ ] **VPN off**, rebuild and run `bash scripts/smoke_test.sh` to a clean pass
- [ ] Re-measure the image size after the rebuild and correct the six 6–9 GB claims if
      needed *(deferred deliberately: the plan says measure, not assume)*
