# Plan — the R package layer fails on a missing `libuv1-dev`

**Status:** awaiting approval. No code has been changed.
**Date:** 2026-08-30
**Trigger:** `bash scripts/smoke_test.sh`, log at `docker/docker_build_logs_260829.txt`.
Failed at step 6/15 after **2133.7 s (35.6 min)**; the whole layer was discarded.

---

## Overview

The TLS fix worked. Steps 3 and 4 — the preflight probe and the R download that failed
yesterday — passed in 1.0 s and 35.7 s. The build got 36 minutes further and then died in
the R package layer with:

```
Error: Failed to install 1 required package(s): qsmooth
```

**`qsmooth` is not the problem. `fs` is.** The R package `fs` 2.1.0 requires a system
`libuv`, and `libuv1-dev` is not in the Dockerfile's apt layer. `fs` said so itself, in the
words of its own configure script:

```
Package libuv was not found in the pkg-config search path.
Configuration failed because libuv was not found. Try installing:
 * deb: libuv1-dev (Debian, Ubuntu, etc)
<stdin>:1:10: fatal error: uv.h: No such file or directory
ERROR: configuration failed for package ‘fs’
```

Everything downstream of `fs` then fell over in a chain that ends at the one package the
build actually verifies:

```
fs  →  sass  →  bslib  →  rmarkdown  →  htmlTable / htmlwidgets  →  Hmisc  →  qsmooth
```

`Hmisc` Imports `rmarkdown`, and `qsmooth` Imports `Hmisc` — these are hard `Imports`
edges, not optional ones, so no amount of dependency trimming avoids `fs`. **`libuv1-dev`
is a genuine, unavoidable requirement of this image.** It is also the only one missing:
`qsmooth` was the sole entry of the 20-name `required` vector left unsatisfied.

Three further defects made a one-package problem cost 36 minutes and read as something
else. Each is addressed below:

- **`dependencies = TRUE` pulls `Suggests`.** The build compiled `V8`, `shiny`, `rgl`,
  `plotly`, `gsl`, `sodium`, `gt`, `kableExtra`, `DHARMa`, `semEff`, `qgam`, `testthat`,
  `reactable`, `sparkline`, `juicyjuice` and `safer` — none of which any registered method,
  imputer or metric group uses, several of which need system libraries the image also does
  not have, and one of which (`V8`) burned about six minutes before failing.
- **The failures were invisible.** `install_or_warn` traps *errors*; `install.packages`
  reports a package that would not build as a **warning**. Three such warnings fired,
  covering 41 package-failure events — and the log contains **not one `WARN:` line**. The
  script's stated promise, that the log "names *all* failures, not just the first," is not
  delivered by the mechanism that claims it.
- **No fail-fast.** `Hmisc` was known dead at **886.7 s**. The build then spent another
  **1247 s** installing Bioconductor and five GitHub packages whose fate was already sealed,
  before the verification block finally reported it.

---

## Background — what the log actually shows

Section-by-section, by the log's own elapsed-second stamps. (The file is stream-interleaved:
R's stdout and stderr flush separately, so line order is not chronological. Timestamps are.)

| Elapsed | Step | Outcome |
|---|---|---|
| 26.6 s | `=== installing compiled prerequisites (quantreg, Hmisc, lme4) ===` | **27 packages failed**, incl. `fs`, `Hmisc`, `V8`, `sodium`, `gsl`, `rgl`, `shiny` |
| 886.7 s | ↳ warning printed | first appearance of the fatal chain — nothing marked it |
| 887.4 s | CRAN packages (`missForest`, `softImpute`, `DWDLargeR`, `huge`) | **6 packages failed**: `fs`, `sass`, `pkgload`, `bslib`, `testthat`, `rmarkdown`. The four *named* packages installed |
| 1358.4 s | `bapred` + Bioconductor deps | `DONE (bapred)` |
| 1602.9 s | Bioconductor packages | **8 packages failed**: `fs`, `sass`, `bslib`, `rmarkdown`, `htmlwidgets`, `htmlTable`, `Hmisc`, **`qsmooth`** |
| 2122.7 s | 5 GitHub packages | all five `DONE` |
| 2133.2 s | consolidated verification | `Error: Failed to install 1 required package(s): qsmooth` |

`fs` is attempted three separate times and fails identically each time, because nothing
between the attempts changes the missing system library.

The `Suggests` tail, timed from the log: `V8` downloaded at 30.1 s and failed its `make` at
390.0 s. `shiny`, `rgl`, `plotly`, `gsl`, `sodium`, `kableExtra`, `DHARMa`, `semEff`,
`qgam`, `reactable`, `sparkline`, `juicyjuice`, `safer` and `testthat` were all fetched
between 30 s and 62 s and compiled through to 886 s. **That step alone ran 860 s**, nearly
all of it on packages nothing in ComboBatch imports.

### Verified fix

`libuv1-dev` in Debian 12 provides exactly the two things `fs`'s configure looked for:

```
pkg-config --exists libuv     -> yes
pkg-config --modversion libuv -> 1.44.2
/usr/include/uv.h present
```

---

## Files to change

### 1. `docker/Dockerfile` — the apt layer

#### 1a. Add `libuv1-dev`

The library list is grouped by what needs it, and the comment above it already explains why
`libpng-dev` and `zlib1g-dev` are load-bearing. `libuv1-dev` joins them for the same reason
and should be documented the same way.

**Before** (lines 32–41, excerpt):

```dockerfile
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libhdf5-dev libbz2-dev liblzma-dev zlib1g-dev \
        libblas-dev liblapack-dev \
```

**After:**

```dockerfile
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libhdf5-dev libbz2-dev liblzma-dev zlib1g-dev libuv1-dev \
        libblas-dev liblapack-dev \
```

#### 1b. Extend the comment above the apt layer

**Before:**

```dockerfile
# libpng-dev and zlib1g-dev are load-bearing: without them Hmisc fails to compile and
# BiocManager reports qsmooth as skipped rather than failed. cmake is needed by the
# build chain of the metrics dependencies.
```

**After:**

```dockerfile
# libpng-dev, zlib1g-dev and libuv1-dev are load-bearing, all three for the same reason:
# qsmooth Imports Hmisc, Hmisc Imports rmarkdown, and that chain runs
# rmarkdown -> bslib -> sass -> fs, where fs 2.x requires a system libuv. Without
# libuv1-dev the 2026-08-29 build spent 36 minutes and died reporting a single missing
# package, qsmooth, five dependency levels away from the actual cause. cmake is needed by
# the build chain of the metrics dependencies.
```

### 2. `docker/install_r_packages.R`

#### 2a. Drop `dependencies = TRUE` from the prerequisites call

R's default is `NA` — `Depends`, `Imports` and `LinkingTo`, which is every edge that can
actually break a load. `TRUE` adds `Suggests`, which is documentation, examples and test
scaffolding.

**Before:**

```r
install_or_warn("compiled prerequisites (quantreg, Hmisc, lme4)", install.packages(
    c("quantreg", "Hmisc", "lme4"),
    dependencies = TRUE
))
```

**After:**

```r
# No `dependencies = TRUE`: that adds Suggests, and on 2026-08-29 it pulled V8, shiny,
# rgl, plotly, gsl, sodium, gt, kableExtra, DHARMa, semEff, qgam, testthat, reactable,
# sparkline, juicyjuice and safer into the image - 860 seconds, five more system
# libraries required, and nothing in the registries imports any of them. The default
# (Depends, Imports, LinkingTo) is every edge that can break a library() call.
install_or_warn(
    "compiled prerequisites (quantreg, Hmisc, lme4)",
    install.packages(c("quantreg", "Hmisc", "lme4"))
)
```

#### 2b. Drop `dependencies = TRUE` from the CRAN call

**Before:**

```r
install_or_warn("CRAN packages", install.packages(
    c(
        "missForest",   # imputer: missforest
        ...
    ),
    dependencies = TRUE
))
```

**After:** the same list, with the `dependencies = TRUE` argument removed.

#### 2c. Report failed installs, which are warnings and not errors

`install.packages` does not stop on a package that fails to build; it emits
`installation of N packages failed: …` as a warning and returns normally. So
`install_or_warn`'s `tryCatch(error = …)` never fires, and the `WARN:` marker it exists to
print has never once appeared in a build log.

**Before:**

```r
install_or_warn <- function(label, expr) {
    cat(sprintf("\n=== installing %s ===\n", label))
    tryCatch(
        force(expr),
        error = function(e) cat(sprintf("WARN: %s failed: %s\n", label, conditionMessage(e)))
    )
    invisible(NULL)
}
```

**After:**

```r
install_or_warn <- function(label, expr) {
    cat(sprintf("\n=== installing %s ===\n", label))
    # Both arms are needed. install.packages() reports a package that would not build as a
    # *warning* and returns normally, so the error handler alone stayed silent through 41
    # package failures on 2026-08-29; BiocManager::install() can also stop() outright.
    withCallingHandlers(
        tryCatch(
            force(expr),
            error = function(e) {
                cat(sprintf("WARN: %s failed: %s\n", label, conditionMessage(e)))
            }
        ),
        warning = function(w) {
            cat(sprintf("WARN: %s: %s\n", label, conditionMessage(w)))
            invokeRestart("muffleWarning")
        }
    )
    invisible(NULL)
}
```

#### 2d. Fail fast on the prerequisites

`Hmisc` is not in the `required` vector — it is a dependency of `qsmooth`, not a method
backend — so its failure is only discovered 1247 seconds later, wearing `qsmooth`'s name.
A checkpoint immediately after the step that installs it turns a 36-minute build into a
15-minute one and names the package that actually failed.

Inserted directly after the prerequisites `install_or_warn` call:

```r
# Checked here rather than in the final block: none of these is a method backend, so none
# is in `required`, yet qsmooth Imports Hmisc and cannot install without it. Failing here
# reports the package that actually broke, 20 minutes before the consolidated check would
# report a symptom five dependency levels downstream.
prerequisites <- c("quantreg", "Hmisc", "lme4")
absent <- prerequisites[!vapply(prerequisites, requireNamespace, logical(1), quietly = TRUE)]
if (length(absent)) {
    stop(
        "prerequisite package(s) failed to install: ", paste(absent, collapse = ", "),
        ". Scroll up for the compiler output - a missing -dev system library in the ",
        "Dockerfile's apt layer is the usual cause."
    )
}
```

#### 2e. Point the final message at the log

The consolidated `stop()` names the missing package but not what to do with it. One clause,
since the answer was five levels up the log and 36 minutes back in time:

```r
    stop(
        "Failed to install ", length(missing), " required package(s): ",
        paste(missing, collapse = ", "),
        ". Search the log above for 'WARN:' and for 'ERROR: configuration failed' - a ",
        "package listed here is often the last casualty of a dependency that failed ",
        "earlier for want of a system library."
    )
```

### 3. `tests/unit/test_docker_assets.py` — three assertions

| Test | Guards against |
|---|---|
| `libuv1-dev` is in the Dockerfile's apt layer | the 36-minute failure returning; it is as load-bearing as `libpng-dev` |
| `dependencies = TRUE` appears nowhere in `install_r_packages.R` | `Suggests` creeping back, with its five extra system libraries |
| the R script checks the prerequisites before the Bioconductor block | the fail-fast checkpoint being removed, restoring the 20-minute delay |

The existing `TestRPackagesMatchTheRegistries` tests are unaffected: `fs`, `Hmisc` and the
rest are transitive dependencies, not entries in the `required` vector, and the registry
comparison is unchanged.

### 4. `docker/BUILD_AND_PUSH.md`

- **Troubleshooting** — add a row for `ERROR: configuration failed for package ‘fs’` →
  `libuv1-dev` missing from the apt layer.
- **Troubleshooting** — the existing row *"Image builds, but `14_qsmooth` skips at run
  time"* describes something that cannot now happen, and its stated remedy sent the reader
  to the wrong library. `qsmooth` is in the `required` vector, so the build goes red instead
  of shipping a silently degraded image — which is exactly what happened on 2026-08-29.
  Rewrite it around the real symptom: `Failed to install N required package(s)`, with
  instructions to search the log for `WARN:` and `configuration failed`.
- **Step 3** — note that the R layer is now materially faster, since `Suggests` no longer
  compiles.

### 5. `docker/CLAUDE.md` — one invariant

*Never install R packages with `dependencies = TRUE`.* It adds `Suggests`, which no
registered method needs, and each one is a new system library the apt layer does not have.
Every package in the image should be traceable to a registry entry or to a hard `Imports`
edge of one.

### 6. `combobatch_tool_plan_260828.md`

One line under Phase 6's unchecked build item recording the second failure and this plan,
next to the reference to the TLS one.

---

## Files that do NOT need to change

| File | Why not |
|---|---|
| the `required <- c(...)` vector | It did its job: it caught a package silently dropped from the image, which is precisely the donor failure mode it exists to prevent. Nothing about it changes. |
| `combobatch/methods/*`, the registries | No method's declared `r_packages` are wrong; `qsmooth` was unbuildable, not misdeclared. |
| `requirements.txt` | The Python layer never ran. Untouched by this failure. |
| `docker/entrypoint.sh`, `scripts/build_and_push_image.sh` | Not reached. |
| `docker/tls_inspection_ca_plan_260829.md` | Its fix worked — steps 3 and 4 passed. Only its final TODO item, the completed smoke test, is still open. |
| `.github/workflows/ci.yml` | Builds the same Dockerfile; it inherits the fix with no edit. |

---

## Side effects and caveats

- **The R layer rebuilds from scratch.** Changing the apt layer invalidates everything
  below it, including the 35-minute R install. Unavoidable — the missing library is *in*
  the apt layer.
- **The image should get smaller and the build shorter.** Sixteen `Suggests` packages and
  their dependency closures no longer compile. The 860-second prerequisites step was
  almost entirely `Suggests`; expect a large part of that back, though the Bioconductor
  closure is unchanged and still dominates.
- **`dependencies = TRUE` was load-bearing for nothing, but verify rather than assume.**
  The `required` vector is checked at the end of the script and
  `tests/integration/test_methods_all_backends.py` asserts zero method skips inside the
  image. If dropping `Suggests` removed something a method needs at run time, those two
  checks are what will say so.
- **`fs` 2.x is the drift here.** `fs` 1.6.x bundled libuv and built anywhere; 2.1.0
  requires the system library, with `USE_BUNDLED_LIBUV=1` as an escape hatch. That
  environment variable is the fallback if `libuv1-dev` ever becomes unavailable — noted in
  the Dockerfile comment, not used, because installing the real library is cleaner.
- **The VPN must still be disconnected**, for the whole build. Nothing in this plan changes
  that; see `tls_inspection_ca_plan_260829.md`.

---

## Verification commands

```bash
# 1. The missing library, confirmed present in Debian 12 (already run, kept for the record)
docker run --rm python:3.11-slim-bookworm bash -c \
    'apt-get update -qq >/dev/null 2>&1; \
     apt-get install -y -qq --no-install-recommends pkg-config libuv1-dev >/dev/null 2>&1; \
     pkg-config --modversion libuv; ls /usr/include/uv.h'
# expect: 1.44.2 and /usr/include/uv.h

# 2. The R script still parses, and the registries still agree with it
Rscript -e 'invisible(parse("docker/install_r_packages.R")); cat("parses\n")'   # if R is local
python -m pytest tests/unit/test_docker_assets.py -v

# 3. Whole laptop suite + formatting
python -m pytest tests/unit -q
black --check .
docker buildx build --check -f docker/Dockerfile .

# 4. The build itself, VPN disconnected
bash scripts/smoke_test.sh

# 5. After it passes, confirm the package that failed is really there
docker run --rm --platform linux/amd64 combobatch:test \
    Rscript -e 'library(qsmooth); library(Hmisc); cat("qsmooth OK\n")'

# 6. And that the Suggests tail is genuinely gone
docker run --rm --platform linux/amd64 combobatch:test \
    Rscript -e 'cat(c("V8","shiny","rgl","plotly") %in% rownames(installed.packages()), "\n")'
# expect: FALSE FALSE FALSE FALSE
```

---

## Two deviations found while implementing

**A stale comment in `install_r_packages.R` said the opposite of what happened.** Directly
above the prerequisites call it claimed that when `Hmisc` fails, "14_qsmooth then vanishes
from the image without the build ever going red." The build *did* go red — `qsmooth` is in
the `required` vector, which is the whole point of that vector. The plan flagged the same
false claim in `BUILD_AND_PUSH.md` but not this copy of it; corrected in both.

**The prerequisites check was reflowed.** `absent <- prerequisites[!vapply(…)]` came to 93
characters on one line, against a file whose comparable line is 84 and a project limit of
88. Split across three lines.

---

## TODO

- [x] `docker/Dockerfile` — add `libuv1-dev` to the apt layer
- [x] `docker/Dockerfile` — extend the load-bearing-libraries comment to explain the
      `fs → … → qsmooth` chain
- [x] `install_r_packages.R` — drop `dependencies = TRUE` from the prerequisites call
- [x] `install_r_packages.R` — drop `dependencies = TRUE` from the CRAN call
- [x] `install_r_packages.R` — `install_or_warn` must trap warnings as well as errors
- [x] `install_r_packages.R` — fail-fast checkpoint after the prerequisites
- [x] `install_r_packages.R` — point the final `stop()` at `WARN:` / `configuration failed`
- [x] `install_r_packages.R` — *(deviation)* correct the stale "never going red" comment
- [x] `tests/unit/test_docker_assets.py` — four assertions (three planned, plus the
      `withCallingHandlers` guard)
- [x] `docker/BUILD_AND_PUSH.md` — troubleshooting row for `configuration failed for ‘fs’`
- [x] `docker/BUILD_AND_PUSH.md` — rewrite the `14_qsmooth` row around the real symptom
- [x] `docker/BUILD_AND_PUSH.md` — note the faster R layer in step 3
- [x] `docker/CLAUDE.md` — the "never `dependencies = TRUE`" invariant
- [x] `combobatch_tool_plan_260828.md` — record the second failure under Phase 6
- [x] Run verification steps 2–3 (laptop, minutes) — R script parses (18 top-level
      expressions) and both `install_or_warn` arms print `WARN:`
- [ ] **VPN disconnected**, run `bash scripts/smoke_test.sh` to completion
- [ ] Run verification steps 5–6 against the built image *(needs the built image)*
