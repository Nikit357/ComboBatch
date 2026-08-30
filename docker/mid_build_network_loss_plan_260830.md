# Plan — the VPN reconnected mid-build, and the R layer had nothing to fall back on

**Status:** awaiting approval. No code has been changed.
**Date:** 2026-08-30
**Trigger:** `bash scripts/smoke_test.sh`, log at `docker/docker_build_logs_260830.txt`.
Failed at step 6/15 after **971.9 s**, having done **857 s of successful work** that was
then discarded.

---

## Overview

**Yesterday's fixes worked. This is a different failure.** `libuv1-dev` did its job — `fs`
compiled at **t=131.3 s**, the prerequisite checkpoint stayed silent, and `Suggests` no
longer drags in V8 and friends. The `WARN:` handler worked too: 273 labelled warnings,
which is why this took minutes to diagnose instead of an hour.

What broke is the network, **partway through the build**:

```
856.7  trying URL 'https://cloud.r-project.org/src/contrib/bapred_1.1.tar.gz'
861.8  Error in download.file(...) : cannot download any files
       …
       URL '…': status was 'SSL peer certificate or SSH remote key was not OK'   × 185
```

That message is libcurl's `CURLE_PEER_FAILED_VERIFICATION` — exit 60, the same class as
Friday's `curl: (60)`. **TLS interception came back at ~857 s.** Everything before it
worked; after it, **227 download attempts, zero successes**, and the run ended reporting 16
of the 21 required packages missing.

The proof is still on the machine as I write this:

```
default   100.65.144.8   UGScg   utun4       ← FortiClient, up again
cloud.r-project.org → issuer /O=BostonGene Corporation/CN=BGSSLInspection
```

Yesterday the tunnel gateway was `100.65.144.3`; it is now `100.65.144.8`. **A different
address means a new session — FortiClient reconnected**, it was not simply left on. Whether
that was an auto-reconnect policy or a manual reconnect while the build ran, the effect is
the same: the build lost verifiable TLS 14 minutes in.

### So why did the preflight probe not catch it?

Because **I built it to answer the wrong question, and this is the flaw in my own fix.**
The probe checks TLS *once*, at the start:

```
[ 3/15] RUN ... curl -sS -o /dev/null https://cdn.posit.co/     1.1s   ✓ passed
```

It passed, honestly, at t=1.1 s. The R layer then needs verifiable TLS continuously for the
next 15–30 minutes. A point-in-time probe cannot promise that, and by passing it created
exactly the false confidence that made this failure surprising. The probe is still worth
keeping — it catches the build you start while already on the VPN, cheaply — but it must be
described as what it is, and the build must stop depending on an uninterrupted network.

### The three things to fix

1. **Nothing is resumable.** This is the third failed build, and the third time every
   completed minute was thrown away. 857 s of correct work — `BiocManager`, `remotes`,
   `quantreg`, `Hmisc`, `lme4`, `fs`, `missForest`, `softImpute`, `DWDLargeR`, `huge` and
   their closures — was discarded because it lived in the same `RUN` as the failure.
2. **The diagnosis misdirects, again.** The final message ends: *"a package listed here is
   often the last casualty of a dependency that failed earlier for want of a system
   library."* I wrote that yesterday, and yesterday it was right. Today the cause is TLS,
   there are 185 `download of package … failed` warnings saying so, and the message still
   points at apt. This is the same defect as the architecture message: **a catch-all that
   names one cause for a failure it has not identified.**
3. **The log is unreadable while it matters.** Every `WARN:` line carries the timestamp
   `971.9` — the moment R exited — because `cat()` writes to stdout, which is block-buffered
   when it is not a terminal. For 16 minutes the operator saw nothing, and afterwards the
   file is interleaved rather than chronological.

---

## Background — the measured timeline

From the log's own elapsed-second stamps on **stderr** lines, which are unbuffered and
therefore truthful. (Stdout lines all read `971.9`; see fix 4.)

| Elapsed | Event |
|---|---|
| 1.1 s | preflight TLS probe passes |
| 90.5 s | R 4.5.3 `.deb` downloads and installs |
| 18–28 s | `BiocManager`, `BiocVersion`, `remotes` |
| 28.8 s | prerequisites start — 58 dependencies, not the previous 100+ |
| **131.3 s** | **`* DONE (fs)`** — yesterday's `libuv1-dev` fix, confirmed working |
| 367.1 s | CRAN packages start (`missForest`, `softImpute`, `DWDLargeR`, `huge`) |
| 848.6 s | `bapred` + Bioconductor dependencies start |
| 856.7 s | last download that had a chance: `bapred_1.1.tar.gz` |
| **861.8 s** | **`cannot download any files`** — TLS gone |
| 862.3 s → 962.9 s | 227 further attempts, **all** failing certificate verification |
| 971.7 s | `Failed to install 16 required package(s)` |

Nothing succeeded after 861.8 s. The cut is total, not intermittent — consistent with a
tunnel coming up and staying up, and inconsistent with a flaky link.

**What the failure cost, by stage:** the first four stages of the script (`BiocManager`,
prerequisites, CRAN, and the start of `bapred`) represent ~14 minutes of compiles that were
complete and correct. All of it is gone.

**Which packages survived:** 5 of 21 — `missForest`, `softImpute`, `DWDLargeR`, `huge`,
`TDM`. The other 16 are everything from `bapred` onward.

---

## Files to change

### 1. `docker/install_r_packages.R` — stages, so a failure costs one stage

#### 1a. Derive `required` from a stage partition

Today `required <- c(...)` is a hand-listed vector at the bottom. It becomes the *union of
what each stage is responsible for*, so there is still exactly one list, and each stage can
verify its own share.

Inserted near the top, replacing the standalone `required` vector at the bottom:

```r
# What each stage installs, and the single source of the final verification vector.
# The install is split into stages because it is not resumable otherwise: on 2026-08-30 a
# VPN reconnected 14 minutes in and every completed compile was discarded with the failure.
# Each stage is its own Docker layer, so a stage that succeeded stays cached and a rebuild
# resumes at the one that did not.
STAGE_PACKAGES <- list(
    prereq = character(0),   # quantreg/Hmisc/lme4: no method's backend, checked separately
    cran   = c("missForest", "softImpute", "DWDLargeR", "huge"),
    bapred = "bapred",
    bioc   = c(
        "DESeq2", "Harman", "NOISeq", "RUVSeq", "batchelor", "edgeR",
        "limma", "qsmooth", "ruv", "sva", "variancePartition"
    ),
    github = c("AMDBNorm", "DBNorm", "FSQN", "HarmonizR", "TDM")
)
required <- unlist(STAGE_PACKAGES, use.names = FALSE)
```

The partition is exhaustive: 4 + 1 + 11 + 5 = **21**, which is the current length of
`required`. `tests/unit/test_docker_assets.py` continues to compare that set against the
registries, so the existing guarantee is unchanged.

#### 1b. Take a `--stage` argument and run only that stage

```r
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 2 && args[[1]] == "--stage") args[[2]] else "all"
valid <- c(names(STAGE_PACKAGES), "verify", "all")
if (!stage %in% valid) {
    stop("unknown --stage ", stage, "; expected one of: ", paste(valid, collapse = ", "))
}
```

Each existing `install_or_warn(...)` block is wrapped in `if (stage %in% c("<name>", "all"))`.
`all` preserves today's behaviour, so the script is still runnable by hand in one go.

`BiocManager` and `remotes` are bootstrapped at the top of every invocation, guarded by
`requireNamespace` — they persist in the image across layers, so this is a no-op after the
first stage.

#### 1c. Each stage verifies its own packages before the layer is cached

This is what makes staging safe. Without it a stage whose installs failed would still exit
0, and BuildKit would cache a layer that is missing packages — worse than the current
behaviour, because the rebuild would skip it.

```r
# Verified per stage, not only at the end: a stage that exits 0 gets cached, so a stage
# must not exit 0 with its own packages missing.
verify_stage <- function(stage_name) {
    expected <- STAGE_PACKAGES[[stage_name]]
    if (!length(expected)) return(invisible(NULL))
    absent <- expected[!vapply(expected, requireNamespace, logical(1), quietly = TRUE)]
    if (length(absent)) stop_with_diagnosis(absent, stage_name)
    cat(sprintf("stage %s: all %d package(s) present\n", stage_name, length(expected)))
}
```

The `verify` stage re-checks the whole `required` vector, unchanged in spirit from today's
consolidated block.

#### 1d. Diagnose the failure that happened, not the one that happened yesterday

`install_or_warn` already sees every warning. Recording them lets the final message tell
a download failure from a build failure instead of asserting one of them.

```r
# Recorded rather than only printed, so the failure message can name the actual cause.
# Yesterday's message blamed a missing system library; today's failure was TLS, and 185
# warnings said so while the message still pointed at apt.
FAILURES <- new.env(parent = emptyenv())
FAILURES$messages <- character(0)

stop_with_diagnosis <- function(absent, stage_name) {
    downloads_failed <- any(grepl(
        "SSL peer certificate|cannot download any files|download of package|unable to access",
        FAILURES$messages
    ))
    if (downloads_failed) {
        stop(
            "stage ", stage_name, " could not install: ", paste(absent, collapse = ", "),
            ". Downloads failed certificate verification partway through this build - a ",
            "VPN or proxy re-established TLS interception while it was running. The ",
            "preflight probe only checks the moment the build starts. Disconnect it, ",
            "confirm it cannot reconnect, and rebuild: the stages that already ",
            "succeeded are cached and will not be repeated."
        )
    }
    stop(
        "stage ", stage_name, " could not install: ", paste(absent, collapse = ", "),
        ". Search the log for 'WARN:' and 'ERROR: configuration failed' - a package ",
        "listed here is often the last casualty of a dependency that failed earlier for ",
        "want of a system library."
    )
}
```

#### 1e. Send `WARN:` to stderr so the log is live and chronological

One-character-class change with a disproportionate effect: stdout is block-buffered when it
is not a terminal, which is why all 273 warnings carry the exit timestamp and the operator
watched a silent terminal for 16 minutes.

```r
        warning = function(w) {
            cat(sprintf("WARN: %s: %s\n", label, conditionMessage(w)), file = stderr())
            FAILURES$messages <- c(FAILURES$messages, conditionMessage(w))
            invokeRestart("muffleWarning")
        }
```

The same for the error arm and for the `=== installing … ===` banners.

#### 1f. Raise the download timeout

R's default is 60 seconds per file. Bioconductor tarballs over a slow link can exceed it,
and the resulting message is indistinguishable from the real network failure above.
Precautionary, not diagnosed here:

```r
options(timeout = 300)
```

### 2. `docker/Dockerfile` — one layer per stage

**Before:**

```dockerfile
COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R && rm -f /tmp/install_r_packages.R
```

**After:**

```dockerfile
# Six layers, not one. The R install takes 15-30 minutes and needs verifiable TLS for all
# of it; on 2026-08-30 a VPN reconnected 14 minutes in and a single RUN meant every
# finished compile was discarded. Each stage verifies its own packages before exiting, so
# a cached layer is a layer that actually installed what it claims.
COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R --stage prereq
RUN Rscript /tmp/install_r_packages.R --stage cran
RUN Rscript /tmp/install_r_packages.R --stage bapred
RUN Rscript /tmp/install_r_packages.R --stage bioc
RUN Rscript /tmp/install_r_packages.R --stage github
RUN Rscript /tmp/install_r_packages.R --stage verify \
    && rm -f /tmp/install_r_packages.R
```

### 3. `tests/unit/test_docker_assets.py`

`verified_packages()` currently parses `^required <- c\(…\)`. That vector becomes derived,
so the helper parses `STAGE_PACKAGES` instead — the same set, from the new single source.

| Test | Guards against |
|---|---|
| `verified_packages()` reads `STAGE_PACKAGES` and still equals the registry union, both ways | the existing guarantee silently lapsing when the vector moves |
| every stage named in `STAGE_PACKAGES` has a `RUN … --stage <name>` line, and `verify` is last | a stage added to the script but never run by the image |
| the Dockerfile runs the stages as separate `RUN`s | the layer being collapsed back into one, restoring the all-or-nothing failure |
| `WARN:` is written to `stderr()` | the log going back to being buffered and unreadable |
| the failure message branches on download failures | the diagnosis hard-coding one cause again |

### 4. `docker/BUILD_AND_PUSH.md`

- **Step 0b** — state plainly that the probe is point-in-time: it checks the moment the
  build starts, and the build needs verifiable TLS for its whole 15–30 minutes. Add the
  check that FortiClient will not reconnect on its own, and a mid-build spot check:
  `netstat -rn -f inet | head -5` should show `en0`, not `utun4`, holding the default route.
- **Troubleshooting** — new row for `status was 'SSL peer certificate or SSH remote key was
  not OK'` → the VPN reconnected mid-build; disconnect, confirm it stays down, rebuild, and
  note that completed stages are cached so the rebuild resumes.
- **Step 3** — record that the R layer is staged, so a failure costs one stage, and that a
  rebuild after a network failure is much shorter than the first.
- **Starting over** — note that `docker builder prune --all` discards the staged cache too,
  which is exactly what you do *not* want after a network failure.

### 5. `docker/CLAUDE.md`

One invariant: **the R install is staged, and must stay staged.** Collapsing it back into a
single `RUN` restores an all-or-nothing 30-minute layer. Each stage verifies its own
packages precisely so that a cached layer can be trusted.

### 6. `combobatch_tool_plan_260828.md`

One clause appended to the Phase 6 build item, naming this plan as the third diagnosis.

---

## Considered and not adopted

**Posit Package Manager binaries for CRAN.** `https://packagemanager.posit.co/cran/__linux__/bookworm/latest`
serves precompiled Debian 12 binaries — verified reachable, HTTP 200. Roughly **820 s of
this build's 970 s was CRAN compiling from source**, so this is the single biggest lever on
build time, and a shorter build is a smaller window for the network to change under it.

Not proposed here for three reasons: PPM does **not** mirror Bioconductor at the equivalent
paths (both candidate URLs return 404, verified), so it addresses only half the work;
`latest` is a moving target, which trades reproducibility for speed unless a dated snapshot
URL is pinned; and it is an orthogonal change that should not ride along with a fix for a
network failure. Worth its own plan if build time becomes the binding constraint.

**Retrying failed downloads in-process.** A bounded retry would survive a brief blip. It
would not have survived this: the interception lasted from 861.8 s to the end and is still
in place now. Staging is the better answer because it makes a *permanent* interruption cost
one stage rather than all of them. If the retry is wanted later it composes cleanly with
stages.

---

## Files that do NOT need to change

| File | Why not |
|---|---|
| the apt layer, `libuv1-dev` | Confirmed working: `* DONE (fs)` at 131.3 s, and the prerequisite checkpoint stayed silent. |
| the preflight probe | It is correct for what it checks and cost 1.1 s. Only its description changes — it is a start-of-build check, not a guarantee. |
| the registries, `MethodSpec.r_packages` | No method is misdeclared; 16 packages were unreachable, not wrong. |
| `requirements.txt`, the Python layers | Never reached. |
| `scripts/smoke_test.sh` | It builds and gates; the resumability lives in the Dockerfile's layering. |
| `docker/tls_inspection_ca_plan_260829.md` | Its Tier 2 remains correctly deferred — an always-on VPN would change that, and this incident is evidence to watch, not yet proof. |

---

## Side effects and caveats

- **This rebuild starts from scratch once.** Editing `install_r_packages.R` invalidates the
  `COPY`, and every stage below it. The staging pays from the *next* failure onward.
- **Editing the R script re-runs all stages.** Correct — a changed script may change any
  stage — but it means script edits and network retries have different costs. Retries after
  a network failure keep the cache; edits do not.
- **Six layers instead of one.** R packages are additive, so total size is essentially
  unchanged; the overhead is layer metadata.
- **A cached stage is only as trustworthy as its verification.** That is why 1c is not
  optional: without per-stage verification, staging would be a regression.
- **The VPN must still be off for the whole build.** Nothing here removes that requirement;
  it only makes violating it cost minutes instead of everything.
- **`docker builder prune --all` after a network failure is now actively harmful** — it
  throws away precisely the stages you want to keep.

---

## Verification commands

```bash
# 1. Is the network clean right now? (utun4 holding the default route = VPN up)
netstat -rn -f inet | head -5
echo | openssl s_client -connect cloud.r-project.org:443 -servername cloud.r-project.org 2>/dev/null \
    | openssl x509 -noout -issuer
# want: an issuer that is not BGSSLInspection

# 2. The R script parses and every stage name is accepted
docker run --rm -v "$PWD/docker/install_r_packages.R:/s.R:ro" python:3.11-slim-bookworm bash -c \
    'apt-get update -qq >/dev/null 2>&1; \
     apt-get install -y -qq --no-install-recommends r-base-core >/dev/null 2>&1; \
     Rscript -e "invisible(parse(\"/s.R\")); cat(\"parses\n\")"'

# 3. Laptop suite and lint
python -m pytest tests/unit -q
black --check .
docker buildx build --check -f docker/Dockerfile .

# 4. The build, VPN disconnected and confirmed unable to reconnect
bash scripts/smoke_test.sh

# 5. Prove the staging works: interrupt a build mid-R-layer, then re-run and watch
#    the earlier stages report CACHED rather than re-running
docker buildx build --platform linux/amd64 --load -t combobatch:test -f docker/Dockerfile . \
    2>&1 | grep -E "CACHED .*--stage|--stage"
```

---

## Three deviations found while implementing

**Two existing tests broke, both correctly.** `test_the_prerequisites_are_checked_before_the_long_installs`
anchored on the literal `stop("prerequisite package(s) failed to install…")`, which is now
`stop_with_diagnosis(absent, "prereq")`; the anchor was moved, the intent kept.
`test_every_verified_package_is_also_installed` split the script on
`# --- one consolidated verification`, a marker this restructure removes — and would have
silently passed against the whole file, quietly testing nothing. Rewritten around what can
now actually drift: the `cran` and `bioc` stages pass `STAGE_PACKAGES` straight to their
installer and cannot drift, so it checks that every `github` package has a repository to
clone from and that `bapred` is really in its install call.

**Two troubleshooting rows described a message whose meaning changed.** Both matched
`Failed to install N required package(s)`, and they duplicated each other. That string now
appears only at the `verify` stage, where it means *a cached layer is missing what it
installed* — a different problem with a different fix (`--no-cache`). The per-stage failure
has its own message, `stage <name> could not install: …`, and its own row.

**`BiocParallel` and `matrixStats` are installed but not verified.** They are
`variancePartition`'s backend and dependency, not method backends, so they belong in the
`bioc` install call but not in `STAGE_PACKAGES` — which is the vector the registry
cross-check compares. Listing them would have failed
`test_nothing_is_verified_that_no_registry_needs`.

---

## TODO

- [x] `install_r_packages.R` — add `STAGE_PACKAGES`; derive `required` from it
- [x] `install_r_packages.R` — accept `--stage`, gate each install block on it
- [x] `install_r_packages.R` — `verify_stage()`, called at the end of every stage
- [x] `install_r_packages.R` — `FAILURES` accumulator and `stop_with_diagnosis()`
- [x] `install_r_packages.R` — `WARN:` and the banners to `stderr()`
- [x] `install_r_packages.R` — `options(timeout = 300)`
- [x] `docker/Dockerfile` — six staged `RUN`s in place of the single one
- [x] `tests/unit/test_docker_assets.py` — `verified_packages()` reads `STAGE_PACKAGES`
- [x] `tests/unit/test_docker_assets.py` — the new assertions (six, not four: `WARN:` to
      stderr and the branching diagnosis got one each)
- [x] `tests/unit/test_docker_assets.py` — *(deviation)* repair the two tests this
      restructure broke
- [x] `docker/BUILD_AND_PUSH.md` — step 0b: the probe is point-in-time; the mid-build check
- [x] `docker/BUILD_AND_PUSH.md` — troubleshooting row for the SSL-peer-certificate message
- [x] `docker/BUILD_AND_PUSH.md` — step 3 and "Starting over": staged layers, and do not
      prune the cache after a network failure
- [x] `docker/CLAUDE.md` — the "the R install stays staged" invariant
- [x] `combobatch_tool_plan_260828.md` — name this plan under Phase 6
- [x] Run verification steps 1–3 — script parses; `--stage verify` and an unknown stage both
      produce the right error; 21 required, 21 unique across the five stages; 852 unit tests
      pass; Black clean; `buildx --check` clean
- [ ] **VPN disconnected and confirmed unable to reconnect**, run `bash scripts/smoke_test.sh`
- [ ] Confirm on a re-run that completed stages report `CACHED`
