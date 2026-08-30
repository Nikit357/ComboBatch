# Building and publishing the ComboBatch image

A step-by-step runbook for producing `ghcr.io/nikit357/combobatch` from this repository.

Read this once before the first build. After that, steps 3 and 6 are the whole routine:

```bash
bash scripts/smoke_test.sh              # build + verify
bash scripts/build_and_push_image.sh    # build + push
```

Everything below is what those two scripts do, why each step exists, and what to do when
one of them fails.

---

## What you are building

One image, **5.07 GB** measured 2026-08-30, carrying Python 3.11, R 4.5.3, Octave and every dependency of all 34
harmonization methods, 4 imputers and 14 metric groups. The Karpenter pod consumes it
directly, which is what replaces the donors' 30–40 minute startup install with an image
pull.

**Two facts drive this entire document:**

1. **The build must target `linux/amd64`.** This Mac is arm64; the `c6a` nodes are amd64.
   A native build produces a pod that dies with `exec format error`.
2. **The image takes a long time to build.** Roughly 30 R packages compile from source,
   several of them large Bioconductor trees. Budget **45–90 minutes** for a cold build and
   plan to leave it running. Rebuilds after a Python-only edit take about a minute, because
   the Dockerfile is layered apt → R → R packages → Octave → pip → source.

---

## Step 0 — Prerequisites, once

| Requirement | Check | Notes |
|---|---|---|
| Docker Desktop running | `docker version` | 29.7.2 verified on this Mac |
| buildx builder | `docker buildx ls` | `desktop-linux` must be `running` and list `linux/amd64` |
| Free disk | `df -h .` | **≥ 25 GB.** The image is ~5 GB, but the build cache holds every R stage and is several times larger |
| Network | step 0b | Public internet only — a 67 MB R package, all of Bioconductor's dependency closure, and `git+https://github.com/BorgwardtLab/reComBat`. Nothing corporate is needed, and **a TLS-inspecting VPN must be disconnected first** |
| GitHub token | see step 5 | Only needed to push, not to build |

Verified on this machine 2026-08-29: 149 GB free, `desktop-linux` builder running,
`linux/amd64` available.

---

## Step 0b — Behind a TLS-inspecting network

A VPN or proxy that terminates and re-signs HTTPS presents a private root certificate.
Your Mac trusts it, because IT installed it in the System keychain. **The container does
not**, and every HTTPS fetch in the build dies with:

```
curl: (60) SSL certificate problem: self-signed certificate in certificate chain
```

This is why the symptom is confusing: the same URL works perfectly in your terminal and
fails inside Docker. On this machine the cause is the FortiClient VPN, which carries a
**full tunnel** — it owns the default route, so every connection reaches the corporate
FortiGate and is inspected there.

Check it in about a minute:

```bash
docker run --rm python:3.11-slim-bookworm sh -c \
    'apt-get update -qq >/dev/null 2>&1; \
     apt-get install -y -qq --no-install-recommends curl ca-certificates >/dev/null 2>&1; \
     curl -fsSI https://cdn.posit.co/r/debian-12/pkgs/r-4.5.3_1_amd64.deb >/dev/null; echo exit=$?'
```

| Result | Meaning |
|---|---|
| `exit=0` | Nothing is intercepting. Build. |
| `exit=60` | HTTPS is being re-signed. **Disconnect the VPN and re-run this probe** before starting an hour-long build. |

The Dockerfile probes for exactly this immediately after the apt layer, so a build started
on an inspected network fails in two seconds with an explanation rather than 155 seconds in
with a misleading one.

**That probe is point-in-time, and cannot be anything else.** It answers "is TLS verifiable
right now", at second one. The R layer then needs verifiable TLS continuously for the next
15–30 minutes. On 2026-08-30 the probe passed, the build ran cleanly for 14 minutes, and
FortiClient reconnected underneath it: every download after that failed certificate
verification. A passing probe is not a promise about the rest of the build.

So, before you start:

- **Confirm the VPN will not reconnect on its own.** FortiClient may be configured to
  re-establish the tunnel; quitting the client is more reliable than disconnecting it.
- **Spot-check mid-build** if the build is long. `en0`, not `utun4`, should hold the
  default route:

  ```bash
  netstat -rn -f inet | head -5
  ```

Two things to remember:

- **Stay disconnected for the whole build.** `install_r_packages.R` downloads from CRAN,
  Bioconductor and GitHub throughout, not just at the start. If the tunnel does come back,
  the R install is staged across six Docker layers, so only the stage that was running is
  lost — reconnect-safe, not reconnect-proof.
- **Reconnect before steps 7 and 8**, and before any S3 work — `rnd-sandbox` and the
  dissertation bucket are reachable only through the VPN. Pushing to GHCR (step 6) works
  either way: Docker Desktop's own VM already trusts the corporate CA, which is why the
  base image pulls fine even while the container cannot reach CRAN.

If disconnecting is not possible — an always-on policy, for instance — the alternative is
to install the corporate root into the image's trust store. That is deliberately *not*
implemented here, because the resulting image trusts a private MITM anchor and must never
be published. See `tls_inspection_ca_plan_260829.md`, "Tier 2".

---

## Step 1 — Pre-flight, on the laptop (30 seconds)

Never start an hour-long build on a tree whose unit tests are red. In particular,
`tests/unit/test_docker_assets.py` compares `docker/install_r_packages.R` against the three
registries, and a mismatch there is the one build failure that costs a full hour to
discover any other way.

```bash
cd /Users/user890/Desktop/ComboBatch
source .venv/bin/activate
python -m pytest tests/unit -q          # expect: 722 passed
```

> Use `python -m pytest`, not bare `pytest`: several test modules import `tests.conftest`,
> which needs the repository root on `sys.path`.

Then lint the Dockerfile without building it:

```bash
docker buildx build --check -f docker/Dockerfile .
# expect: Check complete, no warnings found.
```

---

## Step 2 — Decide on a tag

The build script always writes `:latest` and one dated or named tag.

```bash
bash scripts/build_and_push_image.sh            # → :latest and :20260829
bash scripts/build_and_push_image.sh v0.1.0     # → :latest and :v0.1.0
```

Use a date tag for routine rebuilds and a version tag for anything a manuscript or a
finished run will need to reproduce. `:latest` moves; a dated tag does not, and the pod
manifest can be pointed at it when a result must stay reproducible.

---

## Step 3 — Build and verify locally

The one command that does both:

```bash
bash scripts/smoke_test.sh
```

It builds `combobatch:test` for `linux/amd64`, loads it into the local daemon, and then
runs the acceptance gate. Each check exists because the corresponding failure has actually
happened:

> The R install runs as **six layers** — `prereq`, `cran`, `bapred`, `bioc`, `github`,
> `verify` — each verifying its own packages before it exits. A failure costs the stage it
> was in, not the whole 15–30 minutes, and a rebuild resumes from the last one that
> succeeded. Watch for `CACHED` against the `--stage` lines to see it working.
>
> The R layer is materially faster since 2026-08-30, when `dependencies = TRUE` came out of
> `install_r_packages.R`. That flag pulled `Suggests` — V8, shiny, rgl, plotly, gsl, sodium
> and ten more — and compiled them for 860 seconds before failing. Nothing in the registries
> imports any of them.

| Check | Guards against |
|---|---|
| `R.version.string` is 4.5.x | rpy2 3.6.x has a C-level ABI incompatibility with R 4.6.0 — the symptom is a segfault around `33_amdbnorm` |
| `library(variancePartition); library(bapred)` | metric group F and `37_fabatch` silently losing their backend |
| `octave --version` | `20_shambhala` having no backend |
| `matlab -nosplash -nodesktop` | Shambhala2 calls `system("matlab …")` internally; the wrapper strips MATLAB-only flags |
| `combobatch list-methods` / `list-metrics` | the registries failing to import |
| `pytest tests/integration/test_methods_all_backends.py` | **zero skips** — a method whose R package failed to install would otherwise just SKIP at run time and quietly vanish from the benchmark |

Expected ending:

```
=====================================
combobatch:test passed all checks (1 skipped)
```

The one skip is `combobatch selftest`, which lands in Phase 9. Anything else non-zero is a
real failure — see step 7.

**To build without the checks** (rarely what you want):

```bash
docker buildx build --platform linux/amd64 -t combobatch:test -f docker/Dockerfile .
```

> **A caveat about verifying on this Mac.** The amd64 image runs here only under emulation.
> That is fine for the checks above, but it is slow, and an occasional crash in a compiled
> numeric library is an emulation artifact rather than a real defect. The authoritative
> verification is `combobatch selftest` inside the pod, on a real amd64 node (Phase 7).

---

## Step 4 — Inspect what you built (optional)

```bash
docker images combobatch:test                       # expect ~5 GB
docker run --rm --platform linux/amd64 combobatch:test \
    du -sm /usr/local/lib/python3.11/site-packages   # expect ~1.3 GB, not ~5.8
docker run --rm --platform linux/amd64 combobatch:test combobatch list-methods | wc -l
docker run --rm --platform linux/amd64 combobatch:test \
    Rscript -e 'cat(length(rownames(installed.packages())), "R packages\n")'
```

The entry point runs any real executable as given and hands everything else to the CLI, so
all four of these work:

```bash
docker run --rm --platform linux/amd64 combobatch:test                       # combobatch --help
docker run --rm --platform linux/amd64 combobatch:test combobatch --version
docker run --rm --platform linux/amd64 combobatch:test run --config c.yaml   # → combobatch run …
docker run --rm --platform linux/amd64 combobatch:test octave --version
```

---

## Step 5 — Authenticate to GHCR

GHCR is the only registry this project uses. There is no ECR, no mirror and no
`imagePullSecrets` anywhere, so there is exactly one place an image can be stale.

**Create a token** (once): GitHub → Settings → Developer settings → **Personal access
tokens (classic)** → Generate new token, with the **`write:packages`** scope, which implies
`read:packages`. Use a classic token; fine-grained tokens have incomplete GHCR support.

```bash
export GITHUB_TOKEN=ghp_xxxxxxxxxxxxxxxxxxxx
echo "$GITHUB_TOKEN" | docker login ghcr.io -u nikit357 --password-stdin
# expect: Login Succeeded
```

Do not paste the token as a command-line argument — it lands in your shell history.
`--password-stdin` avoids that.

> The image path must be lowercase (`ghcr.io/nikit357/combobatch`) even though the GitHub
> account renders as `Nikit357`. The login username is not case-sensitive; the path is.

---

## Step 6 — Build and push

```bash
bash scripts/build_and_push_image.sh
```

This rebuilds for `linux/amd64` and pushes `:latest` plus the dated tag straight from the
builder. Layers cached from step 3 are reused, so if nothing changed in between this is
mostly upload time.

To build without pushing (a dry run of the same path):

```bash
COMBOBATCH_PUSH=0 bash scripts/build_and_push_image.sh
```

---

## Step 7 — Make the package public, once

**This is the step that is easy to forget and produces the most confusing failure.** A
newly created GHCR package is **private by default**, and a pod that cannot pull it fails
with a **403 that reads like a missing image**.

1. Open <https://github.com/users/nikit357/packages/container/combobatch/settings>
2. **Danger Zone → Change visibility → Public**

The package links itself to the repository automatically, via the
`org.opencontainers.image.source` label in the Dockerfile.

Then confirm that an unauthenticated client can pull it — this is the exact thing the pod
does:

```bash
docker logout ghcr.io
docker pull ghcr.io/nikit357/combobatch:latest
docker login ghcr.io -u nikit357 --password-stdin <<< "$GITHUB_TOKEN"   # log back in
```

Only after that pull succeeds is the image actually usable by the cluster.

---

## Step 8 — Use it

```bash
kubectl apply -f k8s/pod-combobatch.yaml -n rnd-sandbox        # Phase 7
kubectl exec -it danya-nikitin-combobatch -n rnd-sandbox -- combobatch selftest
```

Startup should be an image pull of a few minutes, not a 40-minute install. If you see
package installation in the pod logs, the manifest is still running a startup script it
should not need.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| Pod: `exec format error` | Image built for arm64 | Rebuild with `--platform linux/amd64`; always go through `scripts/build_and_push_image.sh` |
| Pod: `403 Forbidden` / `denied` on pull | GHCR package still private | Step 7 |
| Push: `denied: permission_denied` | Token lacks `write:packages`, or you are not logged in | Regenerate a classic token with that scope; `docker login ghcr.io` again |
| Build: `the pinned R build is published for amd64 only, this is arm64…` | Built without `--platform linux/amd64` | The Posit `.deb` is amd64-only; use the script. This message now appears only for a genuine architecture mismatch — it used to be printed for *any* download failure, including TLS ones |
| Build: `R version pin failed` | The Posit CDN stopped serving that exact `.deb` | Fall back to the commented alternative in the Dockerfile: `apt-get install r-base=4.5.3-1~bookwormcran.0` + `apt-mark hold` |
| Build: `status was 'SSL peer certificate or SSH remote key was not OK'`, part-way through the R layer | The VPN reconnected **during** the build; the preflight probe only checks the moment it starts | Quit FortiClient rather than disconnecting it, confirm `netstat -rn -f inet` shows `en0` holding the default route, and rebuild. Completed stages are cached, so the rebuild resumes — **do not** prune the build cache first |
| Build: `ERROR: configuration failed for package ‘fs’`, "libuv was not found" | `libuv1-dev` missing from the apt layer | Restore it. `fs` 2.x needs a system libuv, and `qsmooth → Hmisc → rmarkdown → bslib → sass → fs` makes it a hard requirement of this image |
| Any command: `TypeError: 'NoneType' object is not subscriptable` from `rinterop.py` | An R expression whose value is returned *invisibly*. rpy2 3.6 maps that to `None`, so there is nothing to index | Never read the result of `ro.r("requireNamespace(...)")`, `library(...)` or an assignment. Use `rpy2.robjects.packages.isinstalled`, or wrap the expression in `isTRUE(...)` to force visibility |
| The image is far larger than 9 GB | A Python dependency pulled `torch`, and with it `nvidia/*` and `triton` — ~4.5 GB of CUDA in a CPU-only image | `du -sm /usr/local/lib/python3.11/site-packages/*  \| sort -rn \| head`. `harmonypy` 0.2.0 did exactly this and is pinned `<0.2` |
| Octave prints `error: ignoring const execution_exception& while preparing to exit` | An Octave 7.3 shutdown-path artefact | Ignore it. Exit code is 0, stdout is correct, and it appears from plain `octave --no-gui` too. The Shambhala bridge gates on the return code, never on stderr text |
| Build: `The 'sklearn' PyPI package is deprecated` / `metadata-generation-failed`, at the pip layer | Something is pulling reComBat's declared dependencies. It names `sklearn`, the deprecated stub that now only raises | reComBat has its own Dockerfile layer, installed with `--no-deps --ignore-requires-python`. **Do not** add `sklearn` to `requirements.txt` to satisfy it — the stub installs no module, and the real dependency, `scikit-learn`, is already pinned |
| Build: `reComBat requires a different Python: 3.11 … not in '>=3.8,<3.11'` | The same install, without `--ignore-requires-python` | That bound is stale metadata; reComBat runs correctly on 3.11 (measured). Keep the flag, and keep the import check that follows it |
| `21_harmonizr` returns an empty matrix, or `produced an empty result` | HarmonizR writes a three-byte file instead of raising when its ComBat blocks yield nothing. Modes 1 and 3 (`mean.only = FALSE`) do that on any matrix without missing values | Leave `combat_mode` at its default of 2, or try 4. Modes 1 and 3 are selectable but return nothing on complete data |
| Build: `stage <name> could not install: …` | One stage's packages are missing. The message itself says whether the cause was a failed download or a failed build | Follow what it names. Search the log for `WARN:` and for `ERROR: configuration failed` — the named package is often the symptom, and the first `configuration failed` is the cause, normally a missing `-dev` system library. Only that stage and the ones after it re-run |
| Build: `Failed to install N required package(s): …` at the `verify` stage | Every stage verifies itself, so this means a **cached layer** is missing what it installed | `docker buildx build --no-cache …`, or delete the image and rebuild |
| `test_methods_all_backends.py` skips instead of running | `COMBOBATCH_FULL_ENV=1` not set | It is set in the Dockerfile; if you are running the tests outside the image, that skip is correct |
| Build fails with no space left | Build cache | `docker builder prune -f` (add `--all` to also drop cached base layers) |
| Build: `TLS verification failed inside the container` right after the apt layer | TLS-intercepting VPN or proxy | Disconnect it for the build — step 0b |
| Build: `curl: (60) SSL certificate problem: self-signed certificate in certificate chain`, at the R download or during `install_r_packages.R` | The same, on an image built before the preflight probe existed | Disconnect it for the build — step 0b |

## Starting over

```bash
docker builder prune --all -f       # drop the build cache
docker rmi combobatch:test          # drop the local image
```

A clean rebuild is again 45–90 minutes.

> **Do not do this after a network failure.** The build cache is what holds the R stages
> that already succeeded, and pruning it turns a resumable rebuild back into a full one.
> Reach for it when a layer is *wrong* — a stale cached stage, a changed system library —
> not when a download failed.
