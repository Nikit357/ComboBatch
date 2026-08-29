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

One image, ~6–9 GB, carrying Python 3.11, R 4.5.3, Octave and every dependency of all 34
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
| Free disk | `df -h .` | **≥ 25 GB.** The image is 6–9 GB and the build cache is comparable |
| Network | — | The build downloads a 67 MB R package, all of Bioconductor's dependency closure, and `git+https://github.com/BorgwardtLab/reComBat` |
| GitHub token | see step 5 | Only needed to push, not to build |

Verified on this machine 2026-08-29: 149 GB free, `desktop-linux` builder running,
`linux/amd64` available.

A VPN that intercepts TLS will break the build at the CRAN/Bioconductor or GitHub step.
If you use one, turn it off for the build.

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
docker images combobatch:test                       # expect 6-9 GB
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
| Build: `R 4.5.3 is not published for arm64…` | Built without `--platform linux/amd64` | The Posit `.deb` is amd64-only; use the script |
| Build: `R version pin failed` | The Posit CDN stopped serving that exact `.deb` | Fall back to the commented alternative in the Dockerfile: `apt-get install r-base=4.5.3-1~bookwormcran.0` + `apt-mark hold` |
| Build: `Failed to install N required package(s): …` | One or more R installs failed | Scroll up to the `WARN:` lines — every install is attempted, so the log names *all* failures, not just the first. A network or CRAN mirror hiccup is the usual cause; re-run the build |
| Image builds, but `14_qsmooth` skips at run time | `Hmisc` failed to compile, and BiocManager reports qsmooth as *skipped* rather than failed | Check `libpng-dev` and `zlib1g-dev` are still in the apt layer; that is precisely what they are there for |
| `test_methods_all_backends.py` skips instead of running | `COMBOBATCH_FULL_ENV=1` not set | It is set in the Dockerfile; if you are running the tests outside the image, that skip is correct |
| Build fails with no space left | Build cache | `docker builder prune -f` (add `--all` to also drop cached base layers) |
| Build hangs at CRAN/Bioconductor or GitHub | TLS-intercepting VPN | Disable it for the build |

## Starting over

```bash
docker builder prune --all -f       # drop the build cache
docker rmi combobatch:test          # drop the local image
```

A clean rebuild is again 45–90 minutes.
