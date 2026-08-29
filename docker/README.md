# Docker

One image, `combobatch`, containing Python 3.11, R 4.5.3, Octave, and every dependency.
Nothing is installed at pod startup: the image *is* the environment, which is what turns
a 30–40 minute node startup into an image pull.

**[`BUILD_AND_PUSH.md`](BUILD_AND_PUSH.md) is the step-by-step runbook** — prerequisites,
the GHCR token, the one-time visibility flip, and a troubleshooting table. This page is the
summary.

```bash
# amd64 is mandatory — the dev Mac is arm64, the target nodes are not
docker buildx build --platform linux/amd64 -t combobatch:test -f docker/Dockerfile .

# the whole acceptance gate in one command: R version, Octave, wrapper, zero skips
bash scripts/smoke_test.sh combobatch:test

bash scripts/build_and_push_image.sh    # → ghcr.io/nikit357/combobatch
```

| File | Role |
|---|---|
| `Dockerfile` | `python:3.11-slim-bookworm` + R 4.5.3 (pinned) + Octave + all packages. |
| `install_r_packages.R` | R and Bioconductor installs, run during the build. |
| `entrypoint.sh` | Fills in the `combobatch` prefix; execs anything else as given. |

## What the entry point does

`ENTRYPOINT ["combobatch"]` would have made `docker run IMAGE pytest …`, `Rscript -e …`,
`octave --version` and the pod's `sleep infinity` override impossible, so the entry point
execs a real executable as given and hands everything else to the CLI:

```bash
docker run IMAGE                        # → combobatch --help
docker run IMAGE combobatch selftest    # → combobatch selftest
docker run IMAGE run --config c.yaml    # → combobatch run --config c.yaml
docker run IMAGE pytest tests/unit      # → pytest tests/unit
```

Thread pinning is *not* done there. It is in the Python entry point, because a process
launched over SSH into a running pod never passes through the container entry point at
all — the image's `ENV` and `/etc/environment` cover that case instead.

## Why the versions are what they are

**R is pinned to 4.5.3 on purpose.** rpy2 3.6.x has a C-level ABI incompatibility with
R 4.6.0 whose symptom is a segfault around `33_amdbnorm`, and the donor's unpinned
`r-base` from `bookworm-cran40` drifts there on its own schedule. The pin uses the Posit
standalone build (`cdn.posit.co/r/debian-12/pkgs/r-4.5.3_1_amd64.deb`, re-verified
2026-08-29: HTTP 200, 67 MB). The fallback, kept as a comment in the Dockerfile, is
`apt-get install r-base=4.5.3-1~bookwormcran.0` plus `apt-mark hold`.

**The base is `3.11-slim-bookworm`, not `3.11-slim`.** `pandas<2.0` ships a py3.11 wheel
only, and the R build above is a Debian 12 package — a base image that moves to the next
Debian release would break the pin silently.

**`install_r_packages.R` is checked against the registries**, not maintained by hand:
`tests/unit/test_docker_assets.py` fails in a second if a method declares an R package the
script does not install. Three packages the donor installed *and verified* are gone
because each aborted the build unconditionally — `FAbatch` (exists nowhere; the real
package is `bapred`), R `reComBat` (a Python package) and `exploBATCH` (dependency repo
deleted). `DASC` is gone too, because `35_dasc` is not a registered method.

Expected image size 6–9 GB.

## Publishing

GHCR is the only registry — no ECR, no mirror, no `imagePullSecrets` anywhere.

```bash
echo "$GITHUB_TOKEN" | docker login ghcr.io -u nikit357 --password-stdin   # write:packages
bash scripts/build_and_push_image.sh
```

A new GHCR package is **private by default** and the pod then fails to pull with a 403 that
reads like a missing image. Flip it to public once in the GitHub UI, then confirm
anonymously with `docker logout ghcr.io && docker pull ghcr.io/nikit357/combobatch:latest`.
