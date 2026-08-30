# Plan — `docker build` fails behind the FortiClient VPN

**Status:** awaiting approval. No code has been changed.
**Date:** 2026-08-29, revised 2026-08-30 after measuring the VPN's role.
**Trigger:** `bash scripts/smoke_test.sh` failed at Dockerfile step 3/14.

---

## Overview

The build failed because the FortiClient VPN carries a **full tunnel** — it owns the default
route — so every HTTPS connection is terminated and re-signed by the BostonGene FortiGate
with a private root, `CN=BGSSLInspection`. The Mac accepts that certificate because macOS
ships the CA in its System keychain. A Debian container starts from the public bundle, has
never heard of it, and correctly refuses.

**Disconnecting the VPN removes the problem entirely** — measured, not assumed: the same
probe that fails with `curl: (60)` on the tunnel returns `exit=0` off it. The image build
needs nothing corporate (PyPI, CRAN, Bioconductor, GitHub and the Posit CDN are all public),
so the remedy is operational, not a code change:

> **Disconnect FortiClient → build → reconnect for the K8s and S3 steps.**

That is also the *better* path, not merely the easier one. An off-VPN build produces an
image that trusts nothing extra, so it is publishable as-is. Any fix that injected the
corporate CA would bake a private MITM anchor into a public MIT image and would have to be
guarded against ever being pushed.

Two code fixes remain worth making, because both are why this cost an hour to diagnose
rather than a minute, and both help every future user regardless of VPN:

1. **The error message was actively wrong.** `curl` failed with exit 60 (certificate) and
   the Dockerfile's `|| { … }` catch-all reported *"R 4.5.3 is not published for amd64;
   build with --platform linux/amd64"* — false in every particular, on a build that already
   was `linux/amd64`. Same defect class as the 13 donor bugs §6 exists to fix: a failure
   reported as something it is not.
2. **The failure surfaced 155 seconds in**, after the apt layer, at the first HTTPS fetch.
   A two-second preflight probe naming the actual cause belongs right after `curl` exists.

### A correction to the previous revision of this plan

The first draft asserted this was "network-level inspection with no VPN to disable," and on
that basis proposed a `docker/certs/` CA-injection mechanism, an export script, a push guard
and eight new unit tests — about eighteen items. That premise was an assumption I had not
checked, and it was wrong.

It also led me to flag a **correct** line in `BUILD_AND_PUSH.md` as an error. The existing
troubleshooting row `TLS-intercepting VPN | Disable it for the build` had the cause and the
remedy right. Only its symptom text is wrong: the build does not *hang*, it fails with
`curl: (60)`. That row is now amended rather than replaced.

The CA-injection machinery is kept below as **Tier 2, and is not recommended for now** — see
*Deferred*.

---

## Evidence

Measured 2026-08-29 / 08-30 on this Mac.

### The tunnel owns the default route

```
$ netstat -rn -f inet
default            100.65.144.3       UGScg               utun4      ← FortiClient
default            10.18.57.217       UGScIg                en0

$ route -n get 18.244.87.31          # cdn.posit.co
  interface: utun4
```

`ztnafw` is running and `utun40` is up, but it holds exactly one host route
(`198.19.254.1`), so ZTNA is not capturing general traffic. `en0` is wired into
`10.18.57.0/24` with the local gateway as its DNS resolver — not a BostonGene one.

### The interception, and its absence off-tunnel

Issuer of the leaf certificate actually presented, on the tunnel:

| Host | Needed by | Issuer presented |
|---|---|---|
| `cdn.posit.co` | Dockerfile step 2 (the R `.deb`) | **BGSSLInspection** |
| `cloud.r-project.org` | `install_r_packages.R` (CRAN) | **BGSSLInspection** |
| `bioconductor.org` | `install_r_packages.R` (BiocManager) | **BGSSLInspection** |
| `packagemanager.posit.co` | the documented CRAN fallback | **BGSSLInspection** |
| `ghcr.io`, `registry-1.docker.io`, `auth.docker.io` | push, and the base-image pull | **BGSSLInspection** |
| `pypi.org`, `files.pythonhosted.org` | step 5, `pip install` | GlobalSign (not intercepted) |
| `github.com`, `codeload.github.com` | `install_github`, `git+https://…reComBat` | Sectigo (not intercepted) |
| `s3.amazonaws.com` and regional endpoints | **runtime** S3 I/O | Amazon (not intercepted) |

The chain served for an intercepted host is two certificates, both BostonGene — leaf plus a
self-signed root the proxy supplies, valid `2026-05-16` → `2029-05-16`.

The decisive measurement, same command, same container image, same machine:

| | VPN connected | VPN disconnected |
|---|---|---|
| `curl -fsSI https://cdn.posit.co/…/r-4.5.3_1_amd64.deb` | **exit 60** | **exit 0** |

The partial allowlist above (PyPI, GitHub, AWS pass through) is a snapshot of somebody
else's FortiGate policy, not a contract. It is a reason not to build a fix that special-cases
three hostnames, and another reason to prefer being off the tunnel entirely.

### Why `apt` succeeded for 155 seconds first

Debian's default sources are **`http://deb.debian.org`**. There is no TLS to inspect, so
step 1 ran to completion and step 2 was the first thing to touch HTTPS.

### Three things that follow

- **`docker pull` and `docker push` work on the tunnel.** `registry-1.docker.io` is
  intercepted, yet the base image pulled fine: Docker Desktop's Linux VM already trusts the
  host keychain. The daemon's trust store and the *container's* are different things, and
  only the latter was broken. So pushing to GHCR does not require going off-VPN — only
  building does.
- **Runtime S3 was never affected.** The AWS endpoints are not intercepted, so `boto3`'s
  `certifi` bundle needs nothing.
- **The K8s and S3 steps still need the VPN.** `rnd-sandbox` and the dissertation bucket are
  reachable only through it. Disconnect for the build, reconnect afterwards.

---

## Tier 1 — the changes recommended now

Both are small, both are worth having whether or not anyone is on a VPN, and neither adds a
publishing hazard.

### 1. `docker/Dockerfile` — the misleading message, lines 52–64

**Before** (every failure blamed on architecture):

```dockerfile
RUN set -eu; \
    curl -fsSL -o /tmp/r.deb \
        "https://cdn.posit.co/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb" \
        || { echo "R ${R_VERSION} is not published for $(dpkg --print-architecture); build with --platform linux/amd64" >&2; exit 1; }; \
```

**After** (architecture is *checked*, not guessed at; curl's exit code is classified):

```dockerfile
RUN set -eu; \
    arch="$(dpkg --print-architecture)"; \
    [ "$arch" = "amd64" ] \
        || { echo "the pinned R build is published for amd64 only, this is ${arch}; rebuild with --platform linux/amd64" >&2; exit 1; }; \
    url="https://cdn.posit.co/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb"; \
    rc=0; curl -fsSL -o /tmp/r.deb "$url" || rc=$?; \
    if [ "$rc" -ne 0 ]; then \
        case "$rc" in \
            60) echo "TLS verification failed for cdn.posit.co - a VPN or proxy is re-signing HTTPS; disconnect it and rebuild" ;; \
            22) echo "R ${R_VERSION} is no longer served at ${url}; use the r-base fallback in the comment above" ;; \
            6|7|28) echo "cdn.posit.co unreachable (curl exit ${rc})" ;; \
            *) echo "downloading ${url} failed with curl exit ${rc}" ;; \
        esac >&2; \
        exit 1; \
    fi; \
    apt-get update; \
    …unchanged from here…
```

`rc=0; cmd || rc=$?` rather than `if ! cmd`: under `set -e` the latter loses the real exit
code — inside the `then` branch `$?` is the *negated* status, i.e. always 0.

### 2. `docker/Dockerfile` — preflight probe, its own `RUN` after the apt layer (after line 41)

A separate layer on purpose: folding it into the apt `RUN` would invalidate that cached
155-second layer every time this text is touched.

```dockerfile
# Fail here, in two seconds and with the right diagnosis, rather than 155 seconds in.
# A full-tunnel VPN (FortiClient here) re-signs HTTPS with a private root that the host
# trusts via its OS keychain and the container does not - which is why this only ever
# fails inside Docker, and why the fix is to disconnect rather than to debug the build.
RUN set -eu; \
    rc=0; curl -fsS -o /dev/null --max-time 30 https://cdn.posit.co/ || rc=$?; \
    if [ "$rc" -eq 60 ]; then \
        echo "----------------------------------------------------------------" >&2; \
        echo "TLS verification failed inside the container." >&2; \
        echo "A VPN or corporate proxy is re-signing HTTPS with a private root" >&2; \
        echo "this image does not trust. Disconnect it and rebuild." >&2; \
        echo "See docker/BUILD_AND_PUSH.md, 'Behind a TLS-inspecting network'." >&2; \
        echo "----------------------------------------------------------------" >&2; \
        exit 1; \
    elif [ "$rc" -ne 0 ]; then \
        echo "cdn.posit.co unreachable (curl exit $rc); check network access" >&2; \
        exit 1; \
    fi
```

### 3. `tests/unit/test_docker_assets.py` — two assertions in `TestDockerfileInvariants`

| Test | Guards against |
|---|---|
| the R step tests `dpkg --print-architecture` against `amd64` **and** `60)` appears among the case arms | the misleading-message bug returning |
| the preflight probe precedes the R download in the file | the diagnosis drifting back to 155 seconds in |

### 4. `docker/BUILD_AND_PUSH.md`

- **Step 0, *Network* row** — add: the build needs public internet only, and a full-tunnel
  VPN that inspects TLS must be disconnected first.
- **New Step 0b — "Behind a TLS-inspecting network"**: the one-line container probe, what
  `exit=60` versus `exit=0` means, and the reminder to reconnect before the K8s and S3 steps.
- **Troubleshooting** — amend the existing `TLS-intercepting VPN` row: its cause and remedy
  were right, its symptom is not `hangs` but `curl: (60) SSL certificate problem:
  self-signed certificate in certificate chain`, at Dockerfile step 2 or during
  `install_r_packages.R`.
- **Troubleshooting** — update the `R 4.5.3 is not published for arm64…` row to quote the
  new wording, and note that this message no longer appears for TLS failures.

### 5. `docker/README.md` and `docker/CLAUDE.md`

One short paragraph in the README under *Publishing*: build off-VPN, and why an image built
behind an inspecting proxy should not be published. One invariant in `CLAUDE.md`: *diagnose
the failure you had, not the one you expected* — with this incident as the example.

### 6. `combobatch_tool_plan_260828.md`

One line under Phase 6's unchecked build item pointing here, so the phase keeps its history.

---

## Deferred — Tier 2, the CA-injection mechanism

Kept for the record; **not recommended now**.

The shape would be a committed-but-empty `docker/certs/` directory, gitignored `*.crt`, a
`COPY` + `update-ca-certificates` layer ahead of apt, a `scripts/export_corp_ca.sh` that
pulls the root out of the System keychain, a push guard in `build_and_push_image.sh`, a
marker file at `/etc/ssl/certs/.combobatch-extra-ca`, and eight unit tests.

Reasons to leave it unbuilt:

- It solves a problem that disconnecting the VPN already solves, for free.
- Every image it produces trusts a private MITM root, so it must never be pushed — a hazard
  that has to be guarded, documented and remembered forever.
- The Tier 1 preflight message already tells any user, at any company, exactly what is
  happening and what to do. That covers the public-tool case without the machinery.

It becomes worth building only if the VPN turns out to be mandatory-always-on for this
machine, or if CI ever has to run behind an inspecting proxy.

---

## Files that do NOT need to change

| File | Why not |
|---|---|
| `docker/install_r_packages.R` | CRAN and Bioconductor go through libcurl and the system store; off-VPN both are reachable. Verified: the CRAN `PACKAGES` index downloads intact (7,139,786 bytes). |
| `requirements.txt` | PyPI and GitHub are not intercepted even on the tunnel. |
| `.github/workflows/ci.yml` | GitHub runners have no interception. |
| `k8s/*.yaml` | The cluster pulls from GHCR over AWS's network and runs a clean image. |
| `combobatch/io_utils.py` and anything touching S3 | The AWS endpoints are not intercepted; `certifi` is untouched. |
| `scripts/build_and_push_image.sh` | No push guard is needed if no CA is ever injected. The push itself works on the tunnel. |
| `docker/entrypoint.sh` | Runtime only; no network. |
| `.gitignore`, `.dockerignore` | Nothing new to ignore without `docker/certs/`. |

---

## Side effects and caveats

- **One 155-second rebuild.** The new preflight `RUN` sits after the apt layer, so that
  layer stays cached; the layers below it rebuild once.
- **The VPN must be off for the whole build**, 45–90 minutes, because
  `install_r_packages.R` downloads from CRAN, Bioconductor and GitHub throughout — not just
  at the start. Reconnect before `kubectl` or any S3 work.
- **If IT enforces always-on VPN**, disconnecting may drop connectivity rather than route
  around it. That is the trigger for reconsidering Tier 2.
- **This does not make the tunnel trustworthy or untrustworthy.** It records that TLS to the
  R hosts is terminated by a third party, and that the build should not depend on it.
- **The R version pin is unaffected.** The `.deb` is served and intact: 67,072,216 bytes,
  matching the `content-length` the host sees.

---

## Verification commands

```bash
# 1. The diagnosis, on and off the tunnel (60 -> the VPN is inspecting; 0 -> clear)
docker run --rm python:3.11-slim-bookworm sh -c \
    'apt-get update -qq >/dev/null 2>&1; \
     apt-get install -y -qq --no-install-recommends curl ca-certificates >/dev/null 2>&1; \
     curl -fsSI https://cdn.posit.co/r/debian-12/pkgs/r-4.5.3_1_amd64.deb >/dev/null; echo exit=$?'

# 2. Confirm which interface holds the default route
netstat -rn -f inet | head -5

# 3. The laptop suite, including the two new guards
python -m pytest tests/unit/test_docker_assets.py -v

# 4. The real thing, VPN disconnected (45-90 min)
bash scripts/smoke_test.sh

# 5. Reconnect, then push
bash scripts/build_and_push_image.sh
```

---

## Two deviations found while implementing

**The probe must not use `-f`.** The plan's snippet had `curl -fsS … https://cdn.posit.co/`,
and the CDN root answers **403 by design** — with `-f` that is exit 56, so the `elif` would
have reported "cdn.posit.co unreachable" and failed *every* build, on every network. Dropping
`-f` makes the probe test what it is for: TLS, not availability. A 404 or 403 still exits 0;
only a certificate or transport failure trips it.

**The `22)` arm became `22|56)`, and could not mention `r-base`.** 56 is the other way a
dead URL surfaces. The message originally read "use the `r-base` fallback in the comment
above", which `test_r_is_pinned_to_4_5_3` correctly rejected: that test forbids the literal
`r-base` anywhere outside a comment, because an unpinned `r-base` install is the drift that
segfaults rpy2. Reworded to "the pinned apt fallback".

---

## TODO

- [x] `docker/Dockerfile` — rewrite the R download: real arch check + curl exit-code cases
- [x] `docker/Dockerfile` — add the preflight TLS probe as its own layer after apt
- [x] `tests/unit/test_docker_assets.py` — two assertions in `TestDockerfileInvariants`
- [x] `docker/BUILD_AND_PUSH.md` — Step 0 *Network* row
- [x] `docker/BUILD_AND_PUSH.md` — new Step 0b, "Behind a TLS-inspecting network"
- [x] `docker/BUILD_AND_PUSH.md` — amend the `TLS-intercepting VPN` row (symptom only; cause
      and remedy were already correct)
- [x] `docker/BUILD_AND_PUSH.md` — update the `not published for arm64` row to the new wording
- [x] `docker/README.md` — build off-VPN, and why such an image should not be published
- [x] `docker/CLAUDE.md` — the "diagnose the failure you had" invariant
- [x] `combobatch_tool_plan_260828.md` — reference this plan from Phase 6
- [x] Run verification steps 1–3 (laptop, minutes)
- [ ] **Disconnect FortiClient**, run `bash scripts/smoke_test.sh` to completion (45–90 min)
      — the item the original failure blocked *(requires an hour of build time)*
- [ ] Reconnect FortiClient before any `kubectl` or S3 step
- [ ] *(only if the VPN proves mandatory)* revisit Tier 2
