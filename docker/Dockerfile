# combobatch - one image carrying Python 3.11, R 4.5.3, Octave and every dependency.
#
#   docker buildx build --platform linux/amd64 -t combobatch:test -f docker/Dockerfile .
#
# linux/amd64 is mandatory: the dev machine is arm64, the c6a nodes are not, and the
# pinned R build below is published for amd64 only. Build context is the repository root.
#
# Layer order is apt -> R -> R packages -> Octave -> pip -> package install, so that
# editing the Python source rebuilds only the last two layers and not the hour of R.

FROM python:3.11-slim-bookworm

# slim-bookworm, not plain slim: the R pin below is a Debian 12 build, and `3.11-slim`
# will silently move to the next Debian release.

LABEL org.opencontainers.image.title="combobatch" \
      org.opencontainers.image.description="Benchmark imputation x harmonization combinations for bulk transcriptomics." \
      org.opencontainers.image.source="https://github.com/Nikit357/ComboBatch" \
      org.opencontainers.image.licenses="MIT"

ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1

# --- 1. system libraries -----------------------------------------------------------
# libpng-dev, zlib1g-dev and libuv1-dev are load-bearing, all three for the same reason:
# qsmooth Imports Hmisc, Hmisc Imports rmarkdown, and that chain runs
# rmarkdown -> bslib -> sass -> fs, where fs 2.x requires a system libuv (fs 1.6.x bundled
# it; USE_BUNDLED_LIBUV=1 is the escape hatch if the Debian package ever goes away).
# Without libuv1-dev the 2026-08-29 build spent 36 minutes and died reporting a single
# missing package, qsmooth, five dependency levels away from the actual cause. cmake is
# needed by the build chain of the metrics dependencies.
#
# openssh-server, rsync, tmux and procps are for the interactive pod (k8s/pod-combobatch.yaml):
# without them in the image, the pod would have to apt-get them at startup, which is the
# boot-time-install pattern this whole image exists to remove. ~10 MB in a 6-9 GB image.
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates curl wget gnupg2 git \
        build-essential gfortran cmake pkg-config \
        libcurl4-openssl-dev libssl-dev libxml2-dev \
        libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
        libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
        libhdf5-dev libbz2-dev liblzma-dev zlib1g-dev libuv1-dev \
        libblas-dev liblapack-dev \
        openssh-server rsync tmux procps \
    && rm -rf /var/lib/apt/lists/*

# --- 1a. preflight: is HTTPS usable from inside the container? ----------------------
# apt above uses http://deb.debian.org, so it succeeds even on a network that re-signs
# TLS - and the failure then surfaces 155 seconds later, at the first https fetch, as a
# curl error about the R download. A full-tunnel VPN (FortiClient here) presents a
# private root that the host trusts via its OS keychain and the container does not,
# which is why this only ever fails inside Docker, and why the fix is to disconnect
# rather than to debug the build.
#
# No -f: this probe tests TLS, not availability. The CDN root answers 403 by design,
# which with -f would be exit 56 and would fail every build.
RUN set -eu; \
    rc=0; curl -sS -o /dev/null --max-time 30 https://cdn.posit.co/ || rc=$?; \
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

# --- 2. R, pinned ------------------------------------------------------------------
# rpy2 3.6.x has a C-level ABI incompatibility with R 4.6.0 whose symptom is a segfault
# around 33_amdbnorm. The donor installed `r-base` unpinned from bookworm-cran40, which
# drifts into 4.6 on its own schedule; this installs one exact build instead.
#
# Fallback if the Posit CDN ever stops serving this file:
#   apt-get install -y r-base=4.5.3-1~bookwormcran.0 r-base-dev=4.5.3-1~bookwormcran.0
#   apt-mark hold r-base r-base-dev
ARG R_VERSION=4.5.3
# The architecture is checked, not inferred from a failed download: attributing every
# curl failure to the platform is what made a TLS error read as "not published for
# amd64" on a build that already was amd64.
RUN set -eu; \
    arch="$(dpkg --print-architecture)"; \
    [ "$arch" = "amd64" ] \
        || { echo "the pinned R build is published for amd64 only, this is ${arch}; rebuild with --platform linux/amd64" >&2; exit 1; }; \
    url="https://cdn.posit.co/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb"; \
    rc=0; curl -fsSL -o /tmp/r.deb "$url" || rc=$?; \
    if [ "$rc" -ne 0 ]; then \
        case "$rc" in \
            60) echo "TLS verification failed for cdn.posit.co - a VPN or proxy is re-signing HTTPS; disconnect it and rebuild" ;; \
            22|56) echo "R ${R_VERSION} is no longer served at ${url}; use the pinned apt fallback in the comment above" ;; \
            6|7|28) echo "cdn.posit.co unreachable (curl exit ${rc})" ;; \
            *) echo "downloading ${url} failed with curl exit ${rc}" ;; \
        esac >&2; \
        exit 1; \
    fi; \
    apt-get update; \
    apt-get install -y --no-install-recommends /tmp/r.deb; \
    rm -f /tmp/r.deb; \
    rm -rf /var/lib/apt/lists/*; \
    ln -sf "/opt/R/${R_VERSION}/bin/R" /usr/local/bin/R; \
    ln -sf "/opt/R/${R_VERSION}/bin/Rscript" /usr/local/bin/Rscript; \
    R --version | head -1; \
    R --version | head -1 | grep -q " ${R_VERSION} " \
        || { echo "R version pin failed - expected ${R_VERSION}" >&2; exit 1; }

# rpy2 loads libR.so at import time and finds R_HOME by asking `R RHOME`, so only the
# shared-library path needs setting here. No trailing colon: an empty entry in this
# variable makes the loader search the working directory.
ENV LD_LIBRARY_PATH="/opt/R/${R_VERSION}/lib/R/lib"

# --- 3. R packages -----------------------------------------------------------------
# Six layers, not one. The R install takes 15-30 minutes and needs verifiable TLS for all
# of it; on 2026-08-30 a VPN reconnected 14 minutes in and a single RUN meant every
# finished compile was discarded with the failure. Each stage verifies its own packages
# before exiting, so a cached layer is a layer that really installed what it claims.
COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R --stage prereq
RUN Rscript /tmp/install_r_packages.R --stage cran
RUN Rscript /tmp/install_r_packages.R --stage bapred
RUN Rscript /tmp/install_r_packages.R --stage bioc
RUN Rscript /tmp/install_r_packages.R --stage github
RUN Rscript /tmp/install_r_packages.R --stage verify \
    && rm -f /tmp/install_r_packages.R

# --- 4. Octave, for 20_shambhala ---------------------------------------------------
# Shambhala2 calls system("matlab ...") internally, so `matlab` must exist on the PATH.
# The wrapper drops the MATLAB-only flags Octave rejects.
RUN apt-get update && apt-get install -y --no-install-recommends \
        octave octave-statistics \
    && rm -rf /var/lib/apt/lists/* \
    && printf '%s\n' \
        '#!/bin/bash' \
        'args=()' \
        'for a in "$@"; do' \
        '  case "$a" in' \
        '    -nodesktop|-nosplash|-nodisplay) ;;' \
        '    *) args+=("$a") ;;' \
        '  esac' \
        'done' \
        'exec octave-cli "${args[@]}"' \
        > /usr/local/bin/matlab \
    && chmod +x /usr/local/bin/matlab \
    && octave --version | head -1

# --- 5. Python dependencies --------------------------------------------------------
# Copied alone so that a source edit does not re-resolve the whole dependency tree.
# rpy2 compiles against the R installed above, which is why network access is required
# at build time.
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm -f /tmp/requirements.txt

# reComBat, alone, past two stale declarations in its own metadata: `sklearn = "^0.0"`
# (the deprecated stub, which now only raises when built) and `python = ">=3.8,<3.11"`
# (it runs correctly on 3.11 - measured). --no-deps is a global pip flag with no
# per-requirement form, which is why this cannot live in requirements.txt. Its real
# imports - numpy, pandas, scikit-learn, tqdm - are all pinned in the layer above; the
# two dependencies dropped are the sklearn stub and `fire`, which only reComBat's unused
# CLI imports.
#
# Pinned to a commit, and deliberately NOT to the PyPI wheel: the published 0.1.4
# defaults to model='linear' where this commit defaults to model='elastic_net', so the
# wheel would silently change every 30_recombat result and break comparison with the
# published dissertation metrics.
ARG RECOMBAT_COMMIT=b02bd025c9671e75f37dc513d7f466b05448364d
RUN pip install --no-cache-dir --no-deps --ignore-requires-python \
        "reComBat @ git+https://github.com/BorgwardtLab/reComBat@${RECOMBAT_COMMIT}" \
    && python -c "from reComBat import reComBat; print('reComBat import OK')"

# --- 6. combobatch itself ----------------------------------------------------------
# --no-deps: requirements.txt above is the single source of installed versions, and
# resolving pyproject.toml again here could quietly upgrade a pinned one.
WORKDIR /app
COPY . /app
RUN pip install --no-cache-dir --no-deps . \
    && combobatch --version

# --- 7. environment ----------------------------------------------------------------
# Pod-level env vars do not reach processes launched over SSH, so the values are written
# into the image as well as into /etc/environment, which PAM does read for SSH logins.
# combobatch.logging_utils.pin_threads() only ever fills in variables that are unset, so
# these remain the effective defaults and --threads still overrides them.
ENV OMP_NUM_THREADS=1 \
    OPENBLAS_NUM_THREADS=1 \
    MKL_NUM_THREADS=1 \
    NUMEXPR_NUM_THREADS=1 \
    VECLIB_MAXIMUM_THREADS=1 \
    COMBOBATCH_FULL_ENV=1 \
    PYTHONUNBUFFERED=1

RUN printf '%s\n' \
        'OMP_NUM_THREADS=1' \
        'OPENBLAS_NUM_THREADS=1' \
        'MKL_NUM_THREADS=1' \
        'NUMEXPR_NUM_THREADS=1' \
        'VECLIB_MAXIMUM_THREADS=1' \
        'COMBOBATCH_FULL_ENV=1' \
        >> /etc/environment

COPY docker/entrypoint.sh /usr/local/bin/entrypoint.sh
RUN chmod +x /usr/local/bin/entrypoint.sh

ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]
CMD ["combobatch", "--help"]
