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
# libpng-dev and zlib1g-dev are load-bearing: without them Hmisc fails to compile and
# BiocManager reports qsmooth as skipped rather than failed. cmake is needed by the
# build chain of the metrics dependencies.
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
        libhdf5-dev libbz2-dev liblzma-dev zlib1g-dev \
        libblas-dev liblapack-dev \
        openssh-server rsync tmux procps \
    && rm -rf /var/lib/apt/lists/*

# --- 2. R, pinned ------------------------------------------------------------------
# rpy2 3.6.x has a C-level ABI incompatibility with R 4.6.0 whose symptom is a segfault
# around 33_amdbnorm. The donor installed `r-base` unpinned from bookworm-cran40, which
# drifts into 4.6 on its own schedule; this installs one exact build instead.
#
# Fallback if the Posit CDN ever stops serving this file:
#   apt-get install -y r-base=4.5.3-1~bookwormcran.0 r-base-dev=4.5.3-1~bookwormcran.0
#   apt-mark hold r-base r-base-dev
ARG R_VERSION=4.5.3
RUN set -eu; \
    curl -fsSL -o /tmp/r.deb \
        "https://cdn.posit.co/r/debian-12/pkgs/r-${R_VERSION}_1_amd64.deb" \
        || { echo "R ${R_VERSION} is not published for $(dpkg --print-architecture); build with --platform linux/amd64" >&2; exit 1; }; \
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
COPY docker/install_r_packages.R /tmp/install_r_packages.R
RUN Rscript /tmp/install_r_packages.R && rm -f /tmp/install_r_packages.R

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
        'exec octave --no-gui "${args[@]}"' \
        > /usr/local/bin/matlab \
    && chmod +x /usr/local/bin/matlab \
    && octave --version | head -1

# --- 5. Python dependencies --------------------------------------------------------
# Copied alone so that a source edit does not re-resolve the whole dependency tree.
# rpy2 compiles against the R installed above, and reComBat is a git install, which is
# why both git and network access are required at build time.
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm -f /tmp/requirements.txt

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
