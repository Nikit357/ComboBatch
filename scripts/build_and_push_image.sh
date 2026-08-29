#!/usr/bin/env bash
#
# Build the combobatch image for linux/amd64 and push it to GHCR.
#
#   echo "$GITHUB_TOKEN" | docker login ghcr.io -u nikit357 --password-stdin
#   bash scripts/build_and_push_image.sh              # tags :latest and :YYYYMMDD
#   bash scripts/build_and_push_image.sh v0.1.0       # tags :latest and :v0.1.0
#   COMBOBATCH_PUSH=0 bash scripts/build_and_push_image.sh   # build locally, no push
#
# GHCR is the only registry. There is no ECR, no mirror and no imagePullSecrets anywhere
# in this project, so there is exactly one place an image can be stale.
#
# --platform linux/amd64 is mandatory: the build host is arm64, the c6a nodes are amd64,
# and a native build gives a pod that fails with `exec format error`.

set -euo pipefail

IMAGE="${COMBOBATCH_IMAGE:-ghcr.io/nikit357/combobatch}"
VERSION="${1:-$(date +%Y%m%d)}"
PUSH="${COMBOBATCH_PUSH:-1}"

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if ! docker buildx version >/dev/null 2>&1; then
    echo "docker buildx is required for a cross-platform build; install it and retry" >&2
    exit 1
fi

if [ "$PUSH" = "1" ]; then
    # --push publishes straight from the builder. --load would import the image into the
    # local daemon instead, where an amd64 image still runs on this arm64 Mac, but only
    # under emulation - fine for a smoke test, wrong for a release.
    OUTPUT_FLAG="--push"
    echo "==> building ${IMAGE}:${VERSION} and ${IMAGE}:latest for linux/amd64, then pushing"
else
    OUTPUT_FLAG="--load"
    echo "==> building ${IMAGE}:${VERSION} for the local daemon (no push)"
fi

docker buildx build \
    --platform linux/amd64 \
    --tag "${IMAGE}:latest" \
    --tag "${IMAGE}:${VERSION}" \
    --file docker/Dockerfile \
    ${OUTPUT_FLAG} \
    .

if [ "$PUSH" = "1" ]; then
    cat <<EOF

==> pushed ${IMAGE}:${VERSION} and ${IMAGE}:latest

A newly created GHCR package is PRIVATE by default, and the pod then fails to pull with a
403 that reads like a missing image. Flip it to public once, in the GitHub UI:

    https://github.com/users/nikit357/packages/container/combobatch/settings
    -> Danger Zone -> Change visibility -> Public

Then verify anonymously:

    docker logout ghcr.io && docker pull ${IMAGE}:latest
EOF
fi
