#!/usr/bin/env bash
#
# Acceptance gate for a built image. Everything here must pass before the image is
# pushed, because each check corresponds to a failure that has actually happened:
#
#   R 4.5.x        rpy2 3.6.x segfaults against R 4.6.0, around 33_amdbnorm
#   octave         20_shambhala silently has no backend without it
#   matlab wrapper Shambhala2 calls system("matlab ...") internally
#   zero skips     a method whose R package failed to install would otherwise SKIP
#
#   bash scripts/smoke_test.sh                 # builds combobatch:test, then checks it
#   bash scripts/smoke_test.sh my-image:tag    # checks an image that already exists

set -uo pipefail

IMAGE="${1:-}"
PLATFORM="--platform linux/amd64"
FAILED=0
SKIPPED=0

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

if [ -z "$IMAGE" ]; then
    IMAGE="combobatch:test"
    echo "==> building ${IMAGE} (linux/amd64)"
    docker buildx build $PLATFORM --load -t "$IMAGE" -f docker/Dockerfile . || exit 1
fi

check() {
    local label="$1"; shift
    printf '\n==> %s\n' "$label"
    if "$@"; then
        printf 'OK: %s\n' "$label"
    else
        printf 'FAIL: %s\n' "$label"
        FAILED=$((FAILED + 1))
    fi
}

in_image() { docker run --rm $PLATFORM "$IMAGE" "$@"; }

check "R reports 4.5.x through rpy2" \
    in_image python -c \
    'import rpy2.robjects as ro; v=str(ro.r("R.version$version.string")[0]); print(v); assert "4.5" in v, v'

check "R dependencies of the metrics and fabatch paths load" \
    in_image Rscript -e 'library(variancePartition); library(bapred); cat("R deps OK\n")'

check "Octave is present" in_image octave --version

check "the matlab wrapper strips MATLAB-only flags" \
    in_image matlab -nosplash -nodesktop --eval "disp('wrapper ok')"

check "the three registries load" in_image combobatch list-methods
check "metric groups load without a panel file" in_image combobatch list-metrics

check "every backend is present - no method or group may skip" \
    in_image pytest tests/integration/test_methods_all_backends.py -q

printf '\n==> combobatch selftest\n'
in_image combobatch selftest
status=$?
if [ "$status" -eq 0 ]; then
    printf 'OK: selftest\n'
elif [ "$status" -eq 2 ]; then
    # main.py exits 2 and names the phase for a subcommand that has not landed yet.
    printf 'SKIPPED: selftest is not implemented yet (Phase 9)\n'
    SKIPPED=$((SKIPPED + 1))
else
    printf 'FAIL: selftest exited %d\n' "$status"
    FAILED=$((FAILED + 1))
fi

printf '\n=====================================\n'
if [ "$FAILED" -eq 0 ]; then
    printf '%s passed all checks (%d skipped)\n' "$IMAGE" "$SKIPPED"
    exit 0
fi
printf '%s FAILED %d check(s)\n' "$IMAGE" "$FAILED"
exit 1
