"""Inside the Docker image, no method may skip.

``tests/unit/test_methods_all.py`` skips a method whose backend is absent, which is the
right behaviour on a laptop but would let a broken image pass silently. This file is the
counterpart: run it in the image and every registered method must be genuinely runnable.

    docker run --rm --platform linux/amd64 combobatch:test \\
        pytest tests/integration/test_methods_all_backends.py -v

Outside a full environment it self-skips, so it is harmless in the laptop suite.
"""

from __future__ import annotations

import os

import pytest

from combobatch.imputation import IMPUTER_REGISTRY
from combobatch.methods import METHOD_REGISTRY

# The image sets this; a laptop does not. Guards the whole module rather than letting it
# fail everywhere R is simply not installed.
IN_FULL_ENVIRONMENT = os.environ.get("COMBOBATCH_FULL_ENV") == "1"

pytestmark = pytest.mark.skipif(
    not IN_FULL_ENVIRONMENT,
    reason="set COMBOBATCH_FULL_ENV=1 to assert every backend is present",
)


def test_every_method_is_available():
    """No registered method may be missing a backend in the image."""
    unavailable = {
        key: spec.missing_backends()
        for key, spec in METHOD_REGISTRY.items()
        if not spec.is_available()
    }
    assert (
        not unavailable
    ), "these methods cannot run in this environment: " + "; ".join(
        f"{key} needs {', '.join(missing)}" for key, missing in unavailable.items()
    )


def test_every_imputer_is_available():
    """No registered imputer may be missing a backend in the image."""
    unavailable = {
        key: spec.missing_backends()
        for key, spec in IMPUTER_REGISTRY.items()
        if not spec.is_available()
    }
    assert (
        not unavailable
    ), "these imputers cannot run in this environment: " + "; ".join(
        f"{key} needs {', '.join(missing)}" for key, missing in unavailable.items()
    )


def test_r_is_the_pinned_version():
    """rpy2 has a C-level ABI incompatibility with R 4.6, so the pin must hold."""
    from combobatch.methods import rinterop

    assert rinterop.r_available(), "R is not reachable from rpy2"

    import rpy2.robjects as ro

    version = str(ro.r("R.version$version.string")[0])
    assert "4.5" in version, f"expected R 4.5.x, found {version!r}"
