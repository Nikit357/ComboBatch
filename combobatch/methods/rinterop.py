"""rpy2 bridge helpers.

Importable without rpy2 installed: every rpy2 import is deferred into a function body,
so the R-free core install and the whole unit suite work on a laptop with no R.

:func:`r_arglist` is the single place where Python values cross into R source text. The
R-backed methods build their scripts with f-strings, so that boundary needs one
enforced, tested gate rather than an ad-hoc quote at each call site.
"""

from __future__ import annotations

import functools
from typing import Any

_R_TRUE = "TRUE"
_R_FALSE = "FALSE"
_R_NULL = "NULL"


class RNotAvailableError(RuntimeError):
    """Raised when an R-backed method is invoked without a working rpy2/R."""


@functools.lru_cache(maxsize=1)
def r_available() -> bool:
    """Report whether rpy2 is importable and an R runtime is reachable."""
    try:
        import rpy2.robjects  # noqa: F401
    except Exception:
        return False
    return True


@functools.lru_cache(maxsize=None)
def r_package_available(name: str) -> bool:
    """Report whether an R package is installed, without attaching it."""
    if not r_available():
        return False
    import rpy2.robjects as ro

    try:
        result = ro.r(f"requireNamespace('{name}', quietly=TRUE)")
    except Exception:
        return False
    return bool(result[0])


def require_r(method_name: str) -> None:
    """Raise a clear error when an R-backed method is called without R."""
    if not r_available():
        raise NotImplementedError(
            f"{method_name} needs R and rpy2, which are not available here. "
            "Install the extra with `pip install 'combobatch[r]'` (R 4.5.x must "
            "already be installed), or use the combobatch Docker image."
        )


def require_r_packages(method_name: str, packages: tuple[str, ...]) -> None:
    """Raise a clear error naming any R package this method needs but cannot find."""
    require_r(method_name)
    missing = [name for name in packages if not r_package_available(name)]
    if missing:
        raise NotImplementedError(
            f"{method_name} needs R package(s) not installed here: "
            f"{', '.join(missing)}. The combobatch Docker image ships them."
        )


def _converter():
    """Return the pandas<->R converter context, importing rpy2 lazily."""
    import rpy2.robjects as ro
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.conversion import localconverter

    return localconverter(ro.default_converter + pandas2ri.converter)


def py2rpy(obj: Any) -> Any:
    """Convert a pandas object to its R equivalent."""
    import rpy2.robjects as ro

    with _converter():
        return ro.conversion.py2rpy(obj)


def rpy2py(obj: Any) -> Any:
    """Convert an R object back to its pandas/numpy equivalent."""
    import rpy2.robjects as ro

    with _converter():
        return ro.conversion.rpy2py(obj)


def as_numpy(obj: Any):
    """
    Return an R result as a numpy array.

    rpy2 sometimes hands back a numpy array already and sometimes an R matrix,
    depending on the R and rpy2 versions in play; both are accepted.
    """
    import numpy as np

    return obj if isinstance(obj, np.ndarray) else rpy2py(obj)


def r_gc() -> None:
    """Trigger R's garbage collector to release its heap after a large operation."""
    if not r_available():
        return
    import rpy2.robjects as ro

    try:
        ro.r("invisible(gc())")
    except Exception:
        # Reclaiming memory is best-effort; never fail a completed harmonization on it.
        pass


def r_literal(value: Any) -> str:
    """
    Render one Python scalar as R source text.

    Raises
    ------
    TypeError
        For anything that is not a scalar or a flat sequence of scalars. Values are
        interpolated into R source, so unsupported types must be refused rather than
        stringified.
    """
    if isinstance(value, bool):
        return _R_TRUE if value else _R_FALSE
    if value is None:
        return _R_NULL
    if isinstance(value, (int, float)):
        return repr(value)
    if isinstance(value, str):
        escaped = value.replace("\\", "\\\\").replace('"', '\\"')
        return f'"{escaped}"'
    if isinstance(value, (list, tuple)):
        if any(isinstance(item, (list, tuple, dict, set)) for item in value):
            raise TypeError(f"nested sequences cannot be passed to R: {value!r}")
        return f"c({', '.join(r_literal(item) for item in value)})"
    raise TypeError(
        f"cannot pass {type(value).__name__} to R: {value!r}. "
        "Only scalars and flat sequences of scalars are supported."
    )


def r_arglist(params: dict[str, Any], *, leading_comma: bool = False) -> str:
    """
    Render a parameter dict as an R argument list.

    This is the injection boundary: R scripts are built with f-strings, so every value
    that crosses into R source goes through here and anything not a scalar (or a flat
    sequence of scalars) is refused.

    Parameters
    ----------
    params
        Argument names and values. Names are validated as R identifiers.
    leading_comma
        Prepend ``", "`` when the result is non-empty, for splicing after a positional
        argument.

    Returns
    -------
    R source text such as ``k = 5, method = "TMM"``. Empty when ``params`` is empty.
    """
    if not params:
        return ""

    parts = []
    for name in sorted(params):
        if not _is_r_identifier(name):
            raise ValueError(f"not a valid R argument name: {name!r}")
        parts.append(f"{name} = {r_literal(params[name])}")

    rendered = ", ".join(parts)
    return f", {rendered}" if leading_comma else rendered


def _is_r_identifier(name: str) -> bool:
    """Report whether a name is safe to use as an R argument name."""
    if not name or not (name[0].isalpha() or name[0] == "."):
        return False
    return all(char.isalnum() or char in "._" for char in name)
