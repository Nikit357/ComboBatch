"""Logging setup and BLAS/OpenMP thread pinning.

This module must stay importable without numpy, pandas or any other numeric library:
``pin_threads()`` only has an effect if it runs *before* those libraries are imported,
so the CLI entry point imports this first and everything else second.
"""

from __future__ import annotations

import logging
import os
import sys
from datetime import datetime, timezone

# Every threading knob the numeric stack reads at import time. Unpinned, each worker
# fans out to one thread per core; with N workers on one node that oversubscribes badly
# and costs more than the parallelism gains. Pod-level environment variables are not
# enough — they do not reach processes launched over SSH — so this has to happen here.
_THREAD_ENV_VARS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def pin_threads(n_threads: int = 1) -> dict[str, str]:
    """
    Pin every BLAS/OpenMP thread-count variable, unless already set by the caller.

    Parameters
    ----------
    n_threads
        Value to assign. One is correct whenever an outer dispatcher provides the
        parallelism.

    Returns
    -------
    The variables this call actually set, mapping name to value. Variables already
    present in the environment are left alone and omitted from the result.
    """
    applied: dict[str, str] = {}
    for var in _THREAD_ENV_VARS:
        if var not in os.environ:
            os.environ[var] = str(n_threads)
            applied[var] = str(n_threads)
    return applied


def utc_timestamp() -> str:
    """Return the current UTC time as ``HH:MM:SS``, for progress lines."""
    return datetime.now(timezone.utc).strftime("%H:%M:%S")


def setup_logging(level: str = "INFO", *, stream=None) -> logging.Logger:
    """
    Configure the ``combobatch`` logger and return it.

    Logs go to stderr so that stdout stays clean for machine-readable output.
    Calling this more than once replaces the handler rather than adding a second one.

    Parameters
    ----------
    level
        Any level name accepted by :mod:`logging`.
    stream
        Destination stream; defaults to ``sys.stderr``.

    Returns
    -------
    The configured ``combobatch`` logger.
    """
    logger = logging.getLogger("combobatch")
    logger.setLevel(getattr(logging, level.upper()))
    logger.handlers.clear()

    handler = logging.StreamHandler(stream if stream is not None else sys.stderr)
    handler.setFormatter(
        logging.Formatter(
            fmt="[%(asctime)s][%(levelname)s] %(message)s", datefmt="%H:%M:%S"
        )
    )
    logger.addHandler(handler)
    logger.propagate = False
    return logger


def get_logger() -> logging.Logger:
    """Return the ``combobatch`` logger without reconfiguring it."""
    return logging.getLogger("combobatch")
