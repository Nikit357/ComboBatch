"""Free-memory checks for the dispatcher and its workers.

Container-aware on purpose: ``psutil`` reports the *node's* free memory, not the
cgroup's, so inside a memory-limited pod a naive check happily launches workers that
are then OOM-killed. The cgroup files are consulted first and psutil is only a fallback.
"""

from __future__ import annotations

import time
from pathlib import Path

# cgroup v2, then v1. Values are bytes; v2 writes the literal "max" when unlimited.
_CGROUP_V2_LIMIT = Path("/sys/fs/cgroup/memory.max")
_CGROUP_V2_USAGE = Path("/sys/fs/cgroup/memory.current")
_CGROUP_V1_LIMIT = Path("/sys/fs/cgroup/memory/memory.limit_in_bytes")
_CGROUP_V1_USAGE = Path("/sys/fs/cgroup/memory/memory.usage_in_bytes")

# Returned when neither cgroups nor psutil can answer. Deliberately large: an unknown
# memory situation should not block a run, only an unfavourable known one.
_UNKNOWN_FREE_GB = 999.0

MEMORY_POLL_SECONDS = 30.0


def _read_int(path: Path) -> int | None:
    """Read a single integer from a cgroup file, or None if unavailable."""
    try:
        text = path.read_text().strip()
    except (OSError, ValueError):
        return None
    if text == "max":
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _cgroup_free_gb() -> float | None:
    """Free memory within the current cgroup, in GB, or None if not containerised."""
    for limit_path, usage_path in (
        (_CGROUP_V2_LIMIT, _CGROUP_V2_USAGE),
        (_CGROUP_V1_LIMIT, _CGROUP_V1_USAGE),
    ):
        limit = _read_int(limit_path)
        usage = _read_int(usage_path)
        if limit is None or usage is None:
            continue
        # cgroup v1 signals "unlimited" with a sentinel near 2**63, which would otherwise
        # look like an enormous amount of free memory.
        if limit > 2**60:
            continue
        return max(0.0, (limit - usage) / 1e9)
    return None


def free_gb() -> float:
    """
    Return free memory in GB: cgroup limit first, then psutil, then a large sentinel.

    Returns
    -------
    Free memory in gigabytes. ``999.0`` means "could not determine".
    """
    cgroup = _cgroup_free_gb()
    if cgroup is not None:
        return cgroup

    try:
        import psutil
    except ImportError:
        return _UNKNOWN_FREE_GB
    return psutil.virtual_memory().available / 1e9


def check_memory(min_free_gb: float) -> bool:
    """
    Report whether at least ``min_free_gb`` gigabytes are currently free.

    Parameters
    ----------
    min_free_gb
        Threshold in gigabytes.

    Returns
    -------
    True when free memory meets or exceeds the threshold.
    """
    return free_gb() >= min_free_gb


def wait_for_memory(
    min_free_gb: float,
    *,
    poll_seconds: float = MEMORY_POLL_SECONDS,
    timeout_s: float | None = None,
    sleep=time.sleep,
) -> bool:
    """
    Block until ``min_free_gb`` gigabytes are free, or until ``timeout_s`` elapses.

    Used by the dispatcher before launching each worker. Running children are never
    disturbed; only the next launch is delayed.

    Parameters
    ----------
    min_free_gb
        Threshold in gigabytes.
    poll_seconds
        Interval between checks.
    timeout_s
        Give up after this many seconds. ``None`` waits indefinitely.
    sleep
        Sleep callable; injectable so tests need not actually wait.

    Returns
    -------
    True if the threshold was met, False if the timeout expired first.
    """
    waited = 0.0
    while not check_memory(min_free_gb):
        if timeout_s is not None and waited >= timeout_s:
            return False
        sleep(poll_seconds)
        waited += poll_seconds
    return True
