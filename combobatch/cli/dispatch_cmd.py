"""``combobatch dispatch`` — the imputation x method cross-product, in parallel.

A ``ThreadPoolExecutor`` supervising ``subprocess.Popen`` workers. This shape is
load-bearing, not incidental: every rpy2 job needs its own R session, and neither
threads-around-R nor a forked pool can give it one. The threads here only pump pipes and
wait; the work happens in the child processes.

Two donor defects are fixed at the source. The job lists are derived from the registries
instead of being retyped by hand in each dispatcher, and an unparseable key is a loud
error rather than a silently dropped job.
"""

from __future__ import annotations

import argparse
import signal
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from typing import Any

from combobatch import storage
from combobatch.cli import run_cmd
from combobatch.config import (
    ConfigError,
    ImputationSelection,
    MethodSelection,
    RunConfig,
)
from combobatch.logging_utils import get_logger, pin_threads, setup_logging

# Worker exit codes, mirrored from run_cmd so the dispatcher can act on them.
EXIT_INSUFFICIENT_MEMORY = run_cmd.EXIT_INSUFFICIENT_MEMORY
EXIT_MISSING_INPUT = run_cmd.EXIT_MISSING_INPUT
EXIT_TIMEOUT = -9

ALL = "all"


@dataclass
class JobResult:
    """What one worker subprocess produced."""

    key: str
    returncode: int
    elapsed_s: float
    rows: list[dict[str, Any]]

    @property
    def succeeded(self) -> bool:
        """True when the worker exited cleanly and reported no failing output."""
        if self.returncode != 0 or not self.rows:
            return False
        statuses = {row.get("status", "failed") for row in self.rows}
        return not (statuses & {"failed", "upload_failed", "post_removal_failed"})


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the ``dispatch`` flags on top of the shared ones."""
    run_cmd.add_common_arguments(parser)

    axes = parser.add_argument_group("cross-product axes")
    axes.add_argument(
        "--imputations",
        metavar="LIST",
        help="Comma-separated imputer keys, or 'all' for every registered imputer.",
    )
    axes.add_argument(
        "--methods",
        metavar="LIST",
        help="Comma-separated method keys, or 'all' for every registered method.",
    )

    control = parser.add_argument_group("dispatch control")
    control.add_argument(
        "--n-workers",
        type=int,
        default=4,
        metavar="N",
        help="Concurrent worker subprocesses. Default 4.",
    )
    control.add_argument(
        "--timeout-s",
        type=int,
        default=3600,
        metavar="N",
        help="Per-job wall-clock limit. A timed-out job is logged as failed.",
    )
    control.add_argument(
        "--memory-limit-gb",
        type=float,
        default=4.0,
        metavar="GB",
        help="Wait before launching until this much is free; workers get half of it.",
    )
    control.add_argument(
        "--memory-wait-s",
        type=float,
        default=600.0,
        metavar="N",
        help="Give up waiting for --memory-limit-gb after this long and launch anyway; "
        "the worker's own guard then reports the shortfall honestly. Default 600.",
    )
    control.add_argument(
        "--dry-run",
        action="store_true",
        help="Resolve and print the full job matrix, then exit without running anything.",
    )
    control.add_argument(
        "--retry-failed",
        action="store_true",
        help="Run only the jobs listed in this run's failed-jobs log.",
    )


def expand_axis(
    requested: str | None, registry: dict[str, Any], axis: str
) -> list[str]:
    """
    Expand a comma-separated axis, resolving ``all`` against the registry.

    Parameters
    ----------
    requested
        The raw flag value, or ``None`` when the config file supplies the axis.
    registry
        ``METHOD_REGISTRY`` or ``IMPUTER_REGISTRY``.
    axis
        Name used in error messages.

    Returns
    -------
    The requested keys, in registry order when ``all`` was given.

    Raises
    ------
    ConfigError
        If a name is not in the registry.
    """
    if not requested:
        return []
    if requested.strip().lower() == ALL:
        return list(registry)

    names = [name.strip() for name in requested.split(",") if name.strip()]
    unknown = [name for name in names if name not in registry]
    if unknown:
        import difflib

        details = []
        for name in unknown:
            close = difflib.get_close_matches(name, sorted(registry), n=1)
            details.append(
                f"{name!r}" + (f" (did you mean {close[0]!r}?)" if close else "")
            )
        raise ConfigError(f"unknown {axis}(s): {', '.join(details)}")
    return names


def preflight(cfg: RunConfig) -> dict[str, list[str]]:
    """
    Probe every selected method's backends once, before any job launches.

    Returns
    -------
    Mapping of method key to the backends it is missing. Empty means everything selected
    can actually run here.
    """
    from combobatch.methods import METHOD_REGISTRY

    unavailable = {}
    for selection in cfg.methods:
        missing = METHOD_REGISTRY[selection.name].missing_backends()
        if missing:
            unavailable[selection.name] = missing
    return unavailable


def warn_about_nesting(cfg: RunConfig, n_workers: int) -> None:
    """Warn when Shambhala's inner workers multiply against the dispatcher's own."""
    from combobatch.methods.shambhala_method import warn_nested_parallelism

    for selection in cfg.methods:
        if selection.name != "20_shambhala":
            continue
        inner = int(selection.params.get("n_workers", 1))
        warn_nested_parallelism(n_workers, inner)


def build_job_matrix(
    cfg: RunConfig,
) -> list[tuple[ImputationSelection, MethodSelection]]:
    """
    Return every (imputation, method) pair this configuration asks for.

    Two pairs that resolve to the same output key would have the second silently
    overwrite the first — or, under ``--skip-if-exists``, be reported as cached. Config
    validation catches the ordinary case, but only here is the *effective* tag known:
    restating a parameter at its default value produces no tag, so two selections that
    look different in YAML can still collide on disk.

    Raises
    ------
    ConfigError
        If two jobs would write to the same key, naming it.
    """
    jobs = cfg.combinations()

    seen: dict[str, tuple[str, str]] = {}
    for imp_sel, method_sel in jobs:
        key = job_key(imp_sel, method_sel)
        if key in seen:
            raise ConfigError(
                f"two jobs resolve to the same output key {key!r}: "
                f"{seen[key]} and {(imp_sel.name, method_sel.name)}. One would overwrite "
                f"the other — give them different parameters or a distinguishing "
                f"param_tag."
            )
        seen[key] = (imp_sel.name, method_sel.name)
    return jobs


def job_key(imp_sel: ImputationSelection, method_sel: MethodSelection) -> str:
    """Return the stable identifier used in logs and the failed-jobs file."""
    from combobatch.harmonize import combined_param_tag
    from combobatch.params import build_output_key

    return build_output_key(
        imp_sel.name, method_sel.name, combined_param_tag(imp_sel, method_sel), None
    )


def failed_log_key(cfg: RunConfig) -> str:
    """Return the storage key of this run's failed-jobs log."""
    return f"failed_jobs_{cfg.run_id}.txt"


def read_failed_log(backend: storage.StorageBackend, key: str) -> set[str]:
    """Read the failed-jobs log, treating an absent one as empty."""
    if not backend.exists(key):
        return set()
    text = backend.read_bytes(key).decode()
    return {line.strip() for line in text.splitlines() if line.strip()}


def write_failed_log(backend: storage.StorageBackend, key: str, keys: set[str]) -> None:
    """Write the failed-jobs log, deleting it when nothing failed."""
    if not keys:
        backend.delete(key)
        return
    backend.write_bytes(key, ("\n".join(sorted(keys)) + "\n").encode())


def _worker_command(
    cfg: RunConfig,
    imp_sel: ImputationSelection,
    method_sel: MethodSelection,
    *,
    config_path: str | None,
    out_json: str,
    memory_limit_gb: float,
    skip_if_exists: bool,
    log_level: str,
) -> list[str]:
    """
    Build the argv for one worker.

    The config file, when there was one, is passed through and the flags layered on top,
    so settings with no flag of their own (the metric column lists) still reach the
    worker.
    """
    from combobatch.harmonize import _resolved_max_na_frac

    command = [sys.executable, "-m", "combobatch", "run"]
    if config_path:
        command += ["--config", config_path]
    command += [
        "--exp",
        cfg.exp_uri,
        "--ann",
        cfg.ann_uri,
        "--out",
        cfg.out_uri,
        "--batch-col",
        cfg.columns.batch,
        "--bio-col",
        cfg.columns.bio,
        "--imputation",
        imp_sel.name,
        "--method",
        method_sel.name,
        "--max-na-frac",
        str(_resolved_max_na_frac(imp_sel)),
        "--random-seed",
        str(cfg.random_seed),
        "--run-id",
        cfg.run_id,
        "--out-json",
        out_json,
        "--memory-limit-gb",
        str(memory_limit_gb / 2),
        "--log-level",
        log_level,
        # The dispatcher writes the manifest once, up front; concurrent workers must not
        # race to read-modify-write the same file.
        "--no-manifest",
    ]
    if cfg.columns.cohort:
        command += ["--cohort-col", cfg.columns.cohort]
    if cfg.reference_batch:
        command += ["--reference-batch", cfg.reference_batch]
    if cfg.subset_query:
        command += ["--subset-query", cfg.subset_query]
    if cfg.endpoint_url:
        command += ["--endpoint-url", cfg.endpoint_url]
    if imp_sel.params:
        command += ["--impute-params", _params_flag(imp_sel.name, imp_sel.params)]
    if method_sel.params:
        command += ["--method-params", _params_flag(method_sel.name, method_sel.params)]
    if method_sel.param_tag:
        command += ["--param-tag", method_sel.param_tag]
    if cfg.post_removal.enabled:
        command += [
            "--post-removal",
            "--post-removal-n",
            str(cfg.post_removal.n_batches),
            "--post-removal-min-batch-size",
            str(cfg.post_removal.min_batch_size),
            "--post-removal-n-pcs",
            str(cfg.post_removal.n_pcs),
        ]
    if cfg.metrics.enabled:
        command += ["--metrics", "--groups", ",".join(sorted(cfg.metrics.groups))]
        if cfg.metrics.marker_panel:
            command += ["--marker-panel", cfg.metrics.marker_panel]
    if skip_if_exists:
        command.append("--skip-if-exists")
    return command


def _params_flag(target: str, params: dict[str, Any]) -> str:
    """Re-encode a parameter dict into the ``--*-params`` grammar the worker parses."""

    def encode(value: Any) -> str:
        if isinstance(value, (list, tuple)):
            return "+".join(str(item) for item in value)
        return str(value)

    body = ",".join(f"{name}={encode(params[name])}" for name in sorted(params))
    return f"{target}:{body}"


def _stream_output(pipe: Any, prefix: str) -> None:
    """Forward a worker's interleaved stdout/stderr to this process's stderr."""
    for line in iter(pipe.readline, ""):
        print(f"[{prefix}] {line.rstrip()}", file=sys.stderr, flush=True)
    pipe.close()


def run_dispatcher(
    cfg: RunConfig,
    *,
    config_path: str | None,
    n_workers: int,
    timeout_s: int,
    memory_limit_gb: float,
    memory_wait_s: float,
    skip_if_exists: bool,
    log_level: str,
    tmp_dir: str,
    only_keys: set[str] | None = None,
) -> tuple[list[dict[str, Any]], set[str], set[str]]:
    """
    Launch every combination, supervising one subprocess per job.

    Parameters
    ----------
    only_keys
        Restrict the matrix to these job keys. Used by ``--retry-failed``, which must
        rerun exactly the pairs that failed — not the cross-product of the axes those
        pairs happen to mention.

    Returns
    -------
    ``(rows, failed_keys, succeeded_keys)`` — the concatenated sidecar rows plus the two
    key sets the failed-jobs log is rebuilt from.
    """
    import os

    from combobatch import dataio
    from combobatch.memory import check_memory, free_gb, wait_for_memory

    log = get_logger()
    jobs = build_job_matrix(cfg)
    if only_keys is not None:
        jobs = [pair for pair in jobs if job_key(*pair) in only_keys]
    total = len(jobs)
    completed = 0
    lock = threading.Lock()

    def launch(imp_sel: ImputationSelection, method_sel: MethodSelection) -> JobResult:
        key = job_key(imp_sel, method_sel)
        if not check_memory(memory_limit_gb):
            log.info(
                "%s waiting for %.1f GB free (%.1f GB now)",
                key,
                memory_limit_gb,
                free_gb(),
            )
            # Bounded, and it launches anyway on expiry. An unbounded wait is how the
            # donor's dispatcher could sit silently forever on a machine that simply
            # never has that much free; the worker's own guard makes the honest call.
            if not wait_for_memory(memory_limit_gb, timeout_s=memory_wait_s):
                log.warning(
                    "%s: still only %.1f GB free after %.0fs — launching anyway",
                    key,
                    free_gb(),
                    memory_wait_s,
                )
        started = time.time()
        log.info("START %s | %.1f GB free", key, free_gb())

        out_json = os.path.join(tmp_dir, f"{key}.json")
        command = _worker_command(
            cfg,
            imp_sel,
            method_sel,
            config_path=config_path,
            out_json=out_json,
            memory_limit_gb=memory_limit_gb,
            skip_if_exists=skip_if_exists,
            log_level=log_level,
        )

        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
        )
        reader = threading.Thread(
            target=_stream_output, args=(process.stdout, key), daemon=True
        )
        reader.start()
        try:
            process.wait(timeout=timeout_s)
        except subprocess.TimeoutExpired:
            process.kill()
            reader.join(timeout=5)
            elapsed = time.time() - started
            log.error("TIMEOUT %s after %.0fs", key, elapsed)
            return JobResult(key, EXIT_TIMEOUT, elapsed, [])
        reader.join(timeout=5)

        elapsed = time.time() - started
        rows: list[dict[str, Any]] = []
        if os.path.exists(out_json):
            backend, name = storage.resolve(out_json)
            rows = dataio.read_json(backend, name)
        elif process.returncode == EXIT_INSUFFICIENT_MEMORY:
            log.error("%s exited for want of memory", key)
        elif process.returncode != 0:
            log.error("%s exited with code %d", key, process.returncode)

        return JobResult(key, process.returncode, elapsed, rows)

    all_rows: list[dict[str, Any]] = []
    failed: set[str] = set()
    succeeded: set[str] = set()

    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        futures = {
            pool.submit(launch, imp, method): (imp, method) for imp, method in jobs
        }
        for future in as_completed(futures):
            result = future.result()
            with lock:
                completed += 1
                position = completed
            if result.succeeded:
                succeeded.add(result.key)
                all_rows.extend(result.rows)
                statuses = [row.get("status") for row in result.rows]
                log.info(
                    "[%d/%d] %s (%.0fs) -> %s",
                    position,
                    total,
                    result.key,
                    result.elapsed_s,
                    statuses,
                )
            else:
                failed.add(result.key)
                all_rows.extend(result.rows)
                log.error(
                    "[%d/%d] FAILED %s (%.0fs)",
                    position,
                    total,
                    result.key,
                    result.elapsed_s,
                )

    return all_rows, failed, succeeded


def print_dry_run(cfg: RunConfig, unavailable: dict[str, list[str]]) -> None:
    """Print the resolved job matrix, its output count, and any missing backend."""
    from combobatch.harmonize import MATRIX_SUFFIX, combined_param_tag
    from combobatch.params import build_output_key

    jobs = build_job_matrix(cfg)
    variants = [False, True] if cfg.post_removal.enabled else [None]

    print(f"Output root: {cfg.out_uri}")
    print(
        f"Jobs: {len(jobs)}  ({len(cfg.imputations)} imputation(s) x {len(cfg.methods)} method(s))"
    )
    print(f"Outputs: {len(jobs) * len(variants)} matrices")
    print()
    for imp_sel, method_sel in jobs:
        tag = combined_param_tag(imp_sel, method_sel)
        for variant in variants:
            key = build_output_key(imp_sel.name, method_sel.name, tag, variant)
            print(f"  exp/{key}{MATRIX_SUFFIX}")

    if unavailable:
        print()
        print("These methods cannot run in this environment and will skip:")
        for name, missing in sorted(unavailable.items()):
            print(f"  {name}: needs {', '.join(missing)}")


def main(argv: list[str] | None = None) -> int:
    """Run the cross-product. Returns the process exit code."""
    import tempfile

    parser = argparse.ArgumentParser(prog="combobatch dispatch")
    add_arguments(parser)
    args = parser.parse_args(argv)

    pin_threads(args.threads)
    setup_logging(args.log_level)
    log = get_logger()

    from combobatch.imputation import IMPUTER_REGISTRY
    from combobatch.methods import METHOD_REGISTRY

    try:
        # `all` is expanded here so that build_config sees a plain comma-separated list
        # and the two commands share one code path for everything else.
        if args.imputations:
            args.imputations = ",".join(
                expand_axis(args.imputations, IMPUTER_REGISTRY, "imputer")
            )
        if args.methods:
            args.methods = ",".join(
                expand_axis(args.methods, METHOD_REGISTRY, "method")
            )
        cfg = run_cmd.build_config(args)
        cfg.validate(known_methods=METHOD_REGISTRY, known_imputers=IMPUTER_REGISTRY)
    except ConfigError as exc:
        log.error("%s", exc)
        return EXIT_MISSING_INPUT

    unavailable = preflight(cfg)
    warn_about_nesting(cfg, args.n_workers)

    if args.dry_run:
        print_dry_run(cfg, unavailable)
        return run_cmd.EXIT_OK

    for name, missing in sorted(unavailable.items()):
        log.warning("%s will skip: needs %s", name, ", ".join(missing))

    from combobatch import dataio
    from combobatch.harmonize import write_manifest

    backend = storage.backend_for_root(cfg.out_uri, endpoint_url=cfg.endpoint_url)
    log_key = failed_log_key(cfg)
    previously_failed = read_failed_log(backend, log_key)

    only_keys: set[str] | None = None
    if args.retry_failed:
        if not previously_failed:
            log.info("nothing in %s to retry", log_key)
            return run_cmd.EXIT_OK
        only_keys = {
            job_key(imp, method)
            for imp, method in cfg.combinations()
            if job_key(imp, method) in previously_failed
        }
        unknown = previously_failed - only_keys
        if unknown:
            # Loudly, because the donor silently dropped any key it could not parse and
            # a dropped key means a job that quietly never reruns.
            log.warning(
                "%d key(s) in %s are not in this configuration's matrix and will not "
                "be retried: %s",
                len(unknown),
                log_key,
                ", ".join(sorted(unknown)),
            )
        log.info("retrying %d previously failed job(s)", len(only_keys))

    write_manifest(cfg, backend)

    failed_so_far: set[str] = set(previously_failed)

    def flush_and_exit(signum, _frame):  # pragma: no cover - signal path
        """Karpenter's two-minute warning: the resume state must survive it."""
        log.warning("signal %d — flushing the failed-jobs log", signum)
        write_failed_log(backend, log_key, failed_so_far)
        sys.exit(EXIT_TIMEOUT)

    signal.signal(signal.SIGTERM, flush_and_exit)
    signal.signal(signal.SIGINT, flush_and_exit)

    started = time.time()
    with tempfile.TemporaryDirectory(prefix="combobatch-") as tmp_dir:
        rows, failed, succeeded = run_dispatcher(
            cfg,
            config_path=args.config,
            n_workers=args.n_workers,
            timeout_s=args.timeout_s,
            memory_limit_gb=args.memory_limit_gb,
            memory_wait_s=args.memory_wait_s,
            skip_if_exists=args.skip_if_exists,
            log_level=args.log_level,
            tmp_dir=tmp_dir,
            only_keys=only_keys,
        )

    failed_so_far = (previously_failed - succeeded) | failed
    write_failed_log(backend, log_key, failed_so_far)
    if rows:
        dataio.write_json(rows, backend, f"dispatch_rows_{cfg.run_id}.json")

    log.info(
        "Done in %.0fs: %d succeeded, %d failed.%s",
        time.time() - started,
        len(succeeded),
        len(failed),
        f" Retry them with --retry-failed --run-id {cfg.run_id}." if failed else "",
    )
    return run_cmd.EXIT_OK if not failed else run_cmd.EXIT_FAILED


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
