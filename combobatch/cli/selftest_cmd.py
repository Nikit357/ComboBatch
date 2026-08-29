"""``combobatch selftest`` — run every method and imputer on synthetic data.

This is the image's acceptance gate. It answers one question: does this environment
actually have everything the registries claim it needs? Exit code 0 means every method
either PASSed or SKIPped for a named, legitimate reason.

The distinction between SKIP and FAIL is the whole point, and it is subtle. A method that
raises ``NotImplementedError`` is declaring "I cannot run on this input", which is honest.
But rpy2 also surfaces *type-conversion bugs* as ``NotImplementedError`` — and that is a
broken installation wearing an honest skip's clothing. ``_is_rpy2_conversion_error``
separates them, and the second case is a FAIL.

Backends are probed from the registry first, so a missing R package produces a named skip
("needs R package sva") rather than an exception whose text has to be interpreted.
"""

from __future__ import annotations

import argparse
import time
import traceback
import warnings
from dataclasses import dataclass
from typing import Any, Callable

import numpy as np
import pandas as pd

# Matches the donor's fixture exactly: 80 samples, 200 genes, 3 batches, 2 diagnoses.
# The +2.0 offset matters - normalize_vst rounds and clips to 0, and a zero count there
# makes DESeq2's size-factor estimation fail for reasons unrelated to the installation.
N_SAMPLES = 80
N_GENES = 200
N_BATCHES = 3
N_DIAGNOSES = 2
NA_FRACTION = 0.10
RANDOM_SEED = 42

BATCH_COL = "batch"
BIO_COL = "bio"

PASS, SKIP, FAIL = "PASS", "SKIP", "FAIL"


@dataclass
class Result:
    """One method's or imputer's outcome."""

    name: str
    status: str
    elapsed: float
    detail: str = ""
    traceback: str = ""


def make_synthetic() -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build the synthetic expression matrix and annotation.

    Returns
    -------
    ``(exp_df, ann_df)`` — 80 samples x 200 genes, raw scale, with a batch and a biology
    column, aligned on the index.
    """
    rng = np.random.default_rng(RANDOM_SEED)
    samples = [f"S{i:04d}" for i in range(N_SAMPLES)]
    exp_df = pd.DataFrame(
        rng.exponential(5.0, size=(N_SAMPLES, N_GENES)) + 2.0,
        index=samples,
        columns=[f"GENE{j:04d}" for j in range(N_GENES)],
    )
    ann_df = pd.DataFrame(
        {
            BATCH_COL: [f"Batch{i % N_BATCHES}" for i in range(N_SAMPLES)],
            BIO_COL: [f"Class{i % N_DIAGNOSES}" for i in range(N_SAMPLES)],
        },
        index=samples,
    )
    return exp_df, ann_df


def punch_holes(exp_df: pd.DataFrame) -> pd.DataFrame:
    """
    Return a copy with NAs confined to a quarter of the genes.

    Deliberately *not* uniform across the matrix. Real missingness is platform-driven and
    therefore gene-wise, and uniform holes at any appreciable rate leave every gene with
    at least one NA — which makes `strict` correctly keep nothing and turns the whole
    imputation axis into a degenerate test.
    """
    rng = np.random.default_rng(RANDOM_SEED)
    holed = exp_df.copy()
    affected = rng.choice(
        exp_df.columns,
        size=max(1, int(NA_FRACTION * 2.5 * exp_df.shape[1])),
        replace=False,
    )
    mask = rng.random((len(holed), len(affected))) < 0.3
    holed.loc[:, affected] = holed[affected].mask(mask)
    return holed


def _is_rpy2_conversion_error(exc: BaseException) -> bool:
    """
    True if a ``NotImplementedError`` is really an rpy2 type-conversion failure.

    rpy2 raises ``NotImplementedError`` for an unconvertible object, which is
    indistinguishable by type from a method's deliberate "not applicable here". Treating
    the former as a skip is how a broken image passes its own acceptance gate.
    """
    message = str(exc)
    return any(
        marker in message
        for marker in ("Conversion", "conversion", "py2rpy", "rpy2py", "converter")
    )


def _run_one(name: str, call: Callable[[], Any], *, require: bool = False) -> Result:
    """Run one callable and classify the outcome."""
    started = time.time()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            out = call()
        if not isinstance(out, pd.DataFrame):
            raise TypeError(f"returned {type(out).__name__}, not a DataFrame")
        if out.shape[0] == 0:
            raise ValueError("returned 0 rows")
        if out.shape[1] == 0:
            # Distinct from all-NaN, and the more likely defect: an over-eager gene
            # filter empties the matrix without ever raising.
            raise ValueError("returned 0 genes")
        if out.isna().all().all():
            raise ValueError("returned an all-NaN matrix")
        return Result(name, PASS, time.time() - started)
    except NotImplementedError as exc:
        if _is_rpy2_conversion_error(exc):
            return Result(
                name,
                FAIL,
                time.time() - started,
                f"rpy2 conversion failure disguised as a skip: {exc}",
                traceback.format_exc(),
            )
        status = FAIL if require else SKIP
        return Result(
            name, status, time.time() - started, str(exc), traceback.format_exc()
        )
    except Exception as exc:  # noqa: BLE001 - the point is to classify anything
        return Result(
            name, FAIL, time.time() - started, str(exc), traceback.format_exc()
        )


def check_imputers(exp_df: pd.DataFrame, *, require_r: bool) -> list[Result]:
    """Run every registered imputer on a matrix with holes punched in it."""
    from combobatch.imputation import IMPUTER_REGISTRY, impute

    holed = punch_holes(exp_df)
    results = []
    for key, spec in IMPUTER_REGISTRY.items():
        missing = spec.missing_backends()
        if missing and not require_r:
            results.append(Result(key, SKIP, 0.0, f"needs {', '.join(missing)}"))
            continue
        if missing:
            results.append(Result(key, FAIL, 0.0, f"needs {', '.join(missing)}"))
            continue
        results.append(
            _run_one(key, lambda k=key: impute(holed, method=k)[0], require=False)
        )
    return results


def check_methods(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    require_r: bool,
    require_octave: bool,
) -> list[Result]:
    """Run every registered method on the strict-prepared synthetic matrix."""
    from combobatch.methods import METHOD_REGISTRY

    results = []
    for key, spec in METHOD_REGISTRY.items():
        missing = spec.missing_backends()
        if missing:
            # A missing backend is normally a legitimate skip. --with-r / --with-octave
            # turn it into a failure, which is what the image should be checked against:
            # there, a missing backend means the build silently lost a method.
            demanded = (require_r and spec.requires_r) or (
                require_octave and spec.requires_octave
            )
            results.append(
                Result(
                    key,
                    FAIL if demanded else SKIP,
                    0.0,
                    f"needs {', '.join(missing)}",
                )
            )
            continue

        params = spec.defaults()
        if spec.uses_reference_batch and params.get("target_group") is None:
            params["target_group"] = ann_df[BATCH_COL].value_counts().idxmax()

        results.append(
            _run_one(
                key,
                lambda s=spec, p=params: s.fn(
                    exp_df, ann_df, batch_col=BATCH_COL, bio_col=BIO_COL, **p
                ),
            )
        )
    return results


def describe_environment() -> list[str]:
    """Collect version lines for the report, never raising if a backend is absent."""
    import sys

    lines = [f"Python        {sys.version.split()[0]}"]

    try:
        import rpy2

        lines.append(f"rpy2          {rpy2.__version__}")
    except Exception:
        lines.append("rpy2          not installed")

    from combobatch.methods import rinterop

    if rinterop.r_available():
        try:
            import rpy2.robjects as ro

            version = str(ro.r("R.version$version.string")[0])
        except Exception as exc:  # pragma: no cover - only when R is half-installed
            version = f"unknown ({exc})"
        lines.append(f"R             {version}")
        # 4.6 has a C-level ABI incompatibility with rpy2 3.6.x whose symptom is a
        # segfault, so a version check here is worth more than it looks.
        if "4.5" not in version:
            lines.append("              WARNING: expected R 4.5.x; rpy2 may segfault")
    else:
        lines.append("R             not reachable")

    from combobatch.methods.base import _octave_available

    lines.append(f"Octave        {'present' if _octave_available() else 'not found'}")
    return lines


def add_arguments(parser: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Register ``selftest`` arguments."""
    parser.add_argument(
        "--with-r",
        action="store_true",
        help="Treat a missing R backend as a failure rather than a skip.",
    )
    parser.add_argument(
        "--with-octave",
        action="store_true",
        help="Treat a missing Octave backend as a failure rather than a skip.",
    )
    parser.add_argument(
        "--json",
        action="store_true",
        help="Emit the results as JSON instead of a table.",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="Print the traceback of every skip, not just every failure.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """Run the smoke test. Returns 0 when nothing FAILed."""
    parser = add_arguments(
        argparse.ArgumentParser(
            prog="combobatch selftest",
            description="Run every method and imputer on synthetic data.",
        )
    )
    args = parser.parse_args(argv)

    from combobatch.dataio import dumps_json
    from combobatch.logging_utils import pin_threads

    pin_threads(1)

    exp_df, ann_df = make_synthetic()
    imputers = check_imputers(exp_df, require_r=args.with_r)
    methods = check_methods(
        exp_df, ann_df, require_r=args.with_r, require_octave=args.with_octave
    )
    everything = imputers + methods
    failures = [r for r in everything if r.status == FAIL]

    if args.json:
        print(
            dumps_json(
                {
                    "imputers": [vars(r) for r in imputers],
                    "methods": [vars(r) for r in methods],
                    "n_failed": len(failures),
                }
            )
        )
        return 1 if failures else 0

    print("=" * 64)
    print("ENVIRONMENT")
    print("=" * 64)
    for line in describe_environment():
        print(f"  {line}")

    for title, results in (("IMPUTATION", imputers), ("HARMONIZATION", methods)):
        print()
        print("=" * 64)
        print(f"{title}  ({len(results)})")
        print("=" * 64)
        for result in results:
            detail = f"  {result.detail[:70]}" if result.status != PASS else ""
            print(
                f"  {result.name:<22} {result.status}  ({result.elapsed:5.1f}s){detail}"
            )

    if failures:
        print()
        print("=" * 64)
        print("FAILURE DETAILS")
        print("=" * 64)
        for result in failures:
            print(f"\n{'-' * 64}\nFAILED: {result.name}\n{'-' * 64}")
            print(result.traceback or result.detail)

    if args.verbose:
        for result in everything:
            if result.status == SKIP and result.traceback:
                print(f"\nSKIP {result.name}: {result.detail}")

    counts = {
        status: sum(1 for r in everything if r.status == status)
        for status in (PASS, SKIP, FAIL)
    }
    print()
    print(
        f"RESULT: {counts[PASS]} passed, {counts[SKIP]} skipped, {counts[FAIL]} failed"
    )
    if failures:
        print(f"FAILED: {', '.join(r.name for r in failures)}")
        return 1
    return 0
