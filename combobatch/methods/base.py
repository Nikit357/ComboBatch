"""Method descriptors and helpers shared by every harmonization method.

``MethodSpec`` is what the donor's bare ``(function, tier)`` tuple could not express:
which backend a method needs, whether it is restricted to RNA-seq, and — most
importantly — which hyperparameters it accepts. That declaration is what lets the CLI
validate ``--method-params`` before a job launches, lets the dispatcher preflight R and
Octave once instead of failing hundreds of jobs in, and lets the docs be generated rather
than hand-maintained.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Callable

import pandas as pd

from combobatch.logging_utils import get_logger

HARSHNESS_TIERS = ("low", "medium", "high")


@dataclass(frozen=True)
class HyperParam:
    """One tunable parameter of a method or imputer."""

    name: str
    type: type
    default: Any
    help: str
    choices: tuple[Any, ...] | None = None
    minimum: float | None = None
    maximum: float | None = None
    # Set when the parameter is forwarded into an R call rather than consumed in
    # Python, so docs can say where it actually lands.
    r_argument: str | None = None

    def validate(self, value: Any) -> None:
        """
        Check one supplied value against this declaration.

        Raises
        ------
        ValueError
            If the value is outside the declared choices or bounds.
        """
        if self.choices is not None and value not in self.choices:
            allowed = ", ".join(repr(choice) for choice in self.choices)
            raise ValueError(f"{self.name}: expected one of {allowed}, got {value!r}")
        if self.minimum is not None and value < self.minimum:
            raise ValueError(f"{self.name}: must be >= {self.minimum}, got {value!r}")
        if self.maximum is not None and value > self.maximum:
            raise ValueError(f"{self.name}: must be <= {self.maximum}, got {value!r}")


@dataclass(frozen=True)
class MethodSpec:
    """Everything ComboBatch knows about one harmonization method."""

    key: str
    fn: Callable[..., pd.DataFrame]
    harshness: str
    citation: str = ""
    requires_r: bool = False
    r_packages: tuple[str, ...] = ()
    python_packages: tuple[str, ...] = ()
    requires_octave: bool = False
    rnaseq_only: bool = False
    uses_reference_batch: bool = False
    uses_bio_col: bool = False
    hyperparams: dict[str, HyperParam] = field(default_factory=dict)
    notes: str = ""

    def param_names(self) -> set[str]:
        """Return the set of hyperparameter names this method declares."""
        return set(self.hyperparams)

    def defaults(self) -> dict[str, Any]:
        """Return every declared hyperparameter at its default value."""
        return {name: spec.default for name, spec in self.hyperparams.items()}

    def non_default(self, params: dict[str, Any]) -> dict[str, Any]:
        """
        Return only the entries of ``params`` that differ from the declared default.

        This is what goes into the output filename tag: a run at pure defaults produces
        an empty dict and therefore no tag segment at all.
        """
        defaults = self.defaults()
        return {
            name: value
            for name, value in params.items()
            if name not in defaults or defaults[name] != value
        }

    def validate_params(self, params: dict[str, Any]) -> None:
        """
        Validate supplied hyperparameters against this method's declarations.

        Raises
        ------
        ValueError
            If a name is undeclared or a value is out of range.
        """
        from combobatch.params import validate_param_names

        validate_param_names(self.key, params, self.param_names())
        for name, value in params.items():
            self.hyperparams[name].validate(value)

    def missing_backends(self) -> list[str]:
        """
        Return human-readable names of backends this method needs but cannot find.

        An empty list means the method can run here. Used by the dispatcher preflight
        and by the test suite, so that an unavailable backend produces an honest skip
        rather than a silent pass.
        """
        from combobatch.methods import rinterop

        missing: list[str] = []
        if self.requires_r and not rinterop.r_available():
            missing.append("R (rpy2)")
        else:
            missing.extend(
                f"R package {name}"
                for name in self.r_packages
                if not rinterop.r_package_available(name)
            )
        missing.extend(
            f"Python package {name}"
            for name in self.python_packages
            if not _python_package_available(name)
        )
        if self.requires_octave and not _octave_available():
            missing.append("Octave")
        return missing

    def is_available(self) -> bool:
        """Report whether every backend this method needs is present."""
        return not self.missing_backends()


def _python_package_available(name: str) -> bool:
    """Report whether an optional Python package can be imported."""
    import importlib.util

    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def _octave_available() -> bool:
    """Report whether an Octave executable is on PATH."""
    import shutil

    return shutil.which("octave") is not None


def resolve_reference_batch(
    ann_df: pd.DataFrame, batch_col: str, requested: str | None = None
) -> str:
    """
    Decide which batch a reference-based method should normalise toward.

    Nine donor methods defaulted this to one dataset's batch name, which silently
    degraded to a fallback path on anyone else's data. Here it is either what the user
    asked for or, failing that, the largest batch — a defensible, dataset-agnostic
    choice that is logged so it never happens invisibly.

    Parameters
    ----------
    ann_df
        Annotation aligned to the expression matrix.
    batch_col
        Column holding batch identity.
    requested
        The user's ``--reference-batch``, if any.

    Returns
    -------
    The chosen batch label.

    Raises
    ------
    ValueError
        If ``requested`` is not a value of ``batch_col``.
    """
    counts = ann_df[batch_col].astype(str).value_counts()
    if requested is not None:
        if requested not in counts.index:
            raise ValueError(
                f"reference batch {requested!r} is not a value of {batch_col!r} "
                f"({len(counts)} distinct values present)"
            )
        return requested

    chosen = str(counts.index[0])
    get_logger().info(
        "No reference batch given; using the largest, %r (n=%d)", chosen, counts.iloc[0]
    )
    return chosen


def drop_na_genes(exp_df: pd.DataFrame) -> tuple[pd.DataFrame, int]:
    """
    Drop genes containing any NaN, reporting how many were removed.

    Most methods cannot accept missing values. The donor applied this guard inline and
    discarded the count, so a method could quietly return a much smaller gene space than
    it was given; here the number is returned and recorded in the sidecar.

    Returns
    -------
    ``(clean_expression, n_genes_dropped)``.
    """
    if not exp_df.isnull().any().any():
        return exp_df, 0
    cleaned = exp_df.dropna(axis=1)
    return cleaned, exp_df.shape[1] - cleaned.shape[1]


def assert_rnaseq_only(
    ann_df: pd.DataFrame, batch_col: str, method_name: str, prefix: str = "GPL"
) -> None:
    """
    Raise ``NotImplementedError`` if any batch label looks like a microarray platform.

    ``prefix`` defaults to the GEO platform prefix, which is the near-universal
    convention for microarray batch names, but is configurable because it is a
    heuristic over the user's own labels rather than a fact about the data.
    """
    labels = ann_df[batch_col].astype(str)
    if labels.str.startswith(prefix).any():
        offending = sorted(set(labels[labels.str.startswith(prefix)]))[:3]
        raise NotImplementedError(
            f"{method_name} is RNA-seq only, but batch labels starting with "
            f"{prefix!r} are present (e.g. {', '.join(offending)}). Restrict the run "
            f"with --subset-query, or choose a cross-platform method."
        )
