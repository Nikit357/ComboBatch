"""``IMPUTER_REGISTRY`` — the four imputation strategies behind one signature.

The donor split these across two functions: ``strict`` lived in ``prepare_dataset`` with
no parameters at all, while the other three lived in ``prepare_dataset_imputed`` and
raised ``ValueError`` if handed ``"strict"``. Here all four are registry entries obeying
one contract, which is what lets ``--imputations strict,knn,softimpute,missforest`` be a
homogeneous list.

Three donor behaviours are deliberately not reproduced:

* **The silent fallback.** Its workers wrapped imputation in ``try/except`` and fell back
  to ``strict`` on any failure *while still labelling the output* ``__knn__``. A run
  labelled ``knn`` here either used KNN or fails.
* **The invisible no-op.** If the NA-fraction filter already removed every NA, it
  returned before imputing, so a ``knn`` run could be byte-identical to a ``strict`` one
  with nothing to show for it. That case is now recorded.
* **Unreachable parameters.** ``knn_k`` existed but no caller passed it; softImpute's and
  missForest's arguments were not in the Python signature at all. All are exposed.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Any, Callable

import numpy as np
import pandas as pd

from combobatch.methods.base import HyperParam
from combobatch.methods.rinterop import (
    as_numpy,
    py2rpy,
    r_arglist,
    r_gc,
    require_r_packages,
)

DEFAULT_MAX_NA_FRAC = 0.20


@dataclass(frozen=True)
class ImputerSpec:
    """Everything ComboBatch knows about one imputation strategy."""

    key: str
    fn: Callable[..., pd.DataFrame]
    description: str
    requires_r: bool = False
    r_packages: tuple[str, ...] = ()
    hyperparams: dict[str, HyperParam] = field(default_factory=dict)

    def param_names(self) -> set[str]:
        """Return the hyperparameter names this imputer declares."""
        return set(self.hyperparams)

    def defaults(self) -> dict[str, Any]:
        """Return every declared hyperparameter at its default value."""
        return {name: spec.default for name, spec in self.hyperparams.items()}

    def non_default(self, params: dict[str, Any]) -> dict[str, Any]:
        """Return only the entries of ``params`` that differ from the defaults."""
        defaults = self.defaults()
        return {
            name: value
            for name, value in params.items()
            if name not in defaults or defaults[name] != value
        }

    def validate_params(self, params: dict[str, Any]) -> None:
        """Validate supplied hyperparameters, raising on an unknown or out-of-range one."""
        from combobatch.params import validate_param_names

        validate_param_names(self.key, params, self.param_names())
        for name, value in params.items():
            self.hyperparams[name].validate(value)

    def missing_backends(self) -> list[str]:
        """Return names of backends this imputer needs but cannot find."""
        from combobatch.methods import rinterop

        if not self.requires_r:
            return []
        if not rinterop.r_available():
            return ["R (rpy2)"]
        return [
            f"R package {name}"
            for name in self.r_packages
            if not rinterop.r_package_available(name)
        ]

    def is_available(self) -> bool:
        """Report whether this imputer can run here."""
        return not self.missing_backends()


@dataclass(frozen=True)
class ImputationReport:
    """What imputation did, for the run's sidecar."""

    method: str
    n_genes_in: int
    n_genes_dropped_na_frac: int
    n_genes_out: int
    n_cells_imputed: int
    residual_na: int
    was_noop: bool
    elapsed_s: float

    def as_dict(self) -> dict[str, Any]:
        """Return the report as a flat dict for the JSON sidecar."""
        return {
            "imputation_method": self.method,
            "imputation_n_genes_in": self.n_genes_in,
            "imputation_n_genes_dropped_na_frac": self.n_genes_dropped_na_frac,
            "imputation_n_genes_out": self.n_genes_out,
            "imputation_n_cells_imputed": self.n_cells_imputed,
            "imputation_residual_na": self.residual_na,
            "imputation_was_noop": self.was_noop,
            "imputation_elapsed_s": round(self.elapsed_s, 2),
        }


def _impute_strict(exp_df: pd.DataFrame, **params) -> pd.DataFrame:
    """Return the matrix unchanged; the NA-fraction filter has already done the work."""
    return exp_df


def _impute_knn(
    exp_df: pd.DataFrame,
    *,
    knn_k: int = 5,
    weights: str = "uniform",
    metric: str = "nan_euclidean",
    **params,
) -> pd.DataFrame:
    """Impute with scikit-learn's KNNImputer."""
    from sklearn.impute import KNNImputer

    imputer = KNNImputer(n_neighbors=knn_k, weights=weights, metric=metric)
    values = imputer.fit_transform(exp_df.values)
    return pd.DataFrame(values, index=exp_df.index, columns=exp_df.columns)


def _impute_softimpute(
    exp_df: pd.DataFrame,
    *,
    rank_max: int | None = None,
    lambda_: float | None = None,
    thresh: float | None = None,
    maxit: int | None = None,
    type: str = "svd",
    **params,
) -> pd.DataFrame:
    """Impute with R's softImpute, forwarding whichever arguments were supplied."""
    require_r_packages("softimpute", ("softImpute",))
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    base = importr("base")
    importr("softImpute")

    # softImpute's argument names are not valid Python identifiers ("rank.max",
    # "lambda"), so the Python names are mapped across here rather than in the spec.
    supplied = {
        "rank.max": rank_max,
        "lambda": lambda_,
        "thresh": thresh,
        "maxit": maxit,
    }
    r_args: dict[str, Any] = {"type": type}
    r_args.update(
        {name: value for name, value in supplied.items() if value is not None}
    )

    ro.globalenv["cb_mat"] = base.as_matrix(py2rpy(exp_df))
    ro.r(f"cb_fit <- softImpute::softImpute(cb_mat, {r_arglist(r_args)})")
    completed = as_numpy(ro.r("softImpute::complete(cb_mat, cb_fit)"))
    ro.r("rm(cb_mat, cb_fit)")
    r_gc()
    return pd.DataFrame(completed, index=exp_df.index, columns=exp_df.columns)


def _impute_missforest(
    exp_df: pd.DataFrame,
    *,
    ntree: int | None = None,
    maxiter: int | None = None,
    mtry: int | None = None,
    **params,
) -> pd.DataFrame:
    """Impute with R's missForest, forwarding whichever arguments were supplied."""
    require_r_packages("missforest", ("missForest",))
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    importr("missForest")

    r_args = {
        name: value
        for name, value in (("ntree", ntree), ("maxiter", maxiter), ("mtry", mtry))
        if value is not None
    }
    ro.globalenv["cb_mat"] = py2rpy(exp_df)
    ro.r(
        f"cb_fit <- missForest::missForest(cb_mat{r_arglist(r_args, leading_comma=True)})"
    )
    imputed = as_numpy(ro.r("cb_fit$ximp"))
    ro.r("rm(cb_mat, cb_fit)")
    r_gc()
    return pd.DataFrame(imputed, index=exp_df.index, columns=exp_df.columns)


_MAX_NA_FRAC = HyperParam(
    "max_na_frac",
    float,
    DEFAULT_MAX_NA_FRAC,
    "Drop genes whose NA fraction exceeds this before imputing. 0.0 means drop any "
    "gene with a single missing value, which is what 'strict' does.",
    minimum=0.0,
    maximum=1.0,
)

IMPUTER_REGISTRY: dict[str, ImputerSpec] = {
    "strict": ImputerSpec(
        key="strict",
        fn=_impute_strict,
        description=(
            "No imputation: every gene with a missing value is dropped. The "
            "max_na_frac = 0.0 corner of one uniform rule, not a separate code path."
        ),
        hyperparams={
            "max_na_frac": HyperParam(
                "max_na_frac",
                float,
                0.0,
                "Fixed at 0.0 for strict: a gene with any NA is dropped.",
                minimum=0.0,
                maximum=0.0,
            )
        },
    ),
    "knn": ImputerSpec(
        key="knn",
        fn=_impute_knn,
        description="k-nearest-neighbour imputation (scikit-learn).",
        hyperparams={
            "max_na_frac": _MAX_NA_FRAC,
            "knn_k": HyperParam(
                "knn_k", int, 5, "Neighbours used per missing value.", minimum=1
            ),
            "weights": HyperParam(
                "weights",
                str,
                "uniform",
                "Neighbour weighting.",
                choices=("uniform", "distance"),
            ),
            "metric": HyperParam(
                "metric",
                str,
                "nan_euclidean",
                "Distance metric tolerant of missing values.",
                choices=("nan_euclidean",),
            ),
        },
    ),
    "softimpute": ImputerSpec(
        key="softimpute",
        fn=_impute_softimpute,
        description="Soft-thresholded SVD matrix completion (R softImpute).",
        requires_r=True,
        r_packages=("softImpute",),
        hyperparams={
            "max_na_frac": _MAX_NA_FRAC,
            "rank_max": HyperParam(
                "rank_max",
                int,
                None,
                "Maximum solution rank.",
                minimum=1,
                r_argument="softImpute(rank.max=)",
            ),
            "lambda_": HyperParam(
                "lambda_",
                float,
                None,
                "Nuclear-norm penalty.",
                minimum=0.0,
                r_argument="softImpute(lambda=)",
            ),
            "thresh": HyperParam(
                "thresh",
                float,
                None,
                "Convergence threshold.",
                minimum=0.0,
                r_argument="softImpute(thresh=)",
            ),
            "maxit": HyperParam(
                "maxit",
                int,
                None,
                "Maximum iterations.",
                minimum=1,
                r_argument="softImpute(maxit=)",
            ),
            "type": HyperParam(
                "type",
                str,
                "svd",
                "Algorithm variant.",
                choices=("svd", "als"),
                r_argument="softImpute(type=)",
            ),
        },
    ),
    "missforest": ImputerSpec(
        key="missforest",
        fn=_impute_missforest,
        description="Random-forest iterative imputation (R missForest).",
        requires_r=True,
        r_packages=("missForest",),
        hyperparams={
            "max_na_frac": _MAX_NA_FRAC,
            "ntree": HyperParam(
                "ntree",
                int,
                None,
                "Trees per forest.",
                minimum=1,
                r_argument="missForest(ntree=)",
            ),
            "maxiter": HyperParam(
                "maxiter",
                int,
                None,
                "Maximum imputation iterations.",
                minimum=1,
                r_argument="missForest(maxiter=)",
            ),
            "mtry": HyperParam(
                "mtry",
                int,
                None,
                "Variables sampled per split.",
                minimum=1,
                r_argument="missForest(mtry=)",
            ),
        },
    ),
}


def get_imputer(key: str) -> ImputerSpec:
    """Look up an imputer, raising a helpful error for an unknown key."""
    if key not in IMPUTER_REGISTRY:
        raise KeyError(
            f"unknown imputer {key!r}. Available: {', '.join(sorted(IMPUTER_REGISTRY))}"
        )
    return IMPUTER_REGISTRY[key]


def impute(
    exp_df: pd.DataFrame,
    *,
    method: str = "strict",
    max_na_frac: float | None = None,
    params: dict[str, Any] | None = None,
) -> tuple[pd.DataFrame, ImputationReport]:
    """
    Drop over-missing genes, then impute what remains.

    Parameters
    ----------
    exp_df
        Expression, samples x genes.
    method
        An ``IMPUTER_REGISTRY`` key.
    max_na_frac
        Genes whose NA fraction exceeds this are dropped before imputing. ``None`` uses
        the imputer's default, which is 0.0 for ``strict``.
    params
        Imputer hyperparameters, validated against the registry.

    Returns
    -------
    ``(imputed_expression, report)``. The report carries the gene counts, how many cells
    were filled, any residual NA, and whether the call turned out to be a no-op.

    Raises
    ------
    KeyError, ValueError
        For an unknown imputer or an invalid hyperparameter. Imputation failures
        propagate rather than falling back to a different strategy.
    """
    started = time.time()
    spec = get_imputer(method)
    params = dict(params or {})
    spec.validate_params(params)

    resolved_max_na = (
        max_na_frac
        if max_na_frac is not None
        else params.pop("max_na_frac", spec.defaults().get("max_na_frac", 0.0))
    )
    params.pop("max_na_frac", None)

    n_genes_in = exp_df.shape[1]
    na_frac = exp_df.isna().mean(axis=0)
    filtered = exp_df.loc[:, na_frac <= resolved_max_na]
    n_dropped = n_genes_in - filtered.shape[1]

    n_missing_before = int(filtered.isna().sum().sum())
    if n_missing_before == 0:
        # The donor returned here silently, so a "knn" run could be byte-identical to a
        # "strict" one with nothing recording that fact.
        result = filtered
        was_noop = method != "strict"
    else:
        result = spec.fn(filtered, **params)
        was_noop = False

    residual_na = int(result.isna().sum().sum())
    return result, ImputationReport(
        method=method,
        n_genes_in=n_genes_in,
        n_genes_dropped_na_frac=n_dropped,
        n_genes_out=result.shape[1],
        n_cells_imputed=n_missing_before - residual_na,
        residual_na=residual_na,
        was_noop=was_noop,
        elapsed_s=time.time() - started,
    )
