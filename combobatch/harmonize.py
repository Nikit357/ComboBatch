"""One (imputation x method) combination, start to finish.

This is the worker's whole job, and also the unit the dispatcher parallelizes:

    load + align -> subset -> impute [cached] -> resolve reference batch
      -> resolve + validate hyperparameters -> normalize
      -> for each post-removal variant: post-removal -> metrics -> write matrix + sidecar

Three donor behaviours are deliberately not reproduced. Imputation failure no longer
falls back to ``strict`` while keeping the original label; post-removal failure no longer
writes a ``post1`` output identical to ``post0``; and every sidecar is written with the
NaN-safe encoder, so it is valid strict JSON.
"""

from __future__ import annotations

import gc
import time
import traceback
from dataclasses import dataclass
from typing import Any

import pandas as pd

from combobatch import dataio, storage
from combobatch.config import ImputationSelection, MethodSelection, RunConfig
from combobatch.imputation import IMPUTER_REGISTRY, impute
from combobatch.logging_utils import get_logger
from combobatch.methods import METHOD_REGISTRY
from combobatch.methods.base import resolve_reference_batch
from combobatch.params import (
    ParamParseError,
    build_output_key,
    build_prepared_key,
    encode_param_tag,
)
from combobatch.postremoval import PostRemovalError, apply_post_removal

# Every matrix ComboBatch writes uses this format. Gzipped TSV is what the donor used and
# what `dataio` infers a tab separator from.
MATRIX_SUFFIX = ".tsv.gz"

MANIFEST_KEY = "run_manifest.json"


@dataclass(frozen=True)
class PreparedData:
    """An aligned, subset, imputed matrix plus its annotation and provenance."""

    exp_df: pd.DataFrame
    ann_df: pd.DataFrame
    report: dict[str, Any]
    from_cache: bool


def combined_param_tag(
    imp_sel: ImputationSelection, method_sel: MethodSelection
) -> str | None:
    """
    Build the single filename tag covering both the imputer's and the method's parameters.

    The key grammar has one tag slot but two parameter sources, so the non-default
    parameters of each are merged. A name present in both would silently overwrite the
    other, so it raises instead — rare, and always fixable with an explicit
    ``--param-tag``.

    Parameters
    ----------
    imp_sel, method_sel
        The user's two selections.

    Returns
    -------
    The tag, or ``None`` when both run at pure defaults — which is what keeps a default
    run's filename free of any tag segment.

    Raises
    ------
    ParamParseError
        If the same parameter name is tuned on both the imputer and the method.
    """
    if method_sel.param_tag is not None:
        return method_sel.tag()

    imp_params = _non_default_imputer_params(imp_sel)
    method_spec = METHOD_REGISTRY[method_sel.name]
    method_params = method_spec.non_default(method_sel.params)

    clashing = sorted(set(imp_params) & set(method_params))
    if clashing:
        raise ParamParseError(
            f"{imp_sel.name} and {method_sel.name} both tune {', '.join(clashing)}, so "
            f"one would overwrite the other in the output name. Set --param-tag "
            f"explicitly for this combination."
        )
    return encode_param_tag({**imp_params, **method_params})


def _non_default_imputer_params(imp_sel: ImputationSelection) -> dict[str, Any]:
    """Return the imputer parameters that differ from the registry defaults."""
    spec = IMPUTER_REGISTRY[imp_sel.name]
    non_default = dict(spec.non_default(imp_sel.params))
    # max_na_frac is promoted to a first-class CLI concept and so lives outside `params`,
    # but for naming purposes it is a parameter like any other.
    resolved = _resolved_max_na_frac(imp_sel)
    if resolved != spec.defaults().get("max_na_frac"):
        non_default["max_na_frac"] = resolved
    return non_default


def _resolved_max_na_frac(imp_sel: ImputationSelection) -> float:
    """Return the NA ceiling actually in force: the user's, else the imputer's default."""
    if imp_sel.max_na_frac is not None:
        return imp_sel.max_na_frac
    return float(IMPUTER_REGISTRY[imp_sel.name].defaults().get("max_na_frac", 0.0))


def load_inputs(cfg: RunConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read, align and optionally subset the expression and annotation tables.

    Alignment is by index intersection, never positional; the subset query is applied to
    the annotation and the expression follows it.

    Returns
    -------
    ``(expression, annotation)``, both indexed by the same samples in the same order.
    """
    log = get_logger()
    exp_df = dataio.read_table(cfg.exp_uri, endpoint_url=cfg.endpoint_url)
    ann_df = dataio.read_table(cfg.ann_uri, endpoint_url=cfg.endpoint_url)
    exp_df, ann_df, report = dataio.align(exp_df, ann_df)
    log.info("Loaded %s", report.summary())

    cfg.validate(ann_df, known_methods=METHOD_REGISTRY, known_imputers=IMPUTER_REGISTRY)

    if cfg.subset_query:
        from combobatch.subset import apply_subset_query

        ann_df = apply_subset_query(ann_df, cfg.subset_query)
        exp_df = exp_df.loc[ann_df.index]
        log.info("Subset query kept %d sample(s)", len(ann_df))

    return exp_df, ann_df


def prepare(
    cfg: RunConfig,
    imp_sel: ImputationSelection,
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    backend: storage.StorageBackend,
    use_cache: bool = True,
) -> PreparedData:
    """
    Impute one matrix, reusing a cached result when one exists.

    Imputation is far more expensive than most harmonizations — ``missforest`` in
    particular — and every method sharing an ``(imputer, tag)`` pair wants the same
    matrix, so it is written once under ``prepared/`` and read back thereafter.

    The cache is advisory, not a lock: workers that start together all miss it and all
    compute, writing identical content. That costs at most one redundant imputation per
    concurrent worker and is bounded by ``--n-workers``; correctness does not depend on
    who wins.

    Parameters
    ----------
    cfg
        The run configuration, for the output root and endpoint.
    imp_sel
        The imputation strategy and its parameters.
    exp_df, ann_df
        Aligned input matrices.
    backend
        Storage rooted at the output location.
    use_cache
        Set ``False`` to force recomputation.

    Returns
    -------
    The prepared data and the imputation report that goes into every sidecar.
    """
    log = get_logger()
    tag = encode_param_tag(_non_default_imputer_params(imp_sel))
    exp_key = f"prepared/{build_prepared_key(imp_sel.name, tag, 'exp')}{MATRIX_SUFFIX}"
    ann_key = f"prepared/{build_prepared_key(imp_sel.name, tag, 'ann')}{MATRIX_SUFFIX}"

    if use_cache and backend.exists(exp_key) and backend.exists(ann_key):
        log.info("Reusing prepared matrix %s", backend.uri(exp_key))
        cached_exp = dataio.read_table_from(backend, exp_key)
        cached_ann = dataio.read_table_from(backend, ann_key)
        return PreparedData(
            exp_df=cached_exp,
            ann_df=cached_ann,
            report={"imputation": imp_sel.name, "imputation_cached": True},
            from_cache=True,
        )

    imputed, report = impute(
        exp_df,
        method=imp_sel.name,
        max_na_frac=_resolved_max_na_frac(imp_sel),
        params=imp_sel.params,
    )
    log.info(
        "Imputation %s: %d -> %d genes, %d cell(s) filled%s",
        imp_sel.name,
        report.n_genes_in,
        report.n_genes_out,
        report.n_cells_imputed,
        " (no-op: nothing was missing)" if report.was_noop else "",
    )

    dataio.write_table_to(imputed, backend, exp_key)
    dataio.write_table_to(ann_df, backend, ann_key)

    return PreparedData(
        exp_df=imputed,
        ann_df=ann_df,
        report={
            "imputation": imp_sel.name,
            "imputation_cached": False,
            **report.as_dict(),
        },
        from_cache=False,
    )


def run_one_combination(
    cfg: RunConfig,
    imp_sel: ImputationSelection,
    method_sel: MethodSelection,
    *,
    exp_df: pd.DataFrame | None = None,
    ann_df: pd.DataFrame | None = None,
    skip_if_exists: bool = False,
) -> list[dict[str, Any]]:
    """
    Run one imputation x method combination and write every output it produces.

    Parameters
    ----------
    cfg
        The validated run configuration.
    imp_sel, method_sel
        The two axes of this job.
    exp_df, ann_df
        Pre-loaded inputs. ``None`` loads them from ``cfg``.
    skip_if_exists
        Report an existing output as ``cached`` instead of recomputing it.

    Returns
    -------
    One sidecar row per output variant. A row always exists, even for a skipped or
    failed job, so the dispatcher can tell "did not run" from "was never attempted".
    """
    log = get_logger()
    backend = storage.backend_for_root(cfg.out_uri, endpoint_url=cfg.endpoint_url)
    method_spec = METHOD_REGISTRY[method_sel.name]
    param_tag = combined_param_tag(imp_sel, method_sel)
    variants = _post_removal_variants(cfg)

    if exp_df is None or ann_df is None:
        exp_df, ann_df = load_inputs(cfg)

    if skip_if_exists and all(
        backend.exists(_matrix_key(imp_sel.name, method_sel.name, param_tag, variant))
        for variant in variants
    ):
        log.info("%s x %s already present — skipping", imp_sel.name, method_sel.name)
        return [
            _row(cfg, imp_sel, method_sel, param_tag, variant, status="cached")
            for variant in variants
        ]

    prepared = prepare(cfg, imp_sel, exp_df, ann_df, backend=backend)

    params = method_spec.defaults() | dict(method_sel.params)
    method_spec.validate_params(method_sel.params)
    if method_spec.uses_reference_batch and params.get("target_group") is None:
        params["target_group"] = resolve_reference_batch(
            prepared.ann_df, cfg.columns.batch, cfg.reference_batch
        )

    started = time.time()
    try:
        harmonized = method_spec.fn(
            prepared.exp_df,
            prepared.ann_df,
            batch_col=cfg.columns.batch,
            bio_col=cfg.columns.bio,
            **params,
        )
    except NotImplementedError as exc:
        # An honest skip: the method cannot run on this data or in this environment, and
        # says so. Distinct from a failure, and never written as an output.
        log.warning("%s skipped: %s", method_sel.name, exc)
        return [
            _row(
                cfg,
                imp_sel,
                method_sel,
                param_tag,
                variant,
                status="skipped",
                error=str(exc),
                extra=prepared.report,
            )
            for variant in variants
        ]
    except Exception as exc:
        log.error("%s failed:\n%s", method_sel.name, traceback.format_exc())
        return [
            _row(
                cfg,
                imp_sel,
                method_sel,
                param_tag,
                variant,
                status="failed",
                error=f"{type(exc).__name__}: {exc}",
                extra=prepared.report,
            )
            for variant in variants
        ]
    elapsed = time.time() - started
    log.info("%s harmonized in %.0fs", method_sel.name, elapsed)

    rows: list[dict[str, Any]] = []
    for variant in variants:
        rows.append(
            _write_variant(
                cfg,
                imp_sel,
                method_sel,
                param_tag,
                variant,
                harmonized=harmonized,
                prepared=prepared,
                params=params,
                backend=backend,
                elapsed=elapsed,
                skip_if_exists=skip_if_exists,
            )
        )

    del harmonized, prepared
    gc.collect()
    _release_r()
    return rows


def _write_variant(
    cfg: RunConfig,
    imp_sel: ImputationSelection,
    method_sel: MethodSelection,
    param_tag: str | None,
    variant: bool | None,
    *,
    harmonized: pd.DataFrame,
    prepared: PreparedData,
    params: dict[str, Any],
    backend: storage.StorageBackend,
    elapsed: float,
    skip_if_exists: bool,
) -> dict[str, Any]:
    """Apply one post-removal variant, then write its matrix, genes and sidecar."""
    log = get_logger()
    stem = build_output_key(imp_sel.name, method_sel.name, param_tag, variant)
    matrix_key = f"exp/{stem}{MATRIX_SUFFIX}"

    if skip_if_exists and backend.exists(matrix_key):
        return _row(cfg, imp_sel, method_sel, param_tag, variant, status="cached")

    exp_out, ann_out = harmonized, prepared.ann_df
    removal: dict[str, Any] = {}

    if variant:
        try:
            exp_out, ann_out, result = apply_post_removal(
                harmonized,
                prepared.ann_df,
                batch_col=cfg.columns.batch,
                n_batches=cfg.post_removal.n_batches,
                min_batch_size=cfg.post_removal.min_batch_size,
                n_pcs=cfg.post_removal.n_pcs,
            )
            removal = {
                "post_removed_batches": list(result.removed_batches),
                "post_n_samples_removed": result.n_samples_removed,
            }
        except PostRemovalError as exc:
            # The donor logged this and wrote post1 anyway — a file identical to post0,
            # scored as if a batch had been removed.
            log.error("post-removal failed for %s: %s", stem, exc)
            return _row(
                cfg,
                imp_sel,
                method_sel,
                param_tag,
                variant,
                status="post_removal_failed",
                error=str(exc),
                extra=prepared.report,
            )

    metrics: dict[str, Any] = {}
    if cfg.metrics.enabled:
        metrics = _inline_metrics(cfg, exp_out, ann_out, reference=prepared.exp_df)

    try:
        dataio.write_table_to(exp_out, backend, matrix_key)
        status = "ok"
    except Exception as exc:
        log.error("writing %s failed: %s", backend.uri(matrix_key), exc)
        status = "upload_failed"

    row = _row(
        cfg,
        imp_sel,
        method_sel,
        param_tag,
        variant,
        status=status,
        extra={
            **prepared.report,
            **removal,
            **metrics,
            "params": params,
            "compute_time_s": round(elapsed, 3),
            "n_samples": int(exp_out.shape[0]),
            "n_genes": int(exp_out.shape[1]),
            "output": backend.uri(matrix_key),
        },
    )

    if status == "ok":
        dataio.write_json(row, backend, f"metrics/{stem}_metrics.json")
        dataio.write_json(
            {"key": stem, "genes": list(exp_out.columns)},
            backend,
            f"genes/{stem}_genes.json",
        )

    del exp_out, ann_out
    gc.collect()
    return row


def _inline_metrics(
    cfg: RunConfig,
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    reference: pd.DataFrame,
) -> dict[str, Any]:
    """
    Compute the requested metric groups alongside harmonization.

    Group L's reference is this run's own pre-harmonization matrix, held in memory —
    the donor downloaded a separate ``01_raw`` output, which could silently be a
    different baseline than the one actually harmonized.
    """
    from combobatch.metrics.panels import load_panel
    from combobatch.metrics.runner import compute_metrics

    panel = (
        load_panel(cfg.metrics.marker_panel, endpoint_url=cfg.endpoint_url)
        if cfg.metrics.marker_panel
        else None
    )
    return compute_metrics(
        exp_df,
        ann_df,
        cfg.columns,
        spec=cfg.metrics,
        ref_df=reference.reindex(index=exp_df.index),
        panel=panel,
        predict_classes=cfg.metrics.predict_classes,
    )


def _post_removal_variants(cfg: RunConfig) -> list[bool | None]:
    """
    Return the post-removal variants to produce for one combination.

    ``None`` means "no post-removal segment in the name at all", which is what a run
    with post-removal disabled writes — one output, not two identical ones.
    """
    return [False, True] if cfg.post_removal.enabled else [None]


def _matrix_key(
    imputation: str, method: str, param_tag: str | None, variant: bool | None
) -> str:
    """Return the storage key of one output matrix."""
    return (
        f"exp/{build_output_key(imputation, method, param_tag, variant)}{MATRIX_SUFFIX}"
    )


def _row(
    cfg: RunConfig,
    imp_sel: ImputationSelection,
    method_sel: MethodSelection,
    param_tag: str | None,
    variant: bool | None,
    *,
    status: str,
    error: str | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Build one sidecar row. Every job produces one per variant, whatever happened."""
    row: dict[str, Any] = {
        "imp": imp_sel.name,
        "method": method_sel.name,
        "param_tag": param_tag,
        "post_rm": variant,
        "key": build_output_key(imp_sel.name, method_sel.name, param_tag, variant),
        "harshness": METHOD_REGISTRY[method_sel.name].harshness,
        "status": status,
        "run_id": cfg.run_id,
    }
    if error is not None:
        row["error"] = error
    if extra:
        row.update(extra)
    return row


def _release_r() -> None:
    """Run R's garbage collector, when R is in play at all."""
    from combobatch.methods import rinterop

    if rinterop.r_available():
        rinterop.r_gc()


# --------------------------------------------------------------------------------------
# Run manifest
# --------------------------------------------------------------------------------------


def build_manifest(cfg: RunConfig) -> dict[str, Any]:
    """
    Describe every output this configuration will produce, with fully resolved parameters.

    Filenames carry only the *non-default* parameters, and long tags are truncated with a
    hash. The manifest is what keeps a result self-describing anyway: it maps each output
    key back to the complete parameter dict, defaults included.

    Returns
    -------
    The manifest document, ready for :func:`combobatch.dataio.write_json`.
    """
    from combobatch import __version__
    from datetime import datetime, timezone

    entries: dict[str, Any] = {}
    for imp_sel, method_sel in cfg.combinations():
        method_spec = METHOD_REGISTRY[method_sel.name]
        imputer_spec = IMPUTER_REGISTRY[imp_sel.name]
        param_tag = combined_param_tag(imp_sel, method_sel)
        for variant in _post_removal_variants(cfg):
            key = build_output_key(imp_sel.name, method_sel.name, param_tag, variant)
            entries[key] = {
                "imputation": imp_sel.name,
                "imputation_params": imputer_spec.defaults() | dict(imp_sel.params),
                "max_na_frac": _resolved_max_na_frac(imp_sel),
                "method": method_sel.name,
                "method_params": method_spec.defaults() | dict(method_sel.params),
                "param_tag": param_tag,
                "post_rm": variant,
                "harshness": method_spec.harshness,
                "output": f"exp/{key}{MATRIX_SUFFIX}",
            }

    return {
        "combobatch_version": __version__,
        "run_id": cfg.run_id,
        "created_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "input": {"expression": cfg.exp_uri, "annotation": cfg.ann_uri},
        "columns": {
            "batch": cfg.columns.batch,
            "bio": cfg.columns.bio,
            "cohort": cfg.columns.cohort,
        },
        "reference_batch": cfg.reference_batch,
        "subset_query": cfg.subset_query,
        "random_seed": cfg.random_seed,
        "post_removal": {
            "enabled": cfg.post_removal.enabled,
            "n_batches": cfg.post_removal.n_batches,
            "min_batch_size": cfg.post_removal.min_batch_size,
            "n_pcs": cfg.post_removal.n_pcs,
        },
        "outputs": entries,
    }


def write_manifest(
    cfg: RunConfig, backend: storage.StorageBackend | None = None
) -> str:
    """
    Write ``run_manifest.json``, merging into any manifest already at the output root.

    Merging matters because separate invocations add outputs to the same root; the
    dispatcher writes the whole matrix once up front, and standalone ``run`` calls add
    their single entry.

    Returns
    -------
    The URI written.
    """
    backend = backend or storage.backend_for_root(
        cfg.out_uri, endpoint_url=cfg.endpoint_url
    )
    manifest = build_manifest(cfg)

    if backend.exists(MANIFEST_KEY):
        try:
            existing = dataio.read_json(backend, MANIFEST_KEY)
            merged = dict(existing.get("outputs", {}))
            merged.update(manifest["outputs"])
            manifest["outputs"] = merged
        except (ValueError, KeyError, AttributeError) as exc:
            get_logger().warning(
                "existing %s is unreadable (%s); replacing it", MANIFEST_KEY, exc
            )

    dataio.write_json(manifest, backend, MANIFEST_KEY)
    return backend.uri(MANIFEST_KEY)
