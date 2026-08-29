"""``compute_metrics()`` — run the requested groups over one matrix.

Three properties are load-bearing and all three are inherited from the donor:

* **Per-group isolation.** A group that raises records ``error_<LETTER>`` with its
  traceback and the rest still run. Fourteen groups over hundreds of jobs means something
  will fail, and losing thirteen results to the fourteenth is not acceptable.
* **Partial flushing.** ``on_group_done`` fires after every group so a later crash or a
  pod eviction cannot discard groups that already finished.
* **Incremental sentinels.** A group whose sentinel key is already present is not
  recomputed, so a new group can be added to a finished run without redoing the rest.

One donor bug is fixed here: ``--force-groups`` claimed in a comment to beat
``--skip-slow`` and did not, because the slow gate was applied after the force set.
"""

from __future__ import annotations

import time
import traceback
from typing import Any, Callable, Mapping, Sequence

import numpy as np
import pandas as pd

from combobatch.logging_utils import get_logger
from combobatch.metrics import (
    METRIC_GROUP_REGISTRY,
    MetricColumns,
    ordered_groups,
    required_embeddings,
)
from combobatch.metrics.panels import Panel


def resolve_active_groups(
    requested: frozenset[str] | set[str],
    *,
    skip_slow: bool = False,
    skip_wm: bool = False,
    force: frozenset[str] | set[str] = frozenset(),
    prior: Mapping[str, Any] | None = None,
    columns: MetricColumns | None = None,
) -> set[str]:
    """
    Decide which groups actually run.

    Order matters and is the fix for the donor's bug: exclusions are applied first, then
    the force set is unioned back in, so ``--force-groups F`` genuinely beats
    ``--skip-slow`` the way its comment always claimed.

    Parameters
    ----------
    requested
        Letters the user asked for.
    skip_slow
        Drop the groups the registry marks ``very_slow``.
    skip_wm
        Drop the WaterMelon group.
    force
        Letters that run regardless of every exclusion, including an existing sentinel.
    prior
        A previous metrics record; a group whose sentinel is already in it is skipped.
    columns
        Needed to resolve sentinel keys against the configured column names.

    Returns
    -------
    The letters to compute.
    """
    active = {letter.strip().upper() for letter in requested}
    active &= set(METRIC_GROUP_REGISTRY)

    log = get_logger()

    if skip_slow:
        excluded = {
            letter
            for letter in active
            if METRIC_GROUP_REGISTRY[letter].speed == "very_slow"
        }
        if excluded:
            # Named, not silent: a group the user explicitly asked for must never
            # vanish without a word about why or how to get it back.
            log.info(
                "skip_slow excludes group(s) %s; add --force-groups %s to run them",
                "".join(sorted(excluded)),
                "".join(sorted(excluded)),
            )
        active -= excluded
    if skip_wm and "I" in active:
        log.info("skip_wm excludes group I; add --force-groups I to run it")
        active.discard("I")

    if prior and columns is not None:
        already_done = {
            letter
            for letter in active
            if METRIC_GROUP_REGISTRY[letter].sentinel_key(columns) in prior
        }
        active -= already_done

    forced = {letter.strip().upper() for letter in force} & set(METRIC_GROUP_REGISTRY)
    return active | forced


def compute_metrics(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    column_spec: Any,
    *,
    spec: Any = None,
    groups: frozenset[str] | set[str] | None = None,
    ref_df: pd.DataFrame | None = None,
    panel: Panel | None = None,
    prior: Mapping[str, Any] | None = None,
    n_perm: int | None = None,
    predict_classes: Sequence[str] | None = None,
    collect_detail: bool = False,
    on_group_done: Callable[[dict], None] | None = None,
) -> dict:
    """
    Compute the requested metric groups for one expression/annotation pair.

    Parameters
    ----------
    exp_df, ann_df
        Expression (samples x genes) and its aligned annotation.
    column_spec
        A :class:`combobatch.config.ColumnSpec`; every column name comes from it.
    spec
        A :class:`combobatch.config.MetricsSpec` supplying groups, skip flags and
        permutation counts. Explicit keyword arguments override it.
    groups
        Letters to compute, overriding ``spec.groups``.
    ref_df
        The pre-harmonization matrix, required by group L. **This run's own in-memory
        matrix**, not a downloaded baseline: the donor fetched ``01_raw`` from storage,
        which could silently be a different baseline than the one actually harmonized.
    panel
        Marker panel for L and M. ``None`` uses every shared gene.
    prior
        A previous metrics record, for incremental recomputation.
    n_perm
        Permutations for group N's null.
    predict_classes
        Explicit class labels for group N.
    collect_detail
        Ask group L for its full gene-by-cohort matrix.
    on_group_done
        Called with a copy of the accumulated result after each group, success or
        failure. Exceptions from it are logged and swallowed — a failing persistence
        callback must never abort the metrics themselves.

    Returns
    -------
    One flat dict of every metric computed, plus ``error_<LETTER>`` for any group that
    failed and ``metrics_groups_run`` naming what actually ran.

    Raises
    ------
    ValueError
        If the two frames disagree on row count.
    """
    from combobatch.config import MetricsSpec

    if len(exp_df) != len(ann_df):
        raise ValueError(
            f"expression has {len(exp_df)} rows but annotation has {len(ann_df)}; "
            f"they must be aligned before metrics are computed"
        )

    spec = spec or MetricsSpec()
    log = get_logger()
    columns = MetricColumns.from_spec(column_spec)

    active = resolve_active_groups(
        groups if groups is not None else spec.groups,
        skip_slow=spec.skip_slow,
        skip_wm=spec.skip_wm,
        force=spec.force_groups,
        prior=prior,
        columns=columns,
    )
    n_perm = spec.n_perm if n_perm is None else n_perm
    predict_classes = predict_classes or spec.predict_classes

    result: dict = {}

    def flush() -> None:
        if on_group_done is None:
            return
        try:
            on_group_done(dict(result))
        except Exception:
            log.error("on_group_done callback failed:\n%s", traceback.format_exc())

    def run(letter: str, *args, **kwargs) -> None:
        if letter not in active:
            return
        spec_for = METRIC_GROUP_REGISTRY[letter]
        missing = spec_for.missing_backends()
        if missing:
            # An honest skip, named: the same convention the method registry uses, so a
            # missing backend never looks like a computation that failed.
            reason = f"skipped: needs {', '.join(missing)}"
            log.warning("metric group %s %s", letter, reason)
            result[f"error_{letter}"] = reason
            flush()
            return
        started = time.time()
        try:
            result.update(spec_for.fn(*args, **kwargs))
            log.info(
                "metric group %s (%s) done in %.1fs",
                letter,
                spec_for.name,
                time.time() - started,
            )
        except Exception:
            trace = traceback.format_exc()
            log.error("metric group %s FAILED:\n%s", letter, trace)
            result[f"error_{letter}"] = trace
        flush()

    # --- shared embeddings, computed once ------------------------------------------
    pca_coords: np.ndarray | None = None
    eigenvalues: np.ndarray | None = None
    umap_coords: np.ndarray | None = None
    tsne_coords: np.ndarray | None = None

    needed = required_embeddings(active)
    if "pca" in needed:
        from combobatch.metrics.embeddings import compute_pca

        try:
            pca_coords, eigenvalues = compute_pca(exp_df)
        except Exception:
            trace = traceback.format_exc()
            log.error(
                "PCA failed; every embedding-dependent group is skipped:\n%s", trace
            )
            result["error_embeddings"] = trace
            dependent = {
                letter
                for letter in active
                if "pca" in METRIC_GROUP_REGISTRY[letter].needs_embedding
            }
            for letter in dependent:
                result[f"error_{letter}"] = "skipped: PCA failed"
            active -= dependent
            flush()

    if pca_coords is not None and "umap" in required_embeddings(active):
        from combobatch.metrics.embeddings import compute_umap, compute_tsne

        try:
            umap_coords = compute_umap(pca_coords)
            tsne_coords = compute_tsne(pca_coords)
        except Exception:
            trace = traceback.format_exc()
            log.error("UMAP/t-SNE failed; group C is skipped:\n%s", trace)
            result["error_embeddings_2d"] = trace
            result["error_C"] = "skipped: UMAP or t-SNE failed"
            active.discard("C")
            flush()

    # Group L is the only group comparing against a second matrix. Its absence must not
    # affect anything else.
    if "L" in active and ref_df is None:
        result["error_L"] = "skipped: no pre-harmonization reference supplied"
        active.discard("L")
        flush()

    ran = sorted(active)

    run("E", exp_df, ann_df, columns)
    run("K", exp_df)
    run(
        "A",
        pca_coords,
        eigenvalues,
        ann_df,
        columns,
        n_dsc_permutations=spec.dsc_permutations,
    )
    run("J", eigenvalues)
    run("B", pca_coords, ann_df, columns, asw_max_samples=spec.asw_max_samples)
    run("C", umap_coords, tsne_coords, ann_df, columns)
    run("D", exp_df, ann_df, columns, n_genes_max=spec.ks_genes)
    run("G", pca_coords, ann_df, columns)
    run("H", exp_df, ann_df, columns)
    run(
        "I",
        ann_df,
        pca_coords,
        columns,
        n_permutations=spec.wm_permutations,
        max_samples=spec.wm_max_samples,
    )
    run(
        "L",
        exp_df,
        ann_df,
        columns,
        ref_df=ref_df,
        panel=panel,
        min_cohort_n=spec.min_cohort_n,
        collect_detail=collect_detail or spec.collect_detail,
    )
    run("M", exp_df, ann_df, columns, panel=panel, max_samples=spec.xb_max_samples)
    run(
        "N",
        exp_df,
        ann_df,
        columns,
        n_pcs=spec.n_pcs,
        n_perm=n_perm,
        min_test_n=spec.min_test_n,
        predict_classes=predict_classes,
    )
    run("F", exp_df, ann_df, columns)

    result["metrics_groups_run"] = "".join(ran)
    return result
