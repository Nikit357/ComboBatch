"""Post-harmonization removal of the most PCA-deviant batches.

Applied *after* correction, as a second pass: batches whose centroid still sits far from
the global centroid in PCA space are dropped entirely. Off by default, because it
discards data and only sometimes helps.

The donor wrapped this in a bare ``try/except`` that logged and continued, so a failure
still wrote a ``post1`` output — byte-identical to ``post0`` — which was then scored as
though post-removal had happened. Here a failure raises :class:`PostRemovalError` and the
caller records an explicit ``post_removal_failed`` status instead.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from combobatch.logging_utils import get_logger


class PostRemovalError(RuntimeError):
    """Raised when post-removal cannot be carried out."""


@dataclass(frozen=True)
class PostRemovalResult:
    """What post-removal discarded."""

    removed_batches: tuple[str, ...]
    n_samples_removed: int
    n_samples_remaining: int

    def as_dict(self) -> dict[str, object]:
        """Return the result as a flat dict for the JSON sidecar."""
        return {
            "post_removal_batches": list(self.removed_batches),
            "post_removal_n_samples_removed": self.n_samples_removed,
            "post_removal_n_samples_remaining": self.n_samples_remaining,
        }


def identify_outlier_batches(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    n_batches: int = 1,
    min_batch_size: int = 20,
    n_pcs: int = 2,
) -> list[str]:
    """
    Rank batches by PCA centroid distance and return the most deviant ones.

    Parameters
    ----------
    exp_df
        Expression, samples x genes; normally the harmonized matrix.
    ann_df
        Annotation aligned to ``exp_df``.
    batch_col
        Column holding batch identity.
    n_batches
        How many batches to return. The donor hardcoded 1.
    min_batch_size
        Batches smaller than this are never candidates: a small batch's centroid is
        noisy, and removing it would discard data on weak evidence.
    n_pcs
        Principal components used for the distance.

    Returns
    -------
    Up to ``n_batches`` batch labels, most deviant first. May be shorter, or empty, when
    too few batches meet ``min_batch_size``.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    if n_batches < 1:
        raise PostRemovalError(f"n_batches must be at least 1, got {n_batches}")

    common = exp_df.index.intersection(ann_df.index)
    if len(common) == 0:
        raise PostRemovalError("expression and annotation share no samples")

    values = StandardScaler().fit_transform(exp_df.loc[common].fillna(0))
    n_components = min(n_pcs, min(values.shape) - 1)
    if n_components < 1:
        raise PostRemovalError(
            f"not enough samples or genes for PCA: matrix is {values.shape}"
        )

    coords = pd.DataFrame(
        PCA(n_components=n_components).fit_transform(values), index=common
    )
    batches = ann_df.loc[common, batch_col].fillna("Unknown").astype(str)
    global_centroid = coords.mean(axis=0).values

    distances: dict[str, float] = {}
    for batch in batches.unique():
        idx = batches[batches == batch].index
        if len(idx) < min_batch_size:
            continue
        distances[batch] = float(
            np.linalg.norm(coords.loc[idx].mean(axis=0).values - global_centroid)
        )

    ranked = sorted(distances, key=lambda name: distances[name], reverse=True)
    return ranked[:n_batches]


def apply_post_removal(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    n_batches: int = 1,
    min_batch_size: int = 20,
    n_pcs: int = 2,
) -> tuple[pd.DataFrame, pd.DataFrame, PostRemovalResult]:
    """
    Identify and drop the most PCA-deviant batches.

    Returns
    -------
    ``(expression, annotation, result)`` with the offending samples removed.

    Raises
    ------
    PostRemovalError
        If outlier identification fails, or if removal would leave nothing behind.
        Both are failures rather than silent pass-throughs: an output that says
        ``post1`` must actually have had post-removal applied.
    """
    try:
        outliers = identify_outlier_batches(
            exp_df,
            ann_df,
            batch_col=batch_col,
            n_batches=n_batches,
            min_batch_size=min_batch_size,
            n_pcs=n_pcs,
        )
    except PostRemovalError:
        raise
    except Exception as exc:
        raise PostRemovalError(f"outlier-batch identification failed: {exc}") from exc

    logger = get_logger()
    if not outliers:
        # Nothing qualified — every batch is below min_batch_size. Reported as an empty
        # removal rather than an error: the analysis is still valid, just unchanged.
        logger.info(
            "Post-removal: no batch met min_batch_size=%d; nothing removed",
            min_batch_size,
        )
        return (
            exp_df,
            ann_df,
            PostRemovalResult((), 0, len(exp_df)),
        )

    batches = ann_df[batch_col].fillna("Unknown").astype(str)
    keep = ann_df.index[~batches.isin(outliers)]
    if len(keep) == 0:
        raise PostRemovalError(
            f"post-removal of {outliers} would discard every sample; "
            "lower --post-removal-n or raise --post-removal-min-batch"
        )

    kept_index = exp_df.index.intersection(keep)
    logger.info(
        "Post-removal dropped %d batch(es) %s, %d of %d samples",
        len(outliers),
        outliers,
        len(exp_df) - len(kept_index),
        len(exp_df),
    )
    return (
        exp_df.loc[kept_index],
        ann_df.loc[kept_index],
        PostRemovalResult(
            removed_batches=tuple(outliers),
            n_samples_removed=len(exp_df) - len(kept_index),
            n_samples_remaining=len(kept_index),
        ),
    )
