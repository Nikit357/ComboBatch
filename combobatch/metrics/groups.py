"""The fourteen metric groups, A through N.

Ported from the donor's ``compute_batch_metrics.py``. The mathematics is unchanged; what
changed is that **no column name is written down here**. The donor hardcoded four batch
columns, four biology columns, a cohort column and two diagnosis literals across every
group; each one is now resolved from :class:`MetricColumns`, which comes from the user's
``ColumnSpec``. That is what makes these metrics computable on anyone's dataset.

Every group takes and returns the same shape: a flat ``dict`` of metric name to scalar
(with a handful of nested dicts from L, M and N that ``concat`` splits into long tables).
Metric keys are suffixed with the column they were computed on, so a run configured with
several batch columns produces one key per column and the sentinels stay consistent.
"""

from __future__ import annotations

import itertools
import math
import threading
import traceback
from dataclasses import dataclass
from typing import Any, Sequence

import numpy as np
import pandas as pd

from combobatch.logging_utils import get_logger
from combobatch.metrics.panels import Panel, resolve_panel

RANDOM_SEED = 42

# Donor defaults, promoted to parameters. Each was unreachable upstream.
DEFAULT_N_PERM = 20
DEFAULT_MIN_COHORT_N = 20
DEFAULT_MIN_TEST_N = 20
DEFAULT_N_PCS = 10
DEFAULT_XB_MAX_SAMPLES = 8000
DEFAULT_WM_PERMUTATIONS = 200
DEFAULT_WM_MAX_SAMPLES = 2000
DEFAULT_DSC_PERMUTATIONS = 999
DEFAULT_ASW_MAX_SAMPLES = 2000
DEFAULT_KS_GENES = 1000
DEFAULT_D1_TIMEOUT_S = 3600


@dataclass(frozen=True)
class MetricColumns:
    """The annotation columns every group reads, resolved from the user's config.

    ``cohort`` falls back to the batch column when the user names none, because several
    groups partition by cohort and the batch is the closest generic stand-in.
    """

    batch: str
    bio: str
    cohort: str
    batch_cols: tuple[str, ...]
    bio_cols: tuple[str, ...]
    predict_class_col: str

    @classmethod
    def from_spec(cls, spec: Any) -> "MetricColumns":
        """Build from a :class:`combobatch.config.ColumnSpec`."""
        return cls(
            batch=spec.batch,
            bio=spec.bio,
            cohort=spec.cohort or spec.batch,
            batch_cols=tuple(spec.metric_batch_cols or (spec.batch,)),
            bio_cols=tuple(spec.metric_bio_cols or (spec.bio,)),
            predict_class_col=spec.predict_class_col or spec.bio,
        )

    @property
    def all_cols(self) -> tuple[str, ...]:
        """Every column any group groups by, deduplicated in a stable order."""
        names = [*self.batch_cols, *self.bio_cols, self.cohort]
        return tuple(dict.fromkeys(names))


# --------------------------------------------------------------------------------------
# Shared helpers
# --------------------------------------------------------------------------------------


def _usable_mask(ann_df: pd.DataFrame, col: str) -> np.ndarray | None:
    """Rows where ``col`` is present, or ``None`` if it cannot discriminate at all."""
    if col not in ann_df.columns:
        return None
    mask = ann_df[col].notna().values
    if ann_df.loc[mask, col].astype(str).nunique() < 2:
        return None
    return mask


def _labels_for(ann_df: pd.DataFrame, col: str) -> np.ndarray:
    """String labels for the non-NaN rows of ``col``."""
    return ann_df.loc[ann_df[col].notna(), col].astype(str).values


def _nan_keys(result: dict, keys: Sequence[str]) -> None:
    """Set every named key to NaN, so an absent column still yields a complete record."""
    for key in keys:
        result[key] = np.nan


# --------------------------------------------------------------------------------------
# Group E — data quality
# --------------------------------------------------------------------------------------


def compute_group_e(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, columns: MetricColumns
) -> dict:
    """
    Describe the matrix: shape, batch balance, zero inflation, and per-cohort bimodality.

    Returns
    -------
    Dict of ``n_samples``, ``n_genes``, batch-size summaries, zero and below-one
    fractions overall and per batch, expression percentiles, and the bimodality counts.
    """
    from scipy.stats import kurtosis, skew
    from sklearn.mixture import GaussianMixture

    result: dict = {}
    n_samples, n_genes = exp_df.shape
    result["n_samples"] = int(n_samples)
    result["n_genes"] = int(n_genes)

    batch = columns.batch
    if batch in ann_df.columns:
        counts = ann_df[batch].value_counts()
        result["n_batches"] = int(ann_df[batch].nunique())
        result["n_samples_per_batch_min"] = float(counts.min())
        result["n_samples_per_batch_max"] = float(counts.max())
        result["n_samples_per_batch_median"] = float(counts.median())
        result["n_samples_per_batch_sd"] = float(
            counts.std(ddof=1) if len(counts) > 1 else 0.0
        )
    else:
        _nan_keys(
            result,
            [
                "n_batches",
                "n_samples_per_batch_min",
                "n_samples_per_batch_max",
                "n_samples_per_batch_median",
                "n_samples_per_batch_sd",
            ],
        )

    result["n_cohorts"] = (
        int(ann_df[columns.cohort].nunique())
        if columns.cohort in ann_df.columns
        else np.nan
    )
    result["n_bio_groups"] = (
        int(ann_df[columns.bio].nunique()) if columns.bio in ann_df.columns else np.nan
    )

    values = exp_df.values.astype(float)
    result["zero_fraction_global"] = float(np.mean(values == 0))
    result["fraction_genes_below_1"] = float(np.mean(values < 1))

    if batch in ann_df.columns:
        zero_fracs, below1_fracs = [], []
        for level in ann_df[batch].unique():
            block = values[(ann_df[batch] == level).values]
            zero_fracs.append(float(np.mean(block == 0)))
            below1_fracs.append(float(np.mean(block < 1)))
        result["zero_fraction_by_batch_max"] = float(max(zero_fracs))
        result["zero_fraction_by_batch_min"] = float(min(zero_fracs))
        result["below_1_fraction_by_batch_max"] = float(max(below1_fracs))
        result["below_1_fraction_by_batch_min"] = float(min(below1_fracs))
    else:
        _nan_keys(
            result,
            [
                "zero_fraction_by_batch_max",
                "zero_fraction_by_batch_min",
                "below_1_fraction_by_batch_max",
                "below_1_fraction_by_batch_min",
            ],
        )

    flat = values.flatten()
    result["exp_min"] = float(np.nanmin(flat))
    result["exp_max"] = float(np.nanmax(flat))
    result["exp_median"] = float(np.nanmedian(flat))
    result["exp_std"] = float(np.nanstd(flat))
    for percentile, name in [
        (1, "p01"),
        (5, "p05"),
        (25, "p25"),
        (75, "p75"),
        (95, "p95"),
        (99, "p99"),
    ]:
        result[f"exp_{name}"] = float(np.nanpercentile(flat, percentile))

    if batch in ann_df.columns:
        medians = [
            float(np.nanmedian(values[(ann_df[batch] == level).values]))
            for level in ann_df[batch].unique()
        ]
        finite = [value for value in medians if not math.isnan(value)]
        if len(finite) > 1 and abs(np.mean(finite)) > 1e-9:
            result["per_batch_median_cv"] = float(
                np.std(finite, ddof=1) / abs(np.mean(finite))
            )
        else:
            result["per_batch_median_cv"] = np.nan
    else:
        result["per_batch_median_cv"] = np.nan

    if columns.cohort in ann_df.columns:
        n_bimodal = 0
        n_zero_inflated = 0
        cohorts = ann_df[columns.cohort].unique()
        for cohort in cohorts:
            marginal = values[(ann_df[columns.cohort] == cohort).values].flatten()
            marginal = marginal[np.isfinite(marginal)]
            n_c = len(marginal)
            if n_c < 4:
                continue
            try:
                skewness = float(skew(marginal))
                excess = float(kurtosis(marginal, fisher=True))
                denominator = excess + 3.0 * (n_c - 1) ** 2 / ((n_c - 2) * (n_c - 3))
                coefficient = (
                    (skewness**2 + 1) / denominator
                    if abs(denominator) > 1e-9
                    else np.nan
                )
                if not np.isnan(coefficient) and coefficient > 0.555:
                    n_bimodal += 1
                    try:
                        mixture = GaussianMixture(n_components=2, random_state=42)
                        mixture.fit(marginal.reshape(-1, 1))
                        if float(min(mixture.means_.flatten())) < 1.0:
                            n_zero_inflated += 1
                    except Exception:
                        pass
            except Exception:
                pass
        total = len(cohorts)
        result["n_cohorts_bimodal"] = int(n_bimodal)
        result["fraction_cohorts_bimodal"] = (
            float(n_bimodal / total) if total else np.nan
        )
        result["n_cohorts_zero_inflated_bimodal"] = int(n_zero_inflated)
        result["fraction_cohorts_zero_inflated_bimodal"] = (
            float(n_zero_inflated / total) if total else np.nan
        )
    else:
        _nan_keys(
            result,
            [
                "n_cohorts_bimodal",
                "fraction_cohorts_bimodal",
                "n_cohorts_zero_inflated_bimodal",
                "fraction_cohorts_zero_inflated_bimodal",
            ],
        )
    return result


# --------------------------------------------------------------------------------------
# Group A — PCA variance decomposition
# --------------------------------------------------------------------------------------


def _r2_anova_pc(pc_values: np.ndarray, labels: np.ndarray) -> float:
    """One-way ANOVA R-squared of one principal component against one label set."""
    unique = np.unique(labels)
    if len(unique) < 2:
        return np.nan
    grand_mean = pc_values.mean()
    ss_total = float(((pc_values - grand_mean) ** 2).sum())
    if ss_total < 1e-12:
        return 0.0
    ss_between = sum(
        (labels == level).sum() * (pc_values[labels == level].mean() - grand_mean) ** 2
        for level in unique
    )
    return float(ss_between / ss_total)


def _dsc_score(coords_2d: np.ndarray, labels: np.ndarray) -> float:
    """Dispersion separability criterion: between-group over within-group scatter."""
    unique = np.unique(labels)
    if len(unique) < 2:
        return np.nan
    centre = coords_2d.mean(axis=0)
    between = 0.0
    within = 0.0
    for level in unique:
        mask = labels == level
        if not mask.any():
            continue
        group_centre = coords_2d[mask].mean(axis=0)
        offset = group_centre - centre
        between += float(mask.sum() * np.dot(offset, offset))
        within += float(np.sum((coords_2d[mask] - group_centre) ** 2))
    if within < 1e-12:
        return np.nan
    return between / within


def compute_group_a(
    pca_coords: np.ndarray,
    eigenvalues: np.ndarray,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    n_dsc_permutations: int = DEFAULT_DSC_PERMUTATIONS,
) -> dict:
    """
    Decompose PCA variance by each annotation column.

    Returns
    -------
    ``r2_{col}`` (mean over the first ten PCs), ``r2_pc{i}_{col}`` per batch column,
    ``pcr_{col}`` (variance-weighted, so higher is better), and ``dsc_{col}`` with a
    permutation p-value.
    """
    result: dict = {}
    n_pcs_mean = min(10, pca_coords.shape[1])
    n_pcs_pcr = min(100, pca_coords.shape[1])

    for col in columns.all_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            result[f"r2_{col}"] = np.nan
            continue
        labels = _labels_for(ann_df, col)
        result[f"r2_{col}"] = float(
            np.nanmean(
                [_r2_anova_pc(pca_coords[mask, i], labels) for i in range(n_pcs_mean)]
            )
        )

    for col in columns.batch_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            _nan_keys(result, [f"r2_pc{i}_{col}" for i in range(1, n_pcs_mean + 1)])
            continue
        labels = _labels_for(ann_df, col)
        for i in range(n_pcs_mean):
            result[f"r2_pc{i + 1}_{col}"] = _r2_anova_pc(pca_coords[mask, i], labels)

    total_variance = eigenvalues[:n_pcs_pcr].sum()
    for col in columns.all_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None or total_variance < 1e-12:
            result[f"pcr_{col}"] = np.nan
            continue
        labels = _labels_for(ann_df, col)
        weighted = 0.0
        for i in range(n_pcs_pcr):
            r2 = _r2_anova_pc(pca_coords[mask, i], labels)
            if not np.isnan(r2):
                weighted += eigenvalues[i] * r2
        result[f"pcr_{col}"] = float(1.0 - weighted / total_variance)

    coords_2d = pca_coords[:, :2]
    for col in columns.batch_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            _nan_keys(result, [f"dsc_{col}", f"dsc_pvalue_{col}"])
            continue
        labels = _labels_for(ann_df, col)
        subset = coords_2d[mask]
        observed = _dsc_score(subset, labels)
        if np.isnan(observed):
            _nan_keys(result, [f"dsc_{col}", f"dsc_pvalue_{col}"])
            continue
        rng = np.random.default_rng(RANDOM_SEED)
        n_extreme = sum(
            1
            for _ in range(n_dsc_permutations)
            for score in [_dsc_score(subset, rng.permutation(labels))]
            if not np.isnan(score) and score >= observed
        )
        result[f"dsc_{col}"] = float(observed)
        # +1 to both: the observed value is itself one draw from the null, so the
        # p-value can never be exactly zero.
        result[f"dsc_pvalue_{col}"] = float((n_extreme + 1) / (n_dsc_permutations + 1))
    return result


# --------------------------------------------------------------------------------------
# Group J — per-PC variance explained
# --------------------------------------------------------------------------------------


def compute_group_j(eigenvalues: np.ndarray) -> dict:
    """Report the share of variance carried by each of the first ten components."""
    result: dict = {}
    total = float(eigenvalues.sum())
    if total <= 0:
        _nan_keys(result, [f"pct_var_pc{i + 1}" for i in range(10)])
        result["pct_var_cum_top10"] = np.nan
        return result
    for i in range(min(10, len(eigenvalues))):
        result[f"pct_var_pc{i + 1}"] = float(eigenvalues[i] / total * 100)
    result["pct_var_cum_top10"] = float(eigenvalues[:10].sum() / total * 100)
    return result


# --------------------------------------------------------------------------------------
# Group B — neighbourhood integration
# --------------------------------------------------------------------------------------


def _kbet_acceptance(pca_coords: np.ndarray, labels: np.ndarray) -> float:
    """Fraction of test neighbourhoods whose label mix matches the global mix."""
    from scipy.stats import chi2
    from sklearn.neighbors import NearestNeighbors

    n = len(labels)
    if n < 10:
        return np.nan
    k = min(25, n // 4)
    if k < 2:
        return np.nan
    unique, counts = np.unique(labels, return_counts=True)
    if len(unique) < 2:
        return np.nan
    global_freq = {level: count / n for level, count in zip(unique, counts)}

    rng = np.random.default_rng(RANDOM_SEED)
    n_test = max(1, int(n * 0.1))
    test_idx = rng.choice(n, n_test, replace=False)

    neighbours = NearestNeighbors(n_neighbors=k + 1, metric="euclidean", n_jobs=-1)
    neighbours.fit(pca_coords)
    _, indices = neighbours.kneighbors(pca_coords[test_idx])

    degrees = len(unique) - 1
    n_rejected = 0
    for row in indices:
        neighbour_labels = labels[row[1:]]
        observed = np.array(
            [(neighbour_labels == level).sum() for level in unique], dtype=float
        )
        expected = np.array([global_freq[level] * k for level in unique], dtype=float)
        with np.errstate(divide="ignore", invalid="ignore"):
            statistic = np.nansum(
                np.where(expected > 0, (observed - expected) ** 2 / expected, 0.0)
            )
        if 1.0 - chi2.cdf(statistic, degrees) <= 0.05:
            n_rejected += 1
    return 1.0 - n_rejected / n_test


def _lisi_values(pca_coords: np.ndarray, labels: np.ndarray, k: int = 30) -> np.ndarray:
    """Per-sample inverse Simpson index of the label mix among k nearest neighbours."""
    from sklearn.neighbors import NearestNeighbors

    n = len(labels)
    k_eff = min(k, n - 1)
    neighbours = NearestNeighbors(n_neighbors=k_eff + 1, n_jobs=-1)
    neighbours.fit(pca_coords)
    _, indices = neighbours.kneighbors(pca_coords)
    values = np.empty(n, dtype=float)
    for i in range(n):
        _, counts = np.unique(labels[indices[i, 1:]], return_counts=True)
        proportions = counts / counts.sum()
        values[i] = 1.0 / float(np.dot(proportions, proportions))
    return values


def _centroid_dispersion(coords: np.ndarray, labels: np.ndarray) -> float:
    """How far group centroids spread, relative to the overall spread."""
    unique = np.unique(labels)
    if len(unique) < 2:
        return np.nan
    centroids = np.array([coords[labels == level].mean(axis=0) for level in unique])
    global_std = coords.std(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.where(global_std > 1e-9, centroids.std(axis=0) / global_std, np.nan)
    return float(np.nanmean(ratio))


def _local_entropy(
    coords: np.ndarray, labels: np.ndarray, k: int = 30
) -> tuple[float, float]:
    """Mean and max-normalized Shannon entropy of labels among k nearest neighbours."""
    from sklearn.neighbors import NearestNeighbors

    n = len(labels)
    unique = np.unique(labels)
    if len(unique) < 2:
        return np.nan, np.nan
    k_eff = min(k, n - 1)
    neighbours = NearestNeighbors(n_neighbors=k_eff + 1, n_jobs=-1)
    neighbours.fit(coords)
    _, indices = neighbours.kneighbors(coords)
    entropies = np.empty(n, dtype=float)
    for i in range(n):
        _, counts = np.unique(labels[indices[i, 1:]], return_counts=True)
        proportions = counts / counts.sum()
        entropies[i] = -float(np.sum(proportions * np.log(proportions + 1e-15)))
    mean_entropy = float(np.mean(entropies))
    max_entropy = math.log(len(unique))
    return mean_entropy, (mean_entropy / max_entropy if max_entropy > 1e-9 else np.nan)


def _silhouette(
    coords: np.ndarray, labels: np.ndarray, max_samples: int
) -> tuple[np.ndarray, np.ndarray]:
    """Subsample to bound the O(n^2) silhouette computation."""
    if len(labels) <= max_samples:
        return coords, labels
    rng = np.random.default_rng(RANDOM_SEED)
    idx = rng.choice(len(labels), max_samples, replace=False)
    return coords[idx], labels[idx]


def compute_group_b(
    pca_coords: np.ndarray,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    asw_max_samples: int = DEFAULT_ASW_MAX_SAMPLES,
) -> dict:
    """
    Measure how well batches mix, and how well biology stays separated, in PCA space.

    Returns
    -------
    kBET acceptance, iLISI and cLISI, batch and biology silhouettes with their
    normalized forms, and the cell-mixing statistic — each suffixed by its column.
    """
    from scipy.stats import anderson_ksamp
    from sklearn.metrics import silhouette_score
    from sklearn.neighbors import NearestNeighbors

    result: dict = {}

    for col in columns.batch_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            result[f"kbet_acceptance_rate_{col}"] = np.nan
            continue
        result[f"kbet_acceptance_rate_{col}"] = _kbet_acceptance(
            pca_coords[mask], _labels_for(ann_df, col)
        )

    for col in columns.batch_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            _nan_keys(result, [f"ilisi_mean_{col}", f"ilisi_norm_{col}"])
            continue
        labels = _labels_for(ann_df, col)
        mean_lisi = float(np.mean(_lisi_values(pca_coords[mask], labels)))
        n_levels = len(np.unique(labels))
        result[f"ilisi_mean_{col}"] = mean_lisi
        result[f"ilisi_norm_{col}"] = (
            float((mean_lisi - 1.0) / (n_levels - 1)) if n_levels > 1 else np.nan
        )

    for col in columns.bio_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            result[f"clisi_mean_{col}"] = np.nan
            continue
        result[f"clisi_mean_{col}"] = float(
            np.mean(_lisi_values(pca_coords[mask], _labels_for(ann_df, col)))
        )

    for col in columns.batch_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            _nan_keys(result, [f"asw_batch_{col}", f"asw_batch_norm_{col}"])
            continue
        subset, labels = _silhouette(
            pca_coords[mask], _labels_for(ann_df, col), asw_max_samples
        )
        try:
            asw = float(silhouette_score(subset, labels))
        except Exception:
            asw = np.nan
        result[f"asw_batch_{col}"] = asw
        # Batches should not separate, so zero is best and the normalization is 1-|s|.
        result[f"asw_batch_norm_{col}"] = (
            float(1.0 - abs(asw)) if not np.isnan(asw) else np.nan
        )

    for col in columns.bio_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            _nan_keys(result, [f"asw_bio_{col}", f"asw_bio_norm_{col}"])
            continue
        subset, labels = _silhouette(
            pca_coords[mask], _labels_for(ann_df, col), asw_max_samples
        )
        try:
            asw = float(silhouette_score(subset, labels))
        except Exception:
            asw = np.nan
        result[f"asw_bio_{col}"] = asw
        # Biology should separate, so +1 is best and the normalization is (s+1)/2.
        result[f"asw_bio_norm_{col}"] = (
            float((asw + 1.0) / 2.0) if not np.isnan(asw) else np.nan
        )

    for col in columns.batch_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            _nan_keys(result, [f"cms_mean_{col}", f"cms_fraction_mixed_{col}"])
            continue
        labels = _labels_for(ann_df, col)
        subset = pca_coords[mask]
        n = len(labels)
        k_cms = min(30, n // 4)
        if k_cms < 2:
            _nan_keys(result, [f"cms_mean_{col}", f"cms_fraction_mixed_{col}"])
            continue
        neighbours = NearestNeighbors(n_neighbors=k_cms + 1, n_jobs=-1)
        neighbours.fit(subset)
        distances, indices = neighbours.kneighbors(subset)
        unique = np.unique(labels)
        p_values: list[float] = []
        for i in range(n):
            neighbour_labels = labels[indices[i, 1:]]
            neighbour_distances = distances[i, 1:]
            samples = [
                neighbour_distances[neighbour_labels == level] for level in unique
            ]
            samples = [group for group in samples if len(group) > 0]
            if len(samples) < 2:
                continue
            try:
                outcome = anderson_ksamp(samples)
                try:
                    p_values.append(float(outcome.pvalue))
                except AttributeError:
                    p_values.append(float(outcome.significance_level) / 100.0)
            except Exception:
                pass
        if p_values:
            result[f"cms_mean_{col}"] = float(np.mean(p_values))
            result[f"cms_fraction_mixed_{col}"] = float(
                np.mean([value > 0.05 for value in p_values])
            )
        else:
            _nan_keys(result, [f"cms_mean_{col}", f"cms_fraction_mixed_{col}"])
    return result


# --------------------------------------------------------------------------------------
# Group C — 2-D embedding structure
# --------------------------------------------------------------------------------------


def compute_group_c(
    umap_coords: np.ndarray,
    tsne_coords: np.ndarray,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
) -> dict:
    """Centroid dispersion and local entropy of batch labels in UMAP and t-SNE space."""
    result: dict = {}
    for col in columns.batch_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            _nan_keys(
                result,
                [
                    f"umap_centroid_disp_{col}",
                    f"umap_entropy_mean_{col}",
                    f"umap_entropy_norm_{col}",
                    f"tsne_centroid_disp_{col}",
                    f"tsne_entropy_mean_{col}",
                    f"tsne_entropy_norm_{col}",
                ],
            )
            continue
        labels = _labels_for(ann_df, col)
        for prefix, coords in (("umap", umap_coords), ("tsne", tsne_coords)):
            result[f"{prefix}_centroid_disp_{col}"] = _centroid_dispersion(
                coords[mask], labels
            )
            mean_entropy, norm_entropy = _local_entropy(coords[mask], labels)
            result[f"{prefix}_entropy_mean_{col}"] = mean_entropy
            result[f"{prefix}_entropy_norm_{col}"] = norm_entropy
    return result


# --------------------------------------------------------------------------------------
# Group D — distribution comparison
# --------------------------------------------------------------------------------------


def _ks_pairwise(exp_sub: pd.DataFrame, labels: np.ndarray) -> tuple[float, float]:
    """Mean KS statistic and BH-significant fraction over all group pairs and genes."""
    from scipy.stats import ks_2samp
    from statsmodels.stats.multitest import multipletests

    unique = np.unique(labels)
    if len(unique) < 2:
        return np.nan, np.nan
    statistics: list[float] = []
    p_values: list[float] = []
    for first, second in itertools.combinations(unique, 2):
        mask_a = labels == first
        mask_b = labels == second
        for gene in exp_sub.columns:
            x = exp_sub.loc[mask_a, gene].dropna().values
            y = exp_sub.loc[mask_b, gene].dropna().values
            if len(x) < 3 or len(y) < 3:
                continue
            statistic, p_value = ks_2samp(x, y)
            statistics.append(float(statistic))
            p_values.append(float(p_value))
    if not statistics:
        return np.nan, np.nan
    _, adjusted, _, _ = multipletests(p_values, method="fdr_bh")
    return float(np.mean(statistics)), float(np.mean(adjusted < 0.05))


def compute_group_d(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    n_genes_max: int = DEFAULT_KS_GENES,
    timeout_s: int = DEFAULT_D1_TIMEOUT_S,
) -> dict:
    """
    Compare per-gene distributions between groups.

    D1 is all-pairwise Kolmogorov-Smirnov over a gene sample, bounded by a hard timeout
    because its cost is quadratic in the number of batches. D2 asks whether cohort
    identity still separates samples *within* one batch. D3 is the per-gene coefficient
    of variation across group means.
    """
    from scipy.stats import ks_2samp
    from statsmodels.stats.multitest import multipletests

    log = get_logger()
    result: dict = {}

    rng = np.random.default_rng(RANDOM_SEED)
    if exp_df.shape[1] > n_genes_max:
        exp_sub = exp_df[rng.choice(exp_df.columns, n_genes_max, replace=False)]
    else:
        exp_sub = exp_df

    d1_cols = tuple(dict.fromkeys([*columns.batch_cols, columns.cohort]))
    d1_result: dict = {}
    d1_error: str | None = None

    def _run_d1() -> None:
        nonlocal d1_error
        try:
            for col in d1_cols:
                mask = _usable_mask(ann_df, col)
                if mask is None:
                    _nan_keys(d1_result, [f"ks_mean_D_{col}", f"ks_frac_sig_{col}"])
                    continue
                labels = _labels_for(ann_df, col)
                aligned = exp_sub.loc[ann_df[col].notna()]
                mean_d, frac_sig = _ks_pairwise(aligned, labels)
                d1_result[f"ks_mean_D_{col}"] = mean_d
                d1_result[f"ks_frac_sig_{col}"] = frac_sig
        except Exception:
            d1_error = traceback.format_exc()

    worker = threading.Thread(target=_run_d1, daemon=True)
    worker.start()
    worker.join(timeout=timeout_s)

    if worker.is_alive():
        log.warning("group D: D1 timed out after %ds — recording nulls", timeout_s)
        for col in d1_cols:
            _nan_keys(result, [f"ks_mean_D_{col}", f"ks_frac_sig_{col}"])
        result["error_D1"] = "timeout"
    elif d1_error:
        log.error("group D: D1 failed:\n%s", d1_error)
        for col in d1_cols:
            _nan_keys(result, [f"ks_mean_D_{col}", f"ks_frac_sig_{col}"])
        result["error_D1"] = d1_error
    else:
        result.update(d1_result)

    batch, cohort = columns.batch, columns.cohort
    if batch in ann_df.columns and cohort in ann_df.columns and batch != cohort:
        statistics: list[float] = []
        p_values: list[float] = []
        for level in ann_df[batch].unique():
            in_batch = (ann_df[batch] == level).values
            cohort_labels = ann_df.loc[in_batch, cohort].astype(str).values
            if len(np.unique(cohort_labels)) < 2:
                continue
            pooled = exp_sub.loc[in_batch]
            for cohort_level in np.unique(cohort_labels):
                selected = cohort_labels == cohort_level
                rest = ~selected
                if selected.sum() < 3 or rest.sum() < 3:
                    continue
                for gene in exp_sub.columns:
                    x = pooled.loc[selected, gene].dropna().values
                    y = pooled.loc[rest, gene].dropna().values
                    if len(x) < 3 or len(y) < 3:
                        continue
                    statistic, p_value = ks_2samp(x, y)
                    statistics.append(float(statistic))
                    p_values.append(float(p_value))
        if statistics:
            _, adjusted, _, _ = multipletests(p_values, method="fdr_bh")
            result["ks_cohort_within_batch_mean_D"] = float(np.mean(statistics))
            result["ks_cohort_within_batch_frac_sig"] = float(np.mean(adjusted < 0.05))
        else:
            _nan_keys(
                result,
                ["ks_cohort_within_batch_mean_D", "ks_cohort_within_batch_frac_sig"],
            )
    else:
        _nan_keys(
            result, ["ks_cohort_within_batch_mean_D", "ks_cohort_within_batch_frac_sig"]
        )

    for col in d1_cols:
        if col not in ann_df.columns or ann_df[col].nunique() < 2:
            result[f"per_gene_batch_mean_cv_{col}"] = np.nan
            continue
        group_means = np.array(
            [
                exp_df.loc[(ann_df[col] == level).values].mean(axis=0).values
                for level in ann_df[col].unique()
            ]
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            cv = np.std(group_means, axis=0, ddof=1) / np.abs(
                np.mean(group_means, axis=0)
            )
        finite = cv[np.isfinite(cv)]
        result[f"per_gene_batch_mean_cv_{col}"] = (
            float(np.mean(finite)) if len(finite) else np.nan
        )
    return result


# --------------------------------------------------------------------------------------
# Group F — variancePartition (R)
# --------------------------------------------------------------------------------------


def compute_group_f(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    n_genes: int = 2000,
    timeout_s: int = 3600,
) -> dict:
    """
    Partition per-gene variance across the annotation columns, via R's variancePartition.

    The cohort column enters as a random effect and everything else as a fixed effect,
    matching the donor's model. Runs ``Rscript`` over temporary files rather than rpy2,
    so it needs no R session in this process.

    Returns
    -------
    ``vp_median_{col}`` per column, the batch column's interquartile points, and the
    median residual variance.
    """
    import os
    import subprocess
    import tempfile

    result: dict = {}
    available = [
        col
        for col in columns.all_cols
        if col in ann_df.columns and ann_df[col].nunique() >= 2
    ]
    if not available:
        for col in columns.all_cols:
            result[f"vp_median_{col}"] = np.nan
        _nan_keys(
            result,
            [
                f"vp_p25_{columns.batch}",
                f"vp_p75_{columns.batch}",
                "vp_median_residual",
            ],
        )
        return result

    rng = np.random.default_rng(RANDOM_SEED)
    exp_sub = (
        exp_df[rng.choice(exp_df.columns, n_genes, replace=False)]
        if exp_df.shape[1] > n_genes
        else exp_df
    )

    with tempfile.TemporaryDirectory() as tmpdir:
        exp_path = os.path.join(tmpdir, "exp.tsv")
        ann_path = os.path.join(tmpdir, "ann.tsv")
        out_path = os.path.join(tmpdir, "vp.tsv")
        script_path = os.path.join(tmpdir, "run_vp.R")

        exp_sub.T.to_csv(exp_path, sep="\t")
        ann_df[available].to_csv(ann_path, sep="\t")

        cohort_is_random = columns.cohort in available
        fixed = [
            col for col in available if not (cohort_is_random and col == columns.cohort)
        ]
        parts = ([f"(1|{columns.cohort})"] if cohort_is_random else []) + fixed
        formula = "~ " + " + ".join(parts) if parts else "~ 1"

        script = f"""\
library(variancePartition)
library(BiocParallel)
suppressMessages(register(MulticoreParam(2)))

exp  <- as.matrix(read.table("{exp_path}", header=TRUE, row.names=1, sep="\\t",
                             check.names=FALSE))
info <- read.table("{ann_path}", header=TRUE, row.names=1, sep="\\t",
                   check.names=FALSE)

for (col in colnames(info)) {{
  info[[col]] <- as.factor(info[[col]])
}}

vp <- fitExtractVarPartModel(exp, as.formula("{formula}"), info)
write.table(as.data.frame(vp), "{out_path}", sep="\\t", quote=FALSE)
"""
        with open(script_path, "w") as handle:
            handle.write(script)

        completed = subprocess.run(
            ["Rscript", "--no-save", "--no-restore", script_path],
            capture_output=True,
            text=True,
            timeout=timeout_s,
        )
        if completed.returncode != 0:
            raise RuntimeError(f"variancePartition failed:\n{completed.stderr[-2000:]}")
        vp_df = pd.read_csv(out_path, sep="\t", index_col=0)

    for col in columns.all_cols:
        # R silently rewrites '_' to '.' in data-frame names, so both spellings are tried.
        r_name = col.replace("_", ".")
        if col in vp_df.columns:
            result[f"vp_median_{col}"] = float(vp_df[col].median())
        elif r_name in vp_df.columns:
            result[f"vp_median_{col}"] = float(vp_df[r_name].median())
        else:
            result[f"vp_median_{col}"] = np.nan

    batch_name = next(
        (
            name
            for name in (columns.batch, columns.batch.replace("_", "."))
            if name in vp_df.columns
        ),
        None,
    )
    if batch_name is not None:
        result[f"vp_p25_{columns.batch}"] = float(vp_df[batch_name].quantile(0.25))
        result[f"vp_p75_{columns.batch}"] = float(vp_df[batch_name].quantile(0.75))
    else:
        _nan_keys(result, [f"vp_p25_{columns.batch}", f"vp_p75_{columns.batch}"])

    residual = next(
        (name for name in ("Residuals", "residuals") if name in vp_df.columns), None
    )
    result["vp_median_residual"] = (
        float(vp_df[residual].median()) if residual else np.nan
    )
    return result


# --------------------------------------------------------------------------------------
# Group G — graph connectivity
# --------------------------------------------------------------------------------------


def compute_group_g(
    pca_coords: np.ndarray, ann_df: pd.DataFrame, columns: MetricColumns
) -> dict:
    """Ask whether each biology group forms one connected kNN component, or several."""
    from scipy.sparse.csgraph import connected_components
    from sklearn.neighbors import kneighbors_graph

    result: dict = {}
    for col in columns.bio_cols:
        mask = _usable_mask(ann_df, col)
        if mask is None:
            result[f"graph_connectivity_{col}"] = np.nan
            continue
        labels = _labels_for(ann_df, col)
        subset = pca_coords[mask]
        scores: list[float] = []
        for level in np.unique(labels):
            in_group = labels == level
            n_group = in_group.sum()
            if n_group < 2:
                scores.append(1.0)
                continue
            graph = kneighbors_graph(
                subset[in_group],
                n_neighbors=min(15, n_group - 1),
                mode="connectivity",
                include_self=False,
            )
            n_components, _ = connected_components(graph, directed=False)
            scores.append(1.0 / n_components)
        result[f"graph_connectivity_{col}"] = (
            float(np.mean(scores)) if scores else np.nan
        )
    return result


# --------------------------------------------------------------------------------------
# Group H — pairwise distances
# --------------------------------------------------------------------------------------


def compute_group_h(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    subsample: int = 2000,
) -> dict:
    """
    Contrast average within-group and between-group distances in full gene space.

    Returns
    -------
    ``avg_intra_dist_{col}``, ``avg_inter_dist_{col}`` and their ratio. A ratio near 1
    means the column no longer structures the distances at all.
    """
    from sklearn.metrics import pairwise_distances

    result: dict = {}
    n = len(exp_df)
    if n > subsample:
        rng = np.random.default_rng(RANDOM_SEED)
        idx = rng.choice(n, subsample, replace=False)
        values = exp_df.values[idx].astype(float)
        ann_sub = ann_df.iloc[idx]
    else:
        values = exp_df.values.astype(float)
        ann_sub = ann_df

    distances = pairwise_distances(np.nan_to_num(values, nan=0.0), metric="euclidean")
    n_sub = len(ann_sub)
    upper = np.triu(np.ones((n_sub, n_sub), dtype=bool), k=1)

    for col in columns.all_cols:
        if col not in ann_sub.columns:
            _nan_keys(
                result,
                [f"avg_intra_dist_{col}", f"avg_inter_dist_{col}", f"dist_ratio_{col}"],
            )
            continue
        # np.asarray forces a plain ndarray: a pandas StringArray does not support the
        # [:, None] / [None, :] broadcasting below.
        labels = np.asarray(ann_sub[col].astype(str))
        same = (labels[:, None] == labels[None, :]) & upper
        different = (labels[:, None] != labels[None, :]) & upper
        intra = float(distances[same].mean()) if same.any() else np.nan
        inter = float(distances[different].mean()) if different.any() else np.nan
        result[f"avg_intra_dist_{col}"] = intra
        result[f"avg_inter_dist_{col}"] = inter
        result[f"dist_ratio_{col}"] = (
            float(intra / inter)
            if not np.isnan(intra) and not np.isnan(inter) and inter > 1e-9
            else np.nan
        )
    return result


# --------------------------------------------------------------------------------------
# Group K — NA retention
# --------------------------------------------------------------------------------------


def compute_group_k(exp_df: pd.DataFrame) -> dict:
    """Count what survived: genes and samples free of NA, and how many cells are NA."""
    result: dict = {}
    n_samples, n_genes = exp_df.shape
    missing = exp_df.isna()

    n_genes_clean = int((~missing.any(axis=0)).sum())
    n_samples_clean = int((~missing.any(axis=1)).sum())
    n_na_cells = int(missing.sum().sum())

    result["n_genes_noNA"] = n_genes_clean
    result["pct_genes_noNA"] = (
        float(n_genes_clean / n_genes * 100) if n_genes else np.nan
    )
    result["n_samples_noNA"] = n_samples_clean
    result["pct_samples_noNA"] = (
        float(n_samples_clean / n_samples * 100) if n_samples else np.nan
    )
    result["n_genes_allNA"] = int(missing.all(axis=0).sum())
    result["n_samples_allNA"] = int(missing.all(axis=1).sum())
    result["n_na_cells"] = n_na_cells
    result["pct_na_cells"] = (
        float(n_na_cells / (n_samples * n_genes) * 100)
        if n_samples * n_genes
        else np.nan
    )
    return result


# --------------------------------------------------------------------------------------
# Group I — WaterMelon score
# --------------------------------------------------------------------------------------


def _entropy_from_counts(counts: np.ndarray) -> float:
    """Shannon entropy, in bits, of a count vector."""
    total = float(counts.sum())
    if total <= 0:
        return 0.0
    proportions = counts[counts > 0] / total
    return float(-np.sum(proportions * np.log2(proportions)))


def _wm_trajectory(
    linkage: np.ndarray, label_idx: np.ndarray, n_classes: int, n: int
) -> tuple[np.ndarray, float]:
    """
    Information-gain trajectory along a bottom-up hierarchical merge.

    Conditional entropy is updated incrementally at each merge rather than recomputed,
    which is what makes the permutation null affordable.
    """
    global_counts = np.bincount(label_idx, minlength=n_classes).astype(float)
    h_labels = _entropy_from_counts(global_counts)
    if h_labels < 1e-12:
        return np.zeros(n - 1), 0.0

    counts: dict[int, np.ndarray] = {}
    sizes: dict[int, int] = {}
    for i in range(n):
        vector = np.zeros(n_classes)
        vector[label_idx[i]] = 1.0
        counts[i] = vector
        sizes[i] = 1

    conditional = 0.0
    trajectory = np.empty(n - 1)
    for step in range(n - 1):
        left, right = int(linkage[step, 0]), int(linkage[step, 1])
        counts_left, counts_right = counts.pop(left), counts.pop(right)
        size_left, size_right = sizes.pop(left), sizes.pop(right)

        merged_counts = counts_left + counts_right
        merged_size = size_left + size_right
        conditional = (
            conditional
            - (size_left / n) * _entropy_from_counts(counts_left)
            - (size_right / n) * _entropy_from_counts(counts_right)
            + (merged_size / n) * _entropy_from_counts(merged_counts)
        )
        trajectory[step] = (h_labels - conditional) / h_labels

        counts[n + step] = merged_counts
        sizes[n + step] = merged_size
    return trajectory, h_labels


def _stratified_subsample(
    n_total: int,
    strata: np.ndarray | None,
    max_n: int,
    rng: np.random.Generator,
) -> np.ndarray:
    """Positional indices of a subsample that preserves the strata proportions."""
    indices = np.arange(n_total)
    if strata is None or max_n >= n_total:
        return rng.choice(indices, min(max_n, n_total), replace=False)

    selected: list[np.ndarray] = []
    for level in np.unique(strata):
        level_idx = indices[strata == level]
        take = min(max(1, round(max_n * len(level_idx) / n_total)), len(level_idx))
        selected.append(rng.choice(level_idx, take, replace=False))
    chosen = np.concatenate(selected)

    if len(chosen) > max_n:
        return rng.choice(chosen, max_n, replace=False)
    if len(chosen) < max_n:
        remaining = np.setdiff1d(indices, chosen)
        extra = min(max_n - len(chosen), len(remaining))
        if extra > 0:
            chosen = np.concatenate(
                [chosen, rng.choice(remaining, extra, replace=False)]
            )
    return chosen


def compute_group_i(
    ann_df: pd.DataFrame,
    pca_coords: np.ndarray,
    columns: MetricColumns,
    *,
    n_permutations: int = DEFAULT_WM_PERMUTATIONS,
    max_samples: int = DEFAULT_WM_MAX_SAMPLES,
) -> dict:
    """
    Score how much of each label's structure a Ward dendrogram recovers, against a null.

    The reported value is the mean gap between the observed information-gain trajectory
    and the 95th percentile of a label-permutation null: positive means the label
    genuinely structures the hierarchy.
    """
    from scipy.cluster.hierarchy import linkage as scipy_linkage

    result: dict = {}
    rng = np.random.default_rng(RANDOM_SEED)
    n_total = len(ann_df)
    subsampled = n_total > max_samples
    result["wm_subsampled"] = subsampled

    if subsampled:
        strata = (
            ann_df[columns.batch].values.astype(str)
            if columns.batch in ann_df.columns
            else None
        )
        chosen = _stratified_subsample(n_total, strata, max_samples, rng)
        pca_sub = pca_coords[chosen, :]
        ann_sub = ann_df.iloc[chosen]
    else:
        pca_sub = pca_coords
        ann_sub = ann_df
    n_eff = len(ann_sub)

    n_pcs = min(50, pca_sub.shape[1])
    linkage = scipy_linkage(pca_sub[:, :n_pcs], method="ward", metric="euclidean")

    batch_scores: list[float] = []
    bio_scores: list[float] = []

    for col in [*columns.batch_cols, *columns.bio_cols]:
        if _usable_mask(ann_sub, col) is None:
            result[f"wm_{col}"] = np.nan
            continue
        labels = ann_sub[col].astype(str).fillna("__unknown__").values.astype(str)
        classes, raw_idx = np.unique(labels, return_inverse=True)
        label_idx = raw_idx.astype(np.int32)

        observed, h_labels = _wm_trajectory(linkage, label_idx, len(classes), n_eff)
        if h_labels < 1e-12:
            result[f"wm_{col}"] = np.nan
            continue

        null = np.empty((n_permutations, n_eff - 1))
        for m in range(n_permutations):
            null[m], _ = _wm_trajectory(
                linkage, rng.permutation(label_idx), len(classes), n_eff
            )
        score = float(np.mean(observed - np.percentile(null, 95, axis=0)))
        result[f"wm_{col}"] = score

        if col in columns.batch_cols:
            batch_scores.append(score)
        else:
            bio_scores.append(score)

    valid_batch = [value for value in batch_scores if not np.isnan(value)]
    valid_bio = [value for value in bio_scores if not np.isnan(value)]
    result["wm_mean_batch"] = float(np.mean(valid_batch)) if valid_batch else np.nan
    result["wm_mean_bio"] = float(np.mean(valid_bio)) if valid_bio else np.nan
    result["wm_ratio_bio_batch"] = (
        float(np.mean(valid_bio) / (np.mean(valid_batch) + 0.01))
        if valid_batch and valid_bio
        else np.nan
    )
    return result


# --------------------------------------------------------------------------------------
# Groups L / M / N — shared rank machinery
# --------------------------------------------------------------------------------------


def _rank_along(values: np.ndarray, axis: int) -> np.ndarray:
    """Average ranks along an axis, leaving NaN in place to be masked later."""
    from scipy.stats import rankdata

    return rankdata(values, axis=axis, nan_policy="omit")


def _unit_center(ranks: np.ndarray, axis: int) -> np.ndarray:
    """Centre and L2-normalize, so a dot product along the axis is a correlation."""
    centred = ranks - ranks.mean(axis=axis, keepdims=True)
    norm = np.sqrt((centred**2).sum(axis=axis, keepdims=True))
    norm[norm == 0] = np.nan  # a constant vector has no defined correlation
    return centred / norm


def _spearman_columnwise(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    """
    Per-column Spearman correlation between two same-shaped matrices.

    One vectorized pass instead of a scipy call per column.

    Returns
    -------
    An array of length ``n_columns``, NaN where a column is constant.
    """
    ranks_a = _unit_center(_rank_along(first, axis=0), axis=0)
    ranks_b = _unit_center(_rank_along(second, axis=0), axis=0)
    return np.nansum(ranks_a * ranks_b, axis=0) * np.where(
        np.isnan(ranks_a).all(axis=0) | np.isnan(ranks_b).all(axis=0), np.nan, 1.0
    )


def _stratified_subsample_idx(
    ann_df: pd.DataFrame, max_samples: int, strat_col: str
) -> np.ndarray:
    """Sorted positional indices of a stratified subsample, or all rows if small enough."""
    n = len(ann_df)
    if n <= max_samples:
        return np.arange(n)
    rng = np.random.default_rng(RANDOM_SEED)
    if strat_col not in ann_df.columns:
        return np.sort(rng.choice(n, max_samples, replace=False))
    positions = np.arange(n)
    strata = ann_df[strat_col].astype(str).values
    fraction = max_samples / n
    keep: list[int] = []
    for level in np.unique(strata):
        level_idx = positions[strata == level]
        take = min(max(1, int(round(len(level_idx) * fraction))), len(level_idx))
        keep.extend(rng.choice(level_idx, take, replace=False).tolist())
    return np.sort(np.array(keep))


# --------------------------------------------------------------------------------------
# Group L — expression profile preservation
# --------------------------------------------------------------------------------------


def compute_group_l(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    ref_df: pd.DataFrame,
    panel: Panel | None = None,
    min_cohort_n: int = DEFAULT_MIN_COHORT_N,
    collect_detail: bool = False,
) -> dict:
    """
    Ask whether harmonization preserved each gene's profile within each cohort.

    For every cohort and gene, the Spearman correlation between the gene's expression
    across that cohort's samples before harmonization and after.

    Note the argument order: the donor's signature was ``(exp, ref, ann)``, out of step
    with every other group's ``(exp, ann, ...)``. The reference is a keyword here, so the
    two can never be transposed silently.

    **Known ceiling:** any harmonizer applying a per-batch *monotone* transform leaves
    within-cohort ranks untouched and scores ~1.0 here by construction. This is a
    "nothing broke" guard rail; group M's margin carries the discriminative signal.

    Parameters
    ----------
    exp_df
        Harmonized expression, samples x genes.
    ann_df
        Annotation aligned to ``exp_df``.
    columns
        Resolved column roles.
    ref_df
        The pre-harmonization matrix, aligned to ``exp_df``'s rows.
    panel
        Marker panel, or ``None`` for every shared gene.
    min_cohort_n
        Cohorts smaller than this are skipped as too noisy to correlate.
    collect_detail
        Also return the full gene-by-cohort matrix. At roughly 300 KB per job it is off
        by default; ``concat`` splits it into a long table when present.
    """
    result: dict = {}
    shared = exp_df.columns.intersection(ref_df.columns)
    genes, missing = resolve_panel(shared, panel)

    result["mk_n_panel_genes_used"] = len(genes)
    result["mk_n_genes_null"] = len(missing)
    result["mk_panel_coverage_frac"] = (
        float(len(genes) / len(panel)) if panel and len(panel) else np.nan
    )
    result["mk_is_self_reference"] = bool(exp_df is ref_df)

    summary_keys = (
        "mk_rho_mean_all_genes",
        "mk_rho_median_all_genes",
        "mk_rho_p10_all_genes",
        "mk_rho_min_all_genes",
        "mk_rho_frac_genes_above_0.9",
        "mk_rho_mean_markers",
        "mk_rho_mean_housekeeping",
        "mk_rho_marker_minus_hk",
        "mk_rho_mean_by_cohort_mean",
    )
    if not genes:
        _nan_keys(result, summary_keys)
        result["mk_rho_by_gene"] = {}
        result["mk_rho_n_cohorts_by_gene"] = {}
        result["mk_rho_by_cohort"] = {}
        result["mk_n_cohorts_used"] = 0
        result["mk_n_cohorts_skipped_small"] = 0
        return result

    after = exp_df[genes].values.astype(float)
    before = ref_df[genes].values.astype(float)
    cohorts = ann_df[columns.cohort].astype(str).values

    # Accumulated per gene across cohorts, so a gene that is constant in one cohort
    # still contributes from the others rather than being dropped entirely.
    sums = np.zeros(len(genes))
    counts = np.zeros(len(genes))
    per_cohort: dict[str, float] = {}
    detail: dict[str, dict[str, float | None]] = {}
    n_skipped = 0

    for level in np.unique(cohorts):
        mask = cohorts == level
        if mask.sum() < min_cohort_n:
            n_skipped += 1
            continue
        rho = _spearman_columnwise(after[mask], before[mask])
        valid = ~np.isnan(rho)
        sums[valid] += rho[valid]
        counts[valid] += 1
        if valid.any():
            per_cohort[str(level)] = float(np.mean(rho[valid]))
        if collect_detail:
            detail[str(level)] = {
                gene: (float(value) if not np.isnan(value) else None)
                for gene, value in zip(genes, rho)
            }

    with np.errstate(invalid="ignore", divide="ignore"):
        by_gene = np.where(counts > 0, sums / counts, np.nan)

    result["mk_rho_by_gene"] = {
        gene: (float(value) if not np.isnan(value) else None)
        for gene, value in zip(genes, by_gene)
    }
    result["mk_rho_n_cohorts_by_gene"] = {
        gene: int(count) for gene, count in zip(genes, counts)
    }
    result["mk_rho_by_cohort"] = per_cohort
    result["mk_n_cohorts_used"] = len(per_cohort)
    result["mk_n_cohorts_skipped_small"] = n_skipped

    finite = by_gene[~np.isnan(by_gene)]
    if finite.size:
        result["mk_rho_mean_all_genes"] = float(np.mean(finite))
        result["mk_rho_median_all_genes"] = float(np.median(finite))
        result["mk_rho_p10_all_genes"] = float(np.percentile(finite, 10))
        result["mk_rho_min_all_genes"] = float(np.min(finite))
        result["mk_rho_frac_genes_above_0.9"] = float(np.mean(finite > 0.9))
    else:
        _nan_keys(result, summary_keys[:5])

    housekeeping = panel.housekeeping if panel else frozenset()
    is_hk = np.array([gene in housekeeping for gene in genes])
    marker_values = by_gene[~is_hk & ~np.isnan(by_gene)]
    hk_values = by_gene[is_hk & ~np.isnan(by_gene)]
    result["mk_rho_mean_markers"] = (
        float(np.mean(marker_values)) if marker_values.size else np.nan
    )
    result["mk_rho_mean_housekeeping"] = (
        float(np.mean(hk_values)) if hk_values.size else np.nan
    )
    result["mk_rho_marker_minus_hk"] = (
        float(np.mean(marker_values) - np.mean(hk_values))
        if marker_values.size and hk_values.size
        else np.nan
    )
    result["mk_rho_mean_by_cohort_mean"] = (
        float(np.mean(list(per_cohort.values()))) if per_cohort else np.nan
    )
    if collect_detail:
        result["mk_gene_cohort_detail"] = detail
    return result


# --------------------------------------------------------------------------------------
# Group M — cross-batch rank agreement
# --------------------------------------------------------------------------------------


def compute_group_m(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    panel: Panel | None = None,
    max_samples: int = DEFAULT_XB_MAX_SAMPLES,
) -> dict:
    """
    Measure whether samples of the same biology agree across batches, and only those.

    Ranks the panel genes within each sample, then averages the Spearman correlation
    over sample pairs drawn from different batches — separately for pairs sharing a
    biology label and pairs that do not. The **margin** between the two is the signal:
    a method that merely collapses everything toward a common profile raises both and
    leaves the margin flat.
    """
    result: dict = {}
    genes, _ = resolve_panel(exp_df.columns, panel)
    batch, bio = columns.batch, columns.bio

    if len(genes) < 3 or batch not in ann_df.columns or bio not in ann_df.columns:
        _nan_keys(
            result,
            [
                "xb_rank_agree",
                "xb_rank_agree_median",
                "xb_rank_disagree_diffbio",
                "xb_rank_disagree_diffbio_median",
                "xb_rank_agree_ratio",
            ],
        )
        result["xb_rank_agree_by_bio"] = {}
        result["xb_n_pairs_same_bio"] = 0
        result["xb_n_pairs_diff_bio"] = 0
        result["xb_n_samples_used"] = 0
        result["xb_subsampled"] = False
        result["xb_n_panel_genes_used"] = len(genes)
        return result

    keep = _stratified_subsample_idx(ann_df, max_samples, batch)
    result["xb_subsampled"] = bool(len(keep) < len(ann_df))
    result["xb_n_samples_used"] = int(len(keep))

    values = exp_df[genes].values.astype(float)[keep]
    batches = np.asarray(ann_df[batch].astype(str))[keep]
    bios = np.asarray(ann_df[bio].astype(str))[keep]

    ranks = _unit_center(_rank_along(values, axis=1), axis=1)
    ranks = np.nan_to_num(ranks, nan=0.0)

    n = ranks.shape[0]
    sums = {"same": 0.0, "diff": 0.0}
    counts = {"same": 0, "diff": 0}
    collected: dict[str, list[np.ndarray]] = {"same": [], "diff": []}
    per_bio_sum: dict[str, float] = {}
    per_bio_count: dict[str, int] = {}

    # Block-wise: the full n x n correlation matrix is ~412 MB at n = 7,000, and each
    # worker already holds two expression matrices.
    block = 512
    for start in range(0, n, block):
        stop = min(start + block, n)
        correlations = ranks[start:stop] @ ranks.T
        rows = np.arange(start, stop)
        upper = rows[:, None] < np.arange(n)[None, :]
        cross_batch = batches[rows][:, None] != batches[None, :]
        same_bio = bios[rows][:, None] == bios[None, :]

        for tag, selector in (
            ("same", upper & cross_batch & same_bio),
            ("diff", upper & cross_batch & ~same_bio),
        ):
            if selector.any():
                values_here = correlations[selector]
                sums[tag] += float(values_here.sum())
                counts[tag] += int(values_here.size)
                collected[tag].append(values_here)

        for level in np.unique(bios[rows]):
            selector = upper & cross_batch & same_bio & (bios[rows][:, None] == level)
            if selector.any():
                values_here = correlations[selector]
                per_bio_sum[level] = per_bio_sum.get(level, 0.0) + float(
                    values_here.sum()
                )
                per_bio_count[level] = per_bio_count.get(level, 0) + int(
                    values_here.size
                )

    def mean_of(tag: str) -> float:
        return float(sums[tag] / counts[tag]) if counts[tag] else np.nan

    def median_of(tag: str) -> float:
        return (
            float(np.median(np.concatenate(collected[tag])))
            if collected[tag]
            else np.nan
        )

    result["xb_rank_agree"] = mean_of("same")
    result["xb_rank_agree_median"] = median_of("same")
    result["xb_rank_disagree_diffbio"] = mean_of("diff")
    result["xb_rank_disagree_diffbio_median"] = median_of("diff")
    result["xb_rank_agree_ratio"] = (
        result["xb_rank_agree"] - result["xb_rank_disagree_diffbio"]
        if counts["same"] and counts["diff"]
        else np.nan
    )
    result["xb_rank_agree_by_bio"] = {
        level: float(per_bio_sum[level] / per_bio_count[level]) for level in per_bio_sum
    }
    result["xb_n_pairs_same_bio"] = counts["same"]
    result["xb_n_pairs_diff_bio"] = counts["diff"]
    result["xb_n_panel_genes_used"] = len(genes)
    return result


# --------------------------------------------------------------------------------------
# Group N — predictive validation
# --------------------------------------------------------------------------------------


def predictor_labels(
    ann_df: pd.DataFrame,
    col: str,
    *,
    n_classes: int = 3,
    classes: Sequence[str] | None = None,
) -> pd.Series | None:
    """
    Build the classification target for group N.

    The donor collapsed one dataset's diagnosis strings into three hardcoded classes.
    Here the target is either the classes the user named, or simply the ``n_classes``
    most frequent levels of the configured column — dataset-agnostic, and the result
    reports the class count actually used.

    Parameters
    ----------
    ann_df
        Annotation table.
    col
        The column holding the class label.
    n_classes
        How many of the most frequent levels to keep, when ``classes`` is not given.
    classes
        Explicit class labels; every other sample becomes NaN and is excluded.

    Returns
    -------
    Labels indexed like ``ann_df``, NaN where excluded, or ``None`` if fewer than two
    classes survive.
    """
    if col not in ann_df.columns:
        return None
    values = ann_df[col].astype(str)

    if classes:
        kept = values.where(values.isin(set(classes)))
    else:
        top = values.value_counts().index[:n_classes]
        kept = values.where(values.isin(set(top)))

    if kept.dropna().nunique() < 2:
        return None
    return kept


def _safe_auc(
    y_true: np.ndarray, proba: np.ndarray, classes: np.ndarray
) -> float | None:
    """
    One-vs-rest macro AUC over the classes actually present in ``y_true``.

    sklearn's ``multi_class="ovr"`` requires every declared class to appear in the truth
    vector, which never holds for a held-out batch, so the average is assembled here.
    """
    from sklearn.metrics import roc_auc_score

    present = set(np.unique(y_true))
    if len(present) < 2:
        return None
    scores: list[float] = []
    for index, level in enumerate(classes):
        if level not in present:
            continue
        binary = (y_true == level).astype(int)
        if binary.sum() in (0, len(binary)):
            continue
        try:
            scores.append(float(roc_auc_score(binary, proba[:, index])))
        except ValueError:
            continue
    return float(np.mean(scores)) if scores else None


def _fold_pca_projections(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    batch_col: str,
    min_test_n: int,
    n_pcs: int,
) -> tuple[list[dict], dict]:
    """
    Fit one PCA per leave-one-batch-out fold, on the training batches only.

    The projections are label-independent, so they are computed once and reused by the
    observed run and by every permutation. NA cells are dropped, never imputed: the
    matrix was already imputed upstream, so a residual NA means the harmonizer failed
    and filling it would invent data.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    values = exp_df.values.astype(float)

    # Order matters: an all-NA sample is a failed sample and must be counted as one.
    # Dropping genes first would absorb whole failed samples into the gene count.
    sample_ok = ~np.isnan(values).all(axis=1)
    n_samples_dropped = int((~sample_ok).sum())
    positions = np.arange(len(exp_df))[sample_ok]
    values = values[sample_ok]

    gene_ok = (
        ~np.isnan(values).any(axis=0)
        if values.size
        else np.zeros(values.shape[1], bool)
    )
    n_genes_dropped = int((~gene_ok).sum())
    values = values[:, gene_ok]

    diagnostics = {
        "pv_n_genes_dropped_na": n_genes_dropped,
        "pv_n_samples_dropped_na": n_samples_dropped,
    }
    if values.shape[1] < n_pcs or values.shape[0] < 2 * min_test_n:
        return [], diagnostics

    batches = ann_df[batch_col].astype(str).values[sample_ok]
    folds: list[dict] = []
    for level in np.unique(batches):
        test_mask = batches == level
        if test_mask.sum() < min_test_n or (~test_mask).sum() < n_pcs + 1:
            continue
        scaler = StandardScaler().fit(values[~test_mask])
        pca = PCA(n_components=n_pcs, svd_solver="randomized", random_state=RANDOM_SEED)
        folds.append(
            {
                "batch": str(level),
                "train_pos": positions[~test_mask],
                "test_pos": positions[test_mask],
                "Z_train": pca.fit_transform(scaler.transform(values[~test_mask])),
                "Z_test": pca.transform(scaler.transform(values[test_mask])),
            }
        )
    return folds, diagnostics


def _run_lobo(folds: list[dict], labels: pd.Series, collect_folds: bool) -> dict:
    """
    Evaluate leave-one-batch-out classification over pre-computed PCA folds.

    F1 is computed on every fold, single-class test batches included: a predictor facing
    a one-class batch should still label as many of its samples correctly as it can,
    which is the honest deployment test. AUC is restricted to folds holding at least two
    classes, because it is otherwise undefined.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import (
        balanced_accuracy_score,
        f1_score,
        matthews_corrcoef,
    )

    y_all = labels.values
    out: dict = {
        "f1_macro": [],
        "f1_weighted": [],
        "bal_acc": [],
        "mcc": [],
        "auc_macro": [],
        "per_class": {},
        "n_folds": 0,
        "n_folds_multiclass": 0,
        "folds": {},
    }

    for fold in folds:
        y_train = y_all[fold["train_pos"]]
        y_test = y_all[fold["test_pos"]]
        train_ok = pd.notna(y_train)
        test_ok = pd.notna(y_test)
        if train_ok.sum() < 10 or test_ok.sum() < 1:
            continue
        y_train_clean = y_train[train_ok].astype(str)
        y_test_clean = y_test[test_ok].astype(str)
        if len(set(y_train_clean)) < 2:
            continue

        model = LogisticRegression(
            max_iter=2000, class_weight="balanced", random_state=RANDOM_SEED
        )
        model.fit(fold["Z_train"][train_ok], y_train_clean)
        z_test = fold["Z_test"][test_ok]
        predicted = model.predict(z_test)

        present = sorted(set(y_test_clean))
        f1_macro = float(
            f1_score(
                y_test_clean,
                predicted,
                labels=present,
                average="macro",
                zero_division=0,
            )
        )
        auc = _safe_auc(y_test_clean, model.predict_proba(z_test), model.classes_)

        out["f1_macro"].append(f1_macro)
        out["f1_weighted"].append(
            float(
                f1_score(
                    y_test_clean,
                    predicted,
                    labels=present,
                    average="weighted",
                    zero_division=0,
                )
            )
        )
        out["bal_acc"].append(float(balanced_accuracy_score(y_test_clean, predicted)))
        out["mcc"].append(float(matthews_corrcoef(y_test_clean, predicted)))
        out["n_folds"] += 1
        if len(present) >= 2:
            out["n_folds_multiclass"] += 1
        if auc is not None:
            out["auc_macro"].append(auc)

        for level, value in zip(
            present,
            f1_score(
                y_test_clean, predicted, labels=present, average=None, zero_division=0
            ),
        ):
            out["per_class"].setdefault(level, []).append(float(value))

        if collect_folds:
            out["folds"][fold["batch"]] = {
                "n": int(test_ok.sum()),
                "n_classes": len(present),
                "f1_macro": round(f1_macro, 4),
                "auc": round(auc, 4) if auc is not None else None,
            }
    return out


def compute_group_n(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    columns: MetricColumns,
    *,
    n_pcs: int = DEFAULT_N_PCS,
    n_perm: int = DEFAULT_N_PERM,
    min_test_n: int = DEFAULT_MIN_TEST_N,
    predict_classes: Sequence[str] | None = None,
) -> dict:
    """
    Predict biology from a held-out batch, and check it beats a permuted null.

    A logistic regression on the leading principal components is trained on every batch
    but one and evaluated on the held-out batch. PCA is fitted on the training batches
    alone, so no test-set structure leaks into the features. Two targets are evaluated —
    the three and two most frequent classes — and a label-permutation control is run for
    each, cheaply, because the per-fold PCA is label-independent.

    Returns
    -------
    ``pv_lobo3_*`` and ``pv_lobo2_*`` metrics, their permutation means, the observed
    minus null delta, and a permutation p-value.
    """
    result: dict = {"pv_n_perm": int(n_perm)}
    batch = columns.batch

    if batch not in ann_df.columns:
        result["pv_n_genes_dropped_na"] = 0
        result["pv_n_samples_dropped_na"] = 0
        for prefix in ("pv_lobo3", "pv_lobo2"):
            result[f"{prefix}_n_folds"] = 0
        return result

    folds, diagnostics = _fold_pca_projections(exp_df, ann_df, batch, min_test_n, n_pcs)
    result.update(diagnostics)
    rng = np.random.default_rng(RANDOM_SEED)

    empty_suffixes = (
        "f1_macro_mean",
        "f1_macro_median",
        "f1_macro_std",
        "f1_macro_min",
        "f1_weighted_mean",
        "bal_acc_mean",
        "mcc_mean",
        "auc_macro_mean",
        "f1_macro_perm_mean",
        "f1_macro_perm_std",
        "auc_macro_perm_mean",
        "f1_macro_delta",
        "f1_perm_pvalue",
    )

    for prefix, n_classes in (("pv_lobo3", 3), ("pv_lobo2", 2)):
        labels = predictor_labels(
            ann_df,
            columns.predict_class_col,
            n_classes=n_classes,
            classes=predict_classes,
        )
        if labels is None or not folds:
            result[f"{prefix}_n_folds"] = 0
            result[f"{prefix}_n_folds_multiclass"] = 0
            _nan_keys(result, [f"{prefix}_{suffix}" for suffix in empty_suffixes])
            result[f"{prefix}_folds"] = {}
            continue

        result[f"{prefix}_n_classes_used"] = int(labels.dropna().nunique())
        observed = _run_lobo(folds, labels, collect_folds=True)

        def aggregate(values: list[float], fn) -> float:
            return float(fn(values)) if values else np.nan

        result[f"{prefix}_f1_macro_mean"] = aggregate(observed["f1_macro"], np.mean)
        result[f"{prefix}_f1_macro_median"] = aggregate(observed["f1_macro"], np.median)
        result[f"{prefix}_f1_macro_std"] = aggregate(observed["f1_macro"], np.std)
        result[f"{prefix}_f1_macro_min"] = aggregate(observed["f1_macro"], np.min)
        result[f"{prefix}_f1_weighted_mean"] = aggregate(
            observed["f1_weighted"], np.mean
        )
        result[f"{prefix}_bal_acc_mean"] = aggregate(observed["bal_acc"], np.mean)
        result[f"{prefix}_mcc_mean"] = aggregate(observed["mcc"], np.mean)
        result[f"{prefix}_auc_macro_mean"] = aggregate(observed["auc_macro"], np.mean)
        result[f"{prefix}_n_folds"] = observed["n_folds"]
        result[f"{prefix}_n_folds_multiclass"] = observed["n_folds_multiclass"]
        result[f"{prefix}_folds"] = observed["folds"]
        for level, values in observed["per_class"].items():
            result[f"{prefix}_f1_per_class_{level}"] = float(np.mean(values))

        perm_f1: list[float] = []
        perm_auc: list[float] = []
        for _ in range(max(0, n_perm)):
            shuffled = pd.Series(
                rng.permutation(labels.values), index=labels.index, dtype=object
            )
            outcome = _run_lobo(folds, shuffled, collect_folds=False)
            if outcome["f1_macro"]:
                perm_f1.append(float(np.mean(outcome["f1_macro"])))
            if outcome["auc_macro"]:
                perm_auc.append(float(np.mean(outcome["auc_macro"])))

        result[f"{prefix}_f1_macro_perm_mean"] = aggregate(perm_f1, np.mean)
        result[f"{prefix}_f1_macro_perm_std"] = aggregate(perm_f1, np.std)
        result[f"{prefix}_auc_macro_perm_mean"] = aggregate(perm_auc, np.mean)

        observed_mean = result[f"{prefix}_f1_macro_mean"]
        if perm_f1 and not np.isnan(observed_mean):
            result[f"{prefix}_f1_macro_delta"] = float(observed_mean - np.mean(perm_f1))
            # +1 to both: the observed value is one draw from the null, so the p-value
            # can never be exactly zero.
            n_ge = sum(1 for value in perm_f1 if value >= observed_mean)
            result[f"{prefix}_f1_perm_pvalue"] = float((n_ge + 1) / (len(perm_f1) + 1))
        else:
            _nan_keys(result, [f"{prefix}_f1_macro_delta", f"{prefix}_f1_perm_pvalue"])
    return result
