"""Shared embeddings: PCA, UMAP and t-SNE.

Computed once per job and handed to every group that needs them — PCA to A, B, G, I
and J, the two 2-D embeddings to C alone. Each group re-deriving its own would dominate
the runtime and, worse, let two groups disagree about the same matrix.

Every magic number the donor buried in these calls is a parameter here.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

# Defaults reproduce the donor's behaviour exactly; they are parameters so a caller can
# depart from it deliberately rather than by editing the source.
DEFAULT_N_COMPONENTS = 50
DEFAULT_UMAP_NEIGHBORS = 30
DEFAULT_UMAP_MIN_DIST = 0.3
DEFAULT_TSNE_PERPLEXITY = 30
RANDOM_STATE = 42


def compute_pca(
    exp_df: pd.DataFrame,
    n_components: int = DEFAULT_N_COMPONENTS,
    *,
    random_state: int = RANDOM_STATE,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Standardize and project an expression matrix onto its principal components.

    Parameters
    ----------
    exp_df
        Expression, samples x genes. All-NaN gene columns are dropped and any
        remaining NaN is zero-filled; a harmonizer that failed produces those, and
        PCA cannot see them.
    n_components
        Upper bound on the component count, clamped to what the matrix supports.
    random_state
        Seed, so a rerun on the same matrix gives the same coordinates.

    Returns
    -------
    ``(coordinates, explained_variance)`` — the variances are absolute, not ratios,
    because group A weights R-squared by them.

    Raises
    ------
    ValueError
        If no usable gene column survives, or the matrix is too small for one component.
    """
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    values = exp_df.values.astype(float)
    usable = ~np.all(np.isnan(values), axis=0)
    values = np.nan_to_num(values[:, usable], nan=0.0)
    if values.shape[1] == 0:
        raise ValueError(
            f"no usable gene columns after dropping all-NaN genes (input was "
            f"{exp_df.shape[0]} x {exp_df.shape[1]})"
        )

    n_comp = min(n_components, values.shape[0] - 1, values.shape[1])
    if n_comp < 1:
        raise ValueError(f"matrix too small for PCA: shape {values.shape}")

    scaled = StandardScaler().fit_transform(values)
    pca = PCA(n_components=n_comp, random_state=random_state)
    return pca.fit_transform(scaled), pca.explained_variance_


def compute_umap(
    pca_coords: np.ndarray,
    *,
    n_neighbors: int = DEFAULT_UMAP_NEIGHBORS,
    min_dist: float = DEFAULT_UMAP_MIN_DIST,
    random_state: int = RANDOM_STATE,
) -> np.ndarray:
    """
    Reduce PCA coordinates to two UMAP dimensions.

    Returns
    -------
    An ``(n_samples, 2)`` array.
    """
    import umap as umap_lib

    reducer = umap_lib.UMAP(
        n_neighbors=min(n_neighbors, max(2, pca_coords.shape[0] - 1)),
        min_dist=min_dist,
        n_components=2,
        metric="euclidean",
        random_state=random_state,
        verbose=False,
    )
    return reducer.fit_transform(pca_coords)


def compute_tsne(
    pca_coords: np.ndarray,
    *,
    perplexity: int = DEFAULT_TSNE_PERPLEXITY,
    random_state: int = RANDOM_STATE,
) -> np.ndarray:
    """
    Reduce PCA coordinates to two t-SNE dimensions.

    Returns
    -------
    An ``(n_samples, 2)`` array.
    """
    from sklearn.manifold import TSNE

    # t-SNE requires perplexity < n_samples; the donor's //4 rule keeps it well clear.
    effective = max(2, min(perplexity, pca_coords.shape[0] // 4))
    tsne = TSNE(
        n_components=2,
        perplexity=effective,
        random_state=random_state,
        method="barnes_hut",
        n_jobs=-1,
    )
    return tsne.fit_transform(pca_coords)
