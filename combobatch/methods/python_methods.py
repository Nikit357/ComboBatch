"""The 13 harmonization methods implemented without R.

Ported from the donor benchmark with three systematic changes:

* ``target_group`` no longer defaults to one dataset's batch name; callers pass a
  reference resolved by :func:`combobatch.methods.base.resolve_reference_batch`.
* The NaN guard is explicit and its gene count is returned rather than discarded.
* Optional dependencies raise a ``NotImplementedError`` naming the install command,
  so a missing package is an honest skip rather than an import crash.

Every function takes ``exp_df`` as samples x genes and returns the same orientation.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from combobatch.methods.base import drop_na_genes

# Broadly-used housekeeping genes, the default negative controls for RUVg. Overridable,
# because a non-human dataset needs different ones.
DEFAULT_CONTROL_GENES = (
    "ACTB",
    "GAPDH",
    "B2M",
    "HPRT1",
    "RPL13A",
    "SDHA",
    "UBC",
    "YWHAZ",
    "HMBS",
    "TBP",
)


def _require(
    package: str, method_name: str, extra: str, hint: str | None = None
) -> None:
    """Raise a NotImplementedError naming the install command for a missing package.

    Parameters
    ----------
    hint:
        Replaces the default `pip install 'combobatch[<extra>]'` command. Needed where
        an extra cannot carry the package — naming a command that does not install it
        is worse than naming none.
    """
    import importlib.util

    if importlib.util.find_spec(package) is None:
        command = hint or f"pip install 'combobatch[{extra}]'"
        raise NotImplementedError(
            f"{method_name} needs the Python package {package!r}, which is not "
            f"installed. Install it with `{command}`."
        )


def normalize_raw(exp_df: pd.DataFrame, ann_df=None, **kw) -> pd.DataFrame:
    """Return the input unchanged. The baseline every metric is read against."""
    return exp_df.copy()


def normalize_median_scaling(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, **kw
) -> pd.DataFrame:
    """
    Shift each batch's per-gene median onto the global per-gene median.

    Returns
    -------
    Median-centred expression matrix.
    """
    groups = ann_df.loc[exp_df.index, batch_col].fillna("Unknown")
    out = exp_df.copy()
    global_median = exp_df.median(axis=0)
    for group in groups.unique():
        idx = exp_df.index[groups == group]
        shift = global_median - exp_df.loc[idx].median(axis=0)
        out.loc[idx] = exp_df.loc[idx].add(shift, axis=1)
    return out


def normalize_pycombat(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, **kw
) -> pd.DataFrame:
    """
    ComBat empirical-Bayes correction, Python port (Behdenna et al. 2023).

    Returns
    -------
    ComBat-corrected expression matrix.
    """
    _require("combat", "07_pycombat", "methods")
    from combat.pycombat import pycombat

    batches = ann_df.loc[exp_df.index, batch_col].fillna("Unknown").tolist()
    exp_in, _ = drop_na_genes(exp_df)
    return pycombat(exp_in.T, batches).T


def normalize_inmoose_combat_seq(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, bio_col: str, **kw
) -> pd.DataFrame:
    """
    ComBat-seq via InMoose, on integer-rounded counts.

    Returns
    -------
    ComBat-seq-corrected expression matrix.
    """
    _require("inmoose", "08_inmoose_combatseq", "methods")
    from inmoose.pycombat import pycombat_seq

    batches = ann_df.loc[exp_df.index, batch_col].fillna("Unknown").tolist()
    bio = ann_df.loc[exp_df.index, bio_col].fillna("Unknown").astype(str)
    clean = np.nan_to_num(exp_df.values, nan=0.0, posinf=0.0, neginf=0.0)
    counts = np.round(clean).astype(int).clip(0)

    covar_mod = None
    if bio.nunique() >= 2:
        covar_mod = pd.get_dummies(bio, drop_first=True).astype(float).values

    result = pycombat_seq(counts.T, batch=batches, covar_mod=covar_mod)
    return pd.DataFrame(result.T, index=exp_df.index, columns=exp_df.columns)


def normalize_harmony(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    n_pcs: int = 50,
    max_iter: int = 20,
    **kw,
) -> pd.DataFrame:
    """
    Harmony (Korsunsky et al. 2019): correct in PCA space, project back.

    Returns
    -------
    Harmony-corrected expression matrix.
    """
    _require("harmonypy", "11_harmony", "methods")
    import harmonypy
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    scaled = scaler.fit_transform(exp_df.fillna(0).values)
    pca = PCA(n_components=min(n_pcs, exp_df.shape[1] - 1, exp_df.shape[0] - 1))
    coords = pca.fit_transform(scaled)

    meta = ann_df.loc[exp_df.index, [batch_col]].fillna("Unknown")
    harmony = harmonypy.run_harmony(
        coords, meta, vars_use=batch_col, max_iter_harmony=max_iter
    )

    # harmonypy >= 0.0.9 returns (n_samples, n_pcs); older versions transposed it.
    corrected = harmony.Z_corr
    if corrected.shape[0] != coords.shape[0]:
        corrected = corrected.T

    recovered = scaler.inverse_transform(pca.inverse_transform(corrected))
    return pd.DataFrame(recovered, index=exp_df.index, columns=exp_df.columns)


def normalize_scanorama(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    dimred: int = 50,
    **kw,
) -> pd.DataFrame:
    """
    Scanorama (Hie et al. 2019): panoramic stitching, inverse-projected to gene space.

    Returns
    -------
    Scanorama-corrected expression matrix.
    """
    _require("scanorama", "12_scanorama", "methods")
    import scanorama
    from sklearn.decomposition import PCA

    batches = ann_df.loc[exp_df.index, batch_col].fillna("Unknown")
    exp_in, _ = drop_na_genes(exp_df)

    unique_batches = list(batches.unique())
    datasets = [exp_in.loc[batches == b].values for b in unique_batches]
    genes = [list(exp_in.columns)] * len(datasets)
    integrated, _ = scanorama.integrate(datasets, genes, dimred=dimred)

    coords = np.zeros((len(exp_in), integrated[0].shape[1]))
    for position, batch in enumerate(unique_batches):
        coords[(batches == batch).values] = integrated[position]

    pca = PCA(n_components=integrated[0].shape[1]).fit(exp_in.values)
    recovered = coords @ pca.components_ + pca.mean_
    return pd.DataFrame(recovered, index=exp_in.index, columns=exp_in.columns)


def normalize_fsmvn(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    target_group: str | None = None,
    **kw,
) -> pd.DataFrame:
    """
    Feature-specific mean-variance normalization onto a reference batch.

    Returns
    -------
    FSMVN-corrected expression matrix.
    """
    groups = ann_df.loc[exp_df.index, batch_col].astype(str)
    reference = (
        exp_df.loc[groups == target_group]
        if target_group is not None and target_group in groups.values
        else exp_df
    )
    target_mean = reference.mean(axis=0)
    target_std = reference.std(axis=0).replace(0, 1)

    out = exp_df.copy()
    for group in groups.unique():
        idx = exp_df.index[groups == group]
        block = exp_df.loc[idx]
        block_std = block.std(axis=0).replace(0, 1)
        out.loc[idx] = (
            block - block.mean(axis=0)
        ) / block_std * target_std + target_mean
    return out


def normalize_fsqn_py(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    target_group: str | None = None,
    **kw,
) -> pd.DataFrame:
    """
    Feature-specific quantile normalization, Python implementation.

    Returns
    -------
    FSQN-normalized expression matrix.
    """
    groups = ann_df.loc[exp_df.index, batch_col].astype(str)
    use_reference = target_group is not None and target_group in groups.values
    source = exp_df.loc[groups == target_group] if use_reference else exp_df
    target_dist = np.sort(source.values, axis=1).mean(axis=0)

    out = exp_df.copy()
    for sample in exp_df.index:
        order = np.argsort(exp_df.loc[sample].values)
        out.loc[sample, out.columns[order]] = target_dist
    return out


def normalize_quantile(exp_df: pd.DataFrame, ann_df=None, **kw) -> pd.DataFrame:
    """
    Standard quantile normalization across samples (Bolstad et al. 2003).

    Returns
    -------
    Quantile-normalized expression matrix.
    """
    values = exp_df.values.astype(float)
    ranks = np.argsort(np.argsort(values, axis=1), axis=1)
    target = np.sort(values, axis=1).mean(axis=0)
    return pd.DataFrame(target[ranks], index=exp_df.index, columns=exp_df.columns)


def normalize_rank(exp_df: pd.DataFrame, ann_df=None, **kw) -> pd.DataFrame:
    """
    Replace each sample's values with fractional ranks in (0, 1).

    Returns
    -------
    Rank-normalized expression matrix, platform-independent by construction.
    """
    ranked = exp_df.rank(axis=1, method="average", na_option="keep")
    return ranked / (exp_df.notna().sum(axis=1).values[:, None] + 1)


def _platform_variance(exp_df: pd.DataFrame, batch_series: pd.Series) -> pd.Series:
    """Per-gene one-way ANOVA R^2 against batch, in [0, 1]."""
    groups = batch_series.reindex(exp_df.index)
    values = exp_df.values.astype(float)

    grand_mean = np.nanmean(values, axis=0)
    ss_total = np.nansum((values - grand_mean) ** 2, axis=0)

    ss_between = np.zeros(values.shape[1], dtype=float)
    for group in groups.unique():
        mask = (groups == group).values
        if mask.sum() < 2:
            continue
        block = values[mask, :]
        n_present = np.sum(~np.isnan(block), axis=0).astype(float)
        ss_between += n_present * (np.nanmean(block, axis=0) - grand_mean) ** 2

    return pd.Series(
        np.where(ss_total > 0, ss_between / ss_total, 0.0), index=exp_df.columns
    )


def normalize_angel(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    threshold: float = 0.20,
    **kw,
) -> pd.DataFrame:
    """
    Rank-transform, then discard genes whose variance is mostly explained by batch.

    Unlike every other method here this one *reduces the gene space* rather than
    correcting values, so its metrics are not directly comparable with the others.

    Returns
    -------
    Rank-normalized expression restricted to platform-stable genes.
    """
    exp_in, _ = drop_na_genes(exp_df)
    ranked = exp_in.rank(axis=1, method="average", na_option="keep")
    rank_norm = ranked / (exp_in.notna().sum(axis=1).values[:, None] + 1)

    batch_series = ann_df.loc[exp_in.index, batch_col].astype(str)
    gene_r2 = _platform_variance(rank_norm, batch_series)
    stable = gene_r2[gene_r2 < threshold].index

    if len(stable) == 0:
        raise RuntimeError(
            f"25_angel: no genes passed the platform-variance filter at "
            f"threshold={threshold}; all {len(exp_in.columns)} genes had batch "
            f"R^2 >= {threshold}. Raise the threshold to retain more genes."
        )
    return rank_norm[stable]


def normalize_xpn(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    target_group: str | None = None,
    n_quantiles: int = 50,
    **kw,
) -> pd.DataFrame:
    """
    XPN: per-gene piecewise-linear quantile matching onto a reference batch.

    Returns
    -------
    XPN-normalized expression matrix.
    """
    groups = ann_df.loc[exp_df.index, batch_col].astype(str)
    use_reference = target_group is not None and target_group in groups.values
    ref_mask = (
        groups == target_group if use_reference else pd.Series(True, index=exp_df.index)
    )

    levels = np.linspace(0, 100, n_quantiles + 2)
    ref_quantiles = np.percentile(
        exp_df.loc[ref_mask].values.astype(float), levels, axis=0
    )

    out = exp_df.copy().astype(float)
    for batch in groups.unique():
        if use_reference and batch == target_group:
            continue
        idx = exp_df.index[groups == batch]
        block = exp_df.loc[idx].values.astype(float)
        src_quantiles = np.percentile(block, levels, axis=0)

        transformed = np.empty_like(block)
        for gene in range(block.shape[1]):
            transformed[:, gene] = np.interp(
                block[:, gene], src_quantiles[:, gene], ref_quantiles[:, gene]
            )
        out.loc[idx] = transformed
    return out


def normalize_recombat(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame, *, batch_col: str, **kw
) -> pd.DataFrame:
    """
    reComBat (Adossa et al. 2021): ComBat with ridge-regularized batch estimation.

    Returns
    -------
    reComBat-corrected expression matrix.
    """
    # Not `combobatch[methods]`: reComBat declares the deprecated `sklearn` stub and
    # `python <3.11`, so it cannot be an extra at all. See docker/Dockerfile.
    _require(
        "reComBat",
        "30_recombat",
        "methods",
        hint=(
            "pip install --no-deps --ignore-requires-python 'reComBat @ "
            "git+https://github.com/BorgwardtLab/reComBat"
            "@b02bd025c9671e75f37dc513d7f466b05448364d'"
        ),
    )
    from reComBat import reComBat

    exp_in, _ = drop_na_genes(exp_df)
    batch = ann_df.loc[exp_in.index, batch_col]
    corrected = reComBat().fit_transform(exp_in, batch)
    return corrected.reindex(exp_df.index)
