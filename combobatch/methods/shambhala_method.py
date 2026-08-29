"""``20_shambhala`` — the one glue function the source never had.

Upstream exposes ``harmonize_parallel`` as its reusable unit, but the ~40 lines that turn
a user matrix into its arguments (gene intersection, NA strategy, P alignment, Q
statistics, transposition) exist only as **four copies that have already drifted apart** —
in ``run_shambhala.py``, the job runner, the benchmark shim and the test conftest. One of
those copies is where §6 bug 4 lives: it forgets to forward ``skip_qn``, so the approved
speed-up quantile-normalizes twice. This module is the single copy.

Two contracts are load-bearing and easy to violate:

* **Input must be raw, non-log.** ``readExpressionData.m`` applies ``log2(x + 1)`` itself.
  Output is ``exp(rm + rs·log(x + 1))``, also raw, and can reach ~100 000.
* **``harmonize_parallel`` is genes × samples**, unlike every other Shambhala module and
  every other ComboBatch method. The transposes live here and nowhere else.
"""

from __future__ import annotations

import dataclasses
import functools
import gzip
import io
import os
from typing import Any

import numpy as np
import pandas as pd

from combobatch.logging_utils import get_logger
from combobatch.methods.python_methods import _require
from combobatch.vendor.shambhala import OCTAVE_DIR

# Shipped calibration references, so `--method 20_shambhala` works with no extra
# arguments. Users point `P` / `Q` at any other pair — including the 18 benchmark
# variants — without those being baked into the tool.
DEFAULT_P = "P0_standard.csv.gz"
DEFAULT_Q = "Q0_standard.csv.gz"

# A raw expression matrix of any size reaches well past this; a log2 one tops out near 20.
# Shambhala has no scale detection at all upstream, and log-scale input is double-logged
# into plausible-looking but wrong numbers — the exact failure mode this package refuses
# to reproduce quietly. A warning, not an error: it is a heuristic over user data.
_LOG_SCALE_MAX = 30.0


@functools.lru_cache(maxsize=2)
def load_calibration(name: str) -> pd.DataFrame:
    """
    Load one gzipped calibration matrix shipped inside the package.

    Parameters
    ----------
    name
        File name within ``combobatch/data/calibration/``.

    Returns
    -------
    Samples × genes, raw linear scale. Cached, so treat the result as read-only.
    """
    from importlib.resources import files

    resource = files("combobatch.data").joinpath("calibration").joinpath(name)
    raw = gzip.decompress(resource.read_bytes())
    return pd.read_csv(io.BytesIO(raw), index_col=0, low_memory=False).astype(float)


def warn_nested_parallelism(outer_workers: int, inner_workers: int) -> str | None:
    """
    Warn when dispatcher workers times Shambhala workers oversubscribes the CPUs.

    Octave process count is the *product* of the two, not the larger of them, which is why
    the inner default is 1. Returns the warning text (also logged) or ``None``.
    """
    total = outer_workers * inner_workers
    available = os.cpu_count() or 1
    if total <= available:
        return None
    message = (
        f"Shambhala nested parallelism: {outer_workers} dispatcher workers x "
        f"{inner_workers} Shambhala workers = {total} Octave processes on "
        f"{available} CPUs. Lower one of them; the counts multiply."
    )
    get_logger().warning(message)
    return message


def _warn_if_log_scale(exp_df: pd.DataFrame) -> None:
    """Warn when the input looks log-transformed, which Shambhala cannot detect."""
    observed_max = float(np.nanmax(exp_df.values)) if exp_df.size else 0.0
    if observed_max < _LOG_SCALE_MAX:
        get_logger().warning(
            "Shambhala input has maximum %.3g, which suggests log-scale data. Octave "
            "applies log2(x + 1) itself, so log-scale input is silently double-logged. "
            "Pass raw (non-log) expression.",
            observed_max,
        )


def shambhala_harmonize(
    exp_df: pd.DataFrame,
    p_df: pd.DataFrame | None = None,
    q_df: pd.DataFrame | None = None,
    *,
    k: int = 5,
    n_workers: int = 1,
    na_strategy: str = "drop",
    knn_k: int = 5,
    max_na_frac: float = 0.20,
    q_pseudocount: float = 1e-6,
    random_seed: int | None = None,
    octave_bin: str = "octave",
    octave_scripts_dir: str = OCTAVE_DIR,
    timeout_s: int = 6000,
    disable_progress: bool = True,
    precompute_qn_reference: bool = False,
    synthetic_cublock_p: bool = False,
    max_p_samples: int | None = None,
    precompute_cublock_clusters: bool = False,
    python_cublock: bool = False,
) -> pd.DataFrame:
    """
    Harmonize an expression matrix with Shambhala2 (QN + CuBlock + Q-rescaling).

    Each sample is quantile-normalized jointly with the P calibration reference,
    CuBlock-normalized, and finally rescaled to the per-gene log mean and standard
    deviation of the Q reference. Samples are processed independently, so no batch
    labels are involved.

    Parameters
    ----------
    exp_df
        Samples × genes, **raw (non-log) scale**. May contain NaN.
    p_df, q_df
        Samples × genes calibration references. ``None`` loads the shipped P0 / Q0.
    k
        k-means gene clusters for CuBlock.
    n_workers
        Shambhala worker processes. Defaults to 1 because the dispatcher already runs
        jobs in parallel and the two counts multiply into Octave processes.
    na_strategy, knn_k, max_na_frac
        NaN handling before Octave, which requires NaN-free input. ``drop`` removes
        affected genes and restores them as NaN afterwards; ``knn`` imputes and keeps
        the imputed values.
    q_pseudocount
        Added to Q before the log, so zero counts do not become ``-Inf``.
    random_seed
        Seeds Octave's k-means for reproducibility.
    octave_bin, octave_scripts_dir, timeout_s
        Octave invocation details.
    disable_progress
        Skip the progress renderer and its ``multiprocessing.Manager``.
    precompute_qn_reference
        Apply quantile normalization Python-side against a reference distribution taken
        from P, and tell Octave to skip its own. The forwarding upstream forgets is the
        whole point of this parameter existing here.
    synthetic_cublock_p
        Pass a single synthetic P column to CuBlock instead of all P samples.
    max_p_samples
        Subsample P to its N centroid-closest samples.
    precompute_cublock_clusters
        Cluster genes once from P and reuse the assignment for every sample.
    python_cublock
        Experimental pure-Python CuBlock. Upstream reports ~226 % mean relative
        difference against Octave; not equivalent, and not recommended.

    Returns
    -------
    Samples × genes, raw scale, in the input sample order. The gene set is the
    intersection of the input with P and Q — genes absent from either reference cannot
    be harmonized and are dropped, with the count logged.

    Raises
    ------
    ValueError
        If a calibration reference contains NaN, or the three matrices share no genes.
    """
    from combobatch.vendor.shambhala import na_handling, q_rescale
    from combobatch.vendor.shambhala import parallel as parallel_module

    log = get_logger()

    if p_df is None:
        p_df = load_calibration(DEFAULT_P)
    if q_df is None:
        q_df = load_calibration(DEFAULT_Q)
    p_df = p_df.astype(float)
    q_df = q_df.astype(float)

    # Upstream exits the process here; a library raises instead.
    for label, frame in (("P", p_df), ("Q", q_df)):
        if frame.isna().any().any():
            n_bad = int(frame.isna().any(axis=0).sum())
            raise ValueError(
                f"{label} calibration reference has NaN in {n_bad} gene column(s); "
                f"Shambhala requires NaN-free references."
            )

    if python_cublock:
        # Checked here rather than inside the worker: a missing package must not surface
        # as a pickled exception from a subprocess after the pool has already spun up.
        _require("qnorm", "20_shambhala with python_cublock=True", "shambhala")
        log.warning(
            "20_shambhala python_cublock=True is experimental: upstream measures ~226%% "
            "mean relative difference against Octave. Results are not comparable."
        )

    _warn_if_log_scale(exp_df)

    # ── genes × samples from here to the final transpose ──────────────────────────
    input_T = exp_df.T
    p_T = p_df.T
    q_T = q_df.T

    common_genes = input_T.index.intersection(p_T.index).intersection(q_T.index)
    if len(common_genes) == 0:
        raise ValueError(
            "Shambhala found no genes shared by the input, P and Q. All three must use "
            "the same gene identifiers (HGNC symbols in the shipped references)."
        )
    log.info(
        "Shambhala gene intersection: %d genes (input had %d, P %d, Q %d).",
        len(common_genes),
        input_T.shape[0],
        p_T.shape[0],
        q_T.shape[0],
    )

    input_T = input_T.loc[common_genes]
    p_T = p_T.loc[common_genes]
    q_T = q_T.loc[common_genes]

    if max_p_samples is not None and p_T.shape[1] > max_p_samples:
        p_samples = p_T.T
        distances = ((p_samples - p_samples.mean(axis=0)) ** 2).sum(axis=1)
        p_T = p_T[distances.nsmallest(max_p_samples).index]
        log.info(
            "Shambhala subsampled P to %d centroid-closest samples.", max_p_samples
        )

    clean_df, na_mask = na_handling.apply_na_strategy(
        input_T.T, strategy=na_strategy, knn_k=knn_k, max_na_frac=max_na_frac
    )
    clean_T = clean_df.T
    p_T = p_T.loc[clean_T.index]
    log.info(
        "Shambhala NA strategy %r: %d gene(s) dropped.",
        na_mask.strategy,
        len(na_mask.dropped_genes),
    )

    p_for_octave = p_T
    if precompute_qn_reference or synthetic_cublock_p:
        p_qn_reference = np.sort(p_T.values, axis=0).mean(axis=1)

    if precompute_qn_reference:
        # Only this branch needs qnorm, so it is not a declared backend of the method —
        # a plain Shambhala run must not be reported unavailable for want of it.
        _require("qnorm", "20_shambhala with precompute_qn_reference=True", "shambhala")
        import qnorm

        clean_T = qnorm.quantile_normalize(clean_T, axis=1, target=p_qn_reference)
        log.info("Shambhala applied Python-side QN against the P reference.")

    if synthetic_cublock_p:
        p_for_octave = pd.DataFrame(
            p_qn_reference[:, np.newaxis], index=p_T.index, columns=["P_qn_reference"]
        )

    fixed_clusters = None
    if precompute_cublock_clusters:
        from combobatch.vendor.shambhala.octave_bridge import precompute_clusters_from_p

        fixed_clusters = precompute_clusters_from_p(
            p_T=p_for_octave, k=k, random_seed=random_seed
        )

    rm, rs = q_rescale.compute_q_statistics(
        q_T.loc[clean_T.index].T, q_pseudocount=q_pseudocount
    )

    # With q_pseudocount=0 the Q statistics exclude zero-count genes, and the rescale
    # step then silently drops them from the result — which upstream's gene restoration
    # cannot reconcile (it indexes by the *original* gene list and raises a KeyError).
    # Drop them here instead, so they come back as NaN like any other unusable gene.
    unusable_genes = clean_T.index.difference(rm.index)
    if len(unusable_genes) > 0:
        log.warning(
            "Shambhala: %d gene(s) have no Q statistics (q_pseudocount=%g) and cannot "
            "be rescaled; they are returned as NaN.",
            len(unusable_genes),
            q_pseudocount,
        )
        clean_T = clean_T.drop(index=unusable_genes)
        p_for_octave = p_for_octave.loc[clean_T.index]
        na_mask = dataclasses.replace(
            na_mask,
            dropped_genes=na_mask.dropped_genes + list(unusable_genes),
        )

    harmonized_T = parallel_module.harmonize_parallel(
        input_df=clean_T,
        p_df=p_for_octave,
        rm=rm,
        rs=rs,
        k=k,
        n_workers=n_workers,
        octave_scripts_dir=octave_scripts_dir,
        octave_bin=octave_bin,
        timeout_s=timeout_s,
        random_seed=random_seed,
        disable_progress=disable_progress,
        fixed_clusters=fixed_clusters,
        python_cublock=python_cublock,
        # §6 bug 4: the job-runner copy of this glue omits this forward, so Python QN
        # and Octave quantilenorm both run.
        skip_qn=precompute_qn_reference,
    )

    result = na_handling.restore_na_genes(harmonized_T.T, na_mask)
    return result.reindex(index=exp_df.index)


def normalize_shambhala(
    exp_df: pd.DataFrame,
    ann_df: pd.DataFrame,
    *,
    batch_col: str,
    bio_col: str | None = None,
    P: str | None = None,
    Q: str | None = None,
    **params: Any,
) -> pd.DataFrame:
    """
    Method-registry adapter for Shambhala.

    ``batch_col`` and ``bio_col`` are accepted for the common method signature but unused:
    Shambhala harmonizes each sample independently against fixed calibration references,
    so it needs no batch labels and no biology labels.

    Parameters
    ----------
    exp_df
        Samples × genes, raw (non-log) scale.
    ann_df
        Annotation, unused.
    batch_col, bio_col
        Unused; see above.
    P, Q
        Paths or ``s3://`` URIs to calibration matrices. ``None`` uses the shipped pair.
    **params
        Forwarded to :func:`shambhala_harmonize`.

    Returns
    -------
    The harmonized matrix, samples × genes.
    """
    from combobatch import dataio

    p_df = dataio.read_table(P) if P else None
    q_df = dataio.read_table(Q) if Q else None
    return shambhala_harmonize(exp_df, p_df, q_df, **params)
