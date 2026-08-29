"""Shared synthetic fixtures.

Deterministic by construction: one seed, no network, no files. The matrices are small
enough that the whole unit suite stays under a minute.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

SEED = 42
N_SAMPLES = 60
N_GENES = 80
N_BATCHES = 3


def make_expression(
    n_samples: int = N_SAMPLES,
    n_genes: int = N_GENES,
    n_batches: int = N_BATCHES,
    batch_effect: float = 0.0,
    seed: int = SEED,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build a synthetic (expression, annotation) pair.

    Values are exponential plus an offset, so they stay strictly positive — several
    methods round to counts and clip at zero, and an all-zero column makes them
    degenerate for reasons unrelated to what is being tested.

    Parameters
    ----------
    batch_effect
        Size of a per-gene, per-batch additive shift. ``0.0`` gives batches that differ
        only by noise.

    Returns
    -------
    ``(expression, annotation)``, samples x genes and samples x 4 columns.
    """
    rng = np.random.default_rng(seed)
    samples = [f"S{i:04d}" for i in range(n_samples)]
    genes = [f"GENE{j:04d}" for j in range(n_genes)]

    # Batch labels avoid a "GPL" prefix so the RNA-seq-only guard does not fire; the
    # guard has its own dedicated test.
    batches = [f"Batch{i % n_batches}" for i in range(n_samples)]
    values = rng.exponential(5.0, size=(n_samples, n_genes)) + 2.0

    if batch_effect:
        for index in range(n_batches):
            mask = np.array([b == f"Batch{index}" for b in batches])
            shift = rng.normal(0.0, batch_effect, size=n_genes)
            values[mask, :] += shift

    exp = pd.DataFrame(values, index=samples, columns=genes)
    ann = pd.DataFrame(
        {
            "batch": batches,
            "bio": [["TypeA", "TypeB"][i % 2] for i in range(n_samples)],
            "cohort": [f"Cohort{i % 2}" for i in range(n_samples)],
            "binary": [["Yes", "No"][i % 2] for i in range(n_samples)],
        },
        index=samples,
    )
    return exp, ann


@pytest.fixture
def synthetic():
    """A clean (expression, annotation) pair with no batch effect."""
    return make_expression()


@pytest.fixture
def synthetic_with_batch_effect():
    """An (expression, annotation) pair carrying a strong per-batch shift."""
    return make_expression(batch_effect=5.0)


@pytest.fixture
def synthetic_with_na():
    """An expression matrix with a reproducible 10% scattering of missing values."""
    exp, ann = make_expression()
    rng = np.random.default_rng(SEED)
    mask = rng.random(exp.shape) < 0.10
    exp = exp.mask(mask)
    return exp, ann
