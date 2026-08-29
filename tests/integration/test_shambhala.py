"""End-to-end Shambhala, against a real Octave.

Skips itself when Octave is absent, so it is inert on the dev laptop and live in the
image. The matrices are deliberately tiny — Shambhala runs one Octave subprocess per
sample batch, and CuBlock's k-means dominates the wall clock.

The determinism check is seed-based (same seed twice) rather than a committed numeric
baseline: a baseline can only be generated inside the image, where Octave exists, and
pinning one produced on a different Octave build would fail for the wrong reason.
"""

from __future__ import annotations

import shutil

import numpy as np
import pandas as pd
import pytest

from combobatch.methods.shambhala_method import shambhala_harmonize

from tests.conftest import make_expression

pytestmark = [
    pytest.mark.needs_octave,
    pytest.mark.skipif(
        shutil.which("octave") is None, reason="Octave is not installed"
    ),
    pytest.mark.slow,
]

# CuBlock ignores clusters with 100 or fewer genes, so a smaller gene space would come
# back all-NaN for reasons unrelated to what is being tested.
N_GENES = 400
N_SAMPLES = 4


@pytest.fixture(scope="module")
def tiny_inputs():
    exp, _ = make_expression(n_samples=N_SAMPLES, n_genes=N_GENES, seed=1)
    p_df, _ = make_expression(n_samples=6, n_genes=N_GENES, seed=2)
    q_df, _ = make_expression(n_samples=8, n_genes=N_GENES, seed=3)
    # Raw-scale inputs: Octave applies log2(x + 1) itself.
    return exp * 100.0, p_df * 100.0, q_df * 100.0


def test_octave_round_trip_preserves_shape_and_orientation(tiny_inputs):
    exp, p_df, q_df = tiny_inputs

    result = shambhala_harmonize(exp, p_df, q_df, k=2, random_seed=11)

    assert list(result.index) == list(exp.index)
    assert result.shape == (N_SAMPLES, N_GENES)


def test_output_is_raw_scale(tiny_inputs):
    exp, p_df, q_df = tiny_inputs

    result = shambhala_harmonize(exp, p_df, q_df, k=2, random_seed=11)

    # Output is exp(rm + rs * log(x + 1)) — raw linear, not log.
    finite = result.values[np.isfinite(result.values)]
    assert finite.max() > 30.0


def test_same_seed_gives_the_same_matrix(tiny_inputs):
    exp, p_df, q_df = tiny_inputs

    first = shambhala_harmonize(exp, p_df, q_df, k=2, random_seed=11)
    second = shambhala_harmonize(exp, p_df, q_df, k=2, random_seed=11)

    pd.testing.assert_frame_equal(first, second)


def test_precompute_qn_reference_runs_and_differs(tiny_inputs):
    """The A speed-up: Python-side QN plus the pre-QN Octave script."""
    exp, p_df, q_df = tiny_inputs

    plain = shambhala_harmonize(exp, p_df, q_df, k=2, random_seed=11)
    fast = shambhala_harmonize(
        exp, p_df, q_df, k=2, random_seed=11, precompute_qn_reference=True
    )

    assert fast.shape == plain.shape
    assert not np.allclose(
        np.nan_to_num(fast.values), np.nan_to_num(plain.values), rtol=1e-9
    )


def test_approved_speedup_combination_runs(tiny_inputs):
    """A+E together: the combination that previously double-normalized."""
    exp, p_df, q_df = tiny_inputs

    result = shambhala_harmonize(
        exp,
        p_df,
        q_df,
        k=2,
        random_seed=11,
        precompute_qn_reference=True,
        precompute_cublock_clusters=True,
    )

    assert result.shape == (N_SAMPLES, N_GENES)


def test_dropped_na_genes_return_as_nan_columns(tiny_inputs):
    exp, p_df, q_df = tiny_inputs
    exp = exp.copy()
    exp.iloc[0, 5] = np.nan

    result = shambhala_harmonize(exp, p_df, q_df, k=2, random_seed=11)

    assert result.iloc[:, 5].isna().all()
    assert result.shape == (N_SAMPLES, N_GENES)
