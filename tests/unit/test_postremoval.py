"""Post-removal: outlier identification, the n_batches knob, and the loud failure."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from combobatch.postremoval import (
    PostRemovalError,
    apply_post_removal,
    identify_outlier_batches,
)
from tests.conftest import make_expression


@pytest.fixture
def planted_outlier():
    """A dataset where one batch is shifted far away from the rest."""
    exp, ann = make_expression(n_samples=90, n_genes=40, n_batches=3)
    outlier_mask = (ann["batch"] == "Batch2").values
    exp = exp.copy()
    exp.loc[outlier_mask] = exp.loc[outlier_mask] + 500.0
    return exp, ann


class TestIdentify:
    def test_finds_the_planted_outlier(self, planted_outlier):
        exp, ann = planted_outlier
        found = identify_outlier_batches(exp, ann, batch_col="batch", min_batch_size=5)
        assert found == ["Batch2"]

    def test_n_batches_controls_how_many_are_returned(self, planted_outlier):
        exp, ann = planted_outlier
        assert (
            len(
                identify_outlier_batches(
                    exp, ann, batch_col="batch", n_batches=2, min_batch_size=5
                )
            )
            == 2
        )

    def test_n_batches_cannot_exceed_the_candidates(self, planted_outlier):
        exp, ann = planted_outlier
        found = identify_outlier_batches(
            exp, ann, batch_col="batch", n_batches=99, min_batch_size=5
        )
        assert len(found) == 3

    def test_min_batch_size_excludes_small_batches(self, planted_outlier):
        exp, ann = planted_outlier
        # Every batch has 30 samples, so a threshold above that leaves no candidate.
        assert (
            identify_outlier_batches(exp, ann, batch_col="batch", min_batch_size=1000)
            == []
        )

    def test_zero_n_batches_raises(self, planted_outlier):
        exp, ann = planted_outlier
        with pytest.raises(PostRemovalError, match="at least 1"):
            identify_outlier_batches(exp, ann, batch_col="batch", n_batches=0)

    def test_disjoint_indices_raise(self):
        exp, ann = make_expression(n_samples=10, n_genes=5)
        other = ann.copy()
        other.index = [f"X{i}" for i in range(len(other))]
        with pytest.raises(PostRemovalError, match="share no samples"):
            identify_outlier_batches(exp, other, batch_col="batch")


class TestApply:
    def test_removes_the_outlier_batch(self, planted_outlier):
        exp, ann = planted_outlier
        out_exp, out_ann, result = apply_post_removal(
            exp, ann, batch_col="batch", min_batch_size=5
        )
        assert result.removed_batches == ("Batch2",)
        assert "Batch2" not in set(out_ann["batch"])
        assert len(out_exp) == len(out_ann) == result.n_samples_remaining

    def test_removing_two_batches(self, planted_outlier):
        exp, ann = planted_outlier
        _, out_ann, result = apply_post_removal(
            exp, ann, batch_col="batch", n_batches=2, min_batch_size=5
        )
        assert len(result.removed_batches) == 2
        assert out_ann["batch"].nunique() == 1

    def test_no_candidates_removes_nothing_and_says_so(self, planted_outlier):
        exp, ann = planted_outlier
        out_exp, _, result = apply_post_removal(
            exp, ann, batch_col="batch", min_batch_size=1000
        )
        assert result.removed_batches == ()
        assert result.n_samples_removed == 0
        assert len(out_exp) == len(exp)

    def test_removing_everything_raises(self, planted_outlier):
        exp, ann = planted_outlier
        with pytest.raises(PostRemovalError, match="would discard every sample"):
            apply_post_removal(
                exp, ann, batch_col="batch", n_batches=3, min_batch_size=5
            )

    def test_expression_and_annotation_stay_aligned(self, planted_outlier):
        exp, ann = planted_outlier
        out_exp, out_ann, _ = apply_post_removal(
            exp, ann, batch_col="batch", min_batch_size=5
        )
        assert list(out_exp.index) == list(out_ann.index)


class TestFailsLoudly:
    """The donor logged the exception and carried on, writing a post1 output identical
    to post0 which was then scored as though post-removal had been applied."""

    def test_identification_failure_is_not_swallowed(
        self, planted_outlier, monkeypatch
    ):
        exp, ann = planted_outlier

        def boom(*args, **kwargs):
            raise ValueError("PCA exploded")

        monkeypatch.setattr("sklearn.decomposition.PCA.fit_transform", boom)
        with pytest.raises(PostRemovalError, match="identification failed"):
            apply_post_removal(exp, ann, batch_col="batch", min_batch_size=5)

    def test_result_never_silently_equals_the_input(self, planted_outlier):
        exp, ann = planted_outlier
        out_exp, _, result = apply_post_removal(
            exp, ann, batch_col="batch", min_batch_size=5
        )
        # Either something was removed, or the result says plainly that nothing was.
        assert (len(out_exp) < len(exp)) == (result.n_samples_removed > 0)


class TestResultPayload:
    def test_as_dict_is_sidecar_ready(self, planted_outlier):
        exp, ann = planted_outlier
        _, _, result = apply_post_removal(exp, ann, batch_col="batch", min_batch_size=5)
        payload = result.as_dict()
        assert payload["post_removal_batches"] == ["Batch2"]
        assert payload["post_removal_n_samples_removed"] == 30
