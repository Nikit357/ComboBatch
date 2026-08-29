"""The four imputers, their hyperparameters, and the two donor behaviours removed."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from combobatch.imputation import (
    IMPUTER_REGISTRY,
    ImputerSpec,
    get_imputer,
    impute,
)


def run_impute(method, exp, **kwargs):
    """Impute, skipping when the imputer's backend is unavailable."""
    spec = get_imputer(method)
    missing = spec.missing_backends()
    if missing:
        pytest.skip(f"{method} needs {', '.join(missing)}")
    return impute(exp, method=method, **kwargs)


class TestRegistry:
    def test_all_four_strategies_are_registered(self):
        assert set(IMPUTER_REGISTRY) == {"strict", "knn", "softimpute", "missforest"}

    def test_strict_is_a_first_class_entry(self):
        # In the donor, strict lived in a different function with no parameters, and
        # passing "strict" to the imputing one raised ValueError.
        assert isinstance(IMPUTER_REGISTRY["strict"], ImputerSpec)
        assert "max_na_frac" in IMPUTER_REGISTRY["strict"].hyperparams

    def test_every_imputer_declares_the_na_ceiling(self):
        for key, spec in IMPUTER_REGISTRY.items():
            assert "max_na_frac" in spec.hyperparams, f"{key} lacks max_na_frac"

    def test_r_backed_imputers_declare_their_packages(self):
        assert IMPUTER_REGISTRY["softimpute"].r_packages == ("softImpute",)
        assert IMPUTER_REGISTRY["missforest"].r_packages == ("missForest",)
        assert not IMPUTER_REGISTRY["knn"].requires_r

    def test_unknown_imputer_raises(self):
        with pytest.raises(KeyError, match="Available:"):
            get_imputer("nope")

    def test_unknown_parameter_is_rejected(self):
        with pytest.raises(Exception, match="unknown parameter"):
            IMPUTER_REGISTRY["knn"].validate_params({"kk": 3})


class TestStrict:
    def test_drops_every_gene_with_any_na(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        out, report = run_impute("strict", exp)
        assert out.isna().sum().sum() == 0
        assert out.shape[1] < exp.shape[1]
        assert report.n_genes_dropped_na_frac == exp.shape[1] - out.shape[1]

    def test_leaves_a_complete_matrix_untouched(self, synthetic):
        exp, _ = synthetic
        out, report = run_impute("strict", exp)
        pd.testing.assert_frame_equal(out, exp)
        assert report.n_genes_dropped_na_frac == 0

    def test_is_the_zero_threshold_corner_of_one_rule(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        strict_out, _ = run_impute("strict", exp)
        knn_at_zero, _ = run_impute("knn", exp, max_na_frac=0.0)
        assert list(strict_out.columns) == list(knn_at_zero.columns)


class TestKnn:
    def test_fills_every_missing_value(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        out, report = run_impute("knn", exp, max_na_frac=0.5)
        assert out.isna().sum().sum() == 0
        assert report.n_cells_imputed > 0

    def test_preserves_shape_within_the_ceiling(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        out, _ = run_impute("knn", exp, max_na_frac=1.0)
        assert out.shape == exp.shape

    def test_observed_values_are_left_alone(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        out, _ = run_impute("knn", exp, max_na_frac=1.0)
        observed = exp.notna()
        assert np.allclose(out.values[observed.values], exp.values[observed.values])

    def test_k_reaches_the_imputer(self, synthetic_with_na):
        # The donor declared knn_k but no caller ever passed it.
        exp, _ = synthetic_with_na
        few, _ = run_impute("knn", exp, max_na_frac=1.0, params={"knn_k": 2})
        many, _ = run_impute("knn", exp, max_na_frac=1.0, params={"knn_k": 20})
        assert not few.equals(many)

    def test_weights_reaches_the_imputer(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        uniform, _ = run_impute(
            "knn", exp, max_na_frac=1.0, params={"weights": "uniform"}
        )
        distance, _ = run_impute(
            "knn", exp, max_na_frac=1.0, params={"weights": "distance"}
        )
        assert not uniform.equals(distance)


class TestMaxNaFrac:
    @pytest.mark.parametrize("ceiling,expect_all_genes", [(0.0, False), (1.0, True)])
    def test_ceiling_controls_gene_retention(
        self, synthetic_with_na, ceiling, expect_all_genes
    ):
        exp, _ = synthetic_with_na
        out, _ = run_impute("knn", exp, max_na_frac=ceiling)
        assert (out.shape[1] == exp.shape[1]) is expect_all_genes

    def test_monotonic_in_the_ceiling(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        widths = [
            run_impute("knn", exp, max_na_frac=c)[0].shape[1]
            for c in (0.0, 0.05, 0.5, 1.0)
        ]
        assert widths == sorted(widths)


class TestNoSilentFallback:
    """A run labelled 'knn' must have had KNN applied, or have failed."""

    def test_imputer_failure_propagates(self, synthetic_with_na, monkeypatch):
        exp, _ = synthetic_with_na

        def boom(*args, **kwargs):
            raise RuntimeError("imputer exploded")

        monkeypatch.setitem(
            IMPUTER_REGISTRY,
            "knn",
            ImputerSpec(
                key="knn",
                fn=boom,
                description="test double",
                hyperparams=IMPUTER_REGISTRY["knn"].hyperparams,
            ),
        )
        with pytest.raises(RuntimeError, match="imputer exploded"):
            impute(exp, method="knn", max_na_frac=1.0)


class TestNoOpDetection:
    """The donor returned early when the ceiling had already removed every NA, so a
    'knn' run could be byte-identical to a 'strict' one with nothing recording it."""

    def test_noop_is_recorded_when_nothing_needed_imputing(self, synthetic):
        exp, _ = synthetic
        _, report = run_impute("knn", exp, max_na_frac=1.0)
        assert report.was_noop is True
        assert report.n_cells_imputed == 0

    def test_real_work_is_not_flagged_as_a_noop(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        _, report = run_impute("knn", exp, max_na_frac=1.0)
        assert report.was_noop is False

    def test_strict_is_never_a_noop(self, synthetic):
        exp, _ = synthetic
        _, report = run_impute("strict", exp)
        assert report.was_noop is False


class TestReport:
    def test_counts_are_consistent(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        out, report = run_impute("knn", exp, max_na_frac=0.5)
        assert report.n_genes_in == exp.shape[1]
        assert report.n_genes_out == out.shape[1]
        assert report.n_genes_in - report.n_genes_dropped_na_frac == report.n_genes_out

    def test_residual_na_is_reported(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        out, report = run_impute("knn", exp, max_na_frac=1.0)
        assert report.residual_na == int(out.isna().sum().sum())

    def test_as_dict_is_flat_and_sidecar_ready(self, synthetic_with_na):
        exp, _ = synthetic_with_na
        _, report = run_impute("knn", exp, max_na_frac=0.5)
        payload = report.as_dict()
        assert payload["imputation_method"] == "knn"
        assert all(not isinstance(v, (dict, list)) for v in payload.values())
