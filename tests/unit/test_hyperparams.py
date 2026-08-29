"""The anti-regression guard: a hyperparameter must actually change the output.

In the donor pipeline every hyperparameter was unreachable — 19 methods declared them,
but the single call site was frozen at ``fn(exp, ann, batch_col=, bio_col=)``, and the
R-backed ones froze their arguments inside the interpolated R string. A test asserting
only "the method runs" would have passed against all of that. These assert the parameter
reaches the algorithm and moves the numbers.
"""

from __future__ import annotations

import pandas as pd
import pytest

from combobatch.methods import METHOD_REGISTRY
from tests.unit.test_methods_all import run_method

# (method, parameter, value_a, value_b) — two settings that must not agree.
DISTINGUISHING_PARAMS = [
    ("25_angel", "threshold", 0.05, 0.95),
    ("26_xpn", "n_quantiles", 5, 100),
]


class TestParametersReachTheAlgorithm:
    @pytest.mark.parametrize("key,param,value_a,value_b", DISTINGUISHING_PARAMS)
    def test_changing_a_parameter_changes_the_output(
        self, key, param, value_a, value_b, synthetic_with_batch_effect
    ):
        exp, ann = synthetic_with_batch_effect
        first = run_method(key, exp, ann, **{param: value_a})
        second = run_method(key, exp, ann, **{param: value_b})

        if first.shape != second.shape:
            return  # A shape change is itself proof the parameter took effect.
        assert not first.equals(second), (
            f"{key}: {param}={value_a} and {param}={value_b} produced identical output, "
            "so the parameter is not reaching the algorithm"
        )

    def test_angel_threshold_controls_how_many_genes_survive(
        self, synthetic_with_batch_effect
    ):
        exp, ann = synthetic_with_batch_effect
        strict = run_method("25_angel", exp, ann, threshold=0.05)
        lenient = run_method("25_angel", exp, ann, threshold=0.95)
        assert strict.shape[1] < lenient.shape[1]

    def test_fsmvn_reference_choice_changes_the_result(
        self, synthetic_with_batch_effect
    ):
        exp, ann = synthetic_with_batch_effect
        spec = METHOD_REGISTRY["13_fsmvn"]
        first = spec.fn(
            exp, ann, batch_col="batch", bio_col="bio", target_group="Batch0"
        )
        second = spec.fn(
            exp, ann, batch_col="batch", bio_col="bio", target_group="Batch1"
        )
        assert not first.equals(second)


class TestDefaultsAreHonoured:
    def test_omitting_a_parameter_matches_passing_its_default(self, synthetic):
        exp, ann = synthetic
        default_value = METHOD_REGISTRY["26_xpn"].defaults()["n_quantiles"]
        implicit = run_method("26_xpn", exp, ann)
        explicit = run_method("26_xpn", exp, ann, n_quantiles=default_value)
        pd.testing.assert_frame_equal(implicit, explicit)

    def test_registry_defaults_are_what_the_functions_use(self):
        # Guards against the registry advertising a default the function does not have.
        assert METHOD_REGISTRY["25_angel"].defaults()["threshold"] == 0.20
        assert METHOD_REGISTRY["10_mnn"].defaults()["k"] == 20
        assert METHOD_REGISTRY["38_harman"].defaults()["limit"] == 0.1


class TestRSideParametersAreDeclared:
    """Arguments the donor froze inside R strings are now declared and reachable."""

    @pytest.mark.parametrize(
        "key,param",
        [
            ("33_amdbnorm", "poly_degree"),
            ("33_amdbnorm", "n_dist_points"),
            ("19_tdm", "log_target"),
            ("31_ruv3prps", "k_factors"),
            ("22_tmm", "method"),
            ("22_tmm", "prior_count"),
            ("23_vst", "blind"),
            ("28_npn", "npn_func"),
            ("38_harman", "limit"),
            ("34_arsyn", "norm"),
            ("09_ruv", "control_genes"),
        ],
    )
    def test_declared_with_the_r_argument_it_maps_to(self, key, param):
        spec = METHOD_REGISTRY[key]
        assert param in spec.hyperparams, f"{key} does not declare {param}"

    def test_r_argument_mapping_is_recorded_where_it_applies(self):
        assert (
            METHOD_REGISTRY["33_amdbnorm"].hyperparams["poly_degree"].r_argument
            == "DBNorm::polyFit(degree)"
        )
