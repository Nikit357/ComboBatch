"""The method registry: internal consistency, and agreement with reality.

The donor kept three hand-maintained copies of its method list. These tests exist so
that the one copy here cannot drift from the functions it points at.
"""

from __future__ import annotations

import inspect

import pytest

from combobatch.methods import (
    EXCLUDED_METHODS,
    METHOD_REGISTRY,
    MethodSpec,
    get_method,
)
from combobatch.methods.base import HARSHNESS_TIERS

# 39 in the published benchmark, minus the five that cannot or must not run.
EXPECTED_METHOD_COUNT = 34


class TestRegistryShape:
    def test_expected_number_of_methods(self):
        assert len(METHOD_REGISTRY) == EXPECTED_METHOD_COUNT

    def test_excluded_and_registered_account_for_all_39(self):
        assert len(METHOD_REGISTRY) + len(EXCLUDED_METHODS) == 39

    def test_keys_match_their_specs(self):
        for key, spec in METHOD_REGISTRY.items():
            assert spec.key == key

    def test_every_value_is_a_method_spec(self):
        assert all(isinstance(spec, MethodSpec) for spec in METHOD_REGISTRY.values())

    def test_every_function_is_callable(self):
        assert all(callable(spec.fn) for spec in METHOD_REGISTRY.values())

    def test_harshness_tiers_are_valid(self):
        for spec in METHOD_REGISTRY.values():
            assert spec.harshness in HARSHNESS_TIERS

    def test_every_method_has_a_citation(self):
        for key, spec in METHOD_REGISTRY.items():
            assert spec.citation, f"{key} has no citation"

    def test_no_duplicate_functions(self):
        functions = [spec.fn for spec in METHOD_REGISTRY.values()]
        assert len(functions) == len(set(functions))


class TestExcludedMethods:
    """The five methods that are deliberately absent, and stay absent."""

    @pytest.mark.parametrize(
        "key",
        ["24_peer_k10", "32_deepmnn", "35_dasc", "36_explobatch", "39_procrustes"],
    )
    def test_excluded_key_is_not_registered(self, key):
        assert key not in METHOD_REGISTRY

    @pytest.mark.parametrize(
        "key",
        ["24_peer_k10", "32_deepmnn", "35_dasc", "36_explobatch", "39_procrustes"],
    )
    def test_excluded_key_has_a_documented_reason(self, key):
        assert key in EXCLUDED_METHODS
        assert len(EXCLUDED_METHODS[key]) > 20

    def test_lookup_explains_why_rather_than_just_missing(self):
        with pytest.raises(KeyError, match="licence"):
            get_method("39_procrustes")

    def test_unknown_key_lists_what_is_available(self):
        with pytest.raises(KeyError, match="Available:"):
            get_method("99_nonexistent")


class TestHyperparamsMatchSignatures:
    """Every declared hyperparameter must actually be accepted by its function.

    This is the check that makes the generated docs trustworthy: a parameter the
    registry advertises but the function ignores would be silently unreachable, which is
    the exact defect ComboBatch exists to fix.
    """

    @pytest.mark.parametrize("key", sorted(METHOD_REGISTRY))
    def test_declared_params_exist_in_the_signature(self, key):
        spec = METHOD_REGISTRY[key]
        signature = inspect.signature(spec.fn)
        accepts_kwargs = any(
            p.kind is inspect.Parameter.VAR_KEYWORD
            for p in signature.parameters.values()
        )
        for name in spec.hyperparams:
            assert (
                name in signature.parameters or accepts_kwargs
            ), f"{key} declares {name!r} but its function does not accept it"

    @pytest.mark.parametrize("key", sorted(METHOD_REGISTRY))
    def test_declared_defaults_match_the_signature(self, key):
        spec = METHOD_REGISTRY[key]
        signature = inspect.signature(spec.fn)
        for name, param in spec.hyperparams.items():
            if name in signature.parameters:
                actual = signature.parameters[name].default
                if actual is not inspect.Parameter.empty:
                    assert actual == param.default, (
                        f"{key}.{name}: registry says {param.default!r}, "
                        f"function says {actual!r}"
                    )

    @pytest.mark.parametrize("key", sorted(METHOD_REGISTRY))
    def test_every_hyperparam_has_help_text(self, key):
        for name, param in METHOD_REGISTRY[key].hyperparams.items():
            assert param.help, f"{key}.{name} has no help text"

    @pytest.mark.parametrize("key", sorted(METHOD_REGISTRY))
    def test_hyperparam_names_match_their_keys(self, key):
        for name, param in METHOD_REGISTRY[key].hyperparams.items():
            assert param.name == name


class TestCapabilityFlags:
    def test_r_methods_declare_their_packages(self):
        for key, spec in METHOD_REGISTRY.items():
            if spec.requires_r:
                assert spec.r_packages, f"{key} requires R but names no R package"

    def test_non_r_methods_declare_no_r_packages(self):
        for key, spec in METHOD_REGISTRY.items():
            if not spec.requires_r:
                assert not spec.r_packages, f"{key} names R packages but not requires_r"

    def test_reference_batch_users_declare_target_group(self):
        for key, spec in METHOD_REGISTRY.items():
            if spec.uses_reference_batch:
                assert (
                    "target_group" in spec.hyperparams
                ), f"{key} uses a reference batch but does not declare target_group"

    def test_rnaseq_only_methods_are_the_expected_two(self):
        # 39_procrustes was the third; it is excluded on licensing grounds.
        rnaseq_only = {k for k, s in METHOD_REGISTRY.items() if s.rnaseq_only}
        assert rnaseq_only == {"22_tmm", "23_vst"}

    def test_availability_is_reported_not_guessed(self):
        for spec in METHOD_REGISTRY.values():
            assert spec.is_available() == (not spec.missing_backends())


class TestParamHelpers:
    def test_defaults_round_trip(self):
        spec = METHOD_REGISTRY["10_mnn"]
        assert spec.defaults() == {"k": 20}

    def test_non_default_filters_out_defaults(self):
        spec = METHOD_REGISTRY["10_mnn"]
        assert spec.non_default({"k": 20}) == {}
        assert spec.non_default({"k": 50}) == {"k": 50}

    def test_validate_rejects_unknown_names(self):
        with pytest.raises(Exception, match="unknown parameter"):
            METHOD_REGISTRY["10_mnn"].validate_params({"kk": 5})

    def test_validate_rejects_out_of_range_values(self):
        with pytest.raises(ValueError, match="must be >="):
            METHOD_REGISTRY["10_mnn"].validate_params({"k": 0})

    def test_validate_rejects_invalid_choice(self):
        with pytest.raises(ValueError, match="expected one of"):
            METHOD_REGISTRY["21_harmonizr"].validate_params({"algorithm": "nope"})

    def test_validate_accepts_good_values(self):
        METHOD_REGISTRY["21_harmonizr"].validate_params({"algorithm": "limma"})
