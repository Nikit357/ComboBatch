"""Hyperparameter parsing, coercion and tag encoding."""

from __future__ import annotations

import pytest

from combobatch.params import (
    MAX_PARAM_TAG_LENGTH,
    ParamParseError,
    coerce_value,
    coerce_with_spec,
    encode_param_tag,
    params_hash,
    parse_params_arg,
    validate_param_names,
    validate_param_tag,
)


class TestCoercion:
    @pytest.mark.parametrize(
        "raw,expected",
        [
            ("50", 50),
            ("0.05", 0.05),
            ("-3", -3),
            ("1e-6", 1e-6),
            ("true", True),
            ("False", False),
            ("none", None),
            ("ComBat", "ComBat"),
            ("a+b+c", ["a", "b", "c"]),
            ("1+2", [1, 2]),
        ],
    )
    def test_infers_type(self, raw, expected):
        assert coerce_value(raw) == expected

    def test_numeric_one_and_zero_stay_numeric(self):
        # A bare 1 is far more often an int parameter than a boolean, so it must not be
        # silently promoted to True.
        assert coerce_value("1") == 1
        assert not isinstance(coerce_value("1"), bool)
        assert coerce_value("0") == 0

    def test_declared_type_wins(self):
        assert coerce_with_spec("k", "50", int) == 50
        assert coerce_with_spec("limit", "1", float) == pytest.approx(1.0)
        assert coerce_with_spec("flag", "yes", bool) is True

    def test_declared_type_mismatch_raises(self):
        with pytest.raises(ParamParseError, match="expected int"):
            coerce_with_spec("k", "banana", int)
        with pytest.raises(ParamParseError, match="expected a boolean"):
            coerce_with_spec("flag", "banana", bool)


class TestParseParamsArg:
    def test_single_target_single_param(self):
        assert parse_params_arg("10_mnn:k=50") == {"10_mnn": {"k": 50}}

    def test_multiple_targets_and_params(self):
        parsed = parse_params_arg(
            "10_mnn:k=50;38_harman:limit=0.05;31_ruv3prps:k_factors=8,min_cell_size=3"
        )
        assert parsed == {
            "10_mnn": {"k": 50},
            "38_harman": {"limit": 0.05},
            "31_ruv3prps": {"k_factors": 8, "min_cell_size": 3},
        }

    def test_empty_input_is_empty_mapping(self):
        assert parse_params_arg("") == {}
        assert parse_params_arg("   ") == {}

    def test_missing_colon_raises(self):
        with pytest.raises(ParamParseError, match="no ':'"):
            parse_params_arg("k=50")

    def test_missing_equals_raises(self):
        with pytest.raises(ParamParseError, match="<name>=<value>"):
            parse_params_arg("10_mnn:k")

    def test_repeated_parameter_raises(self):
        with pytest.raises(ParamParseError, match="given twice"):
            parse_params_arg("10_mnn:k=5,k=50")

    def test_empty_target_raises(self):
        with pytest.raises(ParamParseError, match="empty target"):
            parse_params_arg(":k=50")


class TestValidateParamNames:
    def test_declared_names_pass(self):
        validate_param_names("10_mnn", {"k": 50}, {"k", "n_pcs"})

    def test_unknown_name_raises_with_suggestion(self):
        with pytest.raises(ParamParseError, match="did you mean 'k_factors'"):
            validate_param_names("31_ruv3prps", {"kfactors": 8}, {"k_factors"})

    def test_unknown_name_lists_declared(self):
        with pytest.raises(ParamParseError, match="Declared: k, n_pcs"):
            validate_param_names("10_mnn", {"zzz": 1}, {"k", "n_pcs"})


class TestParamTag:
    def test_defaults_produce_no_tag(self):
        # This is what keeps a default run's filename free of a tag segment entirely.
        assert encode_param_tag({}) is None

    def test_readable_encoding(self):
        assert encode_param_tag({"k": 50, "limit": 0.05}) == "k50-limit0.05"

    def test_sorted_by_name_so_order_does_not_matter(self):
        assert encode_param_tag({"limit": 0.05, "k": 50}) == encode_param_tag(
            {"k": 50, "limit": 0.05}
        )

    def test_booleans_compact_to_single_letters(self):
        assert encode_param_tag({"flag": True, "other": False}) == "flagT-otherF"

    def test_lists_join_with_plus(self):
        assert encode_param_tag({"cols": ["a", "b"]}) == "colsa+b"

    def test_floats_keep_their_decimal_point(self):
        assert encode_param_tag({"limit": 0.1}) == "limit0.1"

    @pytest.mark.parametrize(
        "value,expected",
        [
            ("a/b", "xa%2Fb"),
            ("a b", "xa%20b"),
            ("a__b", "xa%5F%5Fb"),
            ("a%b", "xa%25b"),
        ],
    )
    def test_unsafe_characters_are_escaped(self, value, expected):
        # "__" in particular must never survive: it separates fields in the key grammar.
        assert encode_param_tag({"x": value}) == expected

    def test_long_tag_truncates_with_hash_and_respects_the_limit(self):
        params = {f"parameter_number_{i}": i for i in range(20)}
        tag = encode_param_tag(params)
        assert len(tag) <= MAX_PARAM_TAG_LENGTH
        assert tag.endswith(params_hash(params))

    def test_truncated_tags_stay_distinct(self):
        a = {f"parameter_number_{i}": i for i in range(20)}
        b = dict(a, parameter_number_0=999)
        assert encode_param_tag(a) != encode_param_tag(b)

    def test_hash_is_order_independent_and_stable(self):
        assert params_hash({"a": 1, "b": 2}) == params_hash({"b": 2, "a": 1})
        assert params_hash({"a": 1}) == params_hash({"a": 1})


class TestValidateParamTag:
    def test_accepts_a_plain_tag(self):
        assert validate_param_tag("k50") == "k50"

    def test_rejects_empty(self):
        with pytest.raises(ParamParseError, match="must not be empty"):
            validate_param_tag("  ")

    def test_rejects_field_separator(self):
        with pytest.raises(ParamParseError, match="separates key fields"):
            validate_param_tag("a__b")

    def test_rejects_path_separator_and_spaces(self):
        with pytest.raises(ParamParseError):
            validate_param_tag("a/b")
        with pytest.raises(ParamParseError):
            validate_param_tag("a b")

    def test_rejects_a_tag_shaped_like_a_post_flag(self):
        # Otherwise a three-segment key would be genuinely ambiguous.
        with pytest.raises(ParamParseError, match="reserved"):
            validate_param_tag("post1")
