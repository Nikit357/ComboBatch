"""The R bridge — in particular ``r_arglist``, the R-injection boundary.

None of these need R: the point of the module is that it is importable and its
value-rendering is testable without one.
"""

from __future__ import annotations

import pytest

from combobatch.methods import rinterop


class TestImportableWithoutR:
    def test_module_imports(self):
        # The core install has no rpy2; importing the bridge must still work.
        assert rinterop.r_available() in (True, False)

    def test_r_package_available_is_false_without_r(self):
        if rinterop.r_available():
            pytest.skip("R is installed here")
        assert rinterop.r_package_available("limma") is False

    def test_require_r_raises_a_helpful_error(self):
        if rinterop.r_available():
            pytest.skip("R is installed here")
        with pytest.raises(NotImplementedError, match=r"combobatch\[r\]"):
            rinterop.require_r("03_limma")

    def test_r_gc_is_a_no_op_without_r(self):
        rinterop.r_gc()


class TestRLiteral:
    @pytest.mark.parametrize(
        "value,expected",
        [
            (True, "TRUE"),
            (False, "FALSE"),
            (None, "NULL"),
            (5, "5"),
            (0.1, "0.1"),
            (-3, "-3"),
        ],
    )
    def test_scalars(self, value, expected):
        assert rinterop.r_literal(value) == expected

    def test_strings_are_quoted(self):
        assert rinterop.r_literal("TMM") == '"TMM"'

    def test_quotes_and_backslashes_are_escaped(self):
        assert rinterop.r_literal('a"b') == '"a\\"b"'
        assert rinterop.r_literal("a\\b") == '"a\\\\b"'

    def test_flat_sequences_become_c(self):
        assert rinterop.r_literal(["a", "b"]) == 'c("a", "b")'
        assert rinterop.r_literal((1, 2)) == "c(1, 2)"


class TestRLiteralRefusesUnsafeValues:
    """R scripts are built with f-strings, so anything not a scalar must be refused
    rather than stringified into source text."""

    @pytest.mark.parametrize(
        "value", [{"a": 1}, {1, 2}, object(), [[1, 2]], [{"a": 1}]]
    )
    def test_unsupported_types_raise(self, value):
        with pytest.raises(TypeError):
            rinterop.r_literal(value)


class TestRArglist:
    def test_empty_params_render_empty(self):
        assert rinterop.r_arglist({}) == ""
        assert rinterop.r_arglist({}, leading_comma=True) == ""

    def test_single_argument(self):
        assert rinterop.r_arglist({"k": 5}) == "k = 5"

    def test_arguments_are_sorted_for_stability(self):
        assert rinterop.r_arglist({"b": 2, "a": 1}) == "a = 1, b = 2"

    def test_mixed_types(self):
        rendered = rinterop.r_arglist({"method": "TMM", "log": True, "prior": 1.0})
        assert rendered == 'log = TRUE, method = "TMM", prior = 1.0'

    def test_leading_comma_for_splicing(self):
        assert rinterop.r_arglist({"k": 5}, leading_comma=True) == ", k = 5"

    def test_dotted_r_names_are_allowed(self):
        # softImpute really does take "rank.max" and "lambda".
        assert rinterop.r_arglist({"rank.max": 30}) == "rank.max = 30"

    @pytest.mark.parametrize("name", ["1bad", "has space", "has-dash", "", "a;b"])
    def test_invalid_argument_names_raise(self, name):
        with pytest.raises(ValueError, match="valid R argument name"):
            rinterop.r_arglist({name: 1})

    def test_unsafe_value_raises_through_arglist(self):
        with pytest.raises(TypeError):
            rinterop.r_arglist({"x": {"nested": 1}})
