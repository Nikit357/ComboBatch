"""The R bridge — in particular ``r_arglist``, the R-injection boundary.

None of these need R: the point of the module is that it is importable and its
value-rendering is testable without one.
"""

from __future__ import annotations

import sys
import types
from pathlib import Path

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


class TestPackageAvailabilityDoesNotUseAnInvisibleResult:
    """The check crashed inside the image for four builds, and could not crash here.

    ``r_package_available`` returns early when rpy2 is absent, so on a laptop it never
    reached the line that mattered: ``ro.r("requireNamespace(...)")`` returns the R value
    *invisibly*, rpy2 3.6 maps that to ``None``, and ``None[0]`` raised TypeError — for
    the installed package and the missing one alike. Only the image caught it, an hour of
    build at a time. These two guards run anywhere.
    """

    def test_the_expression_that_cannot_work_is_gone(self):
        # Comments explaining why it is gone must not count as its return.
        source = "\n".join(
            line
            for line in Path(rinterop.__file__).read_text().splitlines()
            if not line.lstrip().startswith("#")
        )
        assert "requireNamespace" not in source, (
            "requireNamespace returns invisibly, so rpy2 hands back None and there is no "
            "answer to read — use rpy2.robjects.packages.isinstalled instead"
        )

    @pytest.mark.parametrize("installed", [True, False])
    def test_the_answer_is_a_bool_either_way(self, monkeypatch, installed):
        """Both outcomes must survive the round trip, not just the happy one."""
        packages = types.ModuleType("rpy2.robjects.packages")
        packages.isinstalled = lambda name: installed
        monkeypatch.setitem(sys.modules, "rpy2.robjects.packages", packages)
        monkeypatch.setattr(rinterop, "r_available", lambda: True)
        rinterop.r_package_available.cache_clear()

        try:
            assert rinterop.r_package_available("limma") is installed
        finally:
            rinterop.r_package_available.cache_clear()


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
