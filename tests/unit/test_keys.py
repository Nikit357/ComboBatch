"""The output-key grammar: it must round-trip, and it must fail loudly."""

from __future__ import annotations

import pytest

from combobatch.params import (
    KeyParseError,
    OutputKey,
    build_output_key,
    build_prepared_key,
    parse_output_key,
)


class TestBuild:
    def test_minimal_key(self):
        assert build_output_key("strict", "10_mnn") == "strict__10_mnn"

    def test_with_post_flag(self):
        assert (
            build_output_key("strict", "10_mnn", None, False) == "strict__10_mnn__post0"
        )
        assert (
            build_output_key("strict", "10_mnn", None, True) == "strict__10_mnn__post1"
        )

    def test_with_param_tag(self):
        assert build_output_key("knn", "10_mnn", "k50") == "knn__10_mnn__k50"

    def test_with_both(self):
        assert (
            build_output_key("knn", "10_mnn", "k50", False) == "knn__10_mnn__k50__post0"
        )

    def test_post_rm_none_omits_the_segment(self):
        # A run with post-removal disabled must not carry a misleading "post0".
        assert "post" not in build_output_key("strict", "01_raw", None, None)

    def test_rejects_separator_in_a_field(self):
        with pytest.raises(KeyParseError, match="separates fields"):
            build_output_key("stri__ct", "10_mnn")


class TestRoundTrip:
    @pytest.mark.parametrize(
        "key",
        [
            OutputKey("strict", "10_mnn"),
            OutputKey("strict", "10_mnn", None, False),
            OutputKey("strict", "10_mnn", None, True),
            OutputKey("knn", "38_harman", "limit0.05"),
            OutputKey("softimpute", "04_sva", "k50-limit0.05", True),
        ],
    )
    def test_build_then_parse_is_identity(self, key):
        assert parse_output_key(str(key)) == key


class TestParse:
    def test_two_segments(self):
        parsed = parse_output_key("strict__10_mnn")
        assert (parsed.imputation, parsed.method) == ("strict", "10_mnn")
        assert parsed.param_tag is None
        assert parsed.post_rm is None

    def test_third_segment_read_as_post_when_it_looks_like_one(self):
        parsed = parse_output_key("strict__10_mnn__post1")
        assert parsed.post_rm is True
        assert parsed.param_tag is None

    def test_third_segment_read_as_tag_otherwise(self):
        parsed = parse_output_key("strict__10_mnn__k50")
        assert parsed.param_tag == "k50"
        assert parsed.post_rm is None

    def test_four_segments(self):
        parsed = parse_output_key("knn__10_mnn__k50__post0")
        assert parsed.param_tag == "k50"
        assert parsed.post_rm is False

    def test_method_keeps_its_own_underscores(self):
        # Single underscores are inside field values; only "__" separates fields.
        assert parse_output_key("strict__08_inmoose_combatseq").method == (
            "08_inmoose_combatseq"
        )


class TestFailsLoudly:
    """The donor silently dropped keys that did not split into exactly four parts."""

    @pytest.mark.parametrize(
        "stem", ["", "   ", "onlyone", "a__b__c__d__e", "a____b", "__b", "a__"]
    )
    def test_malformed_keys_raise(self, stem):
        with pytest.raises(KeyParseError):
            parse_output_key(stem)

    def test_error_names_the_expected_grammar(self):
        with pytest.raises(KeyParseError, match=r"\{imp\}__\{method\}"):
            parse_output_key("a__b__c__d__e")


class TestPreparedKeys:
    def test_without_tag(self):
        assert build_prepared_key("strict", None, "exp") == "strict__exp"

    def test_with_tag(self):
        assert build_prepared_key("knn", "knn_k10", "ann") == "knn__knn_k10__ann"

    def test_rejects_unknown_kind(self):
        with pytest.raises(KeyParseError, match="'exp' or 'ann'"):
            build_prepared_key("strict", None, "genes")
