"""The metric-group registry — the one that has no upstream counterpart.

Upstream spread this information over three files: a hardcoded call sequence, a sentinel
dict in a different module, and a "needs the reference" constant in a third. These tests
exist so the single table here cannot drift from the functions it points at.
"""

from __future__ import annotations

import inspect

import pytest

from combobatch.config import DEFAULT_METRIC_GROUPS, METRIC_GROUP_LETTERS
from combobatch.metrics import (
    DEFAULT_GROUPS,
    GROUP_LETTERS,
    METRIC_GROUP_REGISTRY,
    MetricColumns,
    MetricGroupSpec,
    get_group,
    ordered_groups,
    required_embeddings,
)

EXPECTED_GROUP_COUNT = 14


@pytest.fixture
def columns():
    """Resolved columns with deliberately un-FL-like names."""
    from combobatch.config import ColumnSpec

    return MetricColumns.from_spec(
        ColumnSpec(batch="site", bio="celltype", cohort="study")
    )


class TestRegistryShape:
    def test_fourteen_groups_a_through_n(self):
        assert len(METRIC_GROUP_REGISTRY) == EXPECTED_GROUP_COUNT
        assert set(METRIC_GROUP_REGISTRY) == set("ABCDEFGHIJKLMN")

    def test_letters_match_their_specs(self):
        for letter, spec in METRIC_GROUP_REGISTRY.items():
            assert spec.letter == letter

    def test_every_value_is_a_spec_with_a_callable(self):
        for spec in METRIC_GROUP_REGISTRY.values():
            assert isinstance(spec, MetricGroupSpec)
            assert callable(spec.fn)

    def test_every_group_has_a_description(self):
        for letter, spec in METRIC_GROUP_REGISTRY.items():
            assert spec.description, f"group {letter} has no description"

    def test_execution_order_covers_every_group_once(self):
        assert sorted(GROUP_LETTERS) == sorted(METRIC_GROUP_REGISTRY)
        assert len(GROUP_LETTERS) == len(set(GROUP_LETTERS))

    def test_speeds_are_from_the_known_set(self):
        for spec in METRIC_GROUP_REGISTRY.values():
            assert spec.speed in {"fast", "moderate", "slow", "very_slow"}


class TestNoDuplicateSourceOfTruth:
    def test_config_defaults_agree_with_the_registry(self):
        # Two literals for one fact is exactly the defect this project exists to fix.
        assert DEFAULT_METRIC_GROUPS == DEFAULT_GROUPS

    def test_config_letters_agree_with_the_registry(self):
        assert METRIC_GROUP_LETTERS == frozenset(METRIC_GROUP_REGISTRY)

    def test_the_default_set_is_cheap(self):
        for letter in DEFAULT_GROUPS:
            assert METRIC_GROUP_REGISTRY[letter].speed in {"fast", "moderate", "slow"}


class TestSentinels:
    def test_every_sentinel_is_unique(self, columns):
        sentinels = [
            spec.sentinel_key(columns) for spec in METRIC_GROUP_REGISTRY.values()
        ]
        assert len(set(sentinels)) == len(sentinels)

    def test_sentinels_resolve_against_the_configured_columns(self, columns):
        assert METRIC_GROUP_REGISTRY["A"].sentinel_key(columns) == "r2_site"
        assert (
            METRIC_GROUP_REGISTRY["G"].sentinel_key(columns)
            == "graph_connectivity_celltype"
        )

    def test_column_free_sentinels_are_unchanged(self, columns):
        assert METRIC_GROUP_REGISTRY["K"].sentinel_key(columns) == "n_genes_noNA"
        assert METRIC_GROUP_REGISTRY["J"].sentinel_key(columns) == "pct_var_pc1"

    def test_no_sentinel_names_a_donor_column(self, columns):
        rendered = " ".join(
            spec.sentinel_key(columns) for spec in METRIC_GROUP_REGISTRY.values()
        )
        for literal in ("RNA_BATCH", "COHORT_LABEL", "Major_group", "PLATFORM_RNA"):
            assert literal not in rendered


class TestDeclarationsMatchReality:
    def test_only_group_l_needs_a_reference(self):
        needing = {
            letter
            for letter, spec in METRIC_GROUP_REGISTRY.items()
            if spec.needs_reference
        }
        assert needing == {"L"}

    def test_group_l_accepts_the_reference_as_a_keyword(self):
        # The donor's signature was (exp, ref, ann) — out of step with every other
        # group, and transposable without a word.
        signature = inspect.signature(METRIC_GROUP_REGISTRY["L"].fn)
        assert signature.parameters["ref_df"].kind is inspect.Parameter.KEYWORD_ONLY

    def test_panel_users_accept_a_panel_keyword(self):
        for letter, spec in METRIC_GROUP_REGISTRY.items():
            if not spec.needs_panel:
                continue
            assert "panel" in inspect.signature(spec.fn).parameters, letter

    def test_only_group_f_needs_r(self):
        needing = {
            letter for letter, spec in METRIC_GROUP_REGISTRY.items() if spec.needs_r
        }
        assert needing == {"F"}

    def test_embedding_declarations_are_from_the_known_set(self):
        for spec in METRIC_GROUP_REGISTRY.values():
            assert set(spec.needs_embedding) <= {"pca", "umap", "tsne"}

    def test_only_group_c_needs_the_two_dimensional_embeddings(self):
        needing = {
            letter
            for letter, spec in METRIC_GROUP_REGISTRY.items()
            if "umap" in spec.needs_embedding
        }
        assert needing == {"C"}


class TestHelpers:
    def test_ordered_groups_follows_the_execution_order(self):
        letters = [spec.letter for spec in ordered_groups(set("NEAK"))]
        assert letters == ["E", "K", "A", "N"]

    def test_ordered_groups_ignores_unknown_letters(self):
        assert [spec.letter for spec in ordered_groups(set("EZ"))] == ["E"]

    def test_required_embeddings_is_the_union(self):
        assert required_embeddings(set("EK")) == set()
        assert required_embeddings(set("A")) == {"pca"}
        assert required_embeddings(set("AC")) == {"pca", "umap", "tsne"}

    def test_get_group_is_case_insensitive(self):
        assert get_group("a") is METRIC_GROUP_REGISTRY["A"]

    def test_get_group_names_the_valid_letters(self):
        with pytest.raises(KeyError, match="ABCDEFGHIJKLMN"):
            get_group("Z")
