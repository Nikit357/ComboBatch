"""--subset-query evaluation and its pre-flight validation."""

from __future__ import annotations

import pandas as pd
import pytest

from combobatch.subset import (
    SubsetQueryError,
    apply_subset_query,
    referenced_columns,
    validate_subset_query,
)


@pytest.fixture
def ann():
    return pd.DataFrame(
        {
            "RNA_BATCH": [
                "RNASeq_FF_PolyA",
                "RNASeq_FFPE_Exome",
                "GPL570_FF_Unknown",
                "GPL570_FFPE_Unknown",
            ],
            "Diagnosis": ["FL", "DLBCL", "FL", "Normal"],
            "COHORT": ["c1", "c1", "c2", "c2"],
        },
        index=["s1", "s2", "s3", "s4"],
    )


class TestReferencedColumns:
    def test_finds_bare_identifiers(self):
        assert "RNA_BATCH" in referenced_columns("RNA_BATCH == 'x'")

    def test_ignores_string_literals(self):
        found = referenced_columns("RNA_BATCH.str.startswith('RNASeq')")
        assert "RNASeq" not in found

    def test_ignores_query_keywords(self):
        found = referenced_columns("A == 1 and B == 2")
        assert found == {"A", "B"}


class TestValidate:
    def test_accepts_a_valid_query(self, ann):
        validate_subset_query("RNA_BATCH.str.startswith('RNASeq')", ann)

    def test_rejects_empty(self, ann):
        with pytest.raises(SubsetQueryError, match="must not be empty"):
            validate_subset_query("   ", ann)

    def test_rejects_unknown_column_with_a_suggestion(self, ann):
        with pytest.raises(SubsetQueryError, match="did you mean 'RNA_BATCH'"):
            validate_subset_query("RNA_BATCHH == 'x'", ann)


class TestApply:
    def test_prefix_predicate(self, ann):
        # The generic replacement for the donor's hardcoded C_rnaseq_only strategy.
        kept = apply_subset_query(ann, "RNA_BATCH.str.startswith('RNASeq')")
        assert list(kept.index) == ["s1", "s2"]

    def test_contains_predicate(self, ann):
        kept = apply_subset_query(ann, "RNA_BATCH.str.contains('_FFPE_')")
        assert list(kept.index) == ["s2", "s4"]

    def test_membership_predicate(self, ann):
        kept = apply_subset_query(ann, "Diagnosis in ['FL', 'DLBCL']")
        assert list(kept.index) == ["s1", "s2", "s3"]

    def test_conjunction(self, ann):
        kept = apply_subset_query(
            ann, "RNA_BATCH.str.startswith('GPL') and Diagnosis == 'FL'"
        )
        assert list(kept.index) == ["s3"]

    def test_zero_matches_is_an_error_not_an_empty_run(self, ann):
        # Almost always a typo in a batch name; running on nothing would waste hours.
        with pytest.raises(SubsetQueryError, match="matched 0 of 4"):
            apply_subset_query(ann, "Diagnosis == 'NoSuchDiagnosis'")

    def test_unknown_column_raises_before_evaluation(self, ann):
        with pytest.raises(SubsetQueryError, match="unknown column"):
            apply_subset_query(ann, "Nonexistent == 'x'")

    def test_malformed_query_raises(self, ann):
        with pytest.raises(SubsetQueryError):
            apply_subset_query(ann, "Diagnosis ===== 'FL'")

    def test_does_not_mutate_the_input(self, ann):
        before = ann.copy()
        apply_subset_query(ann, "Diagnosis == 'FL'")
        pd.testing.assert_frame_equal(ann, before)
