"""The generalization: no metric group knows a column name.

The donor wrote four batch columns, four biology columns, a cohort column and two
diagnosis string literals into every group. These tests assert that renaming a dataset's
columns renames the metrics with it, and that none of the donor's names survive anywhere
in the package.
"""

from __future__ import annotations

import pathlib

import pytest

from combobatch.config import ColumnSpec, MetricsSpec
from combobatch.metrics import MetricColumns
from combobatch.metrics.runner import compute_metrics

from tests.conftest import make_expression

# Names the donor hardcoded. None may appear in the shipped package.
DONOR_LITERALS = (
    "RNA_BATCH",
    "COHORT_LABEL",
    "PLATFORM_RNA",
    "RNASEQ_SOURCE",
    "TUMOR_NORMAL",
    "Major_group",
    "Diagnosis_cell_type_unified",
    "Follicular_Lymphoma",
    "Diffuse_Large_B_Cell_Lymphoma",
    "Normal_B",
)


@pytest.fixture
def renamed():
    """A dataset whose columns are named nothing like the donor's."""
    exp, ann = make_expression(n_samples=40, n_genes=50)
    ann = ann.rename(columns={"batch": "site", "bio": "celltype", "cohort": "study"})
    return exp, ann


class TestColumnResolution:
    def test_roles_come_from_the_spec(self):
        columns = MetricColumns.from_spec(
            ColumnSpec(batch="site", bio="celltype", cohort="study")
        )
        assert columns.batch == "site"
        assert columns.bio == "celltype"
        assert columns.cohort == "study"

    def test_cohort_falls_back_to_batch(self):
        # Several groups partition by cohort; the batch is the closest generic stand-in
        # when the user names no cohort column at all.
        columns = MetricColumns.from_spec(ColumnSpec(batch="site", bio="celltype"))
        assert columns.cohort == "site"

    def test_metric_column_lists_default_to_the_primaries(self):
        columns = MetricColumns.from_spec(ColumnSpec(batch="site", bio="celltype"))
        assert columns.batch_cols == ("site",)
        assert columns.bio_cols == ("celltype",)

    def test_all_cols_is_deduplicated_and_stable(self):
        columns = MetricColumns.from_spec(
            ColumnSpec(
                batch="site",
                bio="celltype",
                cohort="site",
                metric_batch_cols=("site", "platform"),
            )
        )
        assert columns.all_cols == ("site", "platform", "celltype")


class TestMetricKeysFollowTheColumns:
    def test_group_a_keys_carry_the_configured_names(self, renamed):
        exp, ann = renamed
        result = compute_metrics(
            exp,
            ann,
            ColumnSpec(batch="site", bio="celltype", cohort="study"),
            groups=frozenset("A"),
        )
        assert "r2_site" in result
        assert "r2_celltype" in result
        assert "dsc_site" in result
        assert not any("RNA_BATCH" in key for key in result)

    def test_several_batch_columns_give_one_key_each(self, renamed):
        exp, ann = renamed
        result = compute_metrics(
            exp,
            ann,
            ColumnSpec(
                batch="site",
                bio="celltype",
                cohort="study",
                metric_batch_cols=("site", "study"),
            ),
            groups=frozenset("A"),
        )
        assert "r2_pc1_site" in result
        assert "r2_pc1_study" in result

    def test_group_g_keys_follow_the_biology_columns(self, renamed):
        exp, ann = renamed
        result = compute_metrics(
            exp,
            ann,
            ColumnSpec(
                batch="site",
                bio="celltype",
                metric_bio_cols=("celltype", "binary"),
            ),
            groups=frozenset("G"),
        )
        assert "graph_connectivity_celltype" in result
        assert "graph_connectivity_binary" in result

    def test_an_absent_column_yields_nulls_not_a_crash(self, renamed):
        exp, ann = renamed
        result = compute_metrics(
            exp,
            ann,
            ColumnSpec(batch="site", bio="celltype", metric_batch_cols=("nowhere",)),
            groups=frozenset("A"),
        )
        assert "r2_nowhere" in result
        assert result["r2_nowhere"] != result["r2_nowhere"]  # NaN


class TestGroupNIsGeneric:
    def test_classes_default_to_the_most_frequent_levels(self, renamed):
        from combobatch.metrics.groups import predictor_labels

        _, ann = renamed
        labels = predictor_labels(ann, "celltype", n_classes=2)
        assert labels is not None
        assert set(labels.dropna().unique()) <= {"TypeA", "TypeB"}

    def test_explicit_classes_win(self, renamed):
        from combobatch.metrics.groups import predictor_labels

        _, ann = renamed
        labels = predictor_labels(ann, "celltype", classes=["TypeA"])
        # One surviving class is not a classification problem.
        assert labels is None

    def test_a_column_with_one_level_is_untestable(self, renamed):
        from combobatch.metrics.groups import predictor_labels

        _, ann = renamed
        ann = ann.assign(flat="same")
        assert predictor_labels(ann, "flat") is None

    def test_the_class_count_actually_used_is_reported(self, renamed):
        exp, ann = renamed
        result = compute_metrics(
            exp,
            ann,
            ColumnSpec(batch="site", bio="celltype"),
            spec=MetricsSpec(n_perm=1, min_test_n=5, n_pcs=3, skip_slow=False),
            groups=frozenset("N"),
        )
        # The fixture's biology column has two levels, so a 3-class request resolves to
        # 2 — and the result says so rather than implying three.
        assert result.get("pv_lobo3_n_classes_used") == 2


class TestNoDonorLiteralsSurvive:
    def test_the_package_is_free_of_the_donor_column_names(self):
        package = pathlib.Path("combobatch")
        offenders: list[str] = []
        for path in package.rglob("*.py"):
            # subset.py documents user-supplied query strings as examples.
            if path.name == "subset.py":
                continue
            text = path.read_text()
            for literal in DONOR_LITERALS:
                if literal in text:
                    offenders.append(f"{path}: {literal}")
        assert not offenders, "dataset-specific literals leaked in: " + "; ".join(
            offenders
        )
