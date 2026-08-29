"""The orchestrator: isolation, flushing, sentinels, and the force-beats-skip fix."""

from __future__ import annotations

import dataclasses
import json

import numpy as np
import pytest

from combobatch.config import ColumnSpec, MetricsSpec
from combobatch.metrics import METRIC_GROUP_REGISTRY, MetricColumns
from combobatch.metrics.runner import compute_metrics, resolve_active_groups

from tests.conftest import make_expression


@pytest.fixture
def columns():
    return MetricColumns.from_spec(
        ColumnSpec(batch="batch", bio="bio", cohort="cohort")
    )


@pytest.fixture
def data():
    return make_expression(n_samples=50, n_genes=60)


@pytest.fixture
def column_spec():
    return ColumnSpec(batch="batch", bio="bio", cohort="cohort")


def _patch_group(monkeypatch, letter, fn):
    """Replace one group's implementation; MetricGroupSpec is frozen."""
    monkeypatch.setitem(
        METRIC_GROUP_REGISTRY,
        letter,
        dataclasses.replace(METRIC_GROUP_REGISTRY[letter], fn=fn),
    )


class TestGroupSelection:
    def test_skip_slow_excludes_the_very_slow_groups(self):
        active = resolve_active_groups(frozenset("EN"), skip_slow=True)
        assert active == {"E"}

    def test_force_beats_skip_slow(self):
        # The donor's comment claimed this and the code did the opposite, because the
        # slow gate ran after the force set was unioned in.
        active = resolve_active_groups(
            frozenset("EN"), skip_slow=True, force=frozenset("N")
        )
        assert active == {"E", "N"}

    def test_skip_wm_excludes_group_i(self):
        assert resolve_active_groups(frozenset("EI"), skip_wm=True) == {"E"}

    def test_force_beats_skip_wm(self):
        active = resolve_active_groups(
            frozenset("EI"), skip_wm=True, force=frozenset("I")
        )
        assert active == {"E", "I"}

    def test_unknown_letters_are_dropped(self):
        assert resolve_active_groups(frozenset("EZ")) == {"E"}

    def test_a_present_sentinel_skips_recomputation(self, columns):
        prior = {"n_samples": 50}
        assert resolve_active_groups(frozenset("EK"), prior=prior, columns=columns) == {
            "K"
        }

    def test_force_beats_an_existing_sentinel(self, columns):
        active = resolve_active_groups(
            frozenset("EK"),
            prior={"n_samples": 50},
            columns=columns,
            force=frozenset("E"),
        )
        assert active == {"E", "K"}

    def test_sentinels_resolve_against_the_configured_columns(self):
        renamed = MetricColumns.from_spec(ColumnSpec(batch="site", bio="celltype"))
        # A prior keyed on the donor's column name must not satisfy this run.
        assert "A" in resolve_active_groups(
            frozenset("A"), prior={"r2_RNA_BATCH": 0.1}, columns=renamed
        )
        assert "A" not in resolve_active_groups(
            frozenset("A"), prior={"r2_site": 0.1}, columns=renamed
        )


class TestIsolation:
    def test_one_failing_group_does_not_stop_the_others(
        self, data, column_spec, monkeypatch
    ):
        exp, ann = data

        def explode(*args, **kwargs):
            raise RuntimeError("synthetic group failure")

        _patch_group(monkeypatch, "E", explode)
        result = compute_metrics(exp, ann, column_spec, groups=frozenset("EK"))

        assert "synthetic group failure" in result["error_E"]
        assert result["n_genes_noNA"] == 60, "group K must still have run"

    def test_a_pca_failure_skips_only_the_dependent_groups(
        self, data, column_spec, monkeypatch
    ):
        from combobatch.metrics import embeddings

        exp, ann = data

        def explode(*args, **kwargs):
            raise ValueError("synthetic PCA failure")

        monkeypatch.setattr(embeddings, "compute_pca", explode)
        result = compute_metrics(exp, ann, column_spec, groups=frozenset("AJK"))

        assert result["error_A"] == "skipped: PCA failed"
        assert result["error_J"] == "skipped: PCA failed"
        assert result["n_genes_noNA"] == 60, "group K needs no embedding"

    def test_group_l_without_a_reference_is_an_honest_skip(self, data, column_spec):
        exp, ann = data
        result = compute_metrics(exp, ann, column_spec, groups=frozenset("LK"))
        assert "no pre-harmonization reference" in result["error_L"]
        assert result["n_genes_noNA"] == 60

    def test_a_missing_backend_is_a_named_skip_not_a_failure(self, data, column_spec):
        exp, ann = data
        result = compute_metrics(
            exp,
            ann,
            column_spec,
            groups=frozenset("F"),
            spec=MetricsSpec(skip_slow=False),
        )
        if not METRIC_GROUP_REGISTRY["F"].is_available():
            assert result["error_F"].startswith("skipped: needs")
        else:  # pragma: no cover - only inside the image
            assert "error_F" not in result

    def test_misaligned_frames_raise_before_anything_runs(self, data, column_spec):
        exp, ann = data
        with pytest.raises(ValueError, match="must be aligned"):
            compute_metrics(exp, ann.iloc[:10], column_spec, groups=frozenset("K"))


class TestFlushing:
    def test_the_callback_fires_after_every_group(self, data, column_spec):
        exp, ann = data
        seen: list[set[str]] = []
        compute_metrics(
            exp,
            ann,
            column_spec,
            groups=frozenset("EK"),
            on_group_done=lambda partial: seen.append(set(partial)),
        )
        assert len(seen) == 2
        # Each flush is a superset of the last: partial results accumulate.
        assert seen[0] < seen[1]

    def test_a_failing_callback_never_aborts_the_metrics(self, data, column_spec):
        exp, ann = data

        def explode(_partial):
            raise RuntimeError("synthetic persistence failure")

        result = compute_metrics(
            exp, ann, column_spec, groups=frozenset("EK"), on_group_done=explode
        )
        assert result["n_samples"] == 50
        assert result["n_genes_noNA"] == 60

    def test_a_flush_happens_before_a_later_group_fails(
        self, data, column_spec, monkeypatch
    ):
        exp, ann = data
        flushed: list[dict] = []

        def explode(*args, **kwargs):
            raise RuntimeError("synthetic failure")

        _patch_group(monkeypatch, "K", explode)
        compute_metrics(
            exp,
            ann,
            column_spec,
            groups=frozenset("EK"),
            on_group_done=lambda partial: flushed.append(partial),
        )
        # E finished and was persisted before K blew up.
        assert any("n_samples" in partial for partial in flushed)


class TestReporting:
    def test_the_record_names_which_groups_ran(self, data, column_spec):
        exp, ann = data
        result = compute_metrics(exp, ann, column_spec, groups=frozenset("EKJ"))
        assert set(result["metrics_groups_run"]) == {"E", "J", "K"}

    def test_the_result_is_nan_safe_json(self, data, column_spec):
        from combobatch.dataio import dumps_json

        exp, ann = data
        result = compute_metrics(
            exp,
            ann,
            ColumnSpec(batch="batch", bio="bio", metric_batch_cols=("nowhere",)),
            groups=frozenset("A"),
        )
        assert np.isnan(result["r2_nowhere"])
        assert "NaN" not in dumps_json(result)
        assert json.loads(dumps_json(result))["r2_nowhere"] is None


class TestTuningKnobs:
    def test_the_permutation_count_reaches_group_n(self, column_spec):
        exp, ann = make_expression(n_samples=90, n_genes=40)
        result = compute_metrics(
            exp,
            ann,
            column_spec,
            spec=MetricsSpec(n_perm=3, min_test_n=10, n_pcs=5, skip_slow=False),
            groups=frozenset("N"),
        )
        assert result["pv_n_perm"] == 3

    def test_the_dsc_permutation_count_reaches_group_a(self, data, column_spec):
        exp, ann = data
        result = compute_metrics(
            exp,
            ann,
            column_spec,
            spec=MetricsSpec(dsc_permutations=9),
            groups=frozenset("A"),
        )
        # With n permutations every p-value is a multiple of 1/(n+1), which pins the
        # count actually used rather than merely that some number was used.
        assert (result["dsc_pvalue_batch"] * 10) == pytest.approx(
            round(result["dsc_pvalue_batch"] * 10)
        )
        assert result["dsc_pvalue_batch"] >= 1 / 10
