"""Aggregation: nested detail into long tables, and no silently dropped keys."""

from __future__ import annotations

import pytest

from combobatch import dataio, storage
from combobatch.metrics.concat import (
    concat_metrics,
    drop_remaining_containers,
    split_nested,
    write_tables,
)
from combobatch.params import KeyParseError


def _write_sidecar(root, stem: str, record: dict) -> None:
    backend = storage.backend_for_root(str(root))
    dataio.write_json(record, backend, f"metrics/{stem}_metrics.json")


class TestSplitNested:
    def test_gene_correlations_become_rows(self):
        record = {
            "mk_rho_by_gene": {"GENE1": 0.9, "GENE2": None},
            "mk_rho_n_cohorts_by_gene": {"GENE1": 3, "GENE2": 0},
        }
        genes, _, _ = split_nested(record, "strict__01_raw")
        assert len(genes) == 2
        assert genes[0] == {
            "key": "strict__01_raw",
            "gene": "GENE1",
            "rho": 0.9,
            "n_cohorts": 3,
        }
        assert "mk_rho_by_gene" not in record, "the wide table must never see a dict"

    def test_cohort_correlations_become_rows(self):
        record = {"mk_rho_by_cohort": {"C1": 0.8}}
        _, cohorts, _ = split_nested(record, "k")
        assert cohorts == [{"key": "k", "cohort": "C1", "rho": 0.8}]

    def test_prediction_folds_become_rows_tagged_by_target(self):
        record = {
            "pv_lobo3_folds": {
                "B1": {"n": 20, "n_classes": 2, "f1_macro": 0.7, "auc": 0.8}
            },
            "pv_lobo2_folds": {
                "B1": {"n": 20, "n_classes": 2, "f1_macro": 0.9, "auc": 0.95}
            },
        }
        _, _, folds = split_nested(record, "k")
        assert {row["target"] for row in folds} == {"3class", "2class"}
        assert not any(key.endswith("_folds") for key in record)

    def test_the_per_biology_dict_is_pivoted_into_columns(self):
        record = {"xb_rank_agree_by_bio": {"TypeA": 0.5, "TypeB": 0.6}}
        split_nested(record, "k")
        assert record["xb_rank_agree_TypeA"] == 0.5
        assert record["xb_rank_agree_TypeB"] == 0.6
        assert "xb_rank_agree_by_bio" not in record

    def test_the_gene_cohort_detail_is_folded_into_the_gene_table(self):
        record = {"mk_gene_cohort_detail": {"C1": {"GENE1": 0.9, "GENE2": None}}}
        genes, _, _ = split_nested(record, "k")
        assert len(genes) == 2
        assert all(row["cohort"] == "C1" for row in genes)
        assert "mk_gene_cohort_detail" not in record


class TestDropRemainingContainers:
    def test_surviving_containers_are_stripped_and_reported(self):
        rows = [{"a": 1, "params": {"k": 5}, "notes": ["x"]}]
        offenders = drop_remaining_containers(rows)
        assert offenders == {"params", "notes"}
        assert rows == [{"a": 1}]

    def test_scalar_rows_are_untouched(self):
        rows = [{"a": 1, "b": "text", "c": None}]
        assert drop_remaining_containers(rows) == set()
        assert rows == [{"a": 1, "b": "text", "c": None}]


class TestConcatMetrics:
    def test_builds_a_wide_table_with_metadata_first(self, tmp_path):
        _write_sidecar(tmp_path, "strict__01_raw", {"status": "ok", "n_samples": 10})
        _write_sidecar(tmp_path, "knn__10_mnn", {"status": "ok", "n_samples": 12})

        result = concat_metrics(str(tmp_path))
        assert len(result.summary) == 2
        assert list(result.summary.columns)[:2] == ["key", "status"]
        assert sorted(result.summary["key"]) == ["knn__10_mnn", "strict__01_raw"]

    def test_long_tables_are_empty_unless_the_groups_ran(self, tmp_path):
        _write_sidecar(tmp_path, "strict__01_raw", {"status": "ok"})
        result = concat_metrics(str(tmp_path))
        assert result.gene_correlations.empty
        assert result.prediction_folds.empty
        assert set(result.tables()) == {"metrics_summary"}

    def test_nested_detail_lands_in_its_own_table(self, tmp_path):
        _write_sidecar(
            tmp_path,
            "strict__01_raw",
            {"status": "ok", "mk_rho_by_gene": {"GENE1": 0.9}},
        )
        result = concat_metrics(str(tmp_path))
        assert len(result.gene_correlations) == 1
        assert "mk_rho_by_gene" not in result.summary.columns

    def test_an_unparseable_key_is_a_loud_error(self, tmp_path):
        # The donor dropped these silently, so a malformed key meant a result that
        # vanished from every table without a word.
        _write_sidecar(
            tmp_path, "not__a__valid__output__key__at__all", {"status": "ok"}
        )
        with pytest.raises(KeyParseError, match="will not"):
            concat_metrics(str(tmp_path))

    def test_it_can_be_told_to_continue(self, tmp_path):
        _write_sidecar(
            tmp_path, "not__a__valid__output__key__at__all", {"status": "ok"}
        )
        _write_sidecar(tmp_path, "strict__01_raw", {"status": "ok"})
        result = concat_metrics(str(tmp_path), strict_keys=False)
        assert list(result.summary["key"]) == ["strict__01_raw"]

    def test_an_empty_root_yields_an_empty_summary(self, tmp_path):
        assert concat_metrics(str(tmp_path)).summary.empty


class TestWriteTables:
    def test_writes_to_the_run_root(self, tmp_path):
        _write_sidecar(tmp_path, "strict__01_raw", {"status": "ok", "n_samples": 10})
        result = concat_metrics(str(tmp_path))
        backend = storage.backend_for_root(str(tmp_path))

        written = write_tables(result, backend=backend)
        assert backend.exists("metrics_summary.csv")
        assert any("metrics_summary.csv" in uri for uri in written)

    def test_dated_local_copies_do_not_overwrite_a_snapshot(self, tmp_path):
        _write_sidecar(tmp_path, "strict__01_raw", {"status": "ok"})
        result = concat_metrics(str(tmp_path))
        tables = tmp_path / "tables"

        write_tables(result, out_dir=str(tables), date_tag="260101")
        write_tables(result, out_dir=str(tables), date_tag="260102")
        assert (tables / "metrics_summary_260101.csv").exists()
        assert (tables / "metrics_summary_260102.csv").exists()
