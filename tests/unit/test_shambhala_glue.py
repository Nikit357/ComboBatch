"""The Shambhala glue function, with Octave stubbed out.

The glue is the part the donor never had as a single function — it existed as four
drifted copies, and one of them is where the un-forwarded ``skip_qn`` bug lives. These
tests pin the behaviour those copies disagreed about.
"""

from __future__ import annotations

import logging

import numpy as np
import pandas as pd
import pytest

from combobatch.methods import shambhala_method
from combobatch.methods.shambhala_method import shambhala_harmonize
from combobatch.vendor.shambhala import OCTAVE_DIR

from tests.conftest import make_expression


@pytest.fixture
def calibration():
    """A small (P, Q) pair sharing the synthetic gene space."""
    p_exp, _ = make_expression(n_samples=12, n_genes=80, seed=7)
    q_exp, _ = make_expression(n_samples=20, n_genes=80, seed=8)
    return p_exp, q_exp


@pytest.fixture
def captured_call(monkeypatch):
    """Replace ``harmonize_parallel`` with a recorder that echoes its input back."""
    from combobatch.vendor.shambhala import parallel as parallel_module

    recorded: dict = {}

    def fake_harmonize_parallel(input_df, p_df, rm, rs, **kwargs):
        recorded["input_df"] = input_df
        recorded["p_df"] = p_df
        recorded["rm"] = rm
        recorded["rs"] = rs
        recorded.update(kwargs)
        return input_df

    monkeypatch.setattr(parallel_module, "harmonize_parallel", fake_harmonize_parallel)
    return recorded


class TestSkipQnForwarding:
    """§6 bug 4 — the job-runner copy of the glue never forwarded ``skip_qn``."""

    def test_precompute_qn_reference_sets_skip_qn(self, calibration, captured_call):
        exp, _ = make_expression()
        p_df, q_df = calibration

        shambhala_harmonize(exp, p_df, q_df, precompute_qn_reference=True)

        assert captured_call["skip_qn"] is True, (
            "Python-side QN was applied but Octave was not told to skip its own pass; "
            "the matrix would be quantile-normalized twice."
        )

    def test_default_run_does_not_skip_octave_qn(self, calibration, captured_call):
        exp, _ = make_expression()
        p_df, q_df = calibration

        shambhala_harmonize(exp, p_df, q_df)

        assert captured_call["skip_qn"] is False

    def test_python_side_qn_actually_changes_the_matrix(
        self, calibration, captured_call
    ):
        exp, _ = make_expression()
        p_df, q_df = calibration

        shambhala_harmonize(exp, p_df, q_df)
        without_qn = captured_call["input_df"].copy()
        shambhala_harmonize(exp, p_df, q_df, precompute_qn_reference=True)
        with_qn = captured_call["input_df"]

        assert not np.allclose(without_qn.values, with_qn.values)


class TestDispatcherDefaults:
    def test_inner_workers_default_to_one(self, calibration, captured_call):
        exp, _ = make_expression()
        shambhala_harmonize(exp, *calibration)
        assert captured_call["n_workers"] == 1

    def test_progress_is_disabled_by_default(self, calibration, captured_call):
        exp, _ = make_expression()
        shambhala_harmonize(exp, *calibration)
        assert captured_call["disable_progress"] is True

    def test_octave_scripts_dir_defaults_to_the_vendored_directory(
        self, calibration, captured_call
    ):
        import pathlib

        exp, _ = make_expression()
        shambhala_harmonize(exp, *calibration)

        scripts_dir = pathlib.Path(captured_call["octave_scripts_dir"])
        assert scripts_dir == pathlib.Path(OCTAVE_DIR)
        assert (scripts_dir / "Shambhala2_piped.m").is_file()


class TestOrientationAndAlignment:
    def test_harmonize_parallel_receives_genes_by_samples(
        self, calibration, captured_call
    ):
        exp, _ = make_expression(n_samples=60, n_genes=80)
        shambhala_harmonize(exp, *calibration)

        # The one place in the codebase where the axes are the other way round.
        assert captured_call["input_df"].shape == (80, 60)

    def test_output_is_samples_by_genes_in_input_order(
        self, calibration, captured_call
    ):
        exp, _ = make_expression()
        result = shambhala_harmonize(exp, *calibration)

        assert list(result.index) == list(exp.index)
        assert result.shape[0] == exp.shape[0]

    def test_genes_absent_from_a_reference_are_dropped(self, captured_call):
        exp, _ = make_expression(n_genes=80)
        p_exp, _ = make_expression(n_samples=12, n_genes=80, seed=7)
        q_exp, _ = make_expression(n_samples=20, n_genes=80, seed=8)
        q_exp = q_exp.drop(columns=["GENE0000", "GENE0001"])

        result = shambhala_harmonize(exp, p_exp, q_exp)

        assert "GENE0000" not in result.columns
        assert result.shape[1] == 78

    def test_p_is_aligned_to_the_surviving_gene_set(self, captured_call):
        exp, _ = make_expression(n_genes=80)
        p_exp, _ = make_expression(n_samples=12, n_genes=80, seed=7)
        q_exp, _ = make_expression(n_samples=20, n_genes=80, seed=8)
        q_exp = q_exp.drop(columns=["GENE0000"])

        shambhala_harmonize(exp, p_exp, q_exp)

        assert list(captured_call["p_df"].index) == list(
            captured_call["input_df"].index
        )


class TestNaHandling:
    def test_dropped_genes_come_back_as_nan_columns(self, calibration, captured_call):
        exp, _ = make_expression()
        exp = exp.copy()
        exp.loc[exp.index[0], "GENE0005"] = np.nan

        result = shambhala_harmonize(exp, *calibration, na_strategy="drop")

        assert "GENE0005" in result.columns
        assert result["GENE0005"].isna().all()
        assert "GENE0005" not in captured_call["input_df"].index

    def test_knn_strategy_imputes_instead_of_dropping(self, calibration, captured_call):
        exp, _ = make_expression()
        exp = exp.copy()
        exp.loc[exp.index[0], "GENE0005"] = np.nan

        shambhala_harmonize(exp, *calibration, na_strategy="knn")

        assert "GENE0005" in captured_call["input_df"].index

    def test_zero_count_q_genes_are_returned_as_nan_not_a_keyerror(
        self, calibration, captured_call
    ):
        # q_pseudocount=0 excludes zero-count genes from the Q statistics, after which
        # the rescale step drops them. Upstream's restoration then indexes by the
        # original gene list and raises KeyError.
        exp, _ = make_expression()
        p_df, q_df = calibration
        q_df = q_df.copy()
        q_df["GENE0003"] = 0.0

        result = shambhala_harmonize(exp, p_df, q_df, q_pseudocount=0.0)

        assert result["GENE0003"].isna().all()
        assert "GENE0003" not in captured_call["input_df"].index


class TestFailsLoudly:
    def test_nan_in_p_raises(self, calibration, captured_call):
        exp, _ = make_expression()
        p_df, q_df = calibration
        p_df = p_df.copy()
        p_df.iloc[0, 0] = np.nan

        with pytest.raises(ValueError, match="P calibration reference has NaN"):
            shambhala_harmonize(exp, p_df, q_df)

    def test_nan_in_q_raises(self, calibration, captured_call):
        exp, _ = make_expression()
        p_df, q_df = calibration
        q_df = q_df.copy()
        q_df.iloc[0, 0] = np.nan

        with pytest.raises(ValueError, match="Q calibration reference has NaN"):
            shambhala_harmonize(exp, p_df, q_df)

    def test_empty_gene_intersection_raises(self, calibration, captured_call):
        exp, _ = make_expression()
        exp = exp.rename(columns=lambda gene: f"OTHER_{gene}")
        p_df, q_df = calibration

        with pytest.raises(ValueError, match="no genes shared"):
            shambhala_harmonize(exp, p_df, q_df)


class TestScaleContract:
    def test_log_scale_input_warns(self, calibration, captured_call, caplog):
        exp, _ = make_expression()
        p_df, q_df = calibration

        with caplog.at_level(logging.WARNING):
            shambhala_harmonize(np.log2(exp + 1), p_df, q_df)

        assert any("log-scale" in record.message for record in caplog.records)

    def test_raw_scale_input_does_not_warn(self, calibration, captured_call, caplog):
        exp, _ = make_expression()
        p_df, q_df = calibration

        with caplog.at_level(logging.WARNING):
            shambhala_harmonize(exp * 100.0, p_df, q_df)

        assert not any("log-scale" in record.message for record in caplog.records)


class TestNestedParallelismWarning:
    def test_warns_when_the_product_exceeds_the_cpu_count(self, monkeypatch):
        monkeypatch.setattr(shambhala_method.os, "cpu_count", lambda: 4)
        message = shambhala_method.warn_nested_parallelism(4, 4)
        assert message is not None
        assert "16 Octave processes" in message

    def test_silent_when_the_product_fits(self, monkeypatch):
        monkeypatch.setattr(shambhala_method.os, "cpu_count", lambda: 16)
        assert shambhala_method.warn_nested_parallelism(8, 1) is None


class TestShippedCalibration:
    def test_default_p_and_q_load_at_the_documented_shapes(self):
        p_df = shambhala_method.load_calibration(shambhala_method.DEFAULT_P)
        q_df = shambhala_method.load_calibration(shambhala_method.DEFAULT_Q)

        assert p_df.shape == (39, 11768)
        assert q_df.shape == (100, 11887)

    def test_shipped_references_are_raw_scale_and_nan_free(self):
        p_df = shambhala_method.load_calibration(shambhala_method.DEFAULT_P)
        q_df = shambhala_method.load_calibration(shambhala_method.DEFAULT_Q)

        # Octave log2s them itself; a log-scale reference would be double-logged.
        assert p_df.values.max() > 1000
        assert q_df.values.max() > 1000
        assert not p_df.isna().any().any()
        assert not q_df.isna().any().any()


class TestRegistryAdapter:
    def test_p_and_q_are_read_from_paths(self, calibration, captured_call, tmp_path):
        exp, ann = make_expression()
        p_df, q_df = calibration
        p_path, q_path = tmp_path / "P.csv.gz", tmp_path / "Q.tsv"
        p_df.to_csv(p_path, compression="gzip")
        q_df.to_csv(q_path, sep="\t")

        result = shambhala_method.normalize_shambhala(
            exp, ann, batch_col="batch", bio_col="bio", P=str(p_path), Q=str(q_path)
        )

        pd.testing.assert_frame_equal(result, shambhala_harmonize(exp, p_df, q_df))

    def test_batch_labels_do_not_affect_the_result(
        self, calibration, captured_call, tmp_path
    ):
        exp, ann = make_expression()
        p_df, q_df = calibration
        p_path, q_path = tmp_path / "P.csv", tmp_path / "Q.csv"
        p_df.to_csv(p_path)
        q_df.to_csv(q_path)
        scrambled = ann.copy()
        scrambled["batch"] = "AllOneBatch"

        # Shambhala harmonizes each sample against fixed references, so batch identity
        # is genuinely irrelevant — the parameter exists only for the common signature.
        with_real_batches = shambhala_method.normalize_shambhala(
            exp, ann, batch_col="batch", P=str(p_path), Q=str(q_path)
        )
        with_scrambled = shambhala_method.normalize_shambhala(
            exp, scrambled, batch_col="batch", P=str(p_path), Q=str(q_path)
        )

        pd.testing.assert_frame_equal(with_real_batches, with_scrambled)
