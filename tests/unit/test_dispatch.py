"""The cross-product dispatcher.

The pure helpers are tested directly; one real six-job run at the end proves the
ThreadPoolExecutor-around-Popen shape actually works, since that shape is the reason
every rpy2 job gets its own R session.
"""

from __future__ import annotations

import json

import pytest

from combobatch import storage
from combobatch.cli import dispatch_cmd
from combobatch.config import ConfigError, ImputationSelection, MethodSelection
from combobatch.imputation import IMPUTER_REGISTRY
from combobatch.methods import METHOD_REGISTRY

from tests.conftest import make_expression


@pytest.fixture
def dataset(tmp_path):
    """A dataset on disk, returned as the flag list that names it."""
    exp, ann = make_expression(n_samples=20, n_genes=30)
    exp.to_csv(tmp_path / "exp.tsv.gz", sep="\t")
    ann.to_csv(tmp_path / "ann.csv")
    return [
        "--exp",
        str(tmp_path / "exp.tsv.gz"),
        "--ann",
        str(tmp_path / "ann.csv"),
        "--out",
        str(tmp_path / "out"),
        "--batch-col",
        "batch",
        "--bio-col",
        "bio",
    ]


class TestExpandAxis:
    def test_all_expands_from_the_registry(self):
        assert dispatch_cmd.expand_axis("all", METHOD_REGISTRY, "method") == list(
            METHOD_REGISTRY
        )

    def test_all_is_case_insensitive(self):
        assert len(dispatch_cmd.expand_axis("ALL", IMPUTER_REGISTRY, "imputer")) == 4

    def test_a_list_is_taken_literally(self):
        assert dispatch_cmd.expand_axis("01_raw,10_mnn", METHOD_REGISTRY, "method") == [
            "01_raw",
            "10_mnn",
        ]

    def test_an_unknown_name_is_rejected_with_a_suggestion(self):
        with pytest.raises(ConfigError, match="did you mean '10_mnn'"):
            dispatch_cmd.expand_axis("10_mnnn", METHOD_REGISTRY, "method")

    def test_nothing_requested_yields_nothing(self):
        assert dispatch_cmd.expand_axis(None, METHOD_REGISTRY, "method") == []


class TestWorkerCommand:
    def _command(self, tmp_path, **overrides):
        from combobatch.cli import run_cmd

        import argparse

        parser = argparse.ArgumentParser()
        dispatch_cmd.add_arguments(parser)
        args = parser.parse_args(
            [
                "--exp",
                "e.tsv",
                "--ann",
                "a.csv",
                "--out",
                str(tmp_path),
                "--batch-col",
                "b",
                "--bio-col",
                "d",
                "--imputations",
                "knn",
                "--methods",
                "10_mnn",
                "--method-params",
                "10_mnn:k=50",
                "--impute-params",
                "knn:knn_k=7",
            ]
        )
        cfg = run_cmd.build_config(args)
        return dispatch_cmd._worker_command(
            cfg,
            cfg.imputations[0],
            cfg.methods[0],
            config_path=None,
            out_json=str(tmp_path / "rows.json"),
            memory_limit_gb=8.0,
            skip_if_exists=False,
            log_level="INFO",
            **overrides,
        )

    def test_invokes_the_module_not_the_console_script(self, tmp_path):
        command = self._command(tmp_path)
        assert command[1:4] == ["-m", "combobatch", "run"]

    def test_suppresses_the_worker_manifest(self, tmp_path):
        # The dispatcher writes it once up front; concurrent workers must not race to
        # read-modify-write the same file.
        assert "--no-manifest" in self._command(tmp_path)

    def test_halves_the_memory_limit_for_the_worker(self, tmp_path):
        command = self._command(tmp_path)
        assert command[command.index("--memory-limit-gb") + 1] == "4.0"

    def test_parameters_survive_the_round_trip(self, tmp_path):
        from combobatch.params import parse_params_arg

        command = self._command(tmp_path)
        method_spec = command[command.index("--method-params") + 1]
        impute_spec = command[command.index("--impute-params") + 1]

        assert parse_params_arg(method_spec) == {"10_mnn": {"k": 50}}
        assert parse_params_arg(impute_spec) == {"knn": {"knn_k": 7}}

    def test_list_valued_parameters_survive_the_round_trip(self):
        from combobatch.params import parse_params_arg

        encoded = dispatch_cmd._params_flag("09_ruv", {"control_genes": ["A", "B"]})
        assert parse_params_arg(encoded) == {"09_ruv": {"control_genes": ["A", "B"]}}


class TestFailedLog:
    def test_round_trips(self, tmp_path):
        backend = storage.backend_for_root(str(tmp_path))
        dispatch_cmd.write_failed_log(backend, "failed.txt", {"a__b", "c__d"})
        assert dispatch_cmd.read_failed_log(backend, "failed.txt") == {"a__b", "c__d"}

    def test_an_absent_log_is_empty_not_an_error(self, tmp_path):
        backend = storage.backend_for_root(str(tmp_path))
        assert dispatch_cmd.read_failed_log(backend, "nothing.txt") == set()

    def test_an_empty_set_removes_the_file(self, tmp_path):
        backend = storage.backend_for_root(str(tmp_path))
        dispatch_cmd.write_failed_log(backend, "failed.txt", {"a__b"})
        dispatch_cmd.write_failed_log(backend, "failed.txt", set())
        assert not backend.exists("failed.txt")


class TestJobKey:
    def test_matches_the_output_key_without_the_post_segment(self):
        key = dispatch_cmd.job_key(
            ImputationSelection(name="knn"),
            MethodSelection(name="10_mnn", params={"k": 5}),
        )
        assert key == "knn__10_mnn__k5"


class TestJobMatrix:
    """Config validation compares what the YAML says; this compares what gets written."""

    def _config(self, methods):
        from combobatch.config import ColumnSpec, RunConfig

        return RunConfig(
            exp_uri="exp.tsv.gz",
            ann_uri="ann.csv",
            out_uri="out/",
            columns=ColumnSpec(batch="batch", bio="bio"),
            imputations=(ImputationSelection(name="strict"),),
            methods=tuple(methods),
        )

    def test_a_parameter_sweep_produces_distinct_keys(self):
        jobs = dispatch_cmd.build_job_matrix(
            self._config(
                [
                    MethodSelection(name="10_mnn"),
                    MethodSelection(name="10_mnn", params={"k": 50}),
                ]
            )
        )
        assert [dispatch_cmd.job_key(*pair) for pair in jobs] == [
            "strict__10_mnn",
            "strict__10_mnn__k50",
        ]

    def test_a_parameter_restated_at_its_default_is_caught_as_a_collision(self):
        # This one slips past config validation: `{"k": 20}` produces a tag there, but
        # k=20 *is* the default, so no tag survives into the output key and the second
        # job would overwrite the first.
        assert METHOD_REGISTRY["10_mnn"].hyperparams["k"].default == 20
        with pytest.raises(ConfigError, match="same output key"):
            dispatch_cmd.build_job_matrix(
                self._config(
                    [
                        MethodSelection(name="10_mnn"),
                        MethodSelection(name="10_mnn", params={"k": 20}),
                    ]
                )
            )


class TestDryRun:
    def test_resolves_everything_and_runs_nothing(self, dataset, capsys, tmp_path):
        code = dispatch_cmd.main(
            [
                *dataset,
                "--imputations",
                "strict,knn",
                "--methods",
                "01_raw,17_quantile",
                "--dry-run",
            ]
        )
        out = capsys.readouterr().out

        assert code == 0
        assert "Jobs: 4" in out
        assert "exp/knn__17_quantile.tsv.gz" in out
        assert not (tmp_path / "out").exists(), "--dry-run must write nothing at all"

    def test_counts_both_post_variants(self, dataset, capsys):
        dispatch_cmd.main(
            [
                *dataset,
                "--imputations",
                "strict",
                "--methods",
                "01_raw",
                "--post-removal",
                "--dry-run",
            ]
        )
        out = capsys.readouterr().out
        assert "Outputs: 2 matrices" in out
        assert "exp/strict__01_raw__post0.tsv.gz" in out

    def test_names_methods_that_cannot_run_here(self, dataset, capsys):
        dispatch_cmd.main(
            [*dataset, "--imputations", "strict", "--methods", "03_limma", "--dry-run"]
        )
        assert "needs R" in capsys.readouterr().out


class TestPreflight:
    def test_reports_missing_backends_before_anything_launches(self, dataset):
        from combobatch.cli import run_cmd
        import argparse

        parser = argparse.ArgumentParser()
        dispatch_cmd.add_arguments(parser)
        args = parser.parse_args(
            [*dataset, "--imputations", "strict", "--methods", "01_raw,03_limma"]
        )
        unavailable = dispatch_cmd.preflight(run_cmd.build_config(args))

        assert "01_raw" not in unavailable
        assert unavailable["03_limma"] == ["R (rpy2)"]


class TestNestedParallelismWarning:
    def test_shambhala_inner_workers_multiply(self, dataset, monkeypatch, caplog):
        import logging

        from combobatch.cli import run_cmd
        from combobatch.methods import shambhala_method
        import argparse

        monkeypatch.setattr(shambhala_method.os, "cpu_count", lambda: 4)
        parser = argparse.ArgumentParser()
        dispatch_cmd.add_arguments(parser)
        args = parser.parse_args(
            [
                *dataset,
                "--imputations",
                "strict",
                "--methods",
                "20_shambhala",
                "--method-params",
                "20_shambhala:n_workers=4",
            ]
        )
        cfg = run_cmd.build_config(args)

        with caplog.at_level(logging.WARNING):
            dispatch_cmd.warn_about_nesting(cfg, n_workers=4)

        assert any("16 Octave processes" in record.message for record in caplog.records)


class TestEndToEnd:
    """One real dispatch, with real subprocesses. Slow enough to keep to one case."""

    def test_six_jobs_produce_the_full_output_layout(self, dataset, tmp_path):
        code = dispatch_cmd.main(
            [
                *dataset,
                "--imputations",
                "strict,knn",
                "--methods",
                "01_raw,02_median_scaling,17_quantile",
                "--n-workers",
                "3",
                "--run-id",
                "t",
                "--memory-limit-gb",
                "0.1",
            ]
        )
        out_root = tmp_path / "out"

        assert code == 0
        assert len(list((out_root / "exp").glob("*.tsv.gz"))) == 6
        assert len(list((out_root / "metrics").glob("*.json"))) == 6
        assert len(list((out_root / "genes").glob("*.json"))) == 6
        # Two imputers, three methods each: one prepared pair per imputer, shared.
        assert len(list((out_root / "prepared").glob("*__exp.tsv.gz"))) == 2
        assert len(list((out_root / "prepared").glob("*__ann.tsv.gz"))) == 2

        manifest = json.loads((out_root / "run_manifest.json").read_text())
        assert len(manifest["outputs"]) == 6
        assert not (out_root / "failed_jobs_t.txt").exists()

    def test_rerunning_with_skip_if_exists_reports_everything_cached(
        self, dataset, tmp_path
    ):
        argv = [
            *dataset,
            "--imputations",
            "strict",
            "--methods",
            "01_raw",
            "--run-id",
            "t",
            "--memory-limit-gb",
            "0.1",
        ]
        dispatch_cmd.main(argv)
        dispatch_cmd.main([*argv, "--skip-if-exists"])

        rows = json.loads((tmp_path / "out" / "dispatch_rows_t.json").read_text())
        assert [row["status"] for row in rows] == ["cached"]
