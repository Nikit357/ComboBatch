"""The command-line surface: configuration assembly, listings, exit codes.

The listings exist to prove that no CLI choice list is hand-maintained — every one of
them is read out of a registry — and the configuration tests pin the rule that a flag
overrides a config file rather than the other way round.
"""

from __future__ import annotations

import argparse
import json

import pytest

from combobatch.cli import main as main_mod
from combobatch.cli import run_cmd
from combobatch.config import ConfigError

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


def _parse(argv: list[str]) -> argparse.Namespace:
    """Parse a ``run`` argument list without executing anything."""
    parser = argparse.ArgumentParser()
    run_cmd.add_arguments(parser)
    return parser.parse_args(argv)


class TestBuildConfig:
    def test_flags_alone_are_enough(self, dataset):
        cfg = run_cmd.build_config(
            _parse([*dataset, "--imputation", "strict", "--method", "01_raw"])
        )
        assert cfg.columns.batch == "batch"
        assert [sel.name for sel in cfg.methods] == ["01_raw"]

    def test_missing_input_is_named(self):
        with pytest.raises(ConfigError, match="--exp"):
            run_cmd.build_config(
                _parse(
                    [
                        "--batch-col",
                        "b",
                        "--bio-col",
                        "d",
                        "--method",
                        "01_raw",
                        "--imputation",
                        "strict",
                    ]
                )
            )

    def test_missing_columns_are_named(self, dataset):
        argv = [flag for flag in dataset if flag not in {"--bio-col", "bio"}]
        with pytest.raises(ConfigError, match="--bio-col"):
            run_cmd.build_config(
                _parse([*argv, "--imputation", "strict", "--method", "01_raw"])
            )

    def test_flags_override_a_config_file(self, dataset, tmp_path):
        import yaml

        config_path = tmp_path / "cfg.yaml"
        config_path.write_text(
            yaml.safe_dump(
                {
                    "input": {"expression": "a.tsv", "annotation": "b.csv"},
                    "output": "./out",
                    "columns": {"batch": "from_file", "bio": "bio"},
                    "imputations": ["strict"],
                    "methods": ["01_raw"],
                }
            )
        )
        cfg = run_cmd.build_config(_parse(["--config", str(config_path), *dataset]))
        assert cfg.columns.batch == "batch", "the flag must win"
        assert cfg.exp_uri.endswith("exp.tsv.gz")

    def test_config_file_alone_is_enough(self, dataset, tmp_path):
        import yaml

        config_path = tmp_path / "cfg.yaml"
        config_path.write_text(
            yaml.safe_dump(
                {
                    "input": {
                        "expression": str(tmp_path / "exp.tsv.gz"),
                        "annotation": str(tmp_path / "ann.csv"),
                    },
                    "output": str(tmp_path / "out"),
                    "columns": {"batch": "batch", "bio": "bio"},
                    "imputations": ["strict", "knn"],
                    "methods": ["01_raw"],
                }
            )
        )
        cfg = run_cmd.build_config(_parse(["--config", str(config_path)]))
        assert len(cfg.combinations()) == 2


class TestHyperparameterFlags:
    def test_a_declared_parameter_is_carried_through(self, dataset):
        cfg = run_cmd.build_config(
            _parse(
                [
                    *dataset,
                    "--imputation",
                    "strict",
                    "--method",
                    "10_mnn",
                    "--method-params",
                    "10_mnn:k=50",
                ]
            )
        )
        assert cfg.methods[0].params == {"k": 50}

    def test_an_undeclared_parameter_is_rejected_at_parse_time(self, dataset):
        # Not at job time: a typo must never reach an output filename.
        with pytest.raises(ConfigError, match="unknown parameter"):
            run_cmd.build_config(
                _parse(
                    [
                        *dataset,
                        "--imputation",
                        "strict",
                        "--method",
                        "10_mnn",
                        "--method-params",
                        "10_mnn:kk=50",
                    ]
                )
            )

    def test_a_did_you_mean_suggestion_is_offered(self, dataset):
        with pytest.raises(ConfigError, match="did you mean 'k'"):
            run_cmd.build_config(
                _parse(
                    [
                        *dataset,
                        "--imputation",
                        "strict",
                        "--method",
                        "10_mnn",
                        "--method-params",
                        "10_mnn:kk=50",
                    ]
                )
            )

    def test_declared_type_beats_inference(self, dataset):
        # target_group is a string parameter whose value can look like an integer; only
        # the declaration knows that.
        cfg = run_cmd.build_config(
            _parse(
                [
                    *dataset,
                    "--imputation",
                    "strict",
                    "--method",
                    "13_fsmvn",
                    "--method-params",
                    "13_fsmvn:target_group=2024",
                ]
            )
        )
        assert cfg.methods[0].params["target_group"] == "2024"

    def test_imputer_parameters_are_parsed_and_validated(self, dataset):
        cfg = run_cmd.build_config(
            _parse(
                [
                    *dataset,
                    "--imputation",
                    "knn",
                    "--method",
                    "01_raw",
                    "--impute-params",
                    "knn:knn_k=10",
                ]
            )
        )
        assert cfg.imputations[0].params == {"knn_k": 10}

    def test_max_na_frac_flag_reaches_the_selection(self, dataset):
        cfg = run_cmd.build_config(
            _parse(
                [
                    *dataset,
                    "--imputation",
                    "knn",
                    "--method",
                    "01_raw",
                    "--max-na-frac",
                    "0.35",
                ]
            )
        )
        assert cfg.imputations[0].max_na_frac == pytest.approx(0.35)

    def test_unset_max_na_frac_defers_to_the_imputer(self, dataset):
        cfg = run_cmd.build_config(
            _parse([*dataset, "--imputation", "strict", "--method", "01_raw"])
        )
        assert cfg.imputations[0].max_na_frac is None


class TestListings:
    def test_list_methods_covers_the_whole_registry(self, capsys):
        from combobatch.methods import METHOD_REGISTRY

        assert main_mod.main(["list-methods"]) == 0
        out = capsys.readouterr().out
        assert f"({len(METHOD_REGISTRY)})" in out
        for key in METHOD_REGISTRY:
            assert key in out

    def test_list_methods_json_carries_hyperparameters(self, capsys):
        assert main_mod.main(["list-methods", "--json"]) == 0
        payload = json.loads(capsys.readouterr().out)
        mnn = next(row for row in payload["entries"] if row["key"] == "10_mnn")
        assert "k" in mnn["hyperparams"]

    def test_list_methods_can_show_why_five_are_absent(self, capsys):
        from combobatch.methods import EXCLUDED_METHODS

        assert main_mod.main(["list-methods", "--show-excluded"]) == 0
        out = capsys.readouterr().out
        for key in EXCLUDED_METHODS:
            assert key in out
        assert "licence-restricted" in out

    def test_available_only_hides_what_cannot_run_here(self, capsys):
        assert main_mod.main(["list-methods", "--available-only", "--json"]) == 0
        payload = json.loads(capsys.readouterr().out)
        assert all(row["backends"] == "ready" for row in payload["entries"])

    def test_list_imputers_covers_the_whole_registry(self, capsys):
        from combobatch.imputation import IMPUTER_REGISTRY

        assert main_mod.main(["list-imputers"]) == 0
        out = capsys.readouterr().out
        for key in IMPUTER_REGISTRY:
            assert key in out
        assert "max_na_frac" in out


class TestSelftest:
    """`combobatch selftest` — the image's acceptance gate."""

    def test_nothing_is_pending_any_more(self):
        # A stale entry here would tell a user a working command does not exist.
        assert main_mod.PENDING == {}

    def test_it_exits_zero_when_nothing_fails(self, capsys):
        assert main_mod.main(["selftest"]) == 0
        out = capsys.readouterr().out
        assert "RESULT:" in out
        assert " failed" in out

    def test_every_registered_method_and_imputer_is_reported(self, capsys):
        from combobatch.imputation import IMPUTER_REGISTRY
        from combobatch.methods import METHOD_REGISTRY

        main_mod.main(["selftest"])
        out = capsys.readouterr().out
        for key in list(METHOD_REGISTRY) + list(IMPUTER_REGISTRY):
            assert key in out, f"{key} is missing from the selftest report"

    def test_a_missing_backend_is_a_named_skip(self, capsys):
        main_mod.main(["selftest"])
        out = capsys.readouterr().out
        # Honest about *why*, not a bare SKIP: on a laptop most methods cannot run, and
        # the report has to distinguish that from a method that ran and misbehaved.
        assert "SKIP" in out
        assert "needs" in out

    def test_json_output_is_machine_readable(self, capsys):
        assert main_mod.main(["selftest", "--json"]) == 0
        payload = json.loads(capsys.readouterr().out)
        assert payload["n_failed"] == 0
        assert {row["name"] for row in payload["imputers"]} == {
            "strict",
            "knn",
            "softimpute",
            "missforest",
        }

    def test_requiring_an_absent_backend_turns_skips_into_failures(self, capsys):
        from combobatch.methods import rinterop

        if rinterop.r_available():  # pragma: no cover - only inside the image
            pytest.skip("R is present, so --with-r cannot fail here")
        assert main_mod.main(["selftest", "--with-r"]) == 1

    def test_an_rpy2_conversion_error_is_a_failure_not_a_skip(self):
        from combobatch.cli.selftest_cmd import _is_rpy2_conversion_error

        # rpy2 raises NotImplementedError for an unconvertible object, which is
        # indistinguishable by type from a method's deliberate "not applicable here".
        # Treating it as a skip is how a broken image passes its own gate.
        assert _is_rpy2_conversion_error(
            NotImplementedError("Conversion 'py2rpy' not defined for objects of type")
        )
        assert not _is_rpy2_conversion_error(
            NotImplementedError("requires RNA-seq data; this batch is microarray")
        )

    def test_the_synthetic_fixture_leaves_the_imputers_work_to_do(self):
        from combobatch.cli.selftest_cmd import make_synthetic, punch_holes

        exp_df, ann_df = make_synthetic()
        holed = punch_holes(exp_df)
        assert holed.isna().any().any(), "nothing to impute"
        # Holes are confined to a subset of genes on purpose: spread uniformly at any
        # real rate they leave every gene with an NA, `strict` correctly keeps nothing,
        # and the whole imputation axis becomes a degenerate test.
        assert (~holed.isna().any()).sum() > 0, "strict would keep no genes"
        assert len(ann_df) == len(exp_df)


class TestRunExitCodes:
    def test_a_successful_run_exits_zero(self, dataset):
        assert (
            run_cmd.main([*dataset, "--imputation", "strict", "--method", "01_raw"])
            == 0
        )

    def test_a_bad_configuration_exits_three(self, dataset, capsys):
        assert (
            run_cmd.main([*dataset, "--imputation", "strict", "--method", "nope"]) == 3
        )

    def test_insufficient_memory_exits_two(self, dataset, monkeypatch):
        from combobatch import memory

        monkeypatch.setattr(memory, "check_memory", lambda _limit: False)
        code = run_cmd.main(
            [
                *dataset,
                "--imputation",
                "strict",
                "--method",
                "01_raw",
                "--memory-limit-gb",
                "9999",
            ]
        )
        assert code == run_cmd.EXIT_INSUFFICIENT_MEMORY

    def test_out_json_is_written_for_the_dispatcher(self, dataset, tmp_path):
        target = tmp_path / "rows.json"
        run_cmd.main(
            [
                *dataset,
                "--imputation",
                "strict",
                "--method",
                "01_raw",
                "--out-json",
                str(target),
            ]
        )
        rows = json.loads(target.read_text())
        assert rows[0]["status"] == "ok"

    def test_no_manifest_suppresses_the_manifest(self, dataset, tmp_path):
        run_cmd.main(
            [*dataset, "--imputation", "strict", "--method", "01_raw", "--no-manifest"]
        )
        assert not (tmp_path / "out" / "run_manifest.json").exists()

    def test_a_standalone_run_writes_the_manifest(self, dataset, tmp_path):
        run_cmd.main([*dataset, "--imputation", "strict", "--method", "01_raw"])
        assert (tmp_path / "out" / "run_manifest.json").exists()
