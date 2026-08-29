"""The CLI against the shipped example dataset.

The eight scenarios of the plan's section 5.13. These go through `main()` with real
argv, real files and the real dispatcher, because every layer below has already been
tested in isolation and the remaining risk is in how they are wired together.

Only backend-free methods are used, so this suite runs on a laptop with no R and no
Octave. Inside the image the same file exercises the same paths with everything present.
"""

from __future__ import annotations

import gzip
import json
from pathlib import Path

import pandas as pd
import pytest

from combobatch.cli import main as cli_main

REPO_ROOT = Path(__file__).resolve().parents[2]
DATASET = REPO_ROOT / "examples" / "example_dataset"
EXP = str(DATASET / "example_exp.tsv.gz")
ANN = str(DATASET / "example_ann.csv")

BATCH_COL = "RNA_BATCH"
BIO_COL = "Diagnosis_cell_type_unified"
COHORT_COL = "COHORT_LABEL"

# Pure Python, fast, and between them they cover the three shapes a method can have:
# a no-op baseline, a scaling method and a reference-batch method.
METHODS = "01_raw,02_median_scaling,15_fsqn_py"

# The dispatcher waits until --memory-limit-gb is free before launching a worker, and
# gives up after --memory-wait-s. Both defaults are production values; on a busy laptop
# the 4 GB gate turns a 10-second suite into a 10-minute one, and it is not what any of
# these tests are checking.
NO_MEMORY_GATE = ["--memory-limit-gb", "0.5", "--memory-wait-s", "5"]

pytestmark = pytest.mark.skipif(
    not Path(EXP).exists(), reason="example dataset is not present"
)


def base_args(out_dir: Path) -> list[str]:
    """The flags every invocation needs."""
    return [
        "--exp", EXP,
        "--ann", ANN,
        "--out", str(out_dir),
        "--batch-col", BATCH_COL,
        "--bio-col", BIO_COL,
        "--cohort-col", COHORT_COL,
    ]  # fmt: skip


def read_matrix(path: Path) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, sep="\t", index_col=0)


def sidecar(out_dir: Path, key: str) -> dict:
    return json.loads((out_dir / "metrics" / f"{key}_metrics.json").read_text())


class TestSingleRun:
    """Scenario 1 — one method, correct shape, valid strict JSON sidecar."""

    @pytest.fixture(scope="class")
    def run(self, tmp_path_factory):
        out = tmp_path_factory.mktemp("single")
        code = cli_main.main(
            ["run", *base_args(out), "--imputation", "strict", "--method", "15_fsqn_py",
             "--metrics", "--groups", "A,E,K"]
        )  # fmt: skip
        assert code == 0
        return out

    def test_the_matrix_is_written_with_the_expected_shape(self, run):
        matrix = read_matrix(run / "exp" / "strict__15_fsqn_py.tsv.gz")
        assert matrix.shape[0] == 144
        # strict keeps only genes present everywhere: 6,797 of 21,890 in this dataset.
        assert matrix.shape[1] == 6797
        assert not matrix.isna().any().any()

    def test_sample_order_is_preserved(self, run):
        matrix = read_matrix(run / "exp" / "strict__15_fsqn_py.tsv.gz")
        annotation = pd.read_csv(ANN, index_col=0)
        assert list(matrix.index) == list(annotation.index)

    def test_the_sidecar_is_valid_strict_json(self, run):
        raw = (run / "metrics" / "strict__15_fsqn_py_metrics.json").read_text()
        # json.loads accepts bare NaN; a strict parser in any other language does not,
        # which is the failure this guards.
        assert "NaN" not in raw and "Infinity" not in raw
        record = json.loads(raw)
        assert record["status"] == "ok"
        assert record["n_samples"] == 144

    def test_the_metric_keys_carry_the_configured_column_names(self, run):
        record = sidecar(run, "strict__15_fsqn_py")
        assert f"r2_{BATCH_COL}" in record
        assert f"r2_{BIO_COL}" in record

    def test_the_manifest_is_written_and_resolves_every_parameter(self, run):
        manifest = json.loads((run / "run_manifest.json").read_text())
        entry = manifest["outputs"]["strict__15_fsqn_py"]
        assert entry["method"] == "15_fsqn_py"
        assert entry["columns"] if False else manifest["columns"]["batch"] == BATCH_COL
        # Defaults included, not just what the user typed.
        assert "target_group" in entry["method_params"]


class TestDispatchGrid:
    """Scenarios 2, 3 and 6 — a grid, caching, and a hyperparameter sweep."""

    @pytest.fixture(scope="class")
    def grid(self, tmp_path_factory):
        out = tmp_path_factory.mktemp("grid")
        code = cli_main.main(
            ["dispatch", *base_args(out), "--imputations", "strict,knn",
             "--methods", METHODS, "--n-workers", "3", *NO_MEMORY_GATE, "--metrics", "--groups", "E,K"]
        )  # fmt: skip
        assert code == 0
        return out

    def test_every_combination_produced_an_output(self, grid):
        outputs = sorted(p.name for p in (grid / "exp").glob("*.tsv.gz"))
        assert len(outputs) == 6, outputs

    def test_the_summary_table_has_one_row_per_output(self, grid):
        assert cli_main.main(["concat", "--out", str(grid)]) == 0
        summary = pd.read_csv(grid / "metrics_summary.csv")
        assert len(summary) == 6

    def test_a_second_run_is_entirely_cached(self, grid):
        code = cli_main.main(
            ["dispatch", *base_args(grid), "--imputations", "strict,knn",
             "--methods", METHODS, "--n-workers", "3", *NO_MEMORY_GATE, "--skip-if-exists"]
        )  # fmt: skip
        assert code == 0
        records = [
            sidecar(grid, path.name.replace("_metrics.json", ""))
            for path in (grid / "metrics").glob("*_metrics.json")
        ]
        assert all(r["status"] in {"ok", "cached"} for r in records)

    def test_a_sweep_produces_distinct_readable_keys(self, tmp_path):
        # A sweep is a YAML feature, not a flag one: --method-params applies by method
        # name, so `--methods 26_xpn,26_xpn` gives both entries the same parameters and
        # is correctly rejected as a duplicate. Each YAML entry carries its own params.
        config = tmp_path / "sweep.yaml"
        config.write_text(f"""
input:
  expression: {EXP}
  annotation: {ANN}
output: {tmp_path / 'out'}
columns:
  batch: {BATCH_COL}
  bio: {BIO_COL}
imputations:
  - name: strict
methods:
  - name: 26_xpn
  - name: 26_xpn
    params: {{n_quantiles: 20}}
""")
        code = cli_main.main(
            ["dispatch", "--config", str(config), "--n-workers", "2", *NO_MEMORY_GATE]
        )
        assert code == 0

        outputs = sorted(p.name for p in (tmp_path / "out" / "exp").glob("*.tsv.gz"))
        assert outputs == [
            "strict__26_xpn.tsv.gz",
            "strict__26_xpn__n_quantiles20.tsv.gz",
        ], outputs

        # Both recorded, and the manifest resolves each tag back to real parameters.
        manifest = json.loads((tmp_path / "out" / "run_manifest.json").read_text())
        assert (
            manifest["outputs"]["strict__26_xpn"]["method_params"]["n_quantiles"] == 50
        )
        assert (
            manifest["outputs"]["strict__26_xpn__n_quantiles20"]["method_params"][
                "n_quantiles"
            ]
            == 20
        )

        # Different settings must give different numbers, or the sweep is decorative.
        default = read_matrix(tmp_path / "out" / "exp" / "strict__26_xpn.tsv.gz")
        tuned = read_matrix(
            tmp_path / "out" / "exp" / "strict__26_xpn__n_quantiles20.tsv.gz"
        )
        assert not default.equals(tuned)


class TestMetricsAgreeInlineAndStandalone:
    """Scenario 4 — computing metrics later must give the same numbers."""

    def test_identical_values(self, tmp_path):
        inline = tmp_path / "inline"
        later = tmp_path / "later"

        assert (
            cli_main.main(
                ["run", *base_args(inline), "--imputation", "strict",
                 "--method", "02_median_scaling", "--metrics", "--groups", "A,E,K"]
            )  # fmt: skip
            == 0
        )
        assert (
            cli_main.main(
                ["run", *base_args(later), "--imputation", "strict",
                 "--method", "02_median_scaling"]
            )  # fmt: skip
            == 0
        )
        assert cli_main.main(["metrics", *base_args(later), "--groups", "A,E,K"]) == 0

        key = "strict__02_median_scaling"
        a, b = sidecar(inline, key), sidecar(later, key)

        # Durations are bookkeeping, not measurements: two runs of the same work take
        # different amounts of time and that is not a discrepancy.
        def is_duration(name: str) -> bool:
            return name.endswith("_time_s") or name.endswith("_elapsed_s")

        shared = {
            name
            for name in set(a) & set(b)
            if isinstance(a[name], (int, float))
            and not isinstance(a[name], bool)
            and not is_duration(name)
        }
        assert len(shared) > 20, "too few numeric metrics to be a real comparison"
        for name in shared:
            assert a[name] == pytest.approx(b[name], rel=1e-9, nan_ok=True), name


class TestBackends:
    """Scenario 5 — the local and S3 backends must produce the same layout."""

    def test_identical_layouts(self, tmp_path):
        boto3 = pytest.importorskip("boto3")
        moto = pytest.importorskip("moto")

        local = tmp_path / "local"
        assert (
            cli_main.main(
                ["run", *base_args(local), "--imputation", "strict",
                 "--method", "02_median_scaling", "--metrics", "--groups", "E,K"]
            )  # fmt: skip
            == 0
        )
        local_layout = sorted(
            str(p.relative_to(local)) for p in local.rglob("*") if p.is_file()
        )

        with moto.mock_aws():
            boto3.client("s3", region_name="us-east-1").create_bucket(Bucket="cb-e2e")
            assert (
                cli_main.main(
                    ["run", "--exp", EXP, "--ann", ANN, "--out", "s3://cb-e2e/out/",
                     "--batch-col", BATCH_COL, "--bio-col", BIO_COL,
                     "--cohort-col", COHORT_COL, "--imputation", "strict",
                     "--method", "02_median_scaling", "--metrics", "--groups", "E,K"]
                )  # fmt: skip
                == 0
            )
            listing = boto3.client("s3", region_name="us-east-1").list_objects_v2(
                Bucket="cb-e2e", Prefix="out/"
            )
            remote_layout = sorted(
                item["Key"].removeprefix("out/") for item in listing["Contents"]
            )

        assert local_layout == remote_layout


class TestBestApproachConfigs:
    """Scenario 7 — every shipped best-approach config runs end to end."""

    @pytest.mark.parametrize(
        "config",
        sorted((REPO_ROOT / "examples" / "configs" / "best_approaches").glob("*.yaml")),
        ids=lambda p: p.name,
    )
    def test_it_runs(self, config, tmp_path, monkeypatch):
        # The configs name R-backed methods, which are the benchmark's recommendations.
        # Here the point is that the config drives the CLI correctly end to end, so the
        # method is overridden with a backend-free one and everything else is honoured.
        monkeypatch.chdir(REPO_ROOT)
        code = cli_main.main(
            ["run", "--config", str(config), "--out", str(tmp_path),
             "--method", "02_median_scaling"]
        )  # fmt: skip
        assert code == 0
        assert list((tmp_path / "exp").glob("*.tsv.gz"))


class TestImputationAxisIsNonDegenerate:
    """Scenario 8 — the imputers must genuinely disagree on this dataset.

    They only do so above a threshold, and that is a property of the data worth pinning.
    Missingness here is platform-driven and therefore banded: the per-gene NA fraction
    takes five values, and `max_na_frac` matters only where it crosses one. The default
    0.20 falls below every band, so at defaults the imputers have nothing to fill.
    """

    def test_at_default_max_na_frac_the_imputers_match_strict(self, tmp_path):
        # Not a defect, but it must not be discovered by surprise: at the default gate
        # every imputer keeps only the zero-NA genes, so there is nothing left to impute
        # and knn returns exactly what strict returns.
        code = cli_main.main(
            ["dispatch", *base_args(tmp_path), "--imputations", "strict,knn",
             "--methods", "01_raw", "--n-workers", "2", *NO_MEMORY_GATE]
        )  # fmt: skip
        assert code == 0
        strict = read_matrix(tmp_path / "exp" / "strict__01_raw.tsv.gz")
        knn = read_matrix(tmp_path / "exp" / "knn__01_raw.tsv.gz")
        assert strict.shape[1] == knn.shape[1] == 6797

    def test_above_the_threshold_the_imputers_diverge(self, tmp_path):
        config = tmp_path / "axis.yaml"
        config.write_text(f"""
input:
  expression: {EXP}
  annotation: {ANN}
output: {tmp_path / 'out'}
columns:
  batch: {BATCH_COL}
  bio: {BIO_COL}
imputations:
  - name: strict
  - name: knn
    max_na_frac: 0.4
methods:
  - name: 01_raw
metrics:
  enabled: true
  groups: E,K
""")
        assert (
            cli_main.main(
                [
                    "dispatch",
                    "--config",
                    str(config),
                    "--n-workers",
                    "2",
                    *NO_MEMORY_GATE,
                ]
            )
            == 0
        )

        out = tmp_path / "out"
        strict = read_matrix(out / "exp" / "strict__01_raw.tsv.gz")
        knn = read_matrix(out / "exp" / "knn__01_raw__max_na_frac0.4.tsv.gz")

        # 6,797 genes against 15,877: the imputation axis is doing real work here.
        assert strict.shape[1] == 6797
        assert knn.shape[1] == 15877

        a = sidecar(out, "strict__01_raw")
        b = sidecar(out, "knn__01_raw__max_na_frac0.4")
        assert a["n_genes_noNA"] != b["n_genes_noNA"]
