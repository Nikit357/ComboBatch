"""The shipped example dataset and configurations.

Every config here is documentation that runs, which is the only kind that cannot go quietly
stale: a renamed method or a changed schema breaks these files, and the numbers quoted in
the READMEs are asserted rather than remembered.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from combobatch import dataio
from combobatch.config import RunConfig
from combobatch.metrics.panels import load_panel

REPO_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES = REPO_ROOT / "examples"
DATASET = EXAMPLES / "example_dataset"
EXP_PATH = DATASET / "example_exp.tsv.gz"
ANN_PATH = DATASET / "example_ann.csv"
PANEL_PATH = EXAMPLES / "marker_panel_fl.csv"

CONFIGS = sorted(EXAMPLES.glob("configs/**/*.yaml"))


@pytest.fixture(scope="module")
def dataset():
    return dataio.read_table(str(EXP_PATH)), dataio.read_table(str(ANN_PATH))


class TestDataset:
    def test_the_shape_is_what_the_docs_claim(self, dataset):
        exp, ann = dataset
        assert exp.shape == (144, 21890)
        assert ann.shape == (144, 12)

    def test_expression_and_annotation_are_aligned(self, dataset):
        exp, ann = dataset
        assert exp.index.equals(ann.index)

    def test_it_fits_in_the_repository_without_lfs(self):
        # 13.8 MB: under GitHub's 50 MB warning, so no LFS and no download script.
        assert EXP_PATH.stat().st_size < 20_000_000

    def test_the_five_cohorts_are_all_present(self, dataset):
        _, ann = dataset
        assert ann["COHORT_LABEL"].value_counts().to_dict() == {
            "PUB_DLBCL_GSE64555": 40,
            "PUB_FL_GSE148070": 40,
            "GSE69033": 30,
            "PUB_Health_B_cells_GSE119234": 20,
            "PUB_FL_GSE62241": 14,
        }

    def test_the_gse62241_duplicates_were_dropped(self, dataset):
        # Each of these 14 experiments appeared twice: once suffixed with the cohort and
        # fully annotated, once bare and sparser - and the two copies disagreed about the
        # diagnosis for four of them. The suffixed copy is the one kept.
        #
        # Only these ids: GSE119234's own accessions are also bare SRX, and are not
        # duplicates of anything.
        _, ann = dataset
        duplicated = {f"SRX7306{n:02d}" for n in range(0, 13)} | {"SRX730599"}
        assert not duplicated & set(ann.index), "a bare duplicate survived"
        assert {f"{sid}-PUB_FL_GSE62241" for sid in duplicated} <= set(ann.index)

    def test_four_platforms_and_both_preservation_types(self, dataset):
        _, ann = dataset
        assert ann["RNA_BATCH"].nunique() == 4
        batches = set(ann["RNA_BATCH"])
        assert any("_FFPE_" in b for b in batches)
        assert any("_FF_" in b for b in batches)

    def test_the_coverage_is_ragged_which_is_the_point(self, dataset):
        exp, _ = dataset
        present = exp.notna().sum()
        # Ragged, platform-driven missingness is what makes the four imputers give
        # different answers. A dense matrix would make the imputation axis decorative.
        assert int((present == 0).sum()) == 5877
        assert int((present == len(exp)).sum()) == 6797
        assert 0.4 < exp.isna().to_numpy().mean() < 0.45

    def test_the_missingness_is_banded_and_the_default_gate_falls_below_it(
        self, dataset
    ):
        """`max_na_frac` only matters where it crosses a band, and 0.20 crosses none.

        Platform-driven missingness is not scattered: a gene is measured by a whole
        cohort or by none of it. So the per-gene NA fraction takes five values here, and
        the imputer default of 0.20 sits below all of them — at defaults `knn`,
        `softimpute` and `missforest` keep exactly the genes `strict` keeps and have
        nothing left to fill. Every doc that describes this dataset says so; this is what
        keeps that claim true.
        """
        exp, _ = dataset
        na_fraction = exp.isna().mean().round(4)
        bands = sorted(na_fraction.unique())
        assert bands == [0.0, 0.3472, 0.625, 0.6528, 1.0]

        assert int((na_fraction <= 0.20).sum()) == 6797, "the default gate"
        assert int((na_fraction <= 0.40).sum()) == 15877, "above the first band"

    def test_the_values_are_raw_not_log_scale(self, dataset):
        # 20_shambhala applies its own log2(x+1) and has no scale detection, so raw input
        # is a hard requirement rather than a preference.
        exp, _ = dataset
        assert exp.max().max() > 1000

    def test_provenance_is_shipped(self):
        text = (DATASET / "PROVENANCE.md").read_text()
        for accession in ("GSE64555", "GSE148070", "GSE119234", "GSE69033", "GSE62241"):
            assert accession in text


class TestMarkerPanel:
    def test_the_example_panel_loads(self):
        panel = load_panel(str(PANEL_PATH))
        assert len(panel.genes) > 500
        assert panel.housekeeping

    def test_most_of_it_is_measured_in_the_example_matrix(self, dataset):
        exp, _ = dataset
        panel = load_panel(str(PANEL_PATH))
        assert len(set(panel.genes) & set(exp.columns)) > 600


class TestConfigs:
    def test_there_are_configs_to_check(self):
        assert len(CONFIGS) >= 9

    @pytest.mark.parametrize("path", CONFIGS, ids=lambda p: p.name)
    def test_every_config_parses_and_validates(self, path, dataset):
        _, ann = dataset
        cfg = RunConfig.from_yaml(str(path))
        # Against the annotation too, so a renamed column in the example dataset breaks
        # the configs that reference it.
        cfg.validate(ann_df=ann)

    @pytest.mark.parametrize("path", CONFIGS, ids=lambda p: p.name)
    def test_every_config_points_at_the_shipped_dataset(self, path):
        cfg = RunConfig.from_yaml(str(path))
        assert Path(cfg.exp_uri).exists(), f"{path.name}: {cfg.exp_uri} is missing"
        assert Path(cfg.ann_uri).exists(), f"{path.name}: {cfg.ann_uri} is missing"

    @pytest.mark.parametrize("path", CONFIGS, ids=lambda p: p.name)
    def test_a_named_marker_panel_exists(self, path):
        cfg = RunConfig.from_yaml(str(path))
        if cfg.metrics.marker_panel:
            assert Path(cfg.metrics.marker_panel).exists()

    @pytest.mark.parametrize("path", CONFIGS, ids=lambda p: p.name)
    def test_every_subset_query_selects_something(self, path, dataset):
        _, ann = dataset
        cfg = RunConfig.from_yaml(str(path))
        if cfg.subset_query:
            assert len(ann.query(cfg.subset_query, engine="python")) > 0

    def test_the_minimal_config_needs_no_optional_backend(self):
        # It is the first command a new user runs, on a bare `pip install -e .`, so it
        # must not name a method that needs R or Octave.
        from combobatch.methods import METHOD_REGISTRY

        cfg = RunConfig.from_yaml(str(EXAMPLES / "configs" / "minimal.yaml"))
        for selection in cfg.methods:
            spec = METHOD_REGISTRY[selection.name]
            assert not spec.requires_r and not spec.requires_octave

    def test_the_hyperparameter_config_really_sweeps_one_method(self):
        # The feature the tool exists for: the same method twice at different values.
        cfg = RunConfig.from_yaml(str(EXAMPLES / "configs" / "hyperparameters.yaml"))
        names = [selection.name for selection in cfg.methods]
        assert len(names) != len(set(names))

    def test_the_sweep_resolves_to_distinct_output_keys(self):
        from combobatch.cli.dispatch_cmd import build_job_matrix, job_key

        cfg = RunConfig.from_yaml(str(EXAMPLES / "configs" / "hyperparameters.yaml"))
        keys = [job_key(*pair) for pair in build_job_matrix(cfg)]
        assert len(keys) == len(set(keys))
