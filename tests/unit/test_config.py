"""Configuration loading, CLI overlay precedence, and validation."""

from __future__ import annotations

import pandas as pd
import pytest
import yaml

from combobatch.config import (
    ColumnSpec,
    ConfigError,
    ImputationSelection,
    MethodSelection,
    RunConfig,
)

MINIMAL = {
    "input": {"expression": "./exp.tsv.gz", "annotation": "./ann.csv"},
    "output": "./out/",
    "columns": {"batch": "RNA_BATCH", "bio": "Diagnosis"},
    "imputations": ["strict"],
    "methods": ["01_raw"],
}


@pytest.fixture
def ann():
    return pd.DataFrame(
        {
            "RNA_BATCH": ["b1", "b1", "b2", "b2"],
            "Diagnosis": ["FL", "DLBCL", "FL", "DLBCL"],
            "COHORT": ["c1", "c1", "c2", "c2"],
        },
        index=["s1", "s2", "s3", "s4"],
    )


def config(**overrides):
    data = {**MINIMAL, **overrides}
    return RunConfig.from_mapping(data)


class TestColumnSpec:
    def test_metric_columns_default_to_the_primary_columns(self):
        spec = ColumnSpec(batch="B", bio="D")
        assert spec.metric_batch_cols == ("B",)
        assert spec.metric_bio_cols == ("D",)

    def test_predict_class_col_defaults_to_bio(self):
        assert ColumnSpec(batch="B", bio="D").predict_class_col == "D"

    def test_explicit_metric_columns_are_kept(self):
        spec = ColumnSpec(batch="B", bio="D", metric_batch_cols=("B", "PLATFORM"))
        assert spec.metric_batch_cols == ("B", "PLATFORM")

    def test_all_columns_is_deduplicated(self):
        spec = ColumnSpec(batch="B", bio="D", cohort="C", metric_batch_cols=("B",))
        assert spec.all_columns() == ("B", "D", "C")


class TestFromMapping:
    def test_minimal_config(self):
        cfg = config()
        assert cfg.exp_uri == "./exp.tsv.gz"
        assert cfg.out_uri == "./out/"
        assert cfg.imputations == (ImputationSelection("strict"),)
        assert cfg.methods == (MethodSelection("01_raw"),)

    def test_entries_may_be_strings_or_mappings(self):
        cfg = config(methods=["01_raw", {"name": "10_mnn", "params": {"k": 50}}])
        assert cfg.methods[0] == MethodSelection("01_raw")
        assert cfg.methods[1].params == {"k": 50}

    def test_imputation_params_and_na_ceiling(self):
        cfg = config(
            imputations=[{"name": "knn", "max_na_frac": 0.3, "params": {"knn_k": 10}}]
        )
        assert cfg.imputations[0].max_na_frac == pytest.approx(0.3)
        assert cfg.imputations[0].params == {"knn_k": 10}

    def test_metrics_block(self):
        cfg = config(metrics={"enabled": True, "groups": ["A", "B", "K"]})
        assert cfg.metrics.enabled is True
        assert cfg.metrics.groups == frozenset({"A", "B", "K"})

    def test_groups_accept_a_string(self):
        assert config(metrics={"groups": "ABK"}).metrics.groups == frozenset("ABK")

    def test_post_removal_block(self):
        cfg = config(post_removal={"enabled": True, "n_batches": 3})
        assert cfg.post_removal.enabled is True
        assert cfg.post_removal.n_batches == 3

    def test_post_removal_accepts_a_bare_boolean(self):
        assert config(post_removal=True).post_removal.enabled is True

    def test_defaults_are_conservative(self):
        cfg = config()
        assert cfg.metrics.enabled is False
        assert cfg.post_removal.enabled is False
        assert cfg.reference_batch is None
        assert cfg.random_seed == 42

    def test_missing_required_key_raises(self):
        with pytest.raises(ConfigError, match="missing required"):
            RunConfig.from_mapping({"output": "./out/"})

    def test_entry_without_a_name_raises(self):
        with pytest.raises(ConfigError, match="has no 'name'"):
            config(methods=[{"params": {"k": 1}}])

    def test_non_list_methods_raises(self):
        with pytest.raises(ConfigError, match="expected a list"):
            config(methods="01_raw")


class TestFromYaml:
    def test_round_trip_through_a_file(self, tmp_path):
        path = tmp_path / "cfg.yaml"
        path.write_text(yaml.safe_dump(MINIMAL))
        cfg = RunConfig.from_yaml(str(path))
        assert cfg.columns.batch == "RNA_BATCH"

    def test_non_mapping_raises(self, tmp_path):
        path = tmp_path / "cfg.yaml"
        path.write_text("- just\n- a list\n")
        with pytest.raises(ConfigError, match="mapping at the top level"):
            RunConfig.from_yaml(str(path))


class TestOverrides:
    def test_cli_value_wins(self):
        assert config().with_overrides(out_uri="s3://b/x").out_uri == "s3://b/x"

    def test_none_leaves_the_file_value_alone(self):
        # An unset CLI flag arrives as None and must not clobber the config file.
        assert config().with_overrides(out_uri=None).out_uri == "./out/"

    def test_returns_an_equal_object_when_nothing_overridden(self):
        cfg = config()
        assert cfg.with_overrides() == cfg


class TestStaticValidation:
    def test_minimal_config_validates(self):
        config().validate()

    def test_no_methods_raises(self):
        with pytest.raises(ConfigError, match="no harmonization methods"):
            config(methods=[]).validate()

    def test_no_imputations_raises(self):
        with pytest.raises(ConfigError, match="no imputation strategies"):
            config(imputations=[]).validate()

    def test_duplicate_method_raises(self):
        with pytest.raises(ConfigError, match="more than once"):
            config(methods=["01_raw", "01_raw"]).validate()

    def test_the_same_method_at_different_parameters_is_allowed(self):
        # The sweep this tool exists for: two runs of one method, different k, two
        # output files. Rejecting it would make "would k=50 have done better?"
        # unanswerable in a single run.
        config(
            methods=[{"name": "10_mnn"}, {"name": "10_mnn", "params": {"k": 50}}]
        ).validate()

    def test_the_same_imputer_at_different_parameters_is_allowed(self):
        config(
            imputations=[
                {"name": "knn", "params": {"knn_k": 5}},
                {"name": "knn", "params": {"knn_k": 25}},
            ]
        ).validate()

    def test_the_same_method_at_the_same_parameters_still_raises(self):
        with pytest.raises(ConfigError, match="more than once"):
            config(
                methods=[
                    {"name": "10_mnn", "params": {"k": 50}},
                    {"name": "10_mnn", "params": {"k": 50}},
                ]
            ).validate()

    def test_a_repeated_param_tag_raises(self):
        # Different parameters, but the user forced them onto the same filename.
        with pytest.raises(ConfigError, match="more than once"):
            config(
                methods=[
                    {"name": "10_mnn", "params": {"k": 50}, "param_tag": "tuned"},
                    {"name": "10_mnn", "params": {"k": 30}, "param_tag": "tuned"},
                ]
            ).validate()

    def test_bad_metric_group_letter_raises(self):
        with pytest.raises(ConfigError, match="unknown metric group"):
            config(metrics={"groups": ["A", "Z"]}).validate()

    @pytest.mark.parametrize("value", [-0.1, 1.5])
    def test_out_of_range_max_na_frac_raises(self, value):
        with pytest.raises(ConfigError, match="between 0 and 1"):
            config(imputations=[{"name": "knn", "max_na_frac": value}]).validate()

    def test_negative_n_perm_raises(self):
        with pytest.raises(ConfigError, match="non-negative"):
            config(metrics={"n_perm": -1}).validate()

    def test_post_removal_zero_batches_raises(self):
        with pytest.raises(ConfigError, match="at least 1"):
            config(post_removal={"enabled": True, "n_batches": 0}).validate()

    def test_invalid_param_tag_raises(self):
        with pytest.raises(Exception):
            config(methods=[{"name": "01_raw", "param_tag": "a__b"}]).validate()


class TestRegistryValidation:
    def test_unknown_method_raises_with_a_suggestion(self):
        with pytest.raises(ConfigError, match="did you mean '10_mnn'"):
            config(methods=["10_mn"]).validate(known_methods={"10_mnn", "01_raw"})

    def test_unknown_imputer_raises(self):
        with pytest.raises(ConfigError, match="unknown imputer"):
            config(imputations=["knnn"]).validate(known_imputers={"knn", "strict"})

    def test_known_names_pass(self):
        config().validate(known_methods={"01_raw"}, known_imputers={"strict"})

    def test_registries_are_optional(self):
        # Phase 1 has no registries yet; validation must still do everything else.
        config().validate(known_methods=None)


class TestAnnotationValidation:
    def test_valid_columns_pass(self, ann):
        config().validate(ann)

    def test_missing_column_raises_with_a_suggestion(self, ann):
        cfg = config(columns={"batch": "RNA_BATCHH", "bio": "Diagnosis"})
        with pytest.raises(ConfigError, match="did you mean 'RNA_BATCH'"):
            cfg.validate(ann)

    def test_missing_metric_column_is_also_caught(self, ann):
        cfg = config(
            columns={
                "batch": "RNA_BATCH",
                "bio": "Diagnosis",
                "metric_batch_cols": ["RNA_BATCH", "PLATFORM"],
            }
        )
        with pytest.raises(ConfigError, match="PLATFORM"):
            cfg.validate(ann)

    def test_valid_reference_batch_passes(self, ann):
        config(reference_batch="b1").validate(ann)

    def test_unknown_reference_batch_raises(self, ann):
        with pytest.raises(ConfigError, match="not a value of"):
            config(reference_batch="b99").validate(ann)

    def test_subset_query_is_validated_too(self, ann):
        with pytest.raises(ConfigError, match="unknown column"):
            config(subset_query="Nope == 1").validate(ann)

    def test_valid_subset_query_passes(self, ann):
        config(subset_query="Diagnosis == 'FL'").validate(ann)


class TestCombinations:
    def test_cross_product_in_declaration_order(self):
        cfg = config(imputations=["strict", "knn"], methods=["01_raw", "10_mnn"])
        pairs = [(i.name, m.name) for i, m in cfg.combinations()]
        assert pairs == [
            ("strict", "01_raw"),
            ("strict", "10_mnn"),
            ("knn", "01_raw"),
            ("knn", "10_mnn"),
        ]


class TestSelectionTags:
    def test_defaults_have_no_tag(self):
        assert MethodSelection("01_raw").tag() is None

    def test_params_produce_a_readable_tag(self):
        assert MethodSelection("10_mnn", {"k": 50}).tag() == "k50"

    def test_explicit_tag_wins(self):
        assert MethodSelection("10_mnn", {"k": 50}, param_tag="mine").tag() == "mine"
