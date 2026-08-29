"""One combination, end to end, with storage pointed at a temporary directory.

These tests pin the behaviours the donor got wrong: a failure that still produced an
output, a post-removal variant identical to its input, an imputer label that did not
match the imputer that ran, and a sidecar that was not valid JSON.
"""

from __future__ import annotations

import dataclasses
import json

import pandas as pd
import pytest

from combobatch import storage
from combobatch.config import (
    ColumnSpec,
    ImputationSelection,
    MethodSelection,
    MetricsSpec,
    PostRemovalSpec,
    RunConfig,
)
from combobatch.harmonize import (
    MANIFEST_KEY,
    build_manifest,
    combined_param_tag,
    prepare,
    run_one_combination,
    write_manifest,
)
from combobatch.methods import METHOD_REGISTRY
from combobatch.params import ParamParseError

from tests.conftest import make_expression


def _patch_method_fn(monkeypatch, key: str, fn) -> None:
    """Swap one method's implementation. MethodSpec is frozen, so replace the entry."""
    monkeypatch.setitem(
        METHOD_REGISTRY, key, dataclasses.replace(METHOD_REGISTRY[key], fn=fn)
    )


@pytest.fixture
def dataset(tmp_path):
    """Write a small dataset to disk and return its two URIs."""
    exp, ann = make_expression(n_samples=30, n_genes=40)
    exp_uri = tmp_path / "exp.tsv.gz"
    ann_uri = tmp_path / "ann.csv"
    exp.to_csv(exp_uri, sep="\t")
    ann.to_csv(ann_uri)
    return str(exp_uri), str(ann_uri), exp, ann


@pytest.fixture
def config(dataset, tmp_path):
    """A minimal single-combination configuration writing into tmp_path/out."""
    exp_uri, ann_uri, _, _ = dataset

    def build(**overrides):
        base = dict(
            exp_uri=exp_uri,
            ann_uri=ann_uri,
            out_uri=str(tmp_path / "out"),
            columns=ColumnSpec(batch="batch", bio="bio"),
            imputations=(ImputationSelection(name="strict"),),
            methods=(MethodSelection(name="01_raw"),),
        )
        base.update(overrides)
        return RunConfig(**base)

    return build


def _backend(cfg):
    return storage.backend_for_root(cfg.out_uri)


class TestCombinedParamTag:
    def test_pure_defaults_produce_no_tag(self):
        assert (
            combined_param_tag(
                ImputationSelection(name="strict"), MethodSelection(name="01_raw")
            )
            is None
        )

    def test_method_parameters_appear(self):
        tag = combined_param_tag(
            ImputationSelection(name="strict"),
            MethodSelection(name="10_mnn", params={"k": 50}),
        )
        assert tag == "k50"

    def test_imputer_parameters_appear(self):
        tag = combined_param_tag(
            ImputationSelection(name="knn", params={"knn_k": 10}),
            MethodSelection(name="01_raw"),
        )
        assert tag == "knn_k10"

    def test_max_na_frac_is_part_of_the_tag(self):
        tag = combined_param_tag(
            ImputationSelection(name="knn", max_na_frac=0.3),
            MethodSelection(name="01_raw"),
        )
        assert tag == "max_na_frac0.3"

    def test_strict_at_its_own_default_is_untagged(self):
        # strict's default ceiling is 0.0, not the 0.20 the other imputers use; a shared
        # literal default would both mistag it and stop it being strict.
        assert (
            combined_param_tag(
                ImputationSelection(name="strict"), MethodSelection(name="01_raw")
            )
            is None
        )

    def test_both_sources_merge(self):
        tag = combined_param_tag(
            ImputationSelection(name="knn", params={"knn_k": 10}),
            MethodSelection(name="10_mnn", params={"k": 50}),
        )
        assert tag == "k50-knn_k10"

    def test_a_name_tuned_on_both_sides_raises(self, monkeypatch):
        from combobatch.imputation import IMPUTER_REGISTRY
        from combobatch.methods.base import HyperParam

        # No real pair collides today, so the guard is exercised by giving the imputer a
        # parameter named like one of the method's. The specs are frozen, so the registry
        # entry is replaced rather than mutated.
        knn = IMPUTER_REGISTRY["knn"]
        monkeypatch.setitem(
            IMPUTER_REGISTRY,
            "knn",
            dataclasses.replace(
                knn,
                hyperparams={
                    **knn.hyperparams,
                    "k": HyperParam("k", int, 5, "collides with the method's k"),
                },
            ),
        )
        with pytest.raises(ParamParseError, match="both tune k"):
            combined_param_tag(
                ImputationSelection(name="knn", params={"k": 3}),
                MethodSelection(name="10_mnn", params={"k": 50}),
            )

    def test_explicit_tag_wins(self):
        tag = combined_param_tag(
            ImputationSelection(name="knn", params={"knn_k": 10}),
            MethodSelection(name="10_mnn", params={"k": 50}, param_tag="mysweep"),
        )
        assert tag == "mysweep"


class TestPreparedCaching:
    def test_first_call_writes_the_prepared_pair(self, config, dataset):
        cfg = config()
        _, _, exp, ann = dataset
        result = prepare(cfg, cfg.imputations[0], exp, ann, backend=_backend(cfg))
        assert not result.from_cache
        assert _backend(cfg).exists("prepared/strict__exp.tsv.gz")
        assert _backend(cfg).exists("prepared/strict__ann.tsv.gz")

    def test_second_call_reads_the_cache(self, config, dataset):
        cfg = config()
        _, _, exp, ann = dataset
        backend = _backend(cfg)
        prepare(cfg, cfg.imputations[0], exp, ann, backend=backend)
        again = prepare(cfg, cfg.imputations[0], exp, ann, backend=backend)
        assert again.from_cache
        pd.testing.assert_frame_equal(again.exp_df, exp, check_dtype=False)

    def test_different_parameters_use_different_cache_keys(self, config, dataset):
        cfg = config(imputations=(ImputationSelection(name="knn", max_na_frac=0.5),))
        _, _, exp, ann = dataset
        prepare(cfg, cfg.imputations[0], exp, ann, backend=_backend(cfg))
        assert _backend(cfg).exists("prepared/knn__max_na_frac0.5__exp.tsv.gz")
        assert not _backend(cfg).exists("prepared/knn__exp.tsv.gz")

    def test_strict_uses_a_zero_ceiling_not_the_shared_default(self, config, dataset):
        exp, ann = make_expression(n_samples=20, n_genes=30)
        exp = exp.copy()
        exp.iloc[0, 0] = float("nan")  # 5% NA: under 0.20, over 0.0
        cfg = config()
        result = prepare(cfg, cfg.imputations[0], exp, ann, backend=_backend(cfg))
        assert result.exp_df.shape[1] == 29, "strict must drop any gene with a NaN"


class TestRunOneCombination:
    def test_writes_matrix_sidecar_and_genes(self, config):
        cfg = config()
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])
        backend = _backend(cfg)

        assert [row["status"] for row in rows] == ["ok"]
        assert backend.exists("exp/strict__01_raw.tsv.gz")
        assert backend.exists("metrics/strict__01_raw_metrics.json")
        assert backend.exists("genes/strict__01_raw_genes.json")

    def test_sidecar_is_strict_json(self, config):
        cfg = config()
        run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])
        raw = _backend(cfg).read_bytes("metrics/strict__01_raw_metrics.json").decode()

        assert (
            "NaN" not in raw
        ), "bare NaN is invalid JSON and no strict parser accepts it"
        payload = json.loads(raw)
        assert payload["status"] == "ok"
        assert payload["n_samples"] == 30

    def test_skip_if_exists_reports_cached_without_recomputing(self, config):
        cfg = config()
        run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])
        rows = run_one_combination(
            cfg, cfg.imputations[0], cfg.methods[0], skip_if_exists=True
        )
        assert [row["status"] for row in rows] == ["cached"]

    def test_unavailable_backend_is_an_honest_skip(self, config):
        cfg = config(methods=(MethodSelection(name="03_limma"),))
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        assert [row["status"] for row in rows] == ["skipped"]
        assert "R" in rows[0]["error"]
        assert not _backend(cfg).exists("exp/strict__03_limma.tsv.gz")

    def test_a_failing_method_writes_no_output(self, config, monkeypatch):
        def explode(*args, **kwargs):
            raise RuntimeError("synthetic failure")

        _patch_method_fn(monkeypatch, "01_raw", explode)
        cfg = config()
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        assert [row["status"] for row in rows] == ["failed"]
        assert "synthetic failure" in rows[0]["error"]
        assert not _backend(cfg).exists("exp/strict__01_raw.tsv.gz")

    def test_parameters_reach_the_method(self, config, monkeypatch):
        seen = {}

        def record(exp_df, ann_df, **kwargs):
            seen.update(kwargs)
            return exp_df

        _patch_method_fn(monkeypatch, "10_mnn", record)
        cfg = config(methods=(MethodSelection(name="10_mnn", params={"k": 42}),))
        run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        assert seen["k"] == 42, "the whole point of the tool is that this arrives"
        assert seen["batch_col"] == "batch"

    def test_reference_batch_is_resolved_for_methods_that_need_one(
        self, config, monkeypatch
    ):
        seen = {}

        def record(exp_df, ann_df, **kwargs):
            seen.update(kwargs)
            return exp_df

        _patch_method_fn(monkeypatch, "13_fsmvn", record)
        cfg = config(methods=(MethodSelection(name="13_fsmvn"),))
        run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        assert seen["target_group"] in {"Batch0", "Batch1", "Batch2"}


class TestPostRemoval:
    def test_disabled_produces_one_untagged_output(self, config):
        cfg = config()
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])
        assert len(rows) == 1
        assert rows[0]["post_rm"] is None
        assert _backend(cfg).exists("exp/strict__01_raw.tsv.gz")

    def test_enabled_produces_post0_and_post1(self, config):
        cfg = config(
            post_removal=PostRemovalSpec(enabled=True, min_batch_size=5),
        )
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        assert [row["post_rm"] for row in rows] == [False, True]
        assert _backend(cfg).exists("exp/strict__01_raw__post0.tsv.gz")
        assert _backend(cfg).exists("exp/strict__01_raw__post1.tsv.gz")

    def test_post1_actually_differs_from_post0(self, config):
        from combobatch import dataio

        cfg = config(post_removal=PostRemovalSpec(enabled=True, min_batch_size=5))
        run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])
        backend = _backend(cfg)
        post0 = dataio.read_table_from(backend, "exp/strict__01_raw__post0.tsv.gz")
        post1 = dataio.read_table_from(backend, "exp/strict__01_raw__post1.tsv.gz")

        # The donor wrote post1 even when removal failed, giving a file byte-identical to
        # post0 that was then scored as if a batch had been dropped.
        assert post1.shape[0] < post0.shape[0]

    def test_post_removal_failure_writes_nothing_and_says_so(self, config, monkeypatch):
        from combobatch import harmonize
        from combobatch.postremoval import PostRemovalError

        def explode(*args, **kwargs):
            raise PostRemovalError("synthetic post-removal failure")

        monkeypatch.setattr(harmonize, "apply_post_removal", explode)
        cfg = config(post_removal=PostRemovalSpec(enabled=True, min_batch_size=5))
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        assert [row["status"] for row in rows] == ["ok", "post_removal_failed"]
        assert not _backend(cfg).exists("exp/strict__01_raw__post1.tsv.gz")


class TestInlineMetrics:
    def test_metrics_land_in_the_sidecar(self, config):
        import json

        cfg = config(metrics=MetricsSpec(enabled=True, groups=frozenset("EK")))
        run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        sidecar = json.loads(
            _backend(cfg).read_bytes("metrics/strict__01_raw_metrics.json").decode()
        )
        assert sidecar["n_samples"] == 30
        assert sidecar["n_genes_noNA"] == 40
        assert set(sidecar["metrics_groups_run"]) == {"E", "K"}

    def test_group_l_uses_this_run_own_pre_harmonization_matrix(self, config):
        # The donor downloaded a separate 01_raw output, which could silently be a
        # different baseline than the one actually harmonized.
        cfg = config(
            metrics=MetricsSpec(enabled=True, groups=frozenset("L"), min_cohort_n=5)
        )
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])

        assert rows[0]["status"] == "ok"
        # 01_raw returns its input unchanged, so the reference is that same matrix.
        assert rows[0]["mk_rho_mean_all_genes"] == pytest.approx(1.0)

    def test_disabled_metrics_add_nothing(self, config):
        cfg = config()
        rows = run_one_combination(cfg, cfg.imputations[0], cfg.methods[0])
        assert "metrics_groups_run" not in rows[0]


class TestManifest:
    def test_records_full_parameters_including_defaults(self, config):
        cfg = config(methods=(MethodSelection(name="10_mnn", params={"k": 50}),))
        manifest = build_manifest(cfg)
        entry = manifest["outputs"]["strict__10_mnn__k50"]

        assert entry["method_params"]["k"] == 50
        # The tag carries only the tuned parameter; the manifest carries every declared
        # one, which is what keeps an abbreviated name self-describing.
        assert set(entry["method_params"]) == set(METHOD_REGISTRY["10_mnn"].hyperparams)
        assert entry["param_tag"] == "k50"

    def test_covers_every_output_including_post_variants(self, config):
        cfg = config(
            imputations=(
                ImputationSelection(name="strict"),
                ImputationSelection(name="knn"),
            ),
            post_removal=PostRemovalSpec(enabled=True),
        )
        assert len(build_manifest(cfg)["outputs"]) == 4

    def test_writing_twice_merges_rather_than_replaces(self, config):
        from combobatch import dataio

        first = config()
        write_manifest(first)
        second = config(methods=(MethodSelection(name="17_quantile"),))
        write_manifest(second)

        manifest = dataio.read_json(_backend(first), MANIFEST_KEY)
        assert set(manifest["outputs"]) == {"strict__01_raw", "strict__17_quantile"}

    def test_is_strict_json(self, config):
        cfg = config()
        write_manifest(cfg)
        raw = _backend(cfg).read_bytes(MANIFEST_KEY).decode()
        assert "NaN" not in raw
        json.loads(raw)
