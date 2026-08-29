"""Every method in the registry, against the shared output contract.

Parametrized over the whole registry. A method whose backend is missing is **skipped via
its declared capability flags**, never quietly passed — so on a bare laptop this
exercises the dependency-free methods and honestly reports the rest, while inside the
Docker image `tests/integration/test_methods_all_backends.py` asserts that nothing skips.
"""

from __future__ import annotations

import pandas as pd
import pytest

from combobatch.methods import METHOD_REGISTRY, resolve_reference_batch
from tests.conftest import make_expression

ALL_KEYS = sorted(METHOD_REGISTRY)


def run_method(key: str, exp, ann, **params):
    """Invoke a method the way the harmonizer will, skipping if its backend is absent."""
    spec = METHOD_REGISTRY[key]
    missing = spec.missing_backends()
    if missing:
        pytest.skip(f"{key} needs {', '.join(missing)}")

    call: dict[str, object] = {"batch_col": "batch", "bio_col": "bio"}
    if spec.uses_reference_batch:
        call["target_group"] = resolve_reference_batch(ann, "batch")
    call.update(params)
    return spec.fn(exp, ann, **call)


class TestOutputContract:
    @pytest.mark.parametrize("key", ALL_KEYS)
    def test_returns_a_dataframe(self, key, synthetic):
        exp, ann = synthetic
        assert isinstance(run_method(key, exp, ann), pd.DataFrame)

    @pytest.mark.parametrize("key", ALL_KEYS)
    def test_preserves_the_sample_index_and_its_order(self, key, synthetic):
        exp, ann = synthetic
        out = run_method(key, exp, ann)
        assert list(out.index) == list(exp.index)

    @pytest.mark.parametrize("key", ALL_KEYS)
    def test_output_is_not_entirely_missing(self, key, synthetic):
        exp, ann = synthetic
        out = run_method(key, exp, ann)
        assert out.notna().any().any(), f"{key} returned an all-NaN matrix"

    @pytest.mark.parametrize("key", ALL_KEYS)
    def test_gene_space_is_preserved_or_documented(self, key, synthetic):
        exp, ann = synthetic
        out = run_method(key, exp, ann)
        if key == "25_angel":
            # The one method that reduces the gene space by design.
            assert out.shape[1] <= exp.shape[1]
        else:
            assert list(out.columns) == list(exp.columns)

    @pytest.mark.parametrize("key", ALL_KEYS)
    def test_is_deterministic(self, key, synthetic):
        exp, ann = synthetic
        first = run_method(key, exp, ann)
        second = run_method(key, exp, ann)
        pd.testing.assert_frame_equal(first, second)

    @pytest.mark.parametrize("key", ALL_KEYS)
    def test_does_not_mutate_its_input(self, key, synthetic):
        exp, ann = synthetic
        before = exp.copy()
        run_method(key, exp, ann)
        pd.testing.assert_frame_equal(exp, before)


class TestBaselineBehaviour:
    def test_raw_is_an_exact_passthrough(self, synthetic):
        exp, ann = synthetic
        pd.testing.assert_frame_equal(run_method("01_raw", exp, ann), exp)

    def test_raw_returns_a_copy_not_the_same_object(self, synthetic):
        exp, ann = synthetic
        assert run_method("01_raw", exp, ann) is not exp

    def test_median_scaling_aligns_batch_medians(self, synthetic_with_batch_effect):
        exp, ann = synthetic_with_batch_effect
        out = run_method("02_median_scaling", exp, ann)
        medians = out.groupby(ann["batch"]).median().mean(axis=1)
        assert medians.max() - medians.min() < 0.5

    def test_rank_output_is_bounded(self, synthetic):
        exp, ann = synthetic
        out = run_method("18_rank", exp, ann)
        assert out.min().min() > 0.0
        assert out.max().max() < 1.0

    def test_quantile_gives_every_sample_the_same_distribution(self, synthetic):
        exp, ann = synthetic
        out = run_method("17_quantile", exp, ann)
        sorted_rows = out.apply(lambda row: row.sort_values().values, axis=1)
        first = sorted_rows.iloc[0]
        assert all((row == first).all() for row in sorted_rows)


class TestRnaseqOnlyGuard:
    @pytest.mark.parametrize("key", ["22_tmm", "23_vst"])
    def test_microarray_labels_raise_not_implemented(self, key):
        exp, ann = make_expression()
        ann = ann.copy()
        ann.loc[ann.index[:5], "batch"] = "GPL570_FFPE"
        spec = METHOD_REGISTRY[key]
        # The guard runs before any backend check, so this holds with or without R.
        with pytest.raises(NotImplementedError, match="RNA-seq only"):
            spec.fn(exp, ann, batch_col="batch", bio_col="bio")

    def test_guard_names_an_offending_label(self):
        exp, ann = make_expression()
        ann = ann.copy()
        ann.loc[ann.index[:5], "batch"] = "GPL570_FFPE"
        with pytest.raises(NotImplementedError, match="GPL570_FFPE"):
            METHOD_REGISTRY["22_tmm"].fn(exp, ann, batch_col="batch", bio_col="bio")


class TestResolveReferenceBatch:
    def test_picks_the_largest_batch_by_default(self):
        _, ann = make_expression(n_samples=30, n_batches=3)
        ann = ann.copy()
        ann["batch"] = ["Big"] * 20 + ["Small"] * 10
        assert resolve_reference_batch(ann, "batch") == "Big"

    def test_honours_an_explicit_request(self):
        _, ann = make_expression()
        assert resolve_reference_batch(ann, "batch", "Batch1") == "Batch1"

    def test_unknown_request_raises(self):
        _, ann = make_expression()
        with pytest.raises(ValueError, match="not a value of"):
            resolve_reference_batch(ann, "batch", "NoSuchBatch")
