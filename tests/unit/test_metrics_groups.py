"""Each metric group, checked for the thing it is actually supposed to detect.

A group that runs without raising proves nothing — every one of these would pass against
a function returning constants. Each test below moves the input in a known direction and
asserts the metric follows.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from combobatch.config import ColumnSpec
from combobatch.metrics import METRIC_GROUP_REGISTRY, MetricColumns
from combobatch.metrics import groups as g
from combobatch.metrics.embeddings import compute_pca
from combobatch.metrics.panels import Panel

from tests.conftest import make_expression


@pytest.fixture
def columns():
    return MetricColumns.from_spec(
        ColumnSpec(batch="batch", bio="bio", cohort="cohort")
    )


@pytest.fixture
def clean():
    return make_expression(n_samples=60, n_genes=80)


@pytest.fixture
def batched():
    return make_expression(n_samples=60, n_genes=80, batch_effect=8.0)


def _pca(exp_df):
    return compute_pca(exp_df)


class TestGroupE:
    def test_reports_the_shape(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_e(exp, ann, columns)
        assert result["n_samples"] == 60
        assert result["n_genes"] == 80

    def test_counts_the_batches(self, clean, columns):
        exp, ann = clean
        assert g.compute_group_e(exp, ann, columns)["n_batches"] == 3

    def test_detects_zero_inflation(self, clean, columns):
        exp, ann = clean
        zeroed = exp.copy()
        zeroed.iloc[:, :40] = 0.0
        assert (
            g.compute_group_e(zeroed, ann, columns)["zero_fraction_global"]
            > g.compute_group_e(exp, ann, columns)["zero_fraction_global"]
        )


class TestGroupK:
    def test_a_clean_matrix_has_no_missing_values(self, clean):
        exp, _ = clean
        result = g.compute_group_k(exp)
        assert result["n_na_cells"] == 0
        assert result["pct_genes_noNA"] == 100.0

    def test_counts_what_a_failed_harmonizer_destroyed(self, clean):
        exp, _ = clean
        broken = exp.copy()
        broken.iloc[:, :10] = np.nan
        result = g.compute_group_k(broken)
        assert result["n_genes_allNA"] == 10
        assert result["n_genes_noNA"] == 70


class TestGroupA:
    def test_a_batch_effect_raises_r_squared(self, clean, batched, columns):
        clean_exp, ann = clean
        batched_exp, _ = batched
        clean_r2 = g.compute_group_a(
            *_pca(clean_exp), ann, columns, n_dsc_permutations=9
        )
        dirty_r2 = g.compute_group_a(
            *_pca(batched_exp), ann, columns, n_dsc_permutations=9
        )
        assert dirty_r2["r2_batch"] > clean_r2["r2_batch"]

    def test_pcr_moves_the_other_way(self, clean, batched, columns):
        clean_exp, ann = clean
        batched_exp, _ = batched
        # PCR is 1 minus the weighted R-squared, so a stronger effect lowers it.
        assert (
            g.compute_group_a(*_pca(batched_exp), ann, columns, n_dsc_permutations=9)[
                "pcr_batch"
            ]
            < g.compute_group_a(*_pca(clean_exp), ann, columns, n_dsc_permutations=9)[
                "pcr_batch"
            ]
        )

    def test_permutation_p_value_can_never_be_zero(self, batched, columns):
        exp, ann = batched
        result = g.compute_group_a(*_pca(exp), ann, columns, n_dsc_permutations=9)
        assert result["dsc_pvalue_batch"] >= 1 / 10


class TestGroupJ:
    def test_shares_are_ordered_and_bounded(self, clean):
        exp, _ = clean
        _, eigenvalues = _pca(exp)
        result = g.compute_group_j(eigenvalues)
        assert result["pct_var_pc1"] >= result["pct_var_pc2"]
        assert 0 < result["pct_var_cum_top10"] <= 100


class TestGroupB:
    def test_a_batch_effect_lowers_mixing(self, clean, batched, columns):
        clean_exp, ann = clean
        batched_exp, _ = batched
        clean_lisi = g.compute_group_b(_pca(clean_exp)[0], ann, columns)[
            "ilisi_mean_batch"
        ]
        dirty_lisi = g.compute_group_b(_pca(batched_exp)[0], ann, columns)[
            "ilisi_mean_batch"
        ]
        assert dirty_lisi < clean_lisi

    def test_batch_silhouette_normalization_rewards_zero(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_b(_pca(exp)[0], ann, columns)
        assert result["asw_batch_norm_batch"] == pytest.approx(
            1.0 - abs(result["asw_batch_batch"])
        )

    def test_bio_silhouette_normalization_rewards_separation(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_b(_pca(exp)[0], ann, columns)
        assert result["asw_bio_norm_bio"] == pytest.approx(
            (result["asw_bio_bio"] + 1.0) / 2.0
        )


class TestGroupH:
    def test_a_batch_effect_shrinks_the_distance_ratio(self, clean, batched, columns):
        clean_exp, ann = clean
        batched_exp, _ = batched
        assert (
            g.compute_group_h(batched_exp, ann, columns)["dist_ratio_batch"]
            < g.compute_group_h(clean_exp, ann, columns)["dist_ratio_batch"]
        )


class TestGroupD:
    def test_a_batch_effect_raises_the_ks_statistic(self, clean, batched, columns):
        clean_exp, ann = clean
        batched_exp, _ = batched
        assert (
            g.compute_group_d(batched_exp, ann, columns, n_genes_max=25)[
                "ks_mean_D_batch"
            ]
            > g.compute_group_d(clean_exp, ann, columns, n_genes_max=25)[
                "ks_mean_D_batch"
            ]
        )

    def test_a_timeout_records_nulls_rather_than_hanging(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_d(exp, ann, columns, n_genes_max=25, timeout_s=0)
        assert result["error_D1"] == "timeout"
        assert np.isnan(result["ks_mean_D_batch"])


class TestGroupG:
    def test_connectivity_is_bounded(self, clean, columns):
        exp, ann = clean
        value = g.compute_group_g(_pca(exp)[0], ann, columns)["graph_connectivity_bio"]
        assert 0.0 < value <= 1.0


class TestGroupI:
    def test_a_real_batch_effect_scores_above_the_null(self, batched, columns):
        exp, ann = batched
        result = g.compute_group_i(
            ann, _pca(exp)[0], columns, n_permutations=10, max_samples=60
        )
        assert result["wm_batch"] > 0

    def test_noise_does_not(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_i(
            ann, _pca(exp)[0], columns, n_permutations=10, max_samples=60
        )
        assert result["wm_batch"] < 0.1


class TestGroupL:
    def test_a_self_reference_correlates_perfectly(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_l(exp, ann, columns, ref_df=exp, min_cohort_n=5)
        assert result["mk_is_self_reference"] is True
        assert result["mk_rho_mean_all_genes"] == pytest.approx(1.0)

    def test_a_monotone_per_gene_transform_also_scores_one(self, clean, columns):
        exp, ann = clean
        # The documented ceiling: rank correlation cannot see a monotone transform, so
        # this group is a guard rail rather than a ranking.
        transformed = np.log1p(exp)
        result = g.compute_group_l(
            transformed, ann, columns, ref_df=exp, min_cohort_n=5
        )
        assert result["mk_rho_mean_all_genes"] == pytest.approx(1.0)

    def test_shuffling_destroys_the_correlation(self, clean, columns):
        exp, ann = clean
        rng = np.random.default_rng(0)
        shuffled = pd.DataFrame(
            rng.permutation(exp.values), index=exp.index, columns=exp.columns
        )
        result = g.compute_group_l(shuffled, ann, columns, ref_df=exp, min_cohort_n=5)
        assert abs(result["mk_rho_mean_all_genes"]) < 0.3

    def test_small_cohorts_are_skipped_and_counted(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_l(exp, ann, columns, ref_df=exp, min_cohort_n=1000)
        assert result["mk_n_cohorts_used"] == 0
        assert result["mk_n_cohorts_skipped_small"] == 2

    def test_a_panel_restricts_the_genes(self, clean, columns):
        exp, ann = clean
        panel = Panel(genes=("GENE0000", "GENE0001", "GENE0002"))
        result = g.compute_group_l(
            exp, ann, columns, ref_df=exp, panel=panel, min_cohort_n=5
        )
        assert result["mk_n_panel_genes_used"] == 3

    def test_housekeeping_genes_are_contrasted_separately(self, clean, columns):
        exp, ann = clean
        panel = Panel(
            genes=("GENE0000", "GENE0001", "GENE0002"),
            housekeeping=frozenset({"GENE0002"}),
        )
        result = g.compute_group_l(
            exp, ann, columns, ref_df=exp, panel=panel, min_cohort_n=5
        )
        assert not np.isnan(result["mk_rho_mean_markers"])
        assert not np.isnan(result["mk_rho_mean_housekeeping"])


class TestGroupM:
    def test_same_biology_pairs_agree_more_than_different_ones(self, columns):
        # Construct biology that genuinely drives the profile, so the margin is real.
        rng = np.random.default_rng(0)
        n, n_genes = 60, 60
        signature = {
            "TypeA": rng.normal(0, 3, n_genes),
            "TypeB": rng.normal(0, 3, n_genes),
        }
        bios = ["TypeA" if i % 2 == 0 else "TypeB" for i in range(n)]
        values = np.array(
            [signature[bio] + rng.normal(0, 0.4, n_genes) for bio in bios]
        )
        exp = pd.DataFrame(
            values + 20,
            index=[f"S{i}" for i in range(n)],
            columns=[f"GENE{j:04d}" for j in range(n_genes)],
        )
        ann = pd.DataFrame(
            {"batch": [f"Batch{i % 3}" for i in range(n)], "bio": bios, "cohort": "C"},
            index=exp.index,
        )
        result = g.compute_group_m(exp, ann, columns)
        assert result["xb_rank_agree"] > result["xb_rank_disagree_diffbio"]
        assert result["xb_rank_agree_ratio"] > 0

    def test_pairs_are_counted_once_each(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_m(exp, ann, columns)
        total = result["xb_n_pairs_same_bio"] + result["xb_n_pairs_diff_bio"]
        # Every cross-batch pair, counted once: 60 samples over 3 equal batches.
        assert total == 20 * 20 * 3

    def test_too_few_genes_yields_a_complete_null_record(self, clean, columns):
        exp, ann = clean
        result = g.compute_group_m(exp, ann, columns, panel=Panel(genes=("GENE0000",)))
        assert np.isnan(result["xb_rank_agree"])
        assert result["xb_n_samples_used"] == 0


class TestGroupN:
    def test_real_biology_beats_the_permutation_null(self, columns):
        rng = np.random.default_rng(1)
        n, n_genes = 90, 40
        bios = ["TypeA" if i % 2 == 0 else "TypeB" for i in range(n)]
        signature = {
            "TypeA": rng.normal(0, 4, n_genes),
            "TypeB": rng.normal(0, 4, n_genes),
        }
        values = np.array(
            [signature[bio] + rng.normal(0, 1.0, n_genes) for bio in bios]
        )
        exp = pd.DataFrame(
            values + 30,
            index=[f"S{i}" for i in range(n)],
            columns=[f"GENE{j:04d}" for j in range(n_genes)],
        )
        ann = pd.DataFrame(
            {"batch": [f"Batch{i % 3}" for i in range(n)], "bio": bios, "cohort": "C"},
            index=exp.index,
        )
        result = g.compute_group_n(exp, ann, columns, n_pcs=5, n_perm=5, min_test_n=10)
        assert result["pv_lobo2_n_folds"] == 3
        assert result["pv_lobo2_f1_macro_mean"] > result["pv_lobo2_f1_macro_perm_mean"]
        assert result["pv_lobo2_f1_macro_delta"] > 0

    def test_permutation_p_value_can_never_be_zero(self, columns):
        exp, ann = make_expression(n_samples=90, n_genes=40)
        result = g.compute_group_n(exp, ann, columns, n_pcs=5, n_perm=5, min_test_n=10)
        assert result["pv_lobo2_f1_perm_pvalue"] >= 1 / 6

    def test_all_na_samples_are_counted_as_samples_not_genes(self, columns):
        exp, ann = make_expression(n_samples=90, n_genes=40)
        broken = exp.copy()
        broken.iloc[:2, :] = np.nan
        result = g.compute_group_n(
            broken, ann, columns, n_pcs=5, n_perm=1, min_test_n=10
        )
        # Dropping genes first would absorb the failed samples into the gene count and
        # report zero dropped samples.
        assert result["pv_n_samples_dropped_na"] == 2
        assert result["pv_n_genes_dropped_na"] == 0


class TestEveryGroupEmitsItsSentinel:
    @pytest.mark.parametrize("letter", sorted(set("ABEGHIJKLM")))
    def test_the_declared_sentinel_key_is_actually_produced(
        self, letter, clean, columns
    ):
        exp, ann = clean
        pca_coords, eigenvalues = _pca(exp)
        calls = {
            "A": lambda: g.compute_group_a(
                pca_coords, eigenvalues, ann, columns, n_dsc_permutations=9
            ),
            "B": lambda: g.compute_group_b(pca_coords, ann, columns),
            "E": lambda: g.compute_group_e(exp, ann, columns),
            "G": lambda: g.compute_group_g(pca_coords, ann, columns),
            "H": lambda: g.compute_group_h(exp, ann, columns),
            "I": lambda: g.compute_group_i(
                ann, pca_coords, columns, n_permutations=3, max_samples=60
            ),
            "J": lambda: g.compute_group_j(eigenvalues),
            "K": lambda: g.compute_group_k(exp),
            "L": lambda: g.compute_group_l(
                exp, ann, columns, ref_df=exp, min_cohort_n=5
            ),
            "M": lambda: g.compute_group_m(exp, ann, columns),
        }
        result = calls[letter]()
        sentinel = METRIC_GROUP_REGISTRY[letter].sentinel_key(columns)
        assert sentinel in result, f"group {letter} never emits its sentinel {sentinel}"
