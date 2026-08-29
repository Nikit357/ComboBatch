"""The generated docs must match the registries.

This is the test that stops documentation from drifting. In the code ComboBatch was ported
from, the docs claimed 39 methods in one file, 31 in another, 12 filter strategies in a
third and 14 in a fourth — each written by hand at a different time, each wrong somewhere.
Here three files are generated from the registries and this fails when the committed copy
no longer matches.

`test_generated_docs_are_current` is the whole point; the rest assert that the generated
content actually says the things the plan requires it to say.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from combobatch.docsgen import GENERATED_DOCS
from combobatch.imputation import IMPUTER_REGISTRY
from combobatch.methods import EXCLUDED_METHODS, METHOD_REGISTRY
from combobatch.metrics import METRIC_GROUP_REGISTRY

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCS = REPO_ROOT / "docs"

HAND_WRITTEN = ("OUTPUT_LAYOUT.md", "BEST_APPROACHES.md")


@pytest.mark.parametrize("filename", sorted(GENERATED_DOCS))
def test_generated_docs_are_current(filename):
    """Regenerating must reproduce the committed file byte for byte."""
    rendered = GENERATED_DOCS[filename]()
    committed = (DOCS / filename).read_text()
    assert committed == rendered, (
        f"docs/{filename} is out of date with the registries. "
        f"Run: python3.11 scripts/generate_docs.py"
    )


@pytest.mark.parametrize("filename", sorted(GENERATED_DOCS))
def test_generated_docs_say_they_are_generated(filename):
    # Without this line the next person edits the file by hand and loses the work.
    assert "Do not edit by hand" in (DOCS / filename).read_text()


@pytest.mark.parametrize("filename", HAND_WRITTEN)
def test_the_hand_written_docs_exist(filename):
    assert (DOCS / filename).exists()


@pytest.fixture(scope="module")
def methods_text():
    return (DOCS / "METHODS.md").read_text()


@pytest.fixture(scope="module")
def hyperparameters_text():
    return (DOCS / "HYPERPARAMETERS.md").read_text()


@pytest.fixture(scope="module")
def metrics_text():
    return (DOCS / "METRICS.md").read_text()


class TestMethodsDoc:

    def test_every_registered_method_appears(self, methods_text):
        for key in METHOD_REGISTRY:
            assert f"`{key}`" in methods_text

    def test_every_excluded_method_appears_with_a_reason(self, methods_text):
        for key, reason in EXCLUDED_METHODS.items():
            assert f"`{key}`" in methods_text
            assert reason in methods_text

    def test_the_count_is_not_a_hand_typed_literal(self, methods_text):
        assert f"**{len(METHOD_REGISTRY)} harmonization methods**" in methods_text

    def test_the_procrustes_licensing_rationale_is_stated(self, methods_text):
        # The reason it is absent is legal, not technical, and that distinction has to
        # survive into the public docs.
        assert "proprietary" in methods_text and "non-commercial" in methods_text

    def test_the_shambhala_substitution_is_disclosed(self, methods_text):
        assert "not bit-identical" in methods_text

    def test_the_rnaseq_guard_and_reference_rule_are_documented(self, methods_text):
        assert "RNA-seq-only guard" in methods_text
        assert "largest batch in the data" in methods_text


class TestHyperparametersDoc:

    def test_every_tunable_parameter_appears(self, hyperparameters_text):
        for spec in METHOD_REGISTRY.values():
            for name in spec.hyperparams:
                assert f"`{name}`" in hyperparameters_text
        for spec in IMPUTER_REGISTRY.values():
            for name in spec.hyperparams:
                assert f"`{name}`" in hyperparameters_text

    def test_methods_without_parameters_are_named_rather_than_omitted(
        self, hyperparameters_text
    ):
        untunable = [k for k, s in METHOD_REGISTRY.items() if not s.hyperparams]
        assert "nothing to tune" in hyperparameters_text
        for key in untunable:
            assert f"`{key}`" in hyperparameters_text

    def test_the_param_tag_rules_are_explained(self, hyperparameters_text):
        assert "run_manifest.json" in hyperparameters_text
        assert "strict__10_mnn__k50" in hyperparameters_text


class TestMetricsDoc:

    def test_every_group_has_a_section(self, metrics_text):
        for letter, spec in METRIC_GROUP_REGISTRY.items():
            assert f"### {letter} — {spec.name}" in metrics_text

    def test_every_group_declares_a_polarity_and_citation(self, metrics_text):
        for spec in METRIC_GROUP_REGISTRY.values():
            assert spec.polarity, f"group {spec.letter} has no polarity"
            assert spec.citation, f"group {spec.letter} has no citation"
            assert spec.polarity in metrics_text

    def test_the_emitted_keys_are_observed_not_promised(self, metrics_text):
        # The key lists come from running each group, so this checks a real one landed.
        assert "**Keys emitted**" in metrics_text
        assert "`r2_batch`" in metrics_text

    def test_the_group_l_ceiling_is_stated(self, metrics_text):
        assert "saturates" in metrics_text and "guard rail" in metrics_text

    def test_the_r_backed_group_is_honestly_not_enumerated(self, metrics_text):
        # Probing group F needs R, so the file would otherwise differ between a laptop
        # and the image and this suite would fail depending on where it ran.
        assert "not enumerated here" in metrics_text
