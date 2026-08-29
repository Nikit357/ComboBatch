"""The vendored Shambhala core, at the two points where it was repaired.

Neither test needs Octave: the subprocess call is intercepted and the preamble it would
have sent is inspected instead.
"""

from __future__ import annotations

import pathlib
import queue

import numpy as np
import pandas as pd
import pytest

from combobatch.vendor.shambhala import OCTAVE_DIR, octave_bridge, parallel
from combobatch.vendor.shambhala.octave_bridge import (
    _PIPELINE_SCRIPTS,
    _select_pipeline_script,
)
from combobatch.vendor.shambhala.progress_display import ProgressEvent


def _pool(n_genes: int = 6, nh: int = 2, np_: int = 3) -> pd.DataFrame:
    """A tiny genes x (nh + np_) pool matrix."""
    rng = np.random.default_rng(0)
    columns = [f"S{i}" for i in range(nh)] + [f"P{i}" for i in range(np_)]
    return pd.DataFrame(
        rng.exponential(10.0, size=(n_genes, nh + np_)),
        index=[f"GENE{i}" for i in range(n_genes)],
        columns=columns,
    )


class TestScriptSelection:
    """§6 bug 5 — fixed clusters silently overrode the pre-QN script choice."""

    @pytest.mark.parametrize(
        ("skip_qn", "fixed", "expected"),
        [
            (False, False, "Shambhala2_piped.m"),
            (True, False, "Shambhala2_piped_preqn.m"),
            (False, True, "Shambhala2_piped_fixed.m"),
            (True, True, "Shambhala2_piped_preqn_fixed.m"),
        ],
    )
    def test_every_combination_selects_its_own_script(self, skip_qn, fixed, expected):
        assert _select_pipeline_script(skip_qn, fixed) == expected

    def test_the_four_scripts_are_distinct(self):
        assert len(set(_PIPELINE_SCRIPTS.values())) == 4

    def test_every_selectable_script_exists_on_disk(self):
        for name in _PIPELINE_SCRIPTS.values():
            assert (pathlib.Path(OCTAVE_DIR) / name).is_file(), f"{name} is missing"

    def test_approved_speedup_combination_does_not_quantile_normalize_twice(
        self, monkeypatch
    ):
        # The A+E combination: QN precomputed in Python *and* clusters precomputed from
        # P. Upstream sent Shambhala2_piped_fixed.m, whose quantilenorm call then ran on
        # already-normalized data.
        captured: dict = {}

        class _Completed:
            returncode = 0
            stdout = ""
            stderr = ""

        def fake_run(argv, **kwargs):
            captured["preamble"] = argv[-1]
            raise _Stop()

        class _Stop(Exception):
            pass

        monkeypatch.setattr(octave_bridge.subprocess, "run", fake_run)

        with pytest.raises(_Stop):
            octave_bridge.run_octave_normalize(
                pool_df=_pool(),
                nh=2,
                np_=3,
                k=5,
                octave_scripts_dir=OCTAVE_DIR,
                fixed_clusters=np.ones(6, dtype=np.int32),
                skip_qn=True,
            )

        assert "Shambhala2_piped_preqn_fixed.m" in captured["preamble"]
        assert "FIXED_CLUSTERS" in captured["preamble"]

    def test_streaming_variant_selects_the_same_script(self, monkeypatch):
        captured: dict = {}

        class _Stop(Exception):
            pass

        def fake_popen(argv, **kwargs):
            captured["preamble"] = argv[-1]
            raise _Stop()

        monkeypatch.setattr(octave_bridge.subprocess, "Popen", fake_popen)

        with pytest.raises(_Stop):
            octave_bridge.run_octave_normalize_streaming(
                pool_df=_pool(),
                nh=2,
                np_=3,
                k=5,
                octave_scripts_dir=OCTAVE_DIR,
                fixed_clusters=np.ones(6, dtype=np.int32),
                skip_qn=True,
            )

        # The sync and streaming paths had independent copies of this branch; they now
        # share one table, and this is what stops them drifting apart again.
        assert "Shambhala2_piped_preqn_fixed.m" in captured["preamble"]


class TestPythonCublockProgress:
    """§6 bug 6 — ProgressEvent was built with field names it does not have."""

    def test_progress_events_are_constructible(self):
        pytest.importorskip("qnorm")
        rng = np.random.default_rng(0)
        pool = pd.DataFrame(
            rng.exponential(10.0, size=(120, 3)),
            index=[f"GENE{i}" for i in range(120)],
            columns=["S0", "S1", "P0"],
        )
        events: queue.Queue = queue.Queue()

        result = parallel._run_python_cublock_batch(
            pool_df=pool,
            batch_sample_names=["S0", "S1"],
            nh=2,
            k=2,
            random_seed=1,
            progress_queue=events,
            worker_idx=0,
        )

        assert result.shape == (120, 2)
        emitted = [events.get_nowait() for _ in range(events.qsize())]
        assert len(emitted) == 2
        assert all(isinstance(event, ProgressEvent) for event in emitted)
        assert [event.sample_done for event in emitted] == [1, 2]
        assert [event.sample_total for event in emitted] == [2, 2]
        assert emitted[-1].is_final is True


class TestVendoredOctaveDirectory:
    def test_octave_dir_is_inside_the_package(self):
        # Upstream looked one level above the package, so the default only resolved
        # under an editable install.
        octave_dir = pathlib.Path(OCTAVE_DIR)
        assert octave_dir.parent.name == "shambhala"
        assert octave_dir.parent.parent.name == "vendor"

    def test_every_octave_helper_is_present(self):
        required = {
            "CuBlock.m",
            "CuBlock_fixed.m",
            "kmeans.m",
            "quantilenorm.m",
            "readExpressionData.m",
        }
        present = {path.name for path in pathlib.Path(OCTAVE_DIR).glob("*.m")}
        assert required <= present
