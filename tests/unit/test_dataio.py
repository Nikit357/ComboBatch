"""Table reading/writing and index alignment.

Not in the plan's Phase 1 test list, but `dataio` carries source bug #7 — the donor
writer ignored the file extension and wrote plain CSV to a `.tsv.gz` path — and a bug
fix without a regression test is not fixed.
"""

from __future__ import annotations

import gzip

import pandas as pd
import pytest

from combobatch import dataio


@pytest.fixture
def frame():
    return pd.DataFrame({"GENE1": [1.0, 2.0], "GENE2": [3.0, 4.0]}, index=["s1", "s2"])


class TestExtensionDrivesFormat:
    @pytest.mark.parametrize(
        "name,sep,gzipped",
        [
            ("x.tsv", "\t", False),
            ("x.csv", ",", False),
            ("x.tsv.gz", "\t", True),
            ("x.csv.gz", ",", True),
        ],
    )
    def test_round_trip(self, frame, tmp_path, name, sep, gzipped):
        path = str(tmp_path / name)
        dataio.write_table(frame, path)
        pd.testing.assert_frame_equal(dataio.read_table(path), frame)

    def test_gz_path_really_is_gzipped(self, frame, tmp_path):
        # The donor's writer produced an uncompressed comma-separated file here.
        path = tmp_path / "x.tsv.gz"
        dataio.write_table(frame, str(path))
        assert path.read_bytes()[:2] == b"\x1f\x8b"
        assert b"\t" in gzip.decompress(path.read_bytes())

    def test_tsv_path_really_is_tab_separated(self, frame, tmp_path):
        path = tmp_path / "x.tsv"
        dataio.write_table(frame, str(path))
        assert "\t" in path.read_text()
        assert "," not in path.read_text()

    def test_unknown_extension_raises(self, frame, tmp_path):
        with pytest.raises(dataio.DataFormatError, match="supported extensions"):
            dataio.write_table(frame, str(tmp_path / "x.weird"))


class TestBackendHelpers:
    def test_write_and_read_through_a_backend(self, frame, tmp_path):
        from combobatch import storage

        backend = storage.LocalBackend(tmp_path)
        dataio.write_table_to(frame, backend, "exp/a.tsv.gz")
        assert backend.exists("exp/a.tsv.gz")
        pd.testing.assert_frame_equal(
            dataio.read_table_from(backend, "exp/a.tsv.gz"), frame
        )


class TestAlign:
    @staticmethod
    def _pair(exp_index, ann_index):
        exp = pd.DataFrame({"G1": range(len(exp_index))}, index=exp_index)
        ann = pd.DataFrame({"batch": ["b"] * len(ann_index)}, index=ann_index)
        return exp, ann

    def test_intersection_not_position(self):
        # Deliberately different orders: matching by position would silently mislabel.
        exp, ann = self._pair(["s1", "s2", "s3"], ["s3", "s2", "s1"])
        exp_a, ann_a, _ = dataio.align(exp, ann)
        assert list(exp_a.index) == list(ann_a.index)
        assert set(exp_a.index) == {"s1", "s2", "s3"}

    def test_drops_non_overlapping_samples(self):
        exp, ann = self._pair(["s1", "s2", "s9"], ["s1", "s2", "s8"])
        exp_a, ann_a, report = dataio.align(exp, ann)
        assert len(exp_a) == len(ann_a) == 2
        assert report.n_dropped_expression_only == 1
        assert report.n_dropped_annotation_only == 1

    def test_drops_duplicates_keeping_first(self):
        exp, ann = self._pair(["s1", "s1", "s2"], ["s1", "s2"])
        exp_a, _, report = dataio.align(exp, ann)
        assert list(exp_a.index) == ["s1", "s2"]
        assert report.n_duplicate_samples == 1

    def test_no_overlap_raises_with_a_transposition_hint(self):
        exp, ann = self._pair(["s1", "s2"], ["GENE1", "GENE2"])
        with pytest.raises(ValueError, match="samples x genes"):
            dataio.align(exp, ann)

    def test_report_summary_is_readable(self):
        exp, ann = self._pair(["s1", "s2"], ["s1", "s2"])
        _, _, report = dataio.align(exp, ann)
        assert "2 samples" in report.summary()
