"""Shape guards on the R source embedded in ``r_methods.py``.

The R snippets are strings, so nothing type-checks them and nothing on a laptop runs
them. Both defects below reached a finished image and were found only by running it:

* ``27_dwd`` called ``genDWD`` without two of its four required arguments. Every batch
  raised, a ``tryCatch`` handler printed a message, and the function returned its input
  **unchanged** — passing every downstream check as a successful harmonization.
* Every annotation was read without ``check.names=FALSE``, so R silently renamed any
  column that is not a syntactic name. The lookup then returned an empty vector rather
  than raising, and ``37_fabatch`` failed several steps later complaining about
  something else.

Both are cheap to pin here, and neither needs R.
"""

from __future__ import annotations

import re
from pathlib import Path

import combobatch.methods.r_methods as r_methods

SOURCE = Path(r_methods.__file__).read_text()


class TestAnnotationColumnNamesSurviveR:
    """Column names come from the user's annotation and can be anything."""

    def test_every_annotation_read_disables_name_mangling(self):
        # Line-based, not a bracket match: `r_literal(ann_path)` carries its own `)`,
        # so a non-greedy group closes before reaching the arguments that matter.
        reads = [
            line.strip()
            for line in SOURCE.splitlines()
            if "read.csv(" in line and "ann_path" in line
        ]
        assert reads, "no annotation reads found — has the helper been renamed?"
        unguarded = [r for r in reads if "check.names=FALSE" not in r]
        assert not unguarded, (
            "R renames non-syntactic columns on read — 'Cell type' becomes 'Cell.type', "
            "'1st_batch' becomes 'X1st_batch' — and the lookup then returns an empty "
            f"vector instead of raising: {unguarded}"
        )


class TestHarmonizRIsAskedForAModeThatReturnsSomething:
    """HarmonizR's own default writes an empty file and calls it success.

    Modes 1 and 3 (``mean.only = FALSE``) return nothing on a matrix without missing
    values — measured on exponential and normal data alike — and HarmonizR reports that
    by writing three bytes rather than by raising.
    """

    def test_the_mode_is_passed_at_all(self):
        assert (
            "ComBat_mode=" in SOURCE
        ), "not passing ComBat_mode leaves HarmonizR on mode 1, which returns nothing"

    def test_an_empty_result_is_refused(self):
        assert "produced an empty result" in SOURCE


class TestFailuresInsideRAreNotSwallowed:
    def test_no_error_handler_only_prints(self):
        """A handler that just prints leaves the uncorrected input in place.

        ``27_dwd`` shipped that way: the message went to the log, the matrix came back
        untouched, and it counted as a passing method.
        """
        handlers = re.findall(
            r"error\s*=\s*function\(e\)\s*\{\{(.*?)\}\}", SOURCE, re.S
        )
        assert handlers, "no R error handlers found — has the snippet style changed?"
        for body in handlers:
            records = "<<-" in body  # accumulates into a variable checked afterwards
            reraises = "stop(" in body
            assert records or reraises, (
                "this R error handler only reports; the result it was meant to write is "
                f"left as the input and the run still succeeds:\n{body.strip()[:200]}"
            )

    def test_dwd_checks_that_every_batch_was_corrected(self):
        assert "left " in SOURCE and "batch(es) uncorrected" in SOURCE

    def test_dwd_passes_every_required_gendwd_argument(self):
        """``genDWD``'s X, y, C and expon all lack defaults."""
        call = re.search(r"sol <- genDWD\((.*?)\)", SOURCE, re.S)
        assert call, "the genDWD call is missing"
        for argument in ("X =", "y =", "C =", "expon ="):
            assert argument in call.group(1), f"genDWD needs {argument}"
        assert "penaltyParameter(" in SOURCE, "C comes from DWDLargeR's own helper"
