"""Vendored Shambhala2 compute core (pure Python + Octave).

Ported from ``Follicular_lymphoma_disser/shambhala_adoption/Shambhala_containerized/``
at 2026-08-29. The algorithm modules are carried over as-is; the deviations are:

* imports rewritten to this package,
* ``octave/`` moved *inside* the package so it ships as package data (upstream left it
  outside, so the default script directory only resolved under ``pip install -e .``),
* three defects fixed, each marked ``VENDOR FIX`` at the site — the ``skip_qn`` /
  ``fixed_clusters`` script-selection clash, the ``ProgressEvent`` field mismatch, and
  the script-directory default above,
* ``io_utils`` deliberately not vendored: its writer ignored the file extension, so
  ``foo.tsv.gz`` was silently written as plain CSV. All ComboBatch I/O goes through
  :mod:`combobatch.dataio`.

Nothing outside :mod:`combobatch.methods.shambhala_method` should import from here.
"""

from __future__ import annotations

import pathlib

OCTAVE_DIR: str = str(pathlib.Path(__file__).resolve().parent / "octave")

__all__ = ["OCTAVE_DIR"]
