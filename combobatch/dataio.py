"""Reading and writing expression and annotation tables.

Format is driven by the file extension, in both directions. The donor Shambhala writer
ignored the extension and always emitted uncompressed CSV, so ``--output foo.tsv.gz``
silently produced a comma-separated plain-text file; every write here goes through
:func:`write_table` instead.
"""

from __future__ import annotations

import gzip
import io
import json
import numbers
from dataclasses import dataclass
from typing import Any, Mapping

import pandas as pd

from combobatch import storage
from combobatch.logging_utils import get_logger

# Extension -> (pandas separator, gzip?). Longest suffixes first so ".tsv.gz" is tested
# before ".gz".
_FORMATS: tuple[tuple[str, str, bool], ...] = (
    (".tsv.gz", "\t", True),
    (".csv.gz", ",", True),
    (".txt.gz", "\t", True),
    (".tsv", "\t", False),
    (".csv", ",", False),
    (".txt", "\t", False),
)


class DataFormatError(ValueError):
    """Raised when a path's extension is not a supported table format."""


def _format_for(path: str) -> tuple[str, bool]:
    """Return ``(separator, gzipped)`` for a path, by extension."""
    lowered = path.lower()
    for suffix, separator, gzipped in _FORMATS:
        if lowered.endswith(suffix):
            return separator, gzipped
    supported = ", ".join(suffix for suffix, _, _ in _FORMATS)
    raise DataFormatError(
        f"cannot determine table format for {path!r}; supported extensions: {supported}"
    )


def read_table(uri: str, *, endpoint_url: str | None = None) -> pd.DataFrame:
    """
    Read a delimited table from a local path or an ``s3://`` URI.

    The first column is always used as the index, matching the layout every matrix in
    this project uses.

    Parameters
    ----------
    uri
        Local path or ``s3://bucket/key``. Parquet is handled separately.
    endpoint_url
        Alternative S3 endpoint.

    Returns
    -------
    The table, indexed by its first column.
    """
    if uri.lower().endswith(".parquet"):
        backend, key = storage.resolve(uri, endpoint_url=endpoint_url)
        return pd.read_parquet(io.BytesIO(backend.read_bytes(key)))

    separator, gzipped = _format_for(uri)
    backend, key = storage.resolve(uri, endpoint_url=endpoint_url)
    raw = backend.read_bytes(key)
    if gzipped:
        raw = gzip.decompress(raw)
    return pd.read_csv(io.BytesIO(raw), sep=separator, index_col=0, low_memory=False)


def write_table(
    frame: pd.DataFrame, uri: str, *, endpoint_url: str | None = None
) -> None:
    """
    Write a table to a local path or an ``s3://`` URI, honouring the extension.

    Parameters
    ----------
    frame
        The table to write. Its index is written as the first column.
    uri
        Destination. ``.gz`` compresses; the separator follows the extension.
    endpoint_url
        Alternative S3 endpoint.
    """
    if uri.lower().endswith(".parquet"):
        buffer = io.BytesIO()
        frame.to_parquet(buffer)
        payload = buffer.getvalue()
    else:
        separator, gzipped = _format_for(uri)
        text = frame.to_csv(sep=separator)
        payload = gzip.compress(text.encode()) if gzipped else text.encode()

    backend, key = storage.resolve(uri, endpoint_url=endpoint_url)
    backend.write_bytes(key, payload)


def write_table_to(
    frame: pd.DataFrame, backend: storage.StorageBackend, key: str
) -> None:
    """Write a table to an explicit backend and key, honouring the key's extension."""
    if key.lower().endswith(".parquet"):
        buffer = io.BytesIO()
        frame.to_parquet(buffer)
        backend.write_bytes(key, buffer.getvalue())
        return

    separator, gzipped = _format_for(key)
    text = frame.to_csv(sep=separator)
    backend.write_bytes(key, gzip.compress(text.encode()) if gzipped else text.encode())


def read_table_from(backend: storage.StorageBackend, key: str) -> pd.DataFrame:
    """Read a table from an explicit backend and key, honouring the key's extension."""
    raw = backend.read_bytes(key)
    if key.lower().endswith(".parquet"):
        return pd.read_parquet(io.BytesIO(raw))

    separator, gzipped = _format_for(key)
    if gzipped:
        raw = gzip.decompress(raw)
    return pd.read_csv(io.BytesIO(raw), sep=separator, index_col=0, low_memory=False)


def json_ready(value: Any) -> Any:
    """
    Convert a value into something ``json.dumps`` can serialize strictly.

    NaN and infinity become ``None``, and numpy scalars become their Python
    equivalents. The donor wrote sidecars with a bare ``json.dump``, which emits
    ``NaN`` — accepted by Python's own parser and rejected by every strict one, so the
    files could not be read by anything else.
    """
    import math

    if isinstance(value, Mapping):
        return {str(key): json_ready(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, set, frozenset)):
        return [json_ready(item) for item in value]
    if isinstance(value, (bool, str, type(None))):
        return value
    if isinstance(value, numbers.Integral):
        return int(value)
    if isinstance(value, numbers.Real):
        number = float(value)
        return None if math.isnan(number) or math.isinf(number) else number
    return str(value)


def dumps_json(payload: Any, *, indent: int | None = 2) -> str:
    """Serialize to strict JSON, with NaN and infinity rendered as ``null``."""
    return json.dumps(json_ready(payload), indent=indent, sort_keys=False)


def write_json(payload: Any, backend: storage.StorageBackend, key: str) -> None:
    """Write a JSON document to a backend key, NaN-safely."""
    backend.write_bytes(key, dumps_json(payload).encode())


def read_json(backend: storage.StorageBackend, key: str) -> Any:
    """Read a JSON document from a backend key."""
    return json.loads(backend.read_bytes(key).decode())


@dataclass(frozen=True)
class AlignmentReport:
    """What :func:`align` had to discard to make the two tables agree."""

    n_samples: int
    n_genes: int
    n_dropped_expression_only: int
    n_dropped_annotation_only: int
    n_duplicate_samples: int

    def summary(self) -> str:
        """One-line human-readable summary, for the run log."""
        return (
            f"{self.n_samples} samples x {self.n_genes} genes "
            f"(dropped {self.n_dropped_expression_only} expression-only, "
            f"{self.n_dropped_annotation_only} annotation-only, "
            f"{self.n_duplicate_samples} duplicate)"
        )


def align(
    exp_df: pd.DataFrame, ann_df: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame, AlignmentReport]:
    """
    Align expression and annotation by index intersection.

    Never positional: the two tables are matched on sample identifiers, duplicates are
    dropped keeping the first occurrence, and both are returned in the same order.

    Parameters
    ----------
    exp_df
        Expression, samples x genes.
    ann_df
        Annotation, indexed by sample.

    Returns
    -------
    ``(exp, ann, report)`` with identical indices.

    Raises
    ------
    ValueError
        If the two indices do not intersect at all — almost always a transposed matrix
        or a mismatched pair of files, and worth failing on immediately.
    """
    exp_dupes = int(exp_df.index.duplicated().sum())
    ann_dupes = int(ann_df.index.duplicated().sum())
    exp_df = exp_df[~exp_df.index.duplicated(keep="first")]
    ann_df = ann_df[~ann_df.index.duplicated(keep="first")]

    common = exp_df.index.intersection(ann_df.index)
    if len(common) == 0:
        raise ValueError(
            "expression and annotation share no sample identifiers. Check that the "
            "expression matrix is samples x genes (samples as rows) and that both "
            "files describe the same dataset."
        )

    report = AlignmentReport(
        n_samples=len(common),
        n_genes=exp_df.shape[1],
        n_dropped_expression_only=len(exp_df.index.difference(ann_df.index)),
        n_dropped_annotation_only=len(ann_df.index.difference(exp_df.index)),
        n_duplicate_samples=exp_dupes + ann_dupes,
    )

    exp_aligned = exp_df.loc[common]
    ann_aligned = ann_df.loc[common]
    assert len(exp_aligned) == len(ann_aligned)

    logger = get_logger()
    if report.n_dropped_expression_only or report.n_dropped_annotation_only:
        logger.info("Aligned inputs: %s", report.summary())

    return exp_aligned, ann_aligned, report
