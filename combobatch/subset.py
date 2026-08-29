"""Sample subsetting via a pandas query over the user's own annotation columns.

This replaces the donor pipeline's fourteen hardcoded filter strategies, which encoded
one dataset's batch names and diagnosis labels and were meaningless anywhere else. A
query expresses the same thing in one line:

    C_rnaseq_only  ->  "RNA_BATCH.str.startswith('RNASeq')"
    K_ffpe_only    ->  "RNA_BATCH.str.contains('_FFPE_')"
"""

from __future__ import annotations

import re

import pandas as pd

from combobatch.logging_utils import get_logger

# Bare column references in a query expression. Used only to produce a helpful error
# before pandas raises a much less informative one.
_IDENTIFIER_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\b")

# Words that appear in query syntax rather than naming a column.
_QUERY_KEYWORDS = frozenset(
    {
        "and",
        "or",
        "not",
        "in",
        "True",
        "False",
        "None",
        "str",
        "startswith",
        "endswith",
        "contains",
        "isin",
        "isna",
        "notna",
        "lower",
        "upper",
        "strip",
        "len",
        "index",
    }
)


class SubsetQueryError(ValueError):
    """Raised when a ``--subset-query`` is invalid or selects nothing."""


def referenced_columns(query: str) -> set[str]:
    """
    Return the identifiers in ``query`` that look like column references.

    Best-effort: string literals and known query keywords are excluded, so the result is
    suitable for producing an error message but not for semantic analysis.
    """
    without_strings = re.sub(r"'[^']*'|\"[^\"]*\"", " ", query)
    without_backticks = re.sub(r"`[^`]*`", " ", without_strings)
    found = set(_IDENTIFIER_RE.findall(without_backticks))
    return {name for name in found if name not in _QUERY_KEYWORDS}


def validate_subset_query(query: str, ann_df: pd.DataFrame) -> None:
    """
    Check a query against an annotation table before any expensive work happens.

    Parameters
    ----------
    query
        A pandas ``DataFrame.query`` expression.
    ann_df
        The annotation the query will run against.

    Raises
    ------
    SubsetQueryError
        If the query is empty or names a column the annotation does not have.
    """
    if not query or not query.strip():
        raise SubsetQueryError("--subset-query must not be empty")

    unknown = sorted(referenced_columns(query) - set(ann_df.columns))
    if unknown:
        import difflib

        details = []
        for name in unknown:
            close = difflib.get_close_matches(name, list(ann_df.columns), n=1)
            hint = f" (did you mean {close[0]!r}?)" if close else ""
            details.append(f"{name!r}{hint}")
        raise SubsetQueryError(
            f"--subset-query references unknown column(s): {', '.join(details)}. "
            f"The annotation has {len(ann_df.columns)} columns; the first few are: "
            f"{', '.join(list(ann_df.columns)[:10])}"
        )


def apply_subset_query(ann_df: pd.DataFrame, query: str) -> pd.DataFrame:
    """
    Apply a subset query and return the surviving annotation rows.

    Parameters
    ----------
    ann_df
        Annotation to filter.
    query
        A pandas ``DataFrame.query`` expression. ``engine="python"`` is used so that
        ``.str`` accessors work, which the common batch-name predicates rely on.

    Returns
    -------
    The filtered annotation.

    Raises
    ------
    SubsetQueryError
        If the query is invalid, fails to evaluate, or matches no rows. Selecting
        nothing is treated as an error rather than an empty run, because it is almost
        always a typo in a batch name.
    """
    validate_subset_query(query, ann_df)

    try:
        subset = ann_df.query(query, engine="python")
    except Exception as exc:
        raise SubsetQueryError(
            f"--subset-query failed to evaluate: {exc}. Query was: {query!r}"
        ) from exc

    if len(subset) == 0:
        raise SubsetQueryError(
            f"--subset-query matched 0 of {len(ann_df)} samples. Query was: {query!r}"
        )

    get_logger().info("Subset query kept %d of %d samples", len(subset), len(ann_df))
    return subset
