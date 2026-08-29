"""Aggregate per-output metric sidecars into analysis tables.

Groups L, M and N produce nested dicts — per gene, per cohort, per fold. Those must never
reach the wide table: pandas would stringify them into unusable cells. They are split
into long-format tables instead, and anything nested that survives is dropped loudly
rather than silently stringified.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable

import pandas as pd

from combobatch import dataio, storage
from combobatch.logging_utils import get_logger
from combobatch.params import KeyParseError, parse_output_key

# Metadata columns that lead the wide table, in this order.
METADATA_COLUMNS = [
    "key",
    "imp",
    "method",
    "param_tag",
    "post_rm",
    "harshness",
    "status",
    "compute_time_s",
    "n_samples",
    "n_genes",
]

GENE_DICT_KEY = "mk_rho_by_gene"
GENE_COUNT_KEY = "mk_rho_n_cohorts_by_gene"
COHORT_DICT_KEY = "mk_rho_by_cohort"
DETAIL_DICT_KEY = "mk_gene_cohort_detail"
FOLD_DICT_KEYS = {"pv_lobo3_folds": "3class", "pv_lobo2_folds": "2class"}
# Small enough to pivot into the wide table as one column per level.
BIO_DICT_KEY = "xb_rank_agree_by_bio"


@dataclass
class ConcatResult:
    """The four tables ``concat`` produces."""

    summary: pd.DataFrame
    gene_correlations: pd.DataFrame = field(default_factory=pd.DataFrame)
    cohort_correlations: pd.DataFrame = field(default_factory=pd.DataFrame)
    prediction_folds: pd.DataFrame = field(default_factory=pd.DataFrame)

    def tables(self) -> dict[str, pd.DataFrame]:
        """Return ``{stem: frame}`` for every non-empty table."""
        candidates = {
            "metrics_summary": self.summary,
            "marker_gene_correlations_long": self.gene_correlations,
            "marker_cohort_correlations_long": self.cohort_correlations,
            "prediction_folds_long": self.prediction_folds,
        }
        return {stem: frame for stem, frame in candidates.items() if not frame.empty}


def list_metric_keys(backend: storage.StorageBackend) -> list[str]:
    """Return every metrics sidecar key under ``metrics/``, sorted."""
    return sorted(
        key for key in backend.list_keys("metrics") if key.endswith("_metrics.json")
    )


def split_nested(record: dict, key: str) -> tuple[list[dict], list[dict], list[dict]]:
    """
    Pop the nested-dict values out of one record into long-format rows.

    Mutates ``record`` so the wide table only ever sees scalars.

    Parameters
    ----------
    record
        One parsed sidecar.
    key
        The output key, used as the join column of every long table.

    Returns
    -------
    ``(gene_rows, cohort_rows, fold_rows)``.
    """
    gene_rows: list[dict] = []
    cohort_rows: list[dict] = []
    fold_rows: list[dict] = []

    by_gene = record.pop(GENE_DICT_KEY, None) or {}
    cohort_counts = record.pop(GENE_COUNT_KEY, None) or {}
    for gene, rho in by_gene.items():
        gene_rows.append(
            {
                "key": key,
                "gene": gene,
                "rho": rho,
                "n_cohorts": cohort_counts.get(gene),
            }
        )

    for cohort, rho in (record.pop(COHORT_DICT_KEY, None) or {}).items():
        cohort_rows.append({"key": key, "cohort": cohort, "rho": rho})

    # The per-gene-per-cohort detail is the largest object a sidecar can carry; it is
    # folded into the gene table rather than kept as a column.
    for cohort, per_gene in (record.pop(DETAIL_DICT_KEY, None) or {}).items():
        for gene, rho in per_gene.items():
            gene_rows.append({"key": key, "gene": gene, "rho": rho, "cohort": cohort})

    for dict_key, target in FOLD_DICT_KEYS.items():
        for batch, values in (record.pop(dict_key, None) or {}).items():
            fold_rows.append(
                {
                    "key": key,
                    "target": target,
                    "batch": batch,
                    "n": values.get("n"),
                    "n_classes": values.get("n_classes"),
                    "f1_macro": values.get("f1_macro"),
                    "auc": values.get("auc"),
                }
            )

    for level, value in (record.pop(BIO_DICT_KEY, None) or {}).items():
        record[f"xb_rank_agree_{level}"] = value

    return gene_rows, cohort_rows, fold_rows


def drop_remaining_containers(rows: list[dict]) -> set[str]:
    """
    Strip any dict or list value that survived, so the wide table stays scalar-only.

    Returns
    -------
    The keys dropped, so the caller can report them rather than losing them silently.
    """
    offenders: set[str] = set()
    for row in rows:
        for key in [
            name for name, value in row.items() if isinstance(value, (dict, list))
        ]:
            del row[key]
            offenders.add(key)
    return offenders


def concat_metrics(
    out_uri: str, *, endpoint_url: str | None = None, strict_keys: bool = True
) -> ConcatResult:
    """
    Read every metrics sidecar under an output root and build the analysis tables.

    Parameters
    ----------
    out_uri
        The run's output root, local or ``s3://``.
    endpoint_url
        Alternative S3 endpoint.
    strict_keys
        Raise on a sidecar whose filename does not parse. The donor dropped those
        silently, so a malformed key meant a result that vanished from every table
        without a word.

    Returns
    -------
    The four tables, the last three empty unless L, M or N ran.

    Raises
    ------
    KeyParseError
        On an unparseable sidecar name, when ``strict_keys`` is set.
    """
    log = get_logger()
    backend = storage.backend_for_root(out_uri, endpoint_url=endpoint_url)

    records: list[dict] = []
    gene_rows: list[dict] = []
    cohort_rows: list[dict] = []
    fold_rows: list[dict] = []

    for storage_key in list_metric_keys(backend):
        stem = storage_key.split("/")[-1][: -len("_metrics.json")]
        try:
            parse_output_key(stem)
        except KeyParseError as exc:
            if strict_keys:
                raise KeyParseError(
                    f"{storage_key}: {exc}. Fix or remove the file; ComboBatch will not "
                    f"silently drop it from the summary."
                ) from exc
            log.warning("skipping unparseable sidecar %s: %s", storage_key, exc)
            continue

        record = dataio.read_json(backend, storage_key)
        genes, cohorts, folds = split_nested(record, stem)
        gene_rows.extend(genes)
        cohort_rows.extend(cohorts)
        fold_rows.extend(folds)
        record.setdefault("key", stem)
        records.append(record)

    dropped = drop_remaining_containers(records)
    if dropped:
        log.warning("dropped non-scalar key(s) from the summary: %s", sorted(dropped))

    summary = pd.DataFrame(records)
    if not summary.empty:
        leading = [name for name in METADATA_COLUMNS if name in summary.columns]
        rest = sorted(name for name in summary.columns if name not in leading)
        summary = summary[[*leading, *rest]].sort_values("key").reset_index(drop=True)

    return ConcatResult(
        summary=summary,
        gene_correlations=pd.DataFrame(gene_rows),
        cohort_correlations=pd.DataFrame(cohort_rows),
        prediction_folds=pd.DataFrame(fold_rows),
    )


def write_tables(
    result: ConcatResult,
    *,
    out_dir: str | None = None,
    date_tag: str | None = None,
    backend: storage.StorageBackend | None = None,
) -> list[str]:
    """
    Write the tables, to the run's output root and optionally to a local directory.

    Parameters
    ----------
    result
        What :func:`concat_metrics` returned.
    out_dir
        Extra local directory for dated copies, so a rerun never overwrites a snapshot.
    date_tag
        Suffix for those dated copies. Defaults to today in ``YYMMDD``.
    backend
        The run's output root. ``None`` writes only the local copies.

    Returns
    -------
    Every URI written.
    """
    import time

    written: list[str] = []
    tag = date_tag or time.strftime("%y%m%d")

    for stem, frame in result.tables().items():
        payload = frame.to_csv(index=False).encode()
        if backend is not None:
            backend.write_bytes(f"{stem}.csv", payload)
            written.append(backend.uri(f"{stem}.csv"))
        if out_dir:
            from pathlib import Path

            path = Path(out_dir) / f"{stem}_{tag}.csv"
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(payload)
            written.append(str(path))
    return written
