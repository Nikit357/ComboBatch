"""The optional marker-gene panel used by groups L and M.

Loaded **lazily and only when asked for**. Importing the donor's equivalent required
``marker_gene_annotation.csv`` to exist on disk, so a run computing none of the
panel-dependent groups still failed at import if the file was missing — and the file was
a dataset-specific artefact that has no place in a dataset-agnostic tool.

Here the panel is a user-supplied CSV. With none given, groups L and M fall back to every
gene the two matrices share, which is the honest generic default: no panel means no
opinion about which genes matter.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable

# Legacy CD names still carried by older array platforms. Without these a matrix using
# the CD symbol is silently scored as missing the gene entirely.
GENE_ALIASES: dict[str, tuple[str, ...]] = {
    "MS4A1": ("CD20",),
    "MME": ("CD10",),
    "FCER2": ("CD23",),
    "CR2": ("CD21",),
    "CR1": ("CD35",),
    "PTPRC": ("CD45",),
    "PDCD1": ("CD279",),
    "CD274": ("PDL1", "PD-L1"),
    "FAS": ("TNFRSF6", "CD95"),
    "GCSAM": ("HGAL",),
    "NCAM1": ("CD56",),
    "IL2RA": ("CD25",),
    "SELL": ("CD62L",),
    "THY1": ("CD90",),
    "PECAM1": ("CD31",),
    "TNFRSF17": ("BCMA",),
    "HAVCR2": ("TIM3",),
    "VSIR": ("VISTA", "C10orf54"),
    "EBI3": ("IL27B",),
}


class PanelError(ValueError):
    """Raised when a panel file is unreadable or has no usable gene column."""


@dataclass(frozen=True)
class Panel:
    """A marker panel: the genes of interest, and which of them are controls."""

    genes: tuple[str, ...]
    housekeeping: frozenset[str] = field(default_factory=frozenset)

    def markers(self) -> tuple[str, ...]:
        """Return the non-housekeeping genes."""
        return tuple(gene for gene in self.genes if gene not in self.housekeeping)

    def __len__(self) -> int:
        return len(self.genes)


def load_panel(uri: str, *, endpoint_url: str | None = None) -> Panel:
    """
    Load a marker panel from a CSV or TSV file, locally or on S3.

    The file needs a ``gene`` column (or a single unnamed column of symbols) and may
    carry a boolean ``is_housekeeping`` column. Repeated genes are collapsed, since a
    gene may appear once per signature it belongs to.

    Parameters
    ----------
    uri
        Path or ``s3://`` URI.
    endpoint_url
        Alternative S3 endpoint.

    Returns
    -------
    The parsed :class:`Panel`.

    Raises
    ------
    PanelError
        If no gene column can be identified.
    """
    from combobatch import dataio

    frame = dataio.read_table(uri, endpoint_url=endpoint_url).reset_index()
    lowered = {str(name).lower(): name for name in frame.columns}

    gene_col = lowered.get("gene") or lowered.get("genes") or lowered.get("symbol")
    if gene_col is None:
        if frame.shape[1] != 1:
            raise PanelError(
                f"{uri}: no 'gene' column, and the file has {frame.shape[1]} columns so "
                f"the intended one is ambiguous. Name a column 'gene'."
            )
        gene_col = frame.columns[0]

    genes = [str(value).strip() for value in frame[gene_col] if str(value).strip()]
    if not genes:
        raise PanelError(f"{uri}: the gene column is empty")

    housekeeping: set[str] = set()
    hk_col = lowered.get("is_housekeeping") or lowered.get("housekeeping")
    if hk_col is not None:
        flags = frame[hk_col].astype(str).str.lower().isin({"true", "1", "yes", "y"})
        housekeeping = {
            str(gene).strip()
            for gene in frame.loc[flags, gene_col]
            if str(gene).strip()
        }

    return Panel(
        genes=tuple(dict.fromkeys(genes)), housekeeping=frozenset(housekeeping)
    )


def resolve_panel(
    available: Iterable[str], panel: Panel | None
) -> tuple[list[str], list[str]]:
    """
    Match a panel against the genes actually present in a matrix.

    Falls back to :data:`GENE_ALIASES` when the primary symbol is absent, so a matrix
    carrying a legacy CD name is not penalised for it.

    Parameters
    ----------
    available
        Gene symbols present in the expression matrix.
    panel
        The panel, or ``None`` to use every available gene.

    Returns
    -------
    ``(resolved, missing)`` — the symbols to use as they appear in ``available``, and
    the panel entries that matched nothing under any alias.
    """
    present = set(available)
    if panel is None:
        return sorted(present), []

    resolved: list[str] = []
    missing: list[str] = []
    for gene in panel.genes:
        if gene in present:
            resolved.append(gene)
            continue
        alias = next((a for a in GENE_ALIASES.get(gene, ()) if a in present), None)
        if alias is not None:
            resolved.append(alias)
        else:
            missing.append(gene)
    return sorted(set(resolved)), sorted(set(missing))
