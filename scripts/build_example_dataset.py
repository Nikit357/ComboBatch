#!/usr/bin/env python3
"""
Regenerate ``examples/example_dataset/`` from the dissertation's source tables.

Run once by a maintainer; the output is committed, so an ordinary user never needs S3
access or this script. The source tables are read-only donors and are never written to.

    python3.11 scripts/build_example_dataset.py --dry-run     # selection only, no S3
    python3.11 scripts/build_example_dataset.py               # writes the dataset

Five public GEO cohorts are selected to exercise every axis the tool benchmarks at once:
FL and DLBCL, RNA-seq and three microarray platforms, FF and FFPE preservation.

Two properties of the result are deliberate and easy to mistake for defects:

* **The matrix is ragged.** Only 6,797 of 21,890 genes are measured in all five cohorts,
  and thousands of columns are entirely NA. The missingness is platform-driven, which is
  exactly what makes the dataset a real test of the imputation axis - a dense intersection
  matrix would make strict, knn, softimpute and missforest produce near-identical results.
  Being platform-driven, it is also *banded* rather than scattered: `max_na_frac` only
  matters where it crosses a band, and the default 0.20 falls below every one of them.
  PROVENANCE.md tabulates the bands.
* **The values are raw, not log-scale.** ``20_shambhala`` requires that, and nothing here
  transforms them.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from datetime import date
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT_DIR = REPO_ROOT / "examples" / "example_dataset"

DEFAULT_ANN = str(
    REPO_ROOT.parent / "Follicular_lymphoma_disser" / "comb_ann_unified.csv"
)
DEFAULT_BUCKET = "internal-control-data"
DEFAULT_EXP_KEY = "FL_batch_correction/exp/comb_exp.tsv"

# The exact COHORT_LABEL values, with the GEO series each one comes from. GSE69033 carries
# no PUB_ prefix in the source table; that asymmetry is real, not a typo.
COHORTS: dict[str, str] = {
    "PUB_DLBCL_GSE64555": "GSE64555",
    "PUB_FL_GSE148070": "GSE148070",
    "PUB_Health_B_cells_GSE119234": "GSE119234",
    "GSE69033": "GSE69033",
    "PUB_FL_GSE62241": "GSE62241",
}

# An artefact of a previous round-trip through pandas, not annotation.
DROP_COLUMNS = {"Unnamed: 0.1"}


def log(message: str) -> None:
    """Progress to stderr, so stdout stays clean."""
    print(message, file=sys.stderr, flush=True)


def select_samples(ann: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    """
    Filter the annotation to the five cohorts and drop the duplicated GSE62241 block.

    GSE62241 appears twice: once keyed ``SRX...-PUB_FL_GSE62241`` with full metadata and
    once as a bare ``SRX...`` with fewer fields populated. The two copies are the same
    experiments - and they disagree about the diagnosis for four of them, so keeping both
    would put contradictory labels on identical expression profiles.

    The rule is derived rather than hardcoded: group by the accession before the first
    hyphen and keep whichever row carries more annotation. Anything else that collides
    aborts the build, because a second, unexamined duplicate would be silently resolved by
    a rule that was only ever justified for this one.

    Parameters
    ----------
    ann
        The full source annotation table, indexed by sample id.

    Returns
    -------
    ``(kept, dropped, report)`` - the deduplicated annotation, the rows removed, and a
    dict of findings for PROVENANCE.md.
    """
    selected = ann[ann["COHORT_LABEL"].isin(COHORTS)]
    if len(selected) == 0:
        raise SystemExit(
            "no rows matched the five cohort labels; wrong annotation file?"
        )

    base_ids = selected.index.to_series().str.split("-").str[0]
    collisions = sorted(base_ids[base_ids.duplicated(keep=False)].unique())
    unexpected = [base for base in collisions if not base.startswith("SRX")]
    if unexpected:
        raise SystemExit(
            f"unexpected duplicate accessions {unexpected}; the keep-the-richer-row rule "
            f"was only justified for the GSE62241 SRX block - inspect before proceeding"
        )

    disagreements = {}
    for base, group in selected.groupby(base_ids):
        labels = group["Diagnosis_cell_type_unified"]
        if len(group) > 1 and labels.nunique() > 1:
            disagreements[base] = labels.to_dict()

    ranked = selected.assign(
        _base=base_ids, _filled=selected.notna().sum(axis=1)
    ).sort_values("_filled", ascending=False)
    kept_index = ranked[~ranked["_base"].duplicated(keep="first")].index
    kept = selected.loc[sorted(kept_index)]
    dropped = selected.drop(index=kept_index)

    report = {
        "n_source_rows": len(selected),
        "n_kept": len(kept),
        "duplicate_accessions": collisions,
        "dropped_ids": sorted(dropped.index),
        "diagnosis_disagreements": disagreements,
    }
    return kept, dropped, report


def reduce_columns(ann: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only the annotation columns populated for every selected sample.

    The source table has 483 columns, most of them FL-specific clinical fields (``coo``,
    ``mfp``, ``dhitsig``, treatment and survival) populated for one cohort only. Shipping
    them would make an example dataset that looks like it needs an FL study design.
    """
    complete = [
        column
        for column in ann.columns
        if column not in DROP_COLUMNS and ann[column].notna().all()
    ]
    return ann[complete]


def fetch_expression(
    sample_ids: list[str], out_path: Path, *, bucket: str, key: str
) -> dict:
    """
    Stream the source expression TSV and write only the requested rows, gzipped.

    The source is 1.9 GB of samples x genes. Rows are copied as raw bytes rather than
    parsed into a DataFrame, so the values that land in the example dataset are exactly
    the published ones - no float reformatting, no NA-marker translation.

    Returns
    -------
    A dict of shape and size facts for PROVENANCE.md.

    Raises
    ------
    SystemExit
        If any requested sample is absent from the matrix.
    """
    import boto3

    wanted = set(sample_ids)
    found: dict[str, bytes] = {}
    header: bytes | None = None
    scanned = 0

    log(f"streaming s3://{bucket}/{key} for {len(wanted)} samples ...")
    body = boto3.client("s3").get_object(Bucket=bucket, Key=key)["Body"]
    for line in body.iter_lines(chunk_size=4 * 1024 * 1024):
        if header is None:
            header = line
            continue
        scanned += 1
        sample_id = line.split(b"\t", 1)[0].decode()
        if sample_id in wanted:
            found[sample_id] = line
            if len(found) == len(wanted):
                # Every wanted row is in hand; the rest of the 1.9 GB is not needed.
                log(f"  all samples found after {scanned} rows")
                break
        if scanned % 1000 == 0:
            log(f"  {scanned} rows scanned, {len(found)}/{len(wanted)} found")
    body.close()

    missing = sorted(wanted - set(found))
    if missing:
        raise SystemExit(
            f"{len(missing)} annotated samples are absent from the expression matrix: "
            f"{missing[:5]}{' ...' if len(missing) > 5 else ''}"
        )

    assert header is not None
    n_genes = header.count(b"\t")
    uncompressed = len(header) + 1 + sum(len(found[s]) + 1 for s in sample_ids)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    # mtime=0: the gzip header otherwise embeds the build time, so an unchanged rebuild
    # would produce a different file and show up as a spurious diff.
    with gzip.GzipFile(out_path, "wb", compresslevel=9, mtime=0) as handle:
        handle.write(header + b"\n")
        for sample_id in sample_ids:
            handle.write(found[sample_id] + b"\n")

    return {
        "n_samples": len(sample_ids),
        "n_genes": n_genes,
        "uncompressed_bytes": uncompressed,
        "gzipped_bytes": out_path.stat().st_size,
        # Deliberately not carried into PROVENANCE.md: the committed file is public, and
        # an internal bucket URI does not belong in it.
        "source": f"s3://{bucket}/{key}",
    }


def describe_coverage(exp_path: Path) -> dict:
    """Measure the gene-coverage pattern that makes this dataset an imputation test."""
    exp = pd.read_csv(exp_path, sep="\t", index_col=0)
    per_gene_present = exp.notna().sum()
    na_fraction = exp.isna().mean().round(6)
    # Platform-driven missingness is banded, not scattered: a gene is measured by a whole
    # cohort or by none of it, so the per-gene NA fraction takes a handful of values
    # rather than a continuum. Reporting the bands is what tells a reader where
    # max_na_frac actually has an effect - between bands it does nothing at all.
    bands = [
        (float(frac), int(n))
        for frac, n in na_fraction.value_counts().sort_index().items()
    ]
    interior = [frac for frac, _ in bands if 0.0 < frac < 1.0]
    return {
        "shape": exp.shape,
        "all_na_genes": int((per_gene_present == 0).sum()),
        "complete_genes": int((per_gene_present == len(exp)).sum()),
        "overall_na_fraction": float(exp.isna().to_numpy().mean()),
        "value_max": float(exp.max().max()),
        "na_bands": bands,
        "first_gap": interior[0] if interior else 0.0,
        "genes_at_04": int((na_fraction <= 0.4).sum()),
        "genes_strict": int((na_fraction == 0).sum()),
    }


def write_provenance(
    path: Path, ann: pd.DataFrame, report: dict, exp_stats: dict, coverage: dict
) -> None:
    """Emit PROVENANCE.md: what was taken, from where, and every judgement call made."""
    per_cohort = ann.groupby("COHORT_LABEL", observed=True)
    rows = []
    for label, accession in COHORTS.items():
        group = per_cohort.get_group(label)
        batches = ", ".join(sorted(group["RNA_BATCH"].unique()))
        platform = ", ".join(sorted(group["PLATFORM_RNA"].astype(str).unique()))
        rows.append(
            f"| [{accession}](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession})"
            f" | `{label}` | {len(group)} | `{batches}` | {platform} |"
        )

    diagnoses = ann["Diagnosis_cell_type_unified"].value_counts()
    diagnosis_rows = "\n".join(
        f"| {name} | {count} |" for name, count in diagnoses.items()
    )

    disagreement_rows = "\n".join(
        f"| `{base}` | "
        + " vs ".join(f"`{sid}` = {label}" for sid, label in labels.items())
        + " |"
        for base, labels in sorted(report["diagnosis_disagreements"].items())
    )

    path.write_text(f"""# Provenance of `examples/example_dataset/`

Generated by `scripts/build_example_dataset.py` on {date.today().isoformat()}.
Everything here is derived from public GEO series.

## Contents

| File | Shape | Size |
|---|---|---|
| `example_exp.tsv.gz` | {coverage['shape'][0]} samples x {coverage['shape'][1]} genes | {exp_stats['gzipped_bytes'] / 1e6:.2f} MB gzipped ({exp_stats['uncompressed_bytes'] / 1e6:.1f} MB raw) |
| `example_ann.csv` | {ann.shape[0]} samples x {ann.shape[1]} columns | {(path.parent / 'example_ann.csv').stat().st_size / 1e3:.0f} kB |

Expression is **samples x genes**, tab-separated, gzipped - the orientation the whole tool
uses. Values are **raw, not log-transformed**: `20_shambhala` requires raw input and
applies its own `log2(x+1)` internally. The maximum value here is
{coverage['value_max']:,.0f}.

## Source series

| Accession | `COHORT_LABEL` | n | `RNA_BATCH` | Platform |
|---|---|---|---|---|
{chr(10).join(rows)}

Each series remains the property of its depositors and is subject to the terms of the GEO
repository; cite the original publications when using this dataset for anything beyond
exercising this tool.

Expression values were extracted from the dissertation's combined expression matrix, which
is itself an internal reassembly of these same GEO series. Its location is a default in
`scripts/build_example_dataset.py`; regenerating this dataset needs read access to it.

## Composition

| Diagnosis / cell type | n |
|---|---|
{diagnosis_rows}

Four platforms across two preservation types, so the dataset exercises cross-platform
harmonization, FF-vs-FFPE mixing and malignant-vs-normal biology at once.

## Two judgement calls

**1. GSE62241 was deduplicated: {report['n_source_rows']} rows in, {report['n_kept']} out.**
Its {len(report['duplicate_accessions'])} experiments appear twice in the source table -
once keyed `SRX...-PUB_FL_GSE62241` with full metadata, once as a bare `SRX...` with fewer
fields. The build keeps the richer row of each pair. This is not cosmetic: the two copies
disagree about the diagnosis for four samples, so keeping both would attach contradictory
labels to identical expression profiles.

| Accession | Conflict |
|---|---|
{disagreement_rows}

**2. The ragged gene matrix was kept, not intersected.** Of {coverage['shape'][1]} genes,
{coverage['all_na_genes']} are entirely absent across these samples and only
{coverage['complete_genes']} are measured in every one; overall
{coverage['overall_na_fraction']:.1%} of the matrix is NA. The missingness is
platform-driven - GSE69033 and GSE119234 share one gene subspace, the other three another.
Intersecting to a dense matrix would have made `strict`, `knn`, `softimpute` and
`missforest` produce near-identical results and defeated the imputation axis entirely.

`strict` therefore keeps only the genes present everywhere. **Raise `max_na_frac` above
{coverage['first_gap']:.2f} or the other three imputers do nothing at all.**

Platform-driven missingness is banded rather than scattered — a gene is measured by a
whole cohort or by none of it — so the per-gene NA fraction takes just
{len(coverage['na_bands'])} distinct values:

| NA fraction | Genes | Meaning |
|---|---|---|
{chr(10).join(
    f"| {frac:.4f} | {n:,} | missing in {round(frac * coverage['shape'][0])} of "
    f"{coverage['shape'][0]} samples |"
    for frac, n in coverage['na_bands']
)}

`max_na_frac` therefore only matters where it crosses a band. The default **0.20 falls
below every one of them**, so at defaults `knn`, `softimpute` and `missforest` keep
exactly the {coverage['genes_strict']:,} zero-NA genes that `strict` keeps, have nothing
left to fill, and return an identical matrix. At 0.4 the gate admits the
{coverage['first_gap']:.4f} band and they keep {coverage['genes_at_04']:,} genes instead.
`examples/configs/full_crossproduct.yaml` sets it accordingly.

## Annotation columns

The source table has 483 columns; this ships the {ann.shape[1]} that are populated for
every selected sample. The rest are FL-specific clinical fields (`coo`, `mfp`, `dhitsig`,
treatment, survival) present for one cohort only - shipping them would suggest the tool
needs a study design it does not.

## Regenerating

```bash
python3.11 scripts/build_example_dataset.py --dry-run   # selection only, no S3 access
python3.11 scripts/build_example_dataset.py
```

Requires read access to the source bucket. The output is committed, so this is a
maintainer operation only.
""")


def main(argv: list[str] | None = None) -> int:
    """Build the example dataset. Returns a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ann", default=DEFAULT_ANN, help="source annotation CSV")
    parser.add_argument("--bucket", default=DEFAULT_BUCKET)
    parser.add_argument("--exp-key", default=DEFAULT_EXP_KEY)
    parser.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR), type=Path)
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="report the selection and stop, without touching S3.",
    )
    args = parser.parse_args(argv)

    log(f"reading {args.ann}")
    ann = pd.read_csv(args.ann, index_col=0, low_memory=False)
    kept, dropped, report = select_samples(ann)

    log(f"selected {report['n_source_rows']} rows -> {report['n_kept']} unique samples")
    log(f"  dropped {len(dropped)} duplicate rows: {report['dropped_ids'][:3]} ...")
    for base, labels in sorted(report["diagnosis_disagreements"].items()):
        log(f"  diagnosis conflict on {base}: {labels}")
    for label in COHORTS:
        log(f"  {label}: {(kept['COHORT_LABEL'] == label).sum()}")

    reduced = reduce_columns(kept)
    log(f"annotation columns: {kept.shape[1]} -> {reduced.shape[1]} fully populated")

    if args.dry_run:
        log("dry run: nothing written")
        return 0

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    ann_path = out_dir / "example_ann.csv"
    exp_path = out_dir / "example_exp.tsv.gz"

    reduced.to_csv(ann_path)
    log(f"wrote {ann_path} ({ann_path.stat().st_size / 1e3:.0f} kB)")

    exp_stats = fetch_expression(
        list(reduced.index), exp_path, bucket=args.bucket, key=args.exp_key
    )
    log(
        f"wrote {exp_path} ({exp_stats['n_samples']} x {exp_stats['n_genes']}, "
        f"{exp_stats['gzipped_bytes'] / 1e6:.2f} MB gzipped)"
    )

    coverage = describe_coverage(exp_path)
    log(
        f"coverage: {coverage['all_na_genes']} all-NA genes, "
        f"{coverage['complete_genes']} complete, "
        f"{coverage['overall_na_fraction']:.1%} NA overall"
    )

    write_provenance(out_dir / "PROVENANCE.md", reduced, report, exp_stats, coverage)
    log(f"wrote {out_dir / 'PROVENANCE.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
