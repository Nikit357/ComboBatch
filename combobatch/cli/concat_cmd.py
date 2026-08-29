"""``combobatch concat`` — sidecars into analysis tables.

Writes ``metrics_summary.csv`` beside the outputs, plus the three long-format tables when
groups L, M or N contributed nested detail. ``--out-dir`` with ``--date-tag`` also writes
dated local copies, so a rerun never overwrites an earlier snapshot.
"""

from __future__ import annotations

import argparse

from combobatch import storage
from combobatch.cli import run_cmd
from combobatch.logging_utils import get_logger, setup_logging
from combobatch.params import KeyParseError


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the ``concat`` flags."""
    parser.add_argument(
        "--out",
        metavar="URI",
        required=True,
        help="The run's output root, local or s3://. Sidecars are read from its "
        "metrics/ prefix and the tables are written back to it.",
    )
    parser.add_argument(
        "--out-dir",
        metavar="PATH",
        help="Also write dated local copies of every table here.",
    )
    parser.add_argument(
        "--date-tag",
        metavar="YYMMDD",
        help="Suffix for the dated local copies. Defaults to today.",
    )
    parser.add_argument(
        "--endpoint-url", metavar="URL", help="Alternative S3 endpoint."
    )
    parser.add_argument(
        "--allow-unparseable-keys",
        action="store_true",
        help="Warn and continue on a sidecar whose filename does not parse, instead of "
        "failing. Off by default: the donor dropped those silently, and a dropped key "
        "is a result that vanishes from every table without a word.",
    )
    parser.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )


def main(argv: list[str] | None = None) -> int:
    """Aggregate the sidecars. Returns the process exit code."""
    parser = argparse.ArgumentParser(prog="combobatch concat")
    add_arguments(parser)
    args = parser.parse_args(argv)

    setup_logging(args.log_level)
    log = get_logger()

    from combobatch.metrics.concat import concat_metrics, write_tables

    try:
        result = concat_metrics(
            args.out,
            endpoint_url=args.endpoint_url,
            strict_keys=not args.allow_unparseable_keys,
        )
    except KeyParseError as exc:
        log.error("%s", exc)
        return run_cmd.EXIT_FAILED

    if result.summary.empty:
        log.error("no metrics sidecars found under %s/metrics/", args.out)
        return run_cmd.EXIT_MISSING_INPUT

    backend = storage.backend_for_root(args.out, endpoint_url=args.endpoint_url)
    written = write_tables(
        result, out_dir=args.out_dir, date_tag=args.date_tag, backend=backend
    )

    log.info(
        "Summarized %d output(s) into %d column(s).",
        len(result.summary),
        len(result.summary.columns),
    )
    for uri in written:
        log.info("wrote %s", uri)
    return run_cmd.EXIT_OK


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
