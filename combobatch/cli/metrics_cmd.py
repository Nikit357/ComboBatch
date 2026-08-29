"""``combobatch metrics`` — compute metrics over outputs that already exist.

This is the "later" half of the requirement: metrics can run inline with harmonization,
or here, against a finished run. It enumerates ``exp/`` under the output root, reads the
matching sidecar, and computes only the groups whose sentinel key is absent — so a group
added months later is filled in without recomputing the thirteen that already ran.

An unparseable output key is a loud error. The donor's enumerator silently dropped any
key that did not split into exactly four parts, which meant a job that quietly never got
metrics and never appeared in any summary.
"""

from __future__ import annotations

import argparse
import sys

from combobatch import dataio, storage
from combobatch.cli import run_cmd
from combobatch.config import ConfigError, RunConfig
from combobatch.harmonize import MATRIX_SUFFIX
from combobatch.logging_utils import get_logger, pin_threads, setup_logging
from combobatch.params import KeyParseError, parse_output_key


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the ``metrics`` flags on top of the shared ones."""
    run_cmd.add_common_arguments(parser)
    group = parser.add_argument_group("metrics over existing outputs")
    group.add_argument(
        "--from-outputs",
        metavar="URI",
        help="Output root to scan. Defaults to the configured --out.",
    )
    group.add_argument(
        "--force-groups",
        metavar="LETTERS",
        help="Recompute these groups even if their sentinel is already present, and "
        "even if --skip-slow would otherwise exclude them.",
    )
    group.add_argument(
        "--skip-slow",
        action="store_true",
        default=None,
        help="Exclude the groups the registry marks very slow.",
    )
    group.add_argument(
        "--collect-detail",
        action="store_true",
        help="Ask group L for its full gene-by-cohort matrix (large).",
    )
    group.add_argument(
        "--dry-run",
        action="store_true",
        help="List the outputs found and the groups each still needs, then stop.",
    )


def find_outputs(backend: storage.StorageBackend) -> list[str]:
    """
    Return every output key under ``exp/``, parsed and sorted.

    Raises
    ------
    KeyParseError
        On a key that does not match the output grammar.
    """
    stems: list[str] = []
    for storage_key in sorted(backend.list_keys("exp")):
        if not storage_key.endswith(MATRIX_SUFFIX):
            continue
        stem = storage_key.split("/")[-1][: -len(MATRIX_SUFFIX)]
        try:
            parse_output_key(stem)
        except KeyParseError as exc:
            raise KeyParseError(
                f"exp/{stem}{MATRIX_SUFFIX}: {exc}. ComboBatch will not silently skip "
                f"an output it cannot name."
            ) from exc
        stems.append(stem)
    return stems


def pending_groups(prior: dict, cfg: RunConfig, columns) -> list[str]:
    """Return the requested groups whose sentinel is not yet in ``prior``."""
    from combobatch.metrics.runner import resolve_active_groups

    return sorted(
        resolve_active_groups(
            cfg.metrics.groups,
            skip_slow=cfg.metrics.skip_slow,
            skip_wm=cfg.metrics.skip_wm,
            force=cfg.metrics.force_groups,
            prior=prior,
            columns=columns,
        )
    )


def main(argv: list[str] | None = None) -> int:
    """Compute metrics over an existing run. Returns the process exit code."""
    parser = argparse.ArgumentParser(prog="combobatch metrics")
    add_arguments(parser)
    args = parser.parse_args(argv)

    pin_threads(args.threads)
    setup_logging(args.log_level)
    log = get_logger()

    from combobatch.imputation import IMPUTER_REGISTRY
    from combobatch.methods import METHOD_REGISTRY
    from combobatch.metrics import METRIC_GROUP_REGISTRY, MetricColumns

    # `metrics` runs no method and no imputer: it reads matrices that already exist. The
    # shared config builder still requires both axes, so placeholders are supplied here
    # and never used. They are set only when a config file did not provide them.
    args.methods = None
    args.imputations = None
    args.param_tag = None
    args.method = None if args.config else next(iter(METHOD_REGISTRY))
    args.imputation = None if args.config else next(iter(IMPUTER_REGISTRY))

    try:
        cfg = run_cmd.build_config(args)
        if args.skip_slow is not None or args.force_groups:
            from dataclasses import replace

            from combobatch.config import parse_metric_groups

            cfg = replace(
                cfg,
                metrics=replace(
                    cfg.metrics,
                    skip_slow=(
                        args.skip_slow
                        if args.skip_slow is not None
                        else cfg.metrics.skip_slow
                    ),
                    force_groups=parse_metric_groups(
                        args.force_groups, cfg.metrics.force_groups
                    ),
                ),
            )
        cfg.validate(
            known_methods=METHOD_REGISTRY,
            known_imputers=IMPUTER_REGISTRY,
            known_groups=METRIC_GROUP_REGISTRY,
        )
    except ConfigError as exc:
        log.error("%s", exc)
        return run_cmd.EXIT_MISSING_INPUT

    root = args.from_outputs or cfg.out_uri
    backend = storage.backend_for_root(root, endpoint_url=cfg.endpoint_url)
    columns = MetricColumns.from_spec(cfg.columns)

    try:
        stems = find_outputs(backend)
    except KeyParseError as exc:
        log.error("%s", exc)
        return run_cmd.EXIT_FAILED
    if not stems:
        log.error("no outputs found under %s/exp/", root)
        return run_cmd.EXIT_MISSING_INPUT

    if args.dry_run:
        print(f"Output root: {root}")
        print(f"Outputs: {len(stems)}")
        for stem in stems:
            sidecar = f"metrics/{stem}_metrics.json"
            prior = (
                dataio.read_json(backend, sidecar) if backend.exists(sidecar) else {}
            )
            todo = pending_groups(prior, cfg, columns)
            print(f"  {stem}: {''.join(todo) or 'nothing to do'}")
        return run_cmd.EXIT_OK

    from combobatch.metrics.panels import load_panel
    from combobatch.metrics.runner import compute_metrics

    panel = (
        load_panel(cfg.metrics.marker_panel, endpoint_url=cfg.endpoint_url)
        if cfg.metrics.marker_panel
        else None
    )
    ann_df = dataio.read_table(cfg.ann_uri, endpoint_url=cfg.endpoint_url)
    n_failed = 0

    for stem in stems:
        sidecar = f"metrics/{stem}_metrics.json"
        prior = dataio.read_json(backend, sidecar) if backend.exists(sidecar) else {}
        todo = pending_groups(prior, cfg, columns)
        if not todo:
            log.info("%s: every requested group is already present", stem)
            continue

        log.info("%s: computing group(s) %s", stem, "".join(todo))
        exp_df = dataio.read_table_from(backend, f"exp/{stem}{MATRIX_SUFFIX}")
        aligned = ann_df.reindex(exp_df.index)

        reference = None
        if "L" in todo:
            key = parse_output_key(stem)
            prepared = f"prepared/{key.imputation}__exp{MATRIX_SUFFIX}"
            if backend.exists(prepared):
                reference = dataio.read_table_from(backend, prepared).reindex(
                    index=exp_df.index
                )
            else:
                log.warning(
                    "%s: group L needs %s, which is absent; it will skip",
                    stem,
                    prepared,
                )

        def flush(partial: dict, sidecar=sidecar, prior=prior) -> None:
            """Persist after every group, so a later crash cannot lose an earlier one."""
            dataio.write_json({**prior, **partial}, backend, sidecar)

        try:
            computed = compute_metrics(
                exp_df,
                aligned,
                cfg.columns,
                spec=cfg.metrics,
                groups=frozenset(todo),
                ref_df=reference,
                panel=panel,
                collect_detail=args.collect_detail,
                on_group_done=flush,
            )
        except Exception as exc:
            log.error("%s: metrics failed: %s", stem, exc, exc_info=True)
            n_failed += 1
            continue

        dataio.write_json({**prior, **computed}, backend, sidecar)

    log.info("Done: %d output(s), %d failed.", len(stems), n_failed)
    return run_cmd.EXIT_OK if not n_failed else run_cmd.EXIT_FAILED


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
