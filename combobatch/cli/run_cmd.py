"""``combobatch run`` — one imputation x method combination.

This is both a user-facing command and the worker subprocess ``dispatch`` launches. A
fresh process per job is not incidental: it is what gives every rpy2 job its own R
session, which neither threads-around-R nor a forked pool can provide.

It also owns the argument surface that ``dispatch`` reuses, so the two commands cannot
disagree about what a configuration means.
"""

from __future__ import annotations

import argparse
import sys
from typing import Any

from combobatch import storage
from combobatch.config import (
    ColumnSpec,
    ConfigError,
    ImputationSelection,
    MethodSelection,
    MetricsSpec,
    PostRemovalSpec,
    RunConfig,
    parse_metric_groups,
)
from combobatch.logging_utils import get_logger, pin_threads, setup_logging
from combobatch.params import ParamParseError, coerce_with_spec, parse_params_arg

# Exit codes are a contract the dispatcher reads; see cli/CLAUDE.md.
EXIT_OK = 0
EXIT_FAILED = 1
EXIT_INSUFFICIENT_MEMORY = 2
EXIT_MISSING_INPUT = 3


def add_common_arguments(parser: argparse.ArgumentParser) -> None:
    """Register every flag shared by ``run`` and ``dispatch``."""
    source = parser.add_argument_group("input and output")
    source.add_argument("--config", metavar="PATH", help="YAML configuration file.")
    source.add_argument(
        "--exp", metavar="URI", help="Expression matrix, samples x genes."
    )
    source.add_argument(
        "--ann", metavar="URI", help="Annotation table, samples x fields."
    )
    source.add_argument(
        "--out", metavar="URI", help="Output root: a directory or s3://."
    )
    source.add_argument(
        "--endpoint-url", metavar="URL", help="Alternative S3 endpoint."
    )

    columns = parser.add_argument_group("columns")
    columns.add_argument("--batch-col", metavar="NAME", help="Batch identity column.")
    columns.add_argument(
        "--bio-col", metavar="NAME", help="Biology column to preserve."
    )
    columns.add_argument("--cohort-col", metavar="NAME", help="Optional cohort column.")
    columns.add_argument(
        "--reference-batch",
        metavar="LABEL",
        help="Batch that reference-based methods normalize toward. Default: the largest.",
    )
    columns.add_argument(
        "--subset-query",
        metavar="EXPR",
        help="pandas query over annotation columns, applied before anything else.",
    )

    tuning = parser.add_argument_group("hyperparameters")
    tuning.add_argument(
        "--method-params",
        metavar="SPEC",
        help="'10_mnn:k=50;38_harman:limit=0.05' — validated against each method.",
    )
    tuning.add_argument(
        "--impute-params",
        metavar="SPEC",
        help="'knn:knn_k=10;softimpute:rank_max=30' — validated against each imputer.",
    )
    tuning.add_argument(
        "--max-na-frac",
        type=float,
        metavar="FRAC",
        help="Maximum per-gene NA fraction eligible for imputation. Default: the "
        "imputer's own (0.0 for strict, 0.20 otherwise).",
    )

    post = parser.add_argument_group("post-removal")
    post.add_argument(
        "--post-removal",
        action="store_true",
        default=None,
        help="Also write a variant with the most PCA-deviant batch(es) dropped.",
    )
    post.add_argument(
        "--post-removal-n", type=int, metavar="N", help="Batches to drop."
    )
    post.add_argument("--post-removal-min-batch-size", type=int, metavar="N")
    post.add_argument("--post-removal-n-pcs", type=int, metavar="N")

    metrics = parser.add_argument_group("metrics")
    metrics.add_argument(
        "--metrics",
        action="store_true",
        default=None,
        help="Compute quality metrics inline, alongside harmonization.",
    )
    metrics.add_argument(
        "--groups", metavar="LETTERS", help="Metric groups, e.g. A,B,E."
    )
    metrics.add_argument("--marker-panel", metavar="PATH", help="Marker panel for L/M.")

    behaviour = parser.add_argument_group("behaviour")
    behaviour.add_argument(
        "--skip-if-exists",
        action="store_true",
        help="Report an existing output as cached instead of recomputing it.",
    )
    behaviour.add_argument("--random-seed", type=int, metavar="N")
    behaviour.add_argument("--run-id", metavar="ID", help="Names the failed-jobs log.")
    behaviour.add_argument(
        "--threads",
        type=int,
        default=1,
        metavar="N",
        help="BLAS/OpenMP threads per process. Default 1: the dispatcher supplies the "
        "parallelism, and unpinned threads oversubscribe the node.",
    )
    behaviour.add_argument(
        "--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"]
    )


def add_arguments(parser: argparse.ArgumentParser) -> None:
    """Register the ``run``-only flags on top of the shared ones."""
    add_common_arguments(parser)
    single = parser.add_argument_group("this combination")
    single.add_argument(
        "--imputation", metavar="NAME", help="One IMPUTER_REGISTRY key."
    )
    single.add_argument("--method", metavar="KEY", help="One METHOD_REGISTRY key.")
    single.add_argument(
        "--param-tag",
        metavar="TAG",
        help="Override the encoded parameter tag in output names.",
    )
    single.add_argument(
        "--out-json",
        metavar="PATH",
        help="Write the sidecar rows here as JSON. The dispatcher reads this.",
    )
    single.add_argument(
        "--memory-limit-gb",
        type=float,
        metavar="GB",
        help="Exit with code 2 rather than starting if less than this is free.",
    )
    single.add_argument(
        "--no-manifest",
        action="store_true",
        help="Do not write run_manifest.json. The dispatcher writes it once, up front, "
        "so concurrent workers do not race on the same file.",
    )


def build_config(args: argparse.Namespace) -> RunConfig:
    """
    Assemble a :class:`RunConfig` from a YAML file, CLI flags, or both.

    A ``--config`` file supplies the baseline and every flag overrides it, so a config
    can be reused across runs that differ in one axis.

    Returns
    -------
    An unvalidated configuration; the caller validates it before spending anything.

    Raises
    ------
    ConfigError
        If neither a config file nor the required flags are present.
    """
    cfg = RunConfig.from_yaml(args.config) if args.config else None

    imputations = _imputation_selections(args, cfg)
    methods = _method_selections(args, cfg)
    columns = _columns(args, cfg)

    exp_uri = args.exp or (cfg.exp_uri if cfg else None)
    ann_uri = args.ann or (cfg.ann_uri if cfg else None)
    out_uri = args.out or (cfg.out_uri if cfg else None)
    missing = [
        name
        for name, value in (("--exp", exp_uri), ("--ann", ann_uri), ("--out", out_uri))
        if not value
    ]
    if missing:
        raise ConfigError(
            f"missing required input(s): {', '.join(missing)}. Give them as flags or in "
            f"a --config file."
        )

    base_metrics = cfg.metrics if cfg else MetricsSpec()
    metrics = MetricsSpec(
        enabled=base_metrics.enabled if args.metrics is None else args.metrics,
        groups=parse_metric_groups(args.groups, base_metrics.groups),
        skip_slow=base_metrics.skip_slow,
        skip_wm=base_metrics.skip_wm,
        n_perm=base_metrics.n_perm,
        marker_panel=args.marker_panel or base_metrics.marker_panel,
        predict_classes=base_metrics.predict_classes,
        force_groups=base_metrics.force_groups,
    )

    base_post = cfg.post_removal if cfg else PostRemovalSpec()
    post_removal = PostRemovalSpec(
        enabled=base_post.enabled if args.post_removal is None else args.post_removal,
        n_batches=args.post_removal_n or base_post.n_batches,
        min_batch_size=args.post_removal_min_batch_size or base_post.min_batch_size,
        n_pcs=args.post_removal_n_pcs or base_post.n_pcs,
    )

    return RunConfig(
        exp_uri=exp_uri,
        ann_uri=ann_uri,
        out_uri=out_uri,
        columns=columns,
        imputations=imputations,
        methods=methods,
        metrics=metrics,
        post_removal=post_removal,
        reference_batch=args.reference_batch or (cfg.reference_batch if cfg else None),
        subset_query=args.subset_query or (cfg.subset_query if cfg else None),
        random_seed=args.random_seed or (cfg.random_seed if cfg else 42),
        endpoint_url=args.endpoint_url or (cfg.endpoint_url if cfg else None),
        work_dir=cfg.work_dir if cfg else None,
        run_id=args.run_id or (cfg.run_id if cfg else "local"),
    )


def _columns(args: argparse.Namespace, cfg: RunConfig | None) -> ColumnSpec:
    """Resolve the column specification from flags over the config file."""
    batch = args.batch_col or (cfg.columns.batch if cfg else None)
    bio = args.bio_col or (cfg.columns.bio if cfg else None)
    if not batch or not bio:
        raise ConfigError(
            "both --batch-col and --bio-col are required (or a columns: block in "
            "--config)"
        )
    if cfg is not None:
        return ColumnSpec(
            batch=batch,
            bio=bio,
            cohort=args.cohort_col or cfg.columns.cohort,
            metric_batch_cols=cfg.columns.metric_batch_cols,
            metric_bio_cols=cfg.columns.metric_bio_cols,
            predict_class_col=cfg.columns.predict_class_col,
        )
    return ColumnSpec(batch=batch, bio=bio, cohort=args.cohort_col)


def _imputation_selections(
    args: argparse.Namespace, cfg: RunConfig | None
) -> tuple[ImputationSelection, ...]:
    """Resolve which imputers to run, with their parameters."""
    from combobatch.imputation import IMPUTER_REGISTRY

    requested = _requested_names(args, "imputation", "imputations")
    if requested is None:
        selections = list(cfg.imputations) if cfg else []
    else:
        selections = [ImputationSelection(name=name) for name in requested]
    if not selections:
        raise ConfigError("no imputation selected: pass --imputation or --imputations")

    overrides = _parsed_params(args.impute_params, "--impute-params")
    return tuple(
        _apply_overrides(
            selection,
            overrides.get(selection.name, {}),
            IMPUTER_REGISTRY,
            max_na_frac=getattr(args, "max_na_frac", None),
        )
        for selection in selections
    )


def _method_selections(
    args: argparse.Namespace, cfg: RunConfig | None
) -> tuple[MethodSelection, ...]:
    """Resolve which methods to run, with their parameters."""
    from combobatch.methods import METHOD_REGISTRY

    requested = _requested_names(args, "method", "methods")
    if requested is None:
        selections = list(cfg.methods) if cfg else []
    else:
        selections = [MethodSelection(name=name) for name in requested]
    if not selections:
        raise ConfigError("no method selected: pass --method or --methods")

    overrides = _parsed_params(args.method_params, "--method-params")
    param_tag = getattr(args, "param_tag", None)
    return tuple(
        _apply_overrides(
            selection, overrides.get(selection.name, {}), METHOD_REGISTRY, tag=param_tag
        )
        for selection in selections
    )


def _requested_names(
    args: argparse.Namespace, single_flag: str, list_flag: str
) -> list[str] | None:
    """
    Return the names requested on the command line, or ``None`` to use the config file.

    ``run`` takes one name; ``dispatch`` takes a comma-separated list that also accepts
    ``all``, expanded by the dispatcher before this point.
    """
    listed = getattr(args, list_flag, None)
    if listed:
        return [name.strip() for name in listed.split(",") if name.strip()]
    single = getattr(args, single_flag, None)
    return [single] if single else None


def _parsed_params(raw: str | None, flag: str) -> dict[str, dict[str, Any]]:
    """Parse a ``--*-params`` flag, re-raising with the flag named."""
    if not raw:
        return {}
    try:
        return parse_params_arg(raw)
    except ParamParseError as exc:
        raise ConfigError(f"{flag}: {exc}") from exc


def _apply_overrides(
    selection: Any,
    overrides: dict[str, Any],
    registry: dict[str, Any],
    *,
    tag: str | None = None,
    max_na_frac: float | None = None,
) -> Any:
    """
    Merge command-line parameters onto a selection, coercing to each declared type.

    Inference alone is not enough: ``target_group=123`` is a string parameter whose value
    looks like an integer, and only the declaration knows that.
    """
    import dataclasses

    spec = registry.get(selection.name)
    merged = dict(selection.params)
    for name, value in overrides.items():
        declared = spec.hyperparams.get(name) if spec else None
        if declared is not None and not isinstance(value, (list, tuple, type(None))):
            value = coerce_with_spec(name, str(value), declared.type)
        merged[name] = value

    # Names are checked here, at parse time, rather than inside the job: a typo must not
    # survive as far as a parameter tag baked into an output filename.
    if spec is not None:
        try:
            spec.validate_params(merged)
        except (ParamParseError, ValueError) as exc:
            raise ConfigError(str(exc)) from exc

    changes: dict[str, Any] = {"params": merged}
    if tag is not None:
        changes["param_tag"] = tag
    if max_na_frac is not None and hasattr(selection, "max_na_frac"):
        changes["max_na_frac"] = max_na_frac
    return dataclasses.replace(selection, **changes)


def main(argv: list[str] | None = None) -> int:
    """Run one combination. Returns the process exit code."""
    parser = argparse.ArgumentParser(prog="combobatch run")
    add_arguments(parser)
    args = parser.parse_args(argv)

    pin_threads(args.threads)
    setup_logging(args.log_level)
    log = get_logger()

    from combobatch.imputation import IMPUTER_REGISTRY
    from combobatch.methods import METHOD_REGISTRY

    try:
        cfg = build_config(args)
        # With the registries, so an unknown key is a named error with a suggestion here
        # rather than a bare KeyError after the input has already been downloaded.
        cfg.validate(known_methods=METHOD_REGISTRY, known_imputers=IMPUTER_REGISTRY)
    except ConfigError as exc:
        log.error("%s", exc)
        return EXIT_MISSING_INPUT

    if args.memory_limit_gb is not None:
        from combobatch.memory import check_memory, free_gb

        if not check_memory(args.memory_limit_gb):
            log.error(
                "only %.1f GB free, below the %.1f GB required — not starting",
                free_gb(),
                args.memory_limit_gb,
            )
            return EXIT_INSUFFICIENT_MEMORY

    # Imported here so that --help and a bad configuration cost nothing.
    from combobatch import dataio
    from combobatch.harmonize import run_one_combination, write_manifest

    imp_sel, method_sel = cfg.combinations()[0]
    try:
        rows = run_one_combination(
            cfg, imp_sel, method_sel, skip_if_exists=args.skip_if_exists
        )
    except Exception as exc:
        log.error("%s x %s: %s", imp_sel.name, method_sel.name, exc, exc_info=True)
        return EXIT_FAILED

    if not args.no_manifest:
        write_manifest(cfg)

    if args.out_json:
        backend, key = storage.resolve(args.out_json, endpoint_url=cfg.endpoint_url)
        dataio.write_json(rows, backend, key)

    statuses = {row["status"] for row in rows}
    for row in rows:
        log.info("%s -> %s", row["key"], row["status"])
    return EXIT_OK if statuses <= {"ok", "cached", "skipped"} else EXIT_FAILED


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
