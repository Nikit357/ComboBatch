"""Root ``combobatch`` command-line entry point.

Subcommands are dispatched to their own modules. The three ``list-*`` commands read the
registries directly and are the reason no CLI choice list is ever hand-maintained: a
method added to ``METHOD_REGISTRY`` appears here, in ``--methods`` validation and in the
generated docs without another line being written anywhere.

Commands whose implementation phase has not landed exit with code 2 naming that phase,
so the surface is inspectable without pretending to work.
"""

from __future__ import annotations

import argparse
import sys
from typing import Any

from combobatch import __version__

# Subcommands still awaiting their implementation phase, with the phase named.
PENDING: dict[str, tuple[str, str]] = {}

# Subcommands that own their own argument parser and are handed the rest of argv verbatim.
DELEGATED = ("run", "dispatch", "metrics", "concat", "selftest")


def build_parser() -> argparse.ArgumentParser:
    """Construct the root parser and register every subcommand."""
    from combobatch.cli import (
        concat_cmd,
        dispatch_cmd,
        metrics_cmd,
        run_cmd,
        selftest_cmd,
    )

    parser = argparse.ArgumentParser(
        prog="combobatch",
        description=(
            "Benchmark imputation x harmonization combinations for bulk "
            "transcriptomics, with batch-effect quality metrics."
        ),
    )
    parser.add_argument(
        "--version", action="version", version=f"combobatch {__version__}"
    )
    subparsers = parser.add_subparsers(dest="command", metavar="<command>")

    run_cmd.add_arguments(
        subparsers.add_parser("run", help="Run one (imputation, method) combination")
    )
    dispatch_cmd.add_arguments(
        subparsers.add_parser(
            "dispatch", help="Run the imputation x method cross-product in parallel"
        )
    )
    metrics_cmd.add_arguments(
        subparsers.add_parser(
            "metrics", help="Compute quality metrics over existing outputs"
        )
    )
    concat_cmd.add_arguments(
        subparsers.add_parser(
            "concat", help="Aggregate metric sidecars into summary tables"
        )
    )
    selftest_cmd.add_arguments(
        subparsers.add_parser(
            "selftest", help="Run every method and imputer on synthetic data"
        )
    )

    listing = subparsers.add_parser("list-metrics", help="List metric groups A-N")
    listing.add_argument("--json", action="store_true", help="Emit JSON.")
    listing.add_argument(
        "--available-only",
        action="store_true",
        help="List only groups whose backends are present here.",
    )

    for name, help_text in (
        ("list-methods", "List harmonization methods and their hyperparameters"),
        ("list-imputers", "List imputation strategies and their hyperparameters"),
    ):
        listing = subparsers.add_parser(name, help=help_text)
        listing.add_argument(
            "--json", action="store_true", help="Emit JSON instead of a table."
        )
        listing.add_argument(
            "--available-only",
            action="store_true",
            help="List only entries whose backends are present here.",
        )
        if name == "list-methods":
            listing.add_argument(
                "--show-excluded",
                action="store_true",
                help="Also list the five method keys that are deliberately absent, and "
                "why.",
            )

    for name, (help_text, _phase) in PENDING.items():
        subparsers.add_parser(name, help=help_text, add_help=False)

    return parser


def _registry_rows(registry: dict[str, Any], *, available_only: bool) -> list[dict]:
    """Flatten a registry into printable rows, sharing one shape across both kinds."""
    rows = []
    for key, spec in registry.items():
        missing = spec.missing_backends()
        if available_only and missing:
            continue
        rows.append(
            {
                "key": key,
                "harshness": getattr(spec, "harshness", ""),
                "backends": ", ".join(missing) if missing else "ready",
                "hyperparams": {
                    name: {
                        "type": getattr(param.type, "__name__", str(param.type)),
                        "default": param.default,
                        "help": param.help,
                        "choices": list(param.choices) if param.choices else None,
                    }
                    for name, param in spec.hyperparams.items()
                },
                "citation": getattr(spec, "citation", "")
                or getattr(spec, "description", ""),
            }
        )
    return rows


def _print_rows(rows: list[dict], *, title: str) -> None:
    """Print a registry listing as an aligned table with its parameters indented."""
    print(f"{title} ({len(rows)})")
    print()
    width = max((len(row["key"]) for row in rows), default=0)
    for row in rows:
        harshness = f"  [{row['harshness']}]" if row["harshness"] else ""
        status = "" if row["backends"] == "ready" else f"  (needs {row['backends']})"
        print(f"  {row['key']:<{width}}{harshness}{status}")
        if row["citation"]:
            print(f"  {'':<{width}}  {row['citation']}")
        for name, param in row["hyperparams"].items():
            choices = f" one of {param['choices']}" if param["choices"] else ""
            print(
                f"  {'':<{width}}    --{name} ({param['type']}, "
                f"default {param['default']!r}){choices}"
            )
        print()


def _list_metrics(args: argparse.Namespace) -> int:
    """Handle ``list-metrics``, reading the group registry."""
    import json

    from combobatch.metrics import DEFAULT_GROUPS, METRIC_GROUP_REGISTRY, GROUP_LETTERS

    rows = []
    for letter in GROUP_LETTERS:
        spec = METRIC_GROUP_REGISTRY[letter]
        missing = spec.missing_backends()
        if args.available_only and missing:
            continue
        rows.append(
            {
                "letter": letter,
                "name": spec.name,
                "speed": spec.speed,
                "default": letter in DEFAULT_GROUPS,
                "sentinel": spec.sentinel_template,
                "column_roles": list(spec.column_roles),
                "needs_reference": spec.needs_reference,
                "needs_panel": spec.needs_panel,
                "needs_embedding": list(spec.needs_embedding),
                "backends": ", ".join(missing) if missing else "ready",
                "description": spec.description,
            }
        )

    if args.json:
        print(json.dumps({"groups": rows}, indent=2))
        return 0

    print(f"Metric groups ({len(rows)})")
    print()
    for row in rows:
        flags = [row["speed"]]
        if row["default"]:
            flags.append("default")
        if row["needs_reference"]:
            flags.append("needs reference")
        if row["needs_panel"]:
            flags.append("uses panel")
        if row["backends"] != "ready":
            flags.append(f"needs {row['backends']}")
        print(f"  {row['letter']}  {row['name']}  [{', '.join(flags)}]")
        print(f"     {row['description']}")
        print(f"     sentinel: {row['sentinel']}")
        print()
    return 0


def _list_command(args: argparse.Namespace) -> int:
    """Handle ``list-methods`` and ``list-imputers``."""
    import json

    if args.command == "list-metrics":
        return _list_metrics(args)

    if args.command == "list-methods":
        from combobatch.methods import EXCLUDED_METHODS, METHOD_REGISTRY

        registry, title = METHOD_REGISTRY, "Harmonization methods"
    else:
        from combobatch.imputation import IMPUTER_REGISTRY

        registry, title = IMPUTER_REGISTRY, "Imputation strategies"

    rows = _registry_rows(registry, available_only=args.available_only)

    if args.json:
        payload: dict[str, Any] = {"entries": rows}
        if args.command == "list-methods" and args.show_excluded:
            payload["excluded"] = EXCLUDED_METHODS
        print(json.dumps(payload, indent=2, default=str))
        return 0

    _print_rows(rows, title=title)
    if args.command == "list-methods" and args.show_excluded:
        print(f"Deliberately absent ({len(EXCLUDED_METHODS)}):")
        print()
        for key, reason in EXCLUDED_METHODS.items():
            print(f"  {key}: {reason}")
        print()
    return 0


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and dispatch. Returns the process exit code."""
    argv = list(sys.argv[1:] if argv is None else argv)

    # The subcommand modules own their own flags, so the root parser only needs to see
    # the command name; everything after it is handed over verbatim.
    if argv and argv[0] in DELEGATED:
        from combobatch.cli import (
            concat_cmd,
            dispatch_cmd,
            metrics_cmd,
            run_cmd,
            selftest_cmd,
        )

        modules = {
            "run": run_cmd,
            "dispatch": dispatch_cmd,
            "metrics": metrics_cmd,
            "concat": concat_cmd,
            "selftest": selftest_cmd,
        }
        return modules[argv[0]].main(argv[1:])

    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        parser.print_help()
        return 0
    if args.command in PENDING:
        _help_text, phase = PENDING[args.command]
        print(
            f"combobatch {args.command}: not implemented yet (scheduled for {phase}).",
            file=sys.stderr,
        )
        return 2
    return _list_command(args)


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
