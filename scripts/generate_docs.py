#!/usr/bin/env python3
"""
Regenerate the registry-derived documentation.

    python3.11 scripts/generate_docs.py            # write docs/
    python3.11 scripts/generate_docs.py --check    # exit 1 if anything is stale

`--check` is what CI runs. `tests/unit/test_docs_in_sync.py` asserts the same thing, so a
registry change that is not reflected in the docs fails the unit suite too.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_DIR = REPO_ROOT / "docs"


def main(argv: list[str] | None = None) -> int:
    """Write or check the generated docs. Returns a process exit code."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--check",
        action="store_true",
        help="Do not write; exit 1 if any file differs from the registries.",
    )
    parser.add_argument("--docs-dir", type=Path, default=DOCS_DIR)
    args = parser.parse_args(argv)

    from combobatch.docsgen import GENERATED_DOCS

    stale = []
    for filename, render in GENERATED_DOCS.items():
        path = args.docs_dir / filename
        rendered = render()
        if args.check:
            current = path.read_text() if path.exists() else ""
            if current != rendered:
                stale.append(filename)
                print(f"STALE: docs/{filename}", file=sys.stderr)
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(rendered)
        print(f"wrote docs/{filename} ({len(rendered.splitlines())} lines)")

    if args.check:
        if stale:
            print(
                f"\n{len(stale)} file(s) out of date. Run:\n"
                f"    python3.11 scripts/generate_docs.py",
                file=sys.stderr,
            )
            return 1
        print("docs are in sync with the registries")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
