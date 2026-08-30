"""`requirements.txt` and `pyproject.toml` must not drift apart.

`requirements.txt` is the image's environment; `pyproject.toml` is what a laptop install
resolves. They state the same constraints twice, and until 2026-08-30 the claim that a
test enforced their agreement was written in `requirements.txt` and true of nothing — no
test read `pyproject.toml` at all. A drift here means the image and the laptop run
different versions of the same pinned library, which is precisely the class of difference
that makes a result irreproducible without ever failing anything.
"""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
REQUIREMENTS = REPO_ROOT / "requirements.txt"
PYPROJECT = REPO_ROOT / "pyproject.toml"

# Distributions deliberately in the image only, each for a stated reason.
EXPECTED_ONLY_IN_REQUIREMENTS = frozenset(
    {
        # A command-line tool for working inside the pod, not a library this package
        # imports. It has no business in an extra a laptop install would resolve.
        "awscli",
    }
)


def _normalize(name: str) -> str:
    """PEP 503 name normalization: `scikit-learn` and `scikit_learn` are one package."""
    return re.sub(r"[-_.]+", "-", name).lower()


def requirement_constraints() -> dict[str, str]:
    """Every pinned distribution in `requirements.txt`, name → version specifier."""
    found: dict[str, str] = {}
    for line in REQUIREMENTS.read_text().splitlines():
        line = line.split("#", 1)[0].strip()
        if not line:
            continue
        match = re.match(r"^([A-Za-z0-9._-]+)(\[[^\]]*\])?\s*(.*)$", line)
        assert match, f"unparsed requirement line: {line!r}"
        found[_normalize(match.group(1))] = match.group(3).replace(" ", "")
    return found


def pyproject_constraints() -> dict[str, str]:
    """Every pinned distribution across `dependencies` and all extras."""
    data = tomllib.loads(PYPROJECT.read_text())
    project = data["project"]
    entries = list(project.get("dependencies", []))
    for extra in project.get("optional-dependencies", {}).values():
        entries.extend(extra)

    found: dict[str, str] = {}
    for entry in entries:
        if "@" in entry:  # a direct URL, not a version constraint
            continue
        match = re.match(r"^([A-Za-z0-9._-]+)(\[[^\]]*\])?\s*(.*)$", entry.strip())
        assert match, f"unparsed pyproject dependency: {entry!r}"
        name = _normalize(match.group(1))
        if name == "combobatch":  # self-referential extra, e.g. combobatch[shambhala]
            continue
        found[name] = match.group(3).replace(" ", "")
    return found


class TestConstraintsAgree:
    def test_shared_packages_have_identical_constraints(self):
        requirements = requirement_constraints()
        pyproject = pyproject_constraints()
        mismatched = {
            name: (requirements[name], pyproject[name])
            for name in requirements.keys() & pyproject.keys()
            if requirements[name] != pyproject[name]
        }
        assert not mismatched, (
            "these packages are constrained differently in requirements.txt and "
            f"pyproject.toml, so the image and a laptop install can diverge: {mismatched}"
        )

    def test_nothing_is_pinned_only_in_requirements(self):
        extra = (
            requirement_constraints().keys()
            - pyproject_constraints().keys()
            - {_normalize(n) for n in EXPECTED_ONLY_IN_REQUIREMENTS}
        )
        assert not extra, (
            "these are in the image but in no pyproject extra, so a laptop install "
            f"silently lacks them: {sorted(extra)}"
        )
