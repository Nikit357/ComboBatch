"""The image's dependency lists must agree with the registries.

`docker/install_r_packages.R` is the one place in this project where a registry is
restated in another language, and there is no way for the build to notice a drift: a
method whose R package was never installed does not fail the build, it merely SKIPs at
run time and quietly disappears from the benchmark. Building the image to find that out
takes about an hour, so the same check runs here in a second.

The rest of these tests pin the Dockerfile invariants that were each established by a
documented failure, listed in the plan's section 7.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from combobatch.imputation import IMPUTER_REGISTRY
from combobatch.methods import METHOD_REGISTRY
from combobatch.metrics import METRIC_GROUP_REGISTRY

REPO_ROOT = Path(__file__).resolve().parents[2]
DOCKERFILE = REPO_ROOT / "docker" / "Dockerfile"
R_SCRIPT = REPO_ROOT / "docker" / "install_r_packages.R"
ENTRYPOINT = REPO_ROOT / "docker" / "entrypoint.sh"
BUILD_SCRIPT = REPO_ROOT / "scripts" / "build_and_push_image.sh"

# Group F shells out to Rscript for variancePartition; MetricGroupSpec records only that
# R is needed, not which package, so this one name is named here.
METRICS_R_PACKAGES = frozenset({"variancePartition"})

# Installed and verified by the donor script, and each one aborted `docker build` at its
# terminal stop(). None may come back.
IMPOSSIBLE_PACKAGES = ("FAbatch", "exploBATCH")


def code_only(path: Path) -> str:
    """Strip comment lines, so a name mentioned only in an explanation is not a hit."""
    return "\n".join(
        line
        for line in path.read_text().splitlines()
        if not line.lstrip().startswith("#")
    )


def registry_r_packages() -> frozenset[str]:
    """Every R package the three registries declare a need for."""
    packages: set[str] = set(METRICS_R_PACKAGES)
    for spec in METHOD_REGISTRY.values():
        packages.update(spec.r_packages)
    for spec in IMPUTER_REGISTRY.values():
        packages.update(spec.r_packages)
    return frozenset(packages)


def verified_packages() -> frozenset[str]:
    """The names in the R script's final `required <- c(...)` vector."""
    text = R_SCRIPT.read_text()
    match = re.search(r"^required <- c\((.*?)^\)", text, re.MULTILINE | re.DOTALL)
    assert match, "the R script must end with a `required <- c(...)` vector"
    return frozenset(re.findall(r'"([^"]+)"', match.group(1)))


class TestRPackagesMatchTheRegistries:
    def test_every_registry_package_is_verified(self):
        missing = registry_r_packages() - verified_packages()
        assert not missing, (
            "these R packages are needed by a registered method, imputer or metric "
            f"group but are not verified by install_r_packages.R: {sorted(missing)}"
        )

    def test_nothing_is_verified_that_no_registry_needs(self):
        extra = verified_packages() - registry_r_packages()
        assert not extra, (
            "install_r_packages.R verifies packages no registry asks for, which is how "
            f"the donor's build came to fail on a package that cannot exist: {sorted(extra)}"
        )

    def test_every_verified_package_is_also_installed(self):
        """A name in the verification vector but in no install call fails the build."""
        text = R_SCRIPT.read_text()
        install_section = text.split("# --- one consolidated verification")[0]
        for package in sorted(verified_packages()):
            assert package in install_section, (
                f"{package} is verified but never installed; the build would fail with "
                f"it reported as missing"
            )

    @pytest.mark.parametrize("package", IMPOSSIBLE_PACKAGES)
    def test_the_impossible_packages_are_gone(self, package):
        assert package not in code_only(R_SCRIPT)

    def test_the_python_recombat_is_not_installed_as_an_r_package(self):
        # reComBat is a Python package; the donor installed it from GitHub as R and then
        # verified it, which no build could ever satisfy. It arrives via requirements.txt.
        code = code_only(R_SCRIPT)
        assert "bioFAM/reComBat" not in code
        assert '"reComBat"' not in code
        assert "reComBat" in (REPO_ROOT / "requirements.txt").read_text()

    def test_dropped_methods_leave_no_dependency_behind(self):
        """A method removed from the registry must not still be installed."""
        assert "DASC" not in code_only(R_SCRIPT), "35_dasc is not a registered method"

    def test_only_one_metric_group_needs_r(self):
        """``MetricGroupSpec`` records that R is needed, not which package.

        So ``METRICS_R_PACKAGES`` above is the one hand-maintained name in this file, and
        a second R-backed group would slip past every other check here.
        """
        r_groups = {
            letter for letter, spec in METRIC_GROUP_REGISTRY.items() if spec.needs_r
        }
        assert r_groups == {"F"}, (
            f"metric group(s) {sorted(r_groups - {'F'})} now need R; add their packages "
            f"to METRICS_R_PACKAGES and to install_r_packages.R"
        )


class TestDockerfileInvariants:
    def test_r_is_pinned_to_4_5_3(self):
        assert "ARG R_VERSION=4.5.3" in DOCKERFILE.read_text()
        # Unpinned r-base from the CRAN apt repo is exactly the drift that segfaults
        # rpy2; it survives in this file only as a commented-out fallback.
        assert "r-base" not in code_only(DOCKERFILE)

    def test_the_base_image_is_python_3_11_on_bookworm(self):
        # pandas<2.0 has a py3.11 wheel only, and the pinned R build is Debian 12.
        assert "FROM python:3.11-slim-bookworm" in DOCKERFILE.read_text()

    def test_octave_and_the_matlab_wrapper_are_installed(self):
        text = DOCKERFILE.read_text()
        assert "octave-statistics" in text
        assert "/usr/local/bin/matlab" in text

    def test_the_full_environment_flag_is_set(self):
        # tests/integration/test_methods_all_backends.py self-skips without it, which
        # would let an image with a broken backend pass its own acceptance gate.
        assert "COMBOBATCH_FULL_ENV=1" in DOCKERFILE.read_text()

    def test_threads_are_pinned_in_the_image_and_in_etc_environment(self):
        text = DOCKERFILE.read_text()
        assert "ENV OMP_NUM_THREADS=1" in text
        # Pod-level env vars never reach a process launched over SSH; PAM reads this file.
        assert "/etc/environment" in text

    def test_no_proprietary_dependency_is_vendored(self):
        for path in (DOCKERFILE, R_SCRIPT):
            assert "procrustes" not in path.read_text().lower()

    def test_the_entrypoint_exists_and_is_wired_up(self):
        assert 'ENTRYPOINT ["/usr/local/bin/entrypoint.sh"]' in DOCKERFILE.read_text()
        assert ENTRYPOINT.exists()

    def test_the_entrypoint_runs_any_command_it_is_given(self):
        # `ENTRYPOINT ["combobatch"]` would break `docker run IMAGE pytest ...`,
        # `Rscript -e ...`, `octave --version` and the pod's `sleep infinity`.
        assert 'exec "$@"' in ENTRYPOINT.read_text()


class TestEntrypointDispatch:
    """The entry point is a shell script, so run it rather than grep it."""

    @pytest.fixture
    def run_entrypoint(self, tmp_path):
        import subprocess

        bin_dir = tmp_path / "bin"
        bin_dir.mkdir()
        for name in ("combobatch", "pytest"):
            stub = bin_dir / name
            stub.write_text(f'#!/bin/sh\necho "{name}: $*"\n')
            stub.chmod(0o755)

        def run(*args):
            env = {"PATH": f"{bin_dir}:/usr/bin:/bin"}
            done = subprocess.run(
                ["sh", str(ENTRYPOINT), *args],
                capture_output=True,
                text=True,
                env=env,
                timeout=30,
            )
            assert done.returncode == 0, done.stderr
            return done.stdout.strip()

        return run

    def test_no_arguments_prints_the_cli_help(self, run_entrypoint):
        assert run_entrypoint() == "combobatch: --help"

    def test_a_leading_flag_goes_to_the_cli(self, run_entrypoint):
        assert run_entrypoint("--version") == "combobatch: --version"

    def test_a_bare_subcommand_goes_to_the_cli(self, run_entrypoint):
        assert (
            run_entrypoint("run", "--config", "c.yaml")
            == "combobatch: run --config c.yaml"
        )

    def test_an_explicit_cli_invocation_is_passed_through(self, run_entrypoint):
        assert run_entrypoint("combobatch", "selftest") == "combobatch: selftest"

    def test_another_executable_runs_as_itself(self, run_entrypoint):
        # `docker run IMAGE pytest ...` is how the image asserts every backend is present.
        assert run_entrypoint("pytest", "tests/unit") == "pytest: tests/unit"

    def test_the_pods_command_override_still_works(self, run_entrypoint):
        assert run_entrypoint("echo", "sleep", "infinity") == "sleep infinity"


class TestBuildScript:
    def test_the_build_is_amd64(self):
        # The dev Mac is arm64; a native build gives the pod an `exec format error`.
        assert "--platform linux/amd64" in BUILD_SCRIPT.read_text()

    def test_ghcr_is_the_only_registry(self):
        text = BUILD_SCRIPT.read_text()
        assert "ghcr.io/nikit357/combobatch" in text
        assert "dkr.ecr" not in text, "ComboBatch touches no AWS container registry"
