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


def stage_packages() -> dict[str, list[str]]:
    """The R script's `STAGE_PACKAGES` list, one entry per build stage.

    `required` is derived from this in the R script itself, so this is the single place
    the image's package set is stated.
    """
    text = R_SCRIPT.read_text()
    match = re.search(
        r"^STAGE_PACKAGES <- list\((.*?)^\)", text, re.MULTILINE | re.DOTALL
    )
    assert match, "the R script must define a `STAGE_PACKAGES <- list(...)`"
    stages: dict[str, list[str]] = {}
    for name, body in re.findall(
        r"(\w+)\s*=\s*(character\(0\)|\"[^\"]+\"|c\([^)]*\))", match.group(1)
    ):
        stages[name] = re.findall(r'"([^"]+)"', body)
    return stages


def verified_packages() -> frozenset[str]:
    """Every package the R script verifies — the union across all stages."""
    return frozenset(p for names in stage_packages().values() for p in names)


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
        """A name that is verified but never installed fails the build.

        The `cran` and `bioc` stages pass ``STAGE_PACKAGES`` straight to their installer,
        so they cannot drift. The other two restate names, and can.
        """
        text = R_SCRIPT.read_text()
        install_section = text.split("# --- verify ---")[0]
        assert '"bapred"' in install_section

        repos = re.search(r"github_repos <- c\((.*?)^\)", text, re.M | re.S)
        assert repos, "the R script must define `github_repos`"
        cloned = set(re.findall(r"^\s*(\w+)\s*=", repos.group(1), re.M))
        assert cloned == set(stage_packages()["github"]), (
            "every GitHub package verified must have a repository to clone from, and "
            "vice versa"
        )

    @pytest.mark.parametrize("package", IMPOSSIBLE_PACKAGES)
    def test_the_impossible_packages_are_gone(self, package):
        assert package not in code_only(R_SCRIPT)

    def test_the_python_recombat_is_not_installed_as_an_r_package(self):
        # reComBat is a Python package; the donor installed it from GitHub as R and then
        # verified it, which no build could ever satisfy. It arrives through pip — in its
        # own Dockerfile layer rather than requirements.txt, see TestReComBatInstall.
        code = code_only(R_SCRIPT)
        assert "bioFAM/reComBat" not in code
        assert '"reComBat"' not in code
        assert "reComBat" in code_only(DOCKERFILE)

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


class TestReComBatInstall:
    """reComBat cannot be installed the way it declares itself.

    Its metadata names `sklearn` — the deprecated stub, which only raises when built —
    and `python >=3.8,<3.11`, which excludes the version this project pins. Both are
    stale: it runs correctly on 3.11. So it is installed alone, past both, at a pin.
    """

    # The whole block, not the lines naming reComBat: the flags sit on the `RUN pip
    # install` line and the package on its continuation, so a line filter splits them.
    RECOMBAT_LINES = re.search(
        r"ARG RECOMBAT_COMMIT.*?(?=\n\n)", DOCKERFILE.read_text(), re.S
    )
    RECOMBAT_LINES = RECOMBAT_LINES.group(0) if RECOMBAT_LINES else ""

    def test_it_is_installed_past_its_own_metadata(self):
        assert "--no-deps" in self.RECOMBAT_LINES
        assert "--ignore-requires-python" in self.RECOMBAT_LINES

    def test_it_is_pinned_to_a_commit(self):
        """An unpinned `master` is how a silent behaviour change would arrive.

        It nearly did: the published 0.1.4 wheel and this commit differ in one line —
        the default `model`, `linear` against `elastic_net`.
        """
        assert re.search(
            r"RECOMBAT_COMMIT=[0-9a-f]{40}", self.RECOMBAT_LINES
        ), "reComBat must be pinned to a full commit SHA, not a moving branch"

    def test_it_is_not_taken_from_pypi(self):
        """The wheel installs in a second and defaults to a different model."""
        assert "git+https://github.com/BorgwardtLab/reComBat" in self.RECOMBAT_LINES

    def test_the_install_is_checked_by_importing_it(self):
        """`--ignore-requires-python` would otherwise hide a real incompatibility."""
        assert "from reComBat import reComBat" in self.RECOMBAT_LINES

    def test_the_stub_is_never_installed_to_satisfy_the_metadata(self):
        text = (REPO_ROOT / "requirements.txt").read_text()
        assert not re.search(r"^sklearn\b", text, re.M), (
            "the `sklearn` PyPI package is a stub that installs no module; the real "
            "dependency is scikit-learn, which is already pinned"
        )

    def test_tqdm_is_declared(self):
        """--no-deps drops it, and reComBat imports it at module level.

        Without this it arrives only as a transitive dependency of umap-learn.
        """
        assert re.search(
            r"^tqdm>=", (REPO_ROOT / "requirements.txt").read_text(), re.M
        ), "reComBat imports tqdm, and --no-deps means this file must say so"


class TestStagedInstall:
    """The R install is staged so a mid-build failure costs one stage, not all of them.

    On 2026-08-30 a VPN reconnected 14 minutes into the build and 857 seconds of finished
    compiles were discarded, because they lived in the same ``RUN`` as the failure.
    """

    def test_the_dockerfile_runs_every_stage(self):
        invoked = re.findall(r"--stage (\w+)", DOCKERFILE.read_text())
        assert set(stage_packages()) <= set(invoked), (
            "a stage defined in the R script is never run by the image, so its packages "
            "would only be missed by the final verification"
        )

    def test_verify_runs_last(self):
        invoked = re.findall(r"--stage (\w+)", DOCKERFILE.read_text())
        assert invoked[-1] == "verify", "the whole-image check must come after the rest"

    def test_the_stages_are_separate_layers(self):
        """One `RUN` per stage is the entire point — collapsing them undoes the fix."""
        runs = [
            line
            for line in code_only(DOCKERFILE).splitlines()
            if "install_r_packages.R" in line and line.lstrip().startswith("RUN")
        ]
        assert len(runs) >= len(stage_packages()) + 1, (
            f"expected one RUN per stage plus verify, found {len(runs)}; a single RUN is "
            f"the all-or-nothing layer this staging exists to remove"
        )

    def test_each_stage_verifies_itself_before_the_layer_is_cached(self):
        """A stage exiting 0 gets cached, so it must not exit 0 with packages missing."""
        assert "verify_stage <- function" in code_only(R_SCRIPT)
        assert code_only(R_SCRIPT).count("verify_stage(") >= len(stage_packages())


class TestRInstallScript:
    """The 2026-08-29 build died 36 minutes in, on omissions this file can catch."""

    def test_suggests_are_not_installed(self):
        """``dependencies = TRUE`` adds Suggests, which nothing here imports.

        It pulled V8, shiny, rgl, plotly, gsl, sodium and ten more into the build —
        860 seconds, and five system libraries the apt layer does not carry.
        """
        assert "dependencies = TRUE" not in code_only(R_SCRIPT)

    def test_the_prerequisites_are_checked_before_the_long_installs(self):
        """Hmisc is in no registry, so only a checkpoint here can name it.

        qsmooth Imports Hmisc; when Hmisc failed, the consolidated check at the end
        reported ``qsmooth`` — a symptom five dependency levels downstream — after
        another 1247 seconds of installs whose outcome was already decided.
        """
        text = R_SCRIPT.read_text()
        checkpoint = text.find('stop_with_diagnosis(absent, "prereq")')
        bioconductor = text.find('install_or_warn("Bioconductor packages"')
        assert checkpoint != -1, "the prerequisite checkpoint is missing"
        assert (
            checkpoint < bioconductor
        ), "the checkpoint must precede the long installs"

    def test_failed_installs_are_reported_even_though_they_are_warnings(self):
        """``install.packages`` downgrades an unbuildable package to a warning.

        With only an error handler, 41 package failures produced no ``WARN:`` line.
        """
        code = code_only(R_SCRIPT)
        assert "withCallingHandlers(" in code
        assert 'invokeRestart("muffleWarning")' in code

    def test_warnings_go_to_stderr(self):
        """stdout is block-buffered off a terminal.

        On 2026-08-30 that put all 273 warnings at the exit timestamp, after 16 minutes
        of a silent screen and in an order that no longer matched events.
        """
        code = code_only(R_SCRIPT)
        warn_lines = [line for line in code.splitlines() if "WARN: %s" in line]
        assert warn_lines, "the WARN: reporting is missing"
        for line in warn_lines:
            assert "stderr()" in line or "file = stderr()" in code

    def test_the_failure_message_distinguishes_a_download_from_a_build_failure(self):
        """Both earlier messages named a cause they had not identified."""
        code = code_only(R_SCRIPT)
        assert "SSL peer certificate" in code
        assert "cannot download any files" in code


class TestDockerfileInvariants:
    def test_r_is_pinned_to_4_5_3(self):
        assert "ARG R_VERSION=4.5.3" in DOCKERFILE.read_text()
        # Unpinned r-base from the CRAN apt repo is exactly the drift that segfaults
        # rpy2; it survives in this file only as a commented-out fallback.
        assert "r-base" not in code_only(DOCKERFILE)

    def test_the_r_download_diagnoses_the_failure_it_actually_had(self):
        """Every curl failure was once reported as an architecture problem.

        A TLS error behind a full-tunnel VPN came out as "R 4.5.3 is not published for
        amd64" on a build that already was amd64, which cost an hour. The architecture
        is now tested, and curl's exit code decides the message.
        """
        code = code_only(DOCKERFILE)
        assert (
            '[ "$arch" = "amd64" ]' in code
        ), "the architecture must be checked, not inferred"
        assert "60)" in code, "curl's certificate failure needs a message of its own"
        assert "is not published for $(dpkg --print-architecture)" not in code

    def test_the_tls_preflight_runs_before_the_r_download(self):
        """apt uses http, so a TLS problem otherwise surfaces 155 seconds in.

        The probe deliberately omits ``-f``: it tests TLS, not availability, and the
        CDN root answers 403 by design — with ``-f`` that is exit 56 and every build
        would fail here.
        """
        text = DOCKERFILE.read_text()
        probe = text.find("curl -sS -o /dev/null --max-time 30 https://cdn.posit.co/")
        download = text.find('curl -fsSL -o /tmp/r.deb "$url"')
        assert probe != -1, "the preflight TLS probe is missing"
        assert probe < download, "the probe must run before the first large download"

    def test_the_load_bearing_system_libraries_are_installed(self):
        """Each of these three costs an hour of build time when it is absent.

        libuv is the least obvious: fs 2.x requires it, and fs is reached through
        qsmooth -> Hmisc -> rmarkdown -> bslib -> sass, so its absence surfaces as a
        missing qsmooth five levels away.
        """
        code = code_only(DOCKERFILE)
        for library in ("libpng-dev", "zlib1g-dev", "libuv1-dev"):
            assert library in code, f"{library} is load-bearing for the R package layer"

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
