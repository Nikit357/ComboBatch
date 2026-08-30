# CLAUDE.md — `docker/`

## Invariants

- **R is pinned to 4.5.3, deliberately.** rpy2 3.6.x has a C-level ABI incompatibility with
  R 4.6.0 whose symptom is a segfault around `33_amdbnorm`. Never replace the pinned Posit
  `.deb` with an unpinned `r-base` from the CRAN apt repository — that is the drift this
  pin exists to stop.
- **`install_r_packages.R` is checked against the registries**, not maintained by hand.
  `tests/unit/test_docker_assets.py` compares its `required` vector against the union of
  every `r_packages` declaration in both directions. A method that gains an R dependency
  fails the unit suite in a second rather than the image build in an hour.
- **Never re-add `FAbatch`, R `reComBat`, `exploBATCH` or `DASC`.** The first three cannot
  be installed at all and each aborted the donor's build; the fourth backs a method that is
  not registered.
- **Every install is attempted before the build fails.** One consolidated check at the end
  reports all missing packages, not just the first — a broken build should cost one rebuild,
  not one per defect.
- **Build for `linux/amd64`.** The dev machine is arm64, the nodes are not, and the pinned
  R build is published for amd64 only. Always go through `scripts/build_and_push_image.sh`.
- **`ENTRYPOINT` is `entrypoint.sh`, not `combobatch`.** A literal `combobatch` entry point
  would make `docker run IMAGE pytest …`, `Rscript -e …`, `octave --version` and the pod's
  `sleep infinity` override impossible.
- **`COMBOBATCH_FULL_ENV=1` must stay set.** Without it
  `tests/integration/test_methods_all_backends.py` self-skips, and an image that silently
  lost a backend would pass its own acceptance gate.
- **No Procrustes.** No clone, no `sys.path` insertion, no optional dependency: the licence
  is incompatible with an MIT tool.
- **Diagnose the failure you had, not the one you expected.** The R download once answered
  every `curl` failure with "R 4.5.3 is not published for `$(dpkg --print-architecture)`",
  so a certificate error behind a TLS-inspecting VPN read as an architecture problem on a
  build that already was amd64. The architecture is now *checked* and curl's exit code
  chooses the message. Never widen a catch-all to cover a failure it cannot identify.
- **No Python dependency may drag in `torch`.** `harmonypy` 0.2.0 did, and with it
  `nvidia/*` and `triton`: 4.5 GB of CUDA in a CPU-only image for `c6a` nodes, taking the
  build to 12.8 GB against a documented 6–9. It is pinned `<0.2`; check
  `du -sm /usr/local/lib/python3.11/site-packages/* | sort -rn | head` after any
  dependency change.
- **reComBat is installed alone, at a pinned commit, past its own metadata.** It declares
  `sklearn` (the deprecated stub, which only raises when built) and `python >=3.8,<3.11`;
  both are stale, and it runs correctly on 3.11. Hence `--no-deps
  --ignore-requires-python` in its own layer, followed by an import check — never in
  `requirements.txt` (`--no-deps` has no per-requirement form), never `sklearn` or `fire`
  added to satisfy it, and **never the PyPI wheel**: the published 0.1.4 defaults to
  `model='linear'` where the pinned commit defaults to `model='elastic_net'`, so switching
  would silently change every `30_recombat` result.
- **The R install is staged, and stays staged.** Six `RUN` layers — `prereq`, `cran`,
  `bapred`, `bioc`, `github`, `verify` — because the install needs verifiable TLS for
  15–30 minutes and does not always get it: on 2026-08-30 a VPN reconnected 14 minutes in
  and a single `RUN` discarded every finished compile. Each stage verifies its own packages
  before exiting, so that a cached layer can be trusted; never collapse them back into one,
  and never let a stage exit 0 without checking what it installed.
- **Never install R packages with `dependencies = TRUE`.** It adds `Suggests`, which no
  registered method imports, and each one is a new system library the apt layer does not
  have. Every package in the image should be traceable to a registry entry or to a hard
  `Imports` edge of one. The corollary: when a package will not build, look for the
  missing `-dev` library before touching the package list — `libpng-dev`, `zlib1g-dev` and
  `libuv1-dev` are each load-bearing for a chain that ends at a registered method.
- **Never trust an extra CA in this image.** Installing a corporate TLS-inspection root
  would make the build work behind such a proxy and would ship a MITM anchor to everyone
  who pulls a public image. Build with the VPN disconnected instead — `BUILD_AND_PUSH.md`
  step 0b.

## Changing the image

Run `bash scripts/smoke_test.sh` before pushing. It builds and then checks the R version,
Octave, the `matlab` wrapper, the registries, and that zero methods skip.
