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

## Changing the image

Run `bash scripts/smoke_test.sh` before pushing. It builds and then checks the R version,
Octave, the `matlab` wrapper, the registries, and that zero methods skip.
