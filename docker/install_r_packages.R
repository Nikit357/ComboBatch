#!/usr/bin/env Rscript
#
# R and Bioconductor packages for the combobatch image.
#
# `required` below is the union of every `r_packages` declaration in the registries:
#
#     python - <<'PY'
#     from combobatch.imputation import IMPUTER_REGISTRY
#     from combobatch.methods import METHOD_REGISTRY
#     pkgs = {p for s in METHOD_REGISTRY.values() for p in s.r_packages}
#     pkgs |= {p for s in IMPUTER_REGISTRY.values() for p in s.r_packages}
#     print(sorted(pkgs))
#     PY
#
# plus `variancePartition` for metric group F. tests/unit/test_docker_assets.py asserts
# this file and the registries agree, so a method added with a new R dependency fails the
# unit suite in a second rather than the image build in an hour.
#
# Four packages the donor script installed *and verified* are deliberately absent, each
# because it made `docker build` fail unconditionally at the terminal stop():
#   FAbatch    - exists on neither CRAN nor Bioconductor; the real package is `bapred`
#   reComBat   - a *Python* package; it arrives through requirements.txt instead
#   exploBATCH - hard dependency `fMM`, whose GitHub repository is deleted
# `DASC` is gone for a different reason: 35_dasc is not in METHOD_REGISTRY, so installing
# it would leave the image carrying a dependency no registered method can reach.

options(
    repos = c(CRAN = "https://cloud.r-project.org"),
    Ncpus = parallel::detectCores(),
    warn  = 1
)

# Every install is attempted even after one fails. The donor aborted on the first error,
# which meant a broken build reported one missing package per rebuild instead of all of
# them at once; the consolidated check at the bottom is what actually gates the image.
install_or_warn <- function(label, expr) {
    cat(sprintf("\n=== installing %s ===\n", label))
    tryCatch(
        force(expr),
        error = function(e) cat(sprintf("WARN: %s failed: %s\n", label, conditionMessage(e)))
    )
    invisible(NULL)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)

# remotes is needed before the GitHub-only packages below.
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# Hmisc needs a compiler and libpng/zlib present; qsmooth Imports it. When Hmisc fails to
# build, BiocManager reports qsmooth as "skipped" rather than failed, and 14_qsmooth then
# vanishes from the image without the build ever going red. Installing the heavy compiled
# dependencies first turns that into a visible failure here.
install_or_warn("compiled prerequisites (quantreg, Hmisc, lme4)", install.packages(
    c("quantreg", "Hmisc", "lme4"),
    dependencies = TRUE
))

install_or_warn("CRAN packages", install.packages(
    c(
        "missForest",   # imputer: missforest
        "softImpute",   # imputer: softimpute
        "DWDLargeR",    # 27_dwd
        "huge"          # 28_npn
    ),
    dependencies = TRUE
))

# bapred is a CRAN package, but installed through BiocManager so that its Bioconductor
# dependencies (affy, affyPLM, Biobase) resolve from the matching Bioc release rather
# than being reported as unavailable.
install_or_warn("bapred and its Bioconductor dependencies", BiocManager::install(
    c("affy", "affyPLM", "Biobase", "bapred"),
    ask = FALSE, update = FALSE
))

install_or_warn("Bioconductor packages", BiocManager::install(
    c(
        "limma",              # 03_limma, 04_sva
        "sva",                # 04_sva, 05_combat, 06_combat_seq, 29_combat_ref
        "RUVSeq",             # 09_ruv
        "batchelor",          # 10_mnn
        "qsmooth",            # 14_qsmooth
        "edgeR",              # 22_tmm
        "DESeq2",             # 23_vst
        "ruv",                # 31_ruv3prps
        "NOISeq",             # 34_arsyn
        "Harman",             # 38_harman
        "variancePartition",  # metric group F
        "BiocParallel",       # variancePartition's parallel backend
        "matrixStats"         # variancePartition dependency
    ),
    ask = FALSE, update = FALSE
))

# GitHub-only, and the fragile half of this file: each is one deleted repository away
# from disappearing. They are still *required* - every one backs a registered method - so
# a failure here fails the image, just not before the other installs have been attempted.
github_packages <- c(
    FSQN      = "jenniferfranks/FSQN",  # 16_fsqn_r; archived from CRAN
    TDM       = "greenelab/TDM",        # 19_tdm
    HarmonizR = "HSU-HPC/HarmonizR",    # 21_harmonizr
    DBNorm    = "mengqinxue/DBNorm",    # 33_amdbnorm
    AMDBNorm  = "JoevVan/AMDBNorm"      # 33_amdbnorm
)
for (pkg in names(github_packages)) {
    install_or_warn(
        sprintf("%s (github: %s)", pkg, github_packages[[pkg]]),
        remotes::install_github(github_packages[[pkg]], upgrade = "never")
    )
}

# --- one consolidated verification -------------------------------------------------
# Mirrors the registries exactly. Nothing unreachable is listed, so a failure here is
# always a real missing dependency and never a package that could not have existed.
required <- c(
    # imputation
    "missForest", "softImpute",
    # harmonization methods
    "AMDBNorm", "DBNorm", "DESeq2", "DWDLargeR", "FSQN", "Harman", "HarmonizR",
    "NOISeq", "RUVSeq", "TDM", "bapred", "batchelor", "edgeR", "huge", "limma",
    "qsmooth", "ruv", "sva",
    # metric group F
    "variancePartition"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
    stop(
        "Failed to install ", length(missing), " required package(s): ",
        paste(missing, collapse = ", ")
    )
}
cat(sprintf("\nAll %d required R packages installed successfully.\n", length(required)))
