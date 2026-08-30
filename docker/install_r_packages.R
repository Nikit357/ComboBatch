#!/usr/bin/env Rscript
#
# R and Bioconductor packages for the combobatch image.
#
# Run one stage at a time:  Rscript install_r_packages.R --stage bioc
# or the whole thing:       Rscript install_r_packages.R
#
# `STAGE_PACKAGES` below is the union of every `r_packages` declaration in the registries:
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
    warn  = 1,
    # R's default is 60 seconds per file, which a large Bioconductor tarball on a slow
    # link can exceed - and the resulting message is indistinguishable from a real
    # network failure.
    timeout = 300
)

# What each stage installs, and the single source of the final verification vector.
# The install is split into stages because it is not resumable otherwise: on 2026-08-30 a
# VPN reconnected 14 minutes in and every completed compile was discarded along with the
# failure. Each stage is its own Docker layer, so a stage that succeeded stays cached and
# a rebuild resumes at the one that did not.
STAGE_PACKAGES <- list(
    prereq = character(0),  # quantreg/Hmisc/lme4: no method's backend, checked separately
    cran = c(
        "missForest",   # imputer: missforest
        "softImpute",   # imputer: softimpute
        "DWDLargeR",    # 27_dwd
        "huge",         # 28_npn
        "binr"          # 19_tdm: TDM's load_it() installs it interactively otherwise
    ),
    bapred = "bapred",      # 37_fabatch
    bioc = c(
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
        "variancePartition"   # metric group F
    ),
    github = c("FSQN", "TDM", "HarmonizR", "DBNorm", "AMDBNorm")
)
required <- unlist(STAGE_PACKAGES, use.names = FALSE)

args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args) >= 2 && args[[1]] == "--stage") args[[2]] else "all"
valid_stages <- c(names(STAGE_PACKAGES), "verify", "all")
if (!stage %in% valid_stages) {
    stop(
        "unknown --stage ", stage,
        "; expected one of: ", paste(valid_stages, collapse = ", ")
    )
}
runs <- function(name) stage %in% c(name, "all")

# Recorded rather than only printed, so the failure message can name the actual cause.
# The 2026-08-29 message blamed a missing system library; the 2026-08-30 failure was TLS,
# and 185 warnings said so while the message still pointed at apt.
FAILURES <- new.env(parent = emptyenv())
FAILURES$messages <- character(0)

# Every install is attempted even after one fails, so one rebuild reports every missing
# package rather than one per attempt; the per-stage check below is what gates the layer.
#
# Everything here goes to stderr, which is unbuffered. On stdout it is block-buffered
# whenever the build is not a terminal, so on 2026-08-30 all 273 warnings appeared at once
# with the exit timestamp, after 16 minutes of a silent screen and in an order that no
# longer matched when anything happened.
install_or_warn <- function(label, expr) {
    cat(sprintf("\n=== installing %s ===\n", label), file = stderr())
    # Both handlers are needed. install.packages() reports a package that would not build
    # as a *warning* and returns normally, so the error arm alone stayed silent through 41
    # package failures on 2026-08-29 and no WARN: line was ever printed; BiocManager can
    # also stop() outright. Muffling keeps each failure on one greppable line.
    withCallingHandlers(
        tryCatch(
            force(expr),
            error = function(e) {
                FAILURES$messages <- c(FAILURES$messages, conditionMessage(e))
                cat(
                    sprintf("WARN: %s failed: %s\n", label, conditionMessage(e)),
                    file = stderr()
                )
            }
        ),
        warning = function(w) {
            FAILURES$messages <- c(FAILURES$messages, conditionMessage(w))
            cat(
                sprintf("WARN: %s: %s\n", label, conditionMessage(w)),
                file = stderr()
            )
            invokeRestart("muffleWarning")
        }
    )
    invisible(NULL)
}

# A download failure and a build failure need different answers, and guessing between them
# is what made both previous failure messages point somewhere other than the cause.
stop_with_diagnosis <- function(absent, stage_name) {
    downloads_failed <- any(grepl(
        paste(
            "SSL peer certificate", "cannot download any files",
            "download of package", "unable to access", "Timeout",
            sep = "|"
        ),
        FAILURES$messages
    ))
    if (downloads_failed) {
        stop(
            "stage ", stage_name, " could not install: ",
            paste(absent, collapse = ", "),
            ". Downloads failed certificate verification or timed out partway through ",
            "this build - a VPN or proxy re-established TLS interception while it was ",
            "running. The Dockerfile's preflight probe only checks the moment the build ",
            "starts. Disconnect it, confirm it cannot reconnect, and rebuild: the stages ",
            "that already succeeded are cached and will not be repeated."
        )
    }
    stop(
        "stage ", stage_name, " could not install: ", paste(absent, collapse = ", "),
        ". Search the log for 'WARN:' and 'ERROR: configuration failed' - a package ",
        "listed here is often the last casualty of a dependency that failed earlier for ",
        "want of a system library."
    )
}

# Verified per stage, not only at the end: a stage that exits 0 gets cached by BuildKit,
# so a stage must never exit 0 with its own packages missing. Without this, staging would
# be a regression - a rebuild would skip the very layer that failed.
verify_stage <- function(stage_name) {
    expected <- STAGE_PACKAGES[[stage_name]]
    if (!length(expected)) {
        return(invisible(NULL))
    }
    absent <- expected[!vapply(expected, requireNamespace, logical(1), quietly = TRUE)]
    if (length(absent)) stop_with_diagnosis(absent, stage_name)
    cat(
        sprintf("stage %s: all %d package(s) present\n", stage_name, length(expected)),
        file = stderr()
    )
}

# Bootstrapped per invocation but installed once: each stage is a separate R process, and
# the library persists in the image from one layer to the next.
if (runs("bapred") || runs("bioc")) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
    BiocManager::install(version = "3.22", ask = FALSE, update = FALSE)
}
if (runs("github")) {
    if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
}

# --- prereq ------------------------------------------------------------------------
# Hmisc needs a compiler, libpng/zlib and - through rmarkdown -> bslib -> sass -> fs -
# libuv present; qsmooth Imports it. Installing the heavy compiled dependencies first is
# what makes that a visible failure here rather than a `qsmooth` that quietly never
# arrives. It does not go unnoticed either way: qsmooth is in STAGE_PACKAGES, so the
# 2026-08-29 build did go red - it just named the casualty instead of the cause.
# No `dependencies = TRUE`: that adds Suggests, and on 2026-08-29 it pulled V8, shiny,
# rgl, plotly, gsl, sodium, gt, kableExtra, DHARMa, semEff, qgam, testthat, reactable,
# sparkline, juicyjuice and safer into the build - 860 seconds, five more system libraries
# required, and nothing in the registries imports any of them. The default (Depends,
# Imports, LinkingTo) is every edge that can break a library() call.
if (runs("prereq")) {
    install_or_warn(
        "compiled prerequisites (quantreg, Hmisc, lme4)",
        install.packages(c("quantreg", "Hmisc", "lme4"))
    )

    # Checked here rather than in the final block: none of these is a method backend, so
    # none is in STAGE_PACKAGES, yet qsmooth Imports Hmisc and cannot install without it.
    # Failing here names the package that actually broke, 20 minutes before the
    # consolidated check would report a symptom five dependency levels downstream.
    prerequisites <- c("quantreg", "Hmisc", "lme4")
    absent <- prerequisites[
        !vapply(prerequisites, requireNamespace, logical(1), quietly = TRUE)
    ]
    if (length(absent)) stop_with_diagnosis(absent, "prereq")
    verify_stage("prereq")
}

# --- cran --------------------------------------------------------------------------
if (runs("cran")) {
    install_or_warn("CRAN packages", install.packages(STAGE_PACKAGES$cran))
    verify_stage("cran")
}

# --- bapred ------------------------------------------------------------------------
# bapred is a CRAN package, but installed through BiocManager so that its Bioconductor
# dependencies (affy, affyPLM, Biobase) resolve from the matching Bioc release rather
# than being reported as unavailable.
if (runs("bapred")) {
    install_or_warn("bapred and its Bioconductor dependencies", BiocManager::install(
        c("affy", "affyPLM", "Biobase", "bapred"),
        ask = FALSE, update = FALSE
    ))
    verify_stage("bapred")
}

# --- bioc --------------------------------------------------------------------------
# BiocParallel and matrixStats are variancePartition's backend and dependency; they are
# installed explicitly but are not method backends, so they are not verified.
if (runs("bioc")) {
    install_or_warn("Bioconductor packages", BiocManager::install(
        c(STAGE_PACKAGES$bioc, "BiocParallel", "matrixStats"),
        ask = FALSE, update = FALSE
    ))
    verify_stage("bioc")
}

# --- github ------------------------------------------------------------------------
# The fragile half of this file: each is one deleted repository away from disappearing.
# They are still *required* - every one backs a registered method.
github_repos <- c(
    FSQN      = "jenniferfranks/FSQN",  # 16_fsqn_r; archived from CRAN
    TDM       = "greenelab/TDM",        # 19_tdm
    HarmonizR = "HSU-HPC/HarmonizR",    # 21_harmonizr
    DBNorm    = "mengqinxue/DBNorm",    # 33_amdbnorm
    AMDBNorm  = "JoevVan/AMDBNorm"      # 33_amdbnorm
)
if (runs("github")) {
    for (pkg in names(github_repos)) {
        install_or_warn(
            sprintf("%s (github: %s)", pkg, github_repos[[pkg]]),
            remotes::install_github(github_repos[[pkg]], upgrade = "never")
        )
    }
    verify_stage("github")
}

# --- verify ------------------------------------------------------------------------
# Every stage has already checked its own packages, so reaching this with something
# missing means a layer was cached without what it claims to install. Rebuild that stage
# with --no-cache rather than looking for a broken dependency.
if (stage %in% c("verify", "all")) {
    missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing)) {
        stop(
            "Failed to install ", length(missing), " required package(s): ",
            paste(missing, collapse = ", "),
            ". Each stage verifies itself, so this means a cached layer is missing what ",
            "it installed; rebuild with --no-cache."
        )
    }
    cat(
        sprintf("\nAll %d required R packages installed successfully.\n", length(required)),
        file = stderr()
    )
}
