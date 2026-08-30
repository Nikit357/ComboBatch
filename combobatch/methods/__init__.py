"""``METHOD_REGISTRY`` — the single source of truth for harmonization methods.

The donor kept a ``METHODS`` dict in its core library and then *retyped the same list by
hand* in two dispatcher scripts, so adding a method meant editing three files and
forgetting one meant that dispatcher silently skipped it. Here every CLI choice list,
validation check, preflight probe and generated doc reads this one dict.

**34 methods are registered here**, including ``20_shambhala``, which is backed by the
vendored pure-Python + Octave implementation in :mod:`combobatch.methods.shambhala_method`
rather than by the donor's ``Shambhala2`` R package.

Original keys are preserved, so the numbering has gaps at 24, 32, 35, 36 and 39. Those
five are excluded — four can never run, and one is licence-restricted; see
``docs/METHODS.md``. Keeping the original keys means a ComboBatch result can still be
traced to the published benchmark it came from.
"""

from __future__ import annotations

from combobatch.methods import python_methods as _py
from combobatch.methods import r_methods as _r
from combobatch.methods import shambhala_method as _shambhala
from combobatch.methods.base import (
    HyperParam,
    MethodSpec,
    assert_rnaseq_only,
    drop_na_genes,
    resolve_reference_batch,
)
from combobatch.methods.python_methods import DEFAULT_CONTROL_GENES

# Method keys deliberately absent from the registry, with the reason. Surfaced by
# `combobatch list-methods --show-excluded` and asserted by the test suite, so the gap
# between "39 in the article" and "34 here" is always explicable.
EXCLUDED_METHODS: dict[str, str] = {
    "24_peer_k10": "PEER is unavailable for R 4.5 (not on CRAN, Bioconductor or GitHub)",
    "32_deepmnn": "scRNA-seq only; needs AnnData + HVG + PCA, not applicable to bulk",
    "35_dasc": "returns cluster assignments, not a corrected expression matrix",
    "36_explobatch": "hard dependency fMM had its GitHub repository deleted",
    "39_procrustes": (
        "licence-restricted: BostonGene proprietary, non-commercial research use only, "
        "and forbids altering copyright notices - incompatible with an MIT tool"
    ),
}

_TARGET_GROUP = HyperParam(
    name="target_group",
    type=str,
    default=None,
    help=(
        "Reference batch to normalize toward. Defaults to the largest batch, resolved "
        "at run time; set --reference-batch to choose explicitly."
    ),
)


def _spec(key: str, fn, harshness: str, **kwargs) -> MethodSpec:
    """Build a MethodSpec, defaulting hyperparams to an empty dict."""
    kwargs.setdefault("hyperparams", {})
    return MethodSpec(key=key, fn=fn, harshness=harshness, **kwargs)


METHOD_REGISTRY: dict[str, MethodSpec] = {
    # ── Pure Python ──────────────────────────────────────────────────────────────
    "01_raw": _spec(
        "01_raw",
        _py.normalize_raw,
        "low",
        citation="No correction; the baseline every metric is read against.",
    ),
    "02_median_scaling": _spec(
        "02_median_scaling",
        _py.normalize_median_scaling,
        "low",
        citation="Per-batch median centring.",
    ),
    "07_pycombat": _spec(
        "07_pycombat",
        _py.normalize_pycombat,
        "medium",
        citation="Behdenna et al. 2023, BMC Bioinformatics (pyComBat).",
        python_packages=("combat",),
    ),
    "08_inmoose_combatseq": _spec(
        "08_inmoose_combatseq",
        _py.normalize_inmoose_combat_seq,
        "medium",
        citation="Zhang et al. 2020, NAR Genom Bioinform (ComBat-seq), via InMoose.",
        python_packages=("inmoose",),
        uses_bio_col=True,
    ),
    "11_harmony": _spec(
        "11_harmony",
        _py.normalize_harmony,
        "medium",
        citation="Korsunsky et al. 2019, Nat Methods.",
        python_packages=("harmonypy",),
        hyperparams={
            "n_pcs": HyperParam(
                "n_pcs", int, 50, "Principal components used for correction.", minimum=2
            ),
            "max_iter": HyperParam(
                "max_iter", int, 20, "Maximum Harmony iterations.", minimum=1
            ),
        },
    ),
    "12_scanorama": _spec(
        "12_scanorama",
        _py.normalize_scanorama,
        "medium",
        citation="Hie et al. 2019, Nat Biotechnol.",
        python_packages=("scanorama",),
        hyperparams={
            "dimred": HyperParam(
                "dimred", int, 50, "Dimensions of the integrated embedding.", minimum=2
            )
        },
    ),
    "13_fsmvn": _spec(
        "13_fsmvn",
        _py.normalize_fsmvn,
        "medium",
        citation="Feature-specific mean-variance normalization.",
        uses_reference_batch=True,
        hyperparams={"target_group": _TARGET_GROUP},
    ),
    "15_fsqn_py": _spec(
        "15_fsqn_py",
        _py.normalize_fsqn_py,
        "high",
        citation="Franks et al. 2018, Biostatistics (FSQN), Python implementation.",
        uses_reference_batch=True,
        hyperparams={"target_group": _TARGET_GROUP},
    ),
    "17_quantile": _spec(
        "17_quantile",
        _py.normalize_quantile,
        "high",
        citation="Bolstad et al. 2003, Bioinformatics.",
    ),
    "18_rank": _spec(
        "18_rank",
        _py.normalize_rank,
        "high",
        citation="Fractional rank transform; platform-independent by construction.",
    ),
    "25_angel": _spec(
        "25_angel",
        _py.normalize_angel,
        "high",
        citation="Angel et al. 2020, PLOS Comput Biol.",
        notes=(
            "Reduces the gene space rather than correcting values, so its metrics are "
            "not directly comparable with the other methods."
        ),
        hyperparams={
            "threshold": HyperParam(
                "threshold",
                float,
                0.20,
                "Drop genes whose batch variance fraction reaches this. "
                "0.05 is strict, 0.50 lenient.",
                minimum=0.0,
                maximum=1.0,
            )
        },
    ),
    "26_xpn": _spec(
        "26_xpn",
        _py.normalize_xpn,
        "high",
        citation="Shabalin et al. 2008, Bioinformatics (XPN).",
        uses_reference_batch=True,
        hyperparams={
            "target_group": _TARGET_GROUP,
            "n_quantiles": HyperParam(
                "n_quantiles",
                int,
                50,
                "Interior quantile breakpoints for the piecewise-linear map.",
                minimum=1,
            ),
        },
    ),
    "30_recombat": _spec(
        "30_recombat",
        _py.normalize_recombat,
        "medium",
        citation="Adossa et al. 2021 (reComBat).",
        python_packages=("reComBat",),
    ),
    # ── R, in-memory ─────────────────────────────────────────────────────────────
    "03_limma": _spec(
        "03_limma",
        _r.normalize_limma,
        "low",
        citation="Ritchie et al. 2015, NAR (limma::removeBatchEffect).",
        requires_r=True,
        r_packages=("limma",),
    ),
    "04_sva": _spec(
        "04_sva",
        _r.normalize_sva,
        "low",
        citation="Leek & Storey 2012, Bioinformatics (SVA).",
        requires_r=True,
        r_packages=("sva", "limma"),
        uses_bio_col=True,
    ),
    "05_combat": _spec(
        "05_combat",
        _r.normalize_combat,
        "medium",
        citation="Johnson et al. 2007, Biostatistics (ComBat).",
        requires_r=True,
        r_packages=("sva",),
        uses_bio_col=True,
    ),
    "06_combat_seq": _spec(
        "06_combat_seq",
        _r.normalize_combat_seq,
        "medium",
        citation="Zhang et al. 2020, NAR Genom Bioinform (ComBat-seq).",
        requires_r=True,
        r_packages=("sva",),
        uses_bio_col=True,
    ),
    "09_ruv": _spec(
        "09_ruv",
        _r.normalize_ruv,
        "medium",
        citation="Risso et al. 2014, Nat Biotechnol (RUVg).",
        requires_r=True,
        r_packages=("RUVSeq",),
        hyperparams={
            "k": HyperParam(
                "k", int, 2, "Factors of unwanted variation to remove.", minimum=1
            ),
            "control_genes": HyperParam(
                "control_genes",
                tuple,
                DEFAULT_CONTROL_GENES,
                "Negative-control genes. The default is a human housekeeping set; "
                "supply your own for other organisms.",
            ),
        },
    ),
    "10_mnn": _spec(
        "10_mnn",
        _r.normalize_mnn,
        "medium",
        citation="Haghverdi et al. 2018, Nat Biotechnol (fastMNN).",
        requires_r=True,
        r_packages=("batchelor",),
        hyperparams={
            "k": HyperParam(
                "k", int, 20, "Mutual nearest neighbours per batch pair.", minimum=1
            )
        },
    ),
    "14_qsmooth": _spec(
        "14_qsmooth",
        _r.normalize_qsmooth,
        "high",
        citation="Hicks et al. 2018, Biostatistics (qsmooth).",
        requires_r=True,
        r_packages=("qsmooth",),
    ),
    "22_tmm": _spec(
        "22_tmm",
        _r.normalize_tmm,
        "low",
        citation="Robinson & Oshlack 2010, Genome Biol (TMM via edgeR).",
        requires_r=True,
        r_packages=("edgeR",),
        rnaseq_only=True,
        hyperparams={
            "method": HyperParam(
                "method",
                str,
                "TMM",
                "edgeR normalization method.",
                choices=("TMM", "TMMwsp", "RLE", "upperquartile", "none"),
                r_argument="calcNormFactors(method=)",
            ),
            "prior_count": HyperParam(
                "prior_count",
                float,
                1.0,
                "Pseudocount added before the log-CPM transform.",
                minimum=0.0,
                r_argument="cpm(prior.count=)",
            ),
        },
    ),
    "23_vst": _spec(
        "23_vst",
        _r.normalize_vst,
        "low",
        citation="Love et al. 2014, Genome Biol (DESeq2 VST).",
        requires_r=True,
        r_packages=("DESeq2",),
        rnaseq_only=True,
        hyperparams={
            "blind": HyperParam(
                "blind",
                bool,
                True,
                "Ignore the design when estimating dispersions.",
                r_argument="vst(blind=)",
            ),
            "min_genes_for_vst": HyperParam(
                "min_genes_for_vst",
                int,
                1000,
                "Below this gene count, use varianceStabilizingTransformation() "
                "instead of vst(), whose subsetting step needs more genes.",
                minimum=1,
            ),
        },
    ),
    # ── R, file round-trip ───────────────────────────────────────────────────────
    "16_fsqn_r": _spec(
        "16_fsqn_r",
        _r.normalize_fsqn_r,
        "high",
        citation="Franks et al. 2018, Biostatistics (FSQN R package).",
        requires_r=True,
        r_packages=("FSQN",),
        uses_reference_batch=True,
        hyperparams={"target_group": _TARGET_GROUP},
    ),
    "19_tdm": _spec(
        "19_tdm",
        _r.normalize_tdm,
        "high",
        citation="Thompson et al. 2016, BMC Bioinformatics (TDM).",
        requires_r=True,
        # binr is TDM's undeclared run-time dependency. TDM's own `load_it()` helper
        # calls BiocManager::install() for anything missing *at call time*, which
        # prompts, reads EOF in a container and dies with "argument is of length zero".
        # Declaring it here installs it at build time and keeps the method offline.
        r_packages=("TDM", "binr"),
        uses_reference_batch=True,
        hyperparams={
            "target_group": _TARGET_GROUP,
            "log_target": HyperParam(
                "log_target",
                bool,
                False,
                "Whether TDM should log-transform the reference distribution.",
                r_argument="tdm_transform(log_target=)",
            ),
        },
    ),
    "20_shambhala": _spec(
        "20_shambhala",
        _shambhala.normalize_shambhala,
        "high",
        citation="Borisov et al. 2021, Bioinformatics (Shambhala2); CuBlock.",
        requires_octave=True,
        hyperparams={
            "P": HyperParam(
                "P",
                str,
                None,
                "Calibration reference P, the quantile-normalization anchor. Path or "
                "s3:// URI; defaults to the shipped P0_standard.",
            ),
            "Q": HyperParam(
                "Q",
                str,
                None,
                "Calibration reference Q, whose per-gene log mean and std define the "
                "target shape. Path or s3:// URI; defaults to the shipped Q0_standard.",
            ),
            "k": HyperParam(
                "k", int, 5, "k-means gene clusters for CuBlock.", minimum=2
            ),
            "n_workers": HyperParam(
                "n_workers",
                int,
                1,
                "Shambhala worker processes. Defaults to 1: this multiplies with the "
                "dispatcher's own workers to give the total Octave process count.",
                minimum=1,
                maximum=10,
            ),
            "na_strategy": HyperParam(
                "na_strategy",
                str,
                "drop",
                "NaN handling before Octave, which cannot accept missing values. "
                "'drop' restores affected genes as NaN afterwards; 'knn' imputes them.",
                choices=("drop", "knn"),
            ),
            "knn_k": HyperParam(
                "knn_k", int, 5, "Neighbours for na_strategy='knn'.", minimum=1
            ),
            "max_na_frac": HyperParam(
                "max_na_frac",
                float,
                0.20,
                "Genes above this NaN fraction are dropped rather than imputed, for "
                "na_strategy='knn'.",
                minimum=0.0,
                maximum=1.0,
            ),
            "q_pseudocount": HyperParam(
                "q_pseudocount",
                float,
                1e-6,
                "Added to Q before the log so zero counts do not become -Inf. 0 "
                "excludes zero-count genes instead.",
                minimum=0.0,
            ),
            "random_seed": HyperParam(
                "random_seed",
                int,
                None,
                "Seeds Octave's k-means for reproducible CuBlock output.",
            ),
            "octave_bin": HyperParam(
                "octave_bin",
                str,
                "octave",
                "Octave executable. The image places it on PATH.",
            ),
            "timeout_s": HyperParam(
                "timeout_s",
                int,
                6000,
                "Per-batch Octave subprocess timeout, in seconds.",
                minimum=1,
            ),
            "precompute_qn_reference": HyperParam(
                "precompute_qn_reference",
                bool,
                False,
                "Quantile-normalize Python-side against a reference distribution taken "
                "from P and skip Octave's own pass. Large speed-up, near-lossless.",
            ),
            "synthetic_cublock_p": HyperParam(
                "synthetic_cublock_p",
                bool,
                False,
                "Give CuBlock one synthetic P column instead of all P samples.",
            ),
            "max_p_samples": HyperParam(
                "max_p_samples",
                int,
                None,
                "Subsample P to its N centroid-closest samples.",
                minimum=1,
            ),
            "precompute_cublock_clusters": HyperParam(
                "precompute_cublock_clusters",
                bool,
                False,
                "Cluster genes once from P and reuse the assignment for every sample, "
                "instead of 30 k-means runs per sample.",
            ),
            "python_cublock": HyperParam(
                "python_cublock",
                bool,
                False,
                "Experimental pure-Python CuBlock. Upstream measures ~226% mean "
                "relative difference against Octave; not equivalent.",
            ),
        },
        notes=(
            "Replaces the donor's Shambhala2 R package with the containerized pure-"
            "Python + Octave implementation, so results are not bit-identical to the "
            "published benchmark rows. Samples are harmonized independently, so no "
            "batch column is used. Output genes are the input's intersection with P "
            "and Q."
        ),
    ),
    "21_harmonizr": _spec(
        "21_harmonizr",
        _r.normalize_harmonizr,
        "medium",
        citation="Voss et al. 2022, Bioinformatics (HarmonizR).",
        requires_r=True,
        r_packages=("HarmonizR",),
        hyperparams={
            "algorithm": HyperParam(
                "algorithm",
                str,
                "ComBat",
                "Correction applied within each NA-consistent block.",
                choices=("ComBat", "limma"),
                r_argument="harmonizR(algorithm=)",
            ),
            "combat_mode": HyperParam(
                "combat_mode",
                int,
                2,
                "HarmonizR's ComBat parameter set: 1 par.prior/scale, 2 par.prior/"
                "mean-only, 3 non-parametric/scale, 4 non-parametric/mean-only. "
                "Modes 1 and 3 write an empty result on data without missing values, "
                "which is why the default is 2 rather than HarmonizR's own 1.",
                choices=(1, 2, 3, 4),
                r_argument="harmonizR(ComBat_mode=)",
            ),
        },
    ),
    "27_dwd": _spec(
        "27_dwd",
        _r.normalize_dwd,
        "medium",
        citation="Qing & Marron 2018 (DWDLargeR).",
        requires_r=True,
        r_packages=("DWDLargeR",),
        uses_reference_batch=True,
        hyperparams={
            "target_group": _TARGET_GROUP,
            "min_batch_size": HyperParam(
                "min_batch_size",
                int,
                5,
                "Batches smaller than this are left uncorrected; the DWD direction is "
                "unstable below it.",
                minimum=2,
            ),
            "expon": HyperParam(
                "expon",
                float,
                1.0,
                "Exponent of the DWD generalized distance. 1.0 is standard DWD; "
                "DWDLargeR requires it explicitly, with no default of its own.",
                minimum=0.1,
            ),
        },
    ),
    "28_npn": _spec(
        "28_npn",
        _r.normalize_npn,
        "high",
        citation="Liu et al. 2009, JMLR (nonparanormal), via huge.",
        requires_r=True,
        r_packages=("huge",),
        hyperparams={
            "npn_func": HyperParam(
                "npn_func",
                str,
                "truncation",
                "Nonparanormal transform variant.",
                choices=("truncation", "shrinkage", "skeptic"),
                r_argument="huge.npn(npn.func=)",
            )
        },
    ),
    "29_combat_ref": _spec(
        "29_combat_ref",
        _r.normalize_combat_ref,
        "medium",
        citation="M-ComBat: ComBat anchored to a reference batch (sva).",
        requires_r=True,
        r_packages=("sva",),
        uses_bio_col=True,
        uses_reference_batch=True,
        hyperparams={"target_group": _TARGET_GROUP},
    ),
    "31_ruv3prps": _spec(
        "31_ruv3prps",
        _r.normalize_ruv3prps,
        "medium",
        citation="Molania et al. 2022, NAR (RUV-III-PRPS).",
        requires_r=True,
        r_packages=("ruv",),
        uses_bio_col=True,
        hyperparams={
            "k_factors": HyperParam(
                "k_factors",
                int,
                5,
                "Unwanted-variation factors; capped at (pseudo-replicate cells - 1).",
                minimum=1,
                r_argument="RUVIII(k=)",
            ),
            "min_cell_size": HyperParam(
                "min_cell_size",
                int,
                2,
                "Minimum samples in a (biology x batch) cell to form a pseudo-replicate.",
                minimum=2,
            ),
        },
    ),
    "33_amdbnorm": _spec(
        "33_amdbnorm",
        _r.normalize_amdbnorm,
        "high",
        citation="AMDBNorm (PMID 34958674), on DBNorm.",
        requires_r=True,
        r_packages=("AMDBNorm", "DBNorm"),
        uses_reference_batch=True,
        hyperparams={
            "target_group": _TARGET_GROUP,
            "poly_degree": HyperParam(
                "poly_degree",
                int,
                9,
                "Degree of the polynomial fitted to the reference distribution.",
                minimum=1,
                r_argument="DBNorm::polyFit(degree)",
            ),
            "n_dist_points": HyperParam(
                "n_dist_points",
                int,
                500,
                "Points used to estimate the reference distribution.",
                minimum=10,
                r_argument="DBNorm::genDistData(n)",
            ),
        },
    ),
    "34_arsyn": _spec(
        "34_arsyn",
        _r.normalize_arsyn,
        "medium",
        citation="Nueda et al. 2012, Bioinformatics (ARSyN), via NOISeq.",
        requires_r=True,
        r_packages=("NOISeq",),
        uses_bio_col=True,
        hyperparams={
            "norm": HyperParam(
                "norm",
                str,
                "n",
                "NOISeq normalization applied before ARSyN.",
                choices=("n", "rpkm", "uqua", "tmm"),
                r_argument="ARSyNseq(norm=)",
            ),
            "logtransf": HyperParam(
                "logtransf",
                bool,
                False,
                "Whether ARSyN should log-transform first.",
                r_argument="ARSyNseq(logtransf=)",
            ),
        },
    ),
    "37_fabatch": _spec(
        "37_fabatch",
        _r.normalize_fabatch,
        "medium",
        citation="Hornung et al. 2016, BMC Bioinformatics (FAbatch), via bapred.",
        requires_r=True,
        r_packages=("bapred",),
        uses_bio_col=True,
        notes=(
            "bapred::fabatch structurally requires a two-level target, so this method "
            "skips honestly on datasets whose target column has any other number of "
            "levels."
        ),
        hyperparams={
            "target_col": HyperParam(
                "target_col",
                str,
                None,
                "Annotation column supplying the binary target. Defaults to the "
                "biology column; must have exactly two levels.",
            )
        },
    ),
    "38_harman": _spec(
        "38_harman",
        _r.normalize_harman,
        "medium",
        citation="Oytam et al. 2016, BMC Bioinformatics (Harman).",
        requires_r=True,
        r_packages=("Harman",),
        uses_bio_col=True,
        hyperparams={
            "limit": HyperParam(
                "limit",
                float,
                0.1,
                "Confidence limit constraining overcorrection; lower is more "
                "conservative.",
                minimum=0.0,
                maximum=1.0,
                r_argument="harman(limit=)",
            )
        },
    ),
}


def get_method(key: str) -> MethodSpec:
    """
    Look up a method, raising a helpful error for an unknown or excluded key.

    Raises
    ------
    KeyError
        If the key is unknown. An excluded key gets an error explaining *why* it is
        absent rather than a bare miss.
    """
    if key in METHOD_REGISTRY:
        return METHOD_REGISTRY[key]
    if key in EXCLUDED_METHODS:
        raise KeyError(
            f"{key!r} is deliberately not included in ComboBatch: "
            f"{EXCLUDED_METHODS[key]}. See docs/METHODS.md."
        )
    raise KeyError(
        f"unknown method {key!r}. Available: {', '.join(sorted(METHOD_REGISTRY))}"
    )


def available_methods() -> list[str]:
    """Return the keys whose backends are all present in this environment."""
    return [key for key, spec in METHOD_REGISTRY.items() if spec.is_available()]


__all__ = [
    "METHOD_REGISTRY",
    "EXCLUDED_METHODS",
    "HyperParam",
    "MethodSpec",
    "assert_rnaseq_only",
    "available_methods",
    "drop_na_genes",
    "get_method",
    "resolve_reference_batch",
]
