"""``METRIC_GROUP_REGISTRY`` — the third registry, and the only one written from scratch.

Upstream has no equivalent. Group dispatch there is a hardcoded sequence of
``_run_group("A", ...)`` calls, the sentinel keys live in a **different file** as a
parallel dict, and "which groups need the reference matrix" is a third constant in a
third place. Adding a group meant editing all three and forgetting one meant a group that
computed but never recorded, or recorded but never recomputed.

Here one table drives everything: ``--groups`` validation, ``--skip-slow``, the
incremental sentinels, which embeddings to compute, ``combobatch list-metrics``, and the
generated ``docs/METRICS.md``.

``column_roles`` is what makes the sentinels dataset-agnostic. A group declares that it
keys off *a batch column*; the concrete name arrives from the user's ``ColumnSpec``, so
the same group emits ``r2_site`` on one dataset and ``r2_sequencing_run`` on another, and
its sentinel follows.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable

from combobatch.metrics.groups import MetricColumns

# Every group letter, in the order the runner executes them: cheap and structural first,
# so a run that dies late still carries the descriptive groups.
GROUP_LETTERS = "EKAJBCDGHILMNF"


@dataclass(frozen=True)
class MetricGroupSpec:
    """Everything ComboBatch knows about one metric group."""

    letter: str
    name: str
    fn: Callable[..., dict]
    sentinel_template: str
    needs_annotation: bool = True
    column_roles: tuple[str, ...] = ()
    needs_reference: bool = False
    needs_panel: bool = False
    needs_embedding: tuple[str, ...] = ()
    needs_r: bool = False
    speed: str = "fast"
    default_active: bool = True
    description: str = ""
    # Which direction is good. Several groups emit both batch-keyed and biology-keyed
    # metrics whose polarities are opposite, hence "mixed" - reading one of those without
    # knowing which column it keyed off is how a metric gets quoted backwards.
    polarity: str = "mixed"
    citation: str = ""

    def sentinel_key(self, columns: MetricColumns) -> str:
        """
        Return the key whose presence means this group already ran.

        Resolved against the configured columns, so the sentinel names the column the
        group actually keyed off rather than a hardcoded one.
        """
        return self.sentinel_template.format(
            batch=columns.batch, bio=columns.bio, cohort=columns.cohort
        )

    def missing_backends(self) -> list[str]:
        """Return backends this group needs but cannot find here."""
        missing: list[str] = []
        if self.needs_r and not _rscript_available():
            missing.append("R (Rscript)")
        for package in self.needs_embedding:
            if package == "umap" and not _module_available("umap"):
                missing.append("Python package umap-learn")
        return missing

    def is_available(self) -> bool:
        """Report whether every backend this group needs is present."""
        return not self.missing_backends()


def _module_available(name: str) -> bool:
    """Report whether an optional module can be imported."""
    import importlib.util

    try:
        return importlib.util.find_spec(name) is not None
    except (ImportError, ValueError):
        return False


def _rscript_available() -> bool:
    """Report whether Rscript is on PATH. Group F shells out rather than using rpy2."""
    import shutil

    return shutil.which("Rscript") is not None


def _spec(letter: str, name: str, fn, sentinel: str, **kwargs) -> MetricGroupSpec:
    """Build a MetricGroupSpec with the common defaults."""
    return MetricGroupSpec(
        letter=letter, name=name, fn=fn, sentinel_template=sentinel, **kwargs
    )


def _build_registry() -> dict[str, MetricGroupSpec]:
    """Assemble the registry, importing the group implementations lazily."""
    from combobatch.metrics import groups as g

    return {
        "A": _spec(
            "A",
            "PCA variance decomposition",
            g.compute_group_a,
            "r2_{batch}",
            column_roles=("batch", "bio", "cohort"),
            needs_embedding=("pca",),
            speed="moderate",
            description=(
                "How much of each principal component's variance each annotation column "
                "explains, plus a permutation-tested separability score."
            ),
            polarity="lower is better for batch, higher for biology",
            citation="Buttner et al. 2019, Nat Methods (PCR/PCReg); Teschendorff & Zhu 2011 (DSC).",
        ),
        "B": _spec(
            "B",
            "Neighbourhood integration",
            g.compute_group_b,
            "kbet_acceptance_rate_{batch}",
            column_roles=("batch", "bio"),
            needs_embedding=("pca",),
            speed="slow",
            description=(
                "kBET, iLISI/cLISI, silhouettes and cell mixing: whether batches blend "
                "locally while biology stays apart."
            ),
            polarity="higher is better for batch mixing, higher for cLISI",
            citation="Buttner et al. 2019, Nat Methods (kBET); Korsunsky et al. 2019, Nat Methods (LISI).",
        ),
        "C": _spec(
            "C",
            "Embedding structure",
            g.compute_group_c,
            "umap_centroid_disp_{batch}",
            column_roles=("batch",),
            needs_embedding=("pca", "umap", "tsne"),
            speed="slow",
            default_active=False,
            description=(
                "Centroid dispersion and local entropy of batch labels in UMAP and "
                "t-SNE space."
            ),
            polarity="lower dispersion and higher entropy are better",
            citation="McInnes et al. 2018, arXiv (UMAP); van der Maaten & Hinton 2008, JMLR (t-SNE).",
        ),
        "D": _spec(
            "D",
            "Distribution comparison",
            g.compute_group_d,
            "ks_mean_D_{batch}",
            column_roles=("batch", "cohort"),
            speed="slow",
            default_active=False,
            description=(
                "Pairwise Kolmogorov-Smirnov distances between groups, cohort effects "
                "within a batch, and per-gene coefficient of variation."
            ),
            polarity="lower is better",
            citation="Massey 1951, JASA (Kolmogorov-Smirnov).",
        ),
        "E": _spec(
            "E",
            "Data quality",
            g.compute_group_e,
            "n_samples",
            column_roles=("batch", "bio", "cohort"),
            speed="fast",
            description=(
                "Shape, batch balance, zero inflation, expression percentiles and "
                "per-cohort bimodality."
            ),
            polarity="descriptive; no direction",
            citation="Descriptive statistics; no single source.",
        ),
        "F": _spec(
            "F",
            "Variance partition",
            g.compute_group_f,
            "vp_median_{batch}",
            column_roles=("batch", "bio", "cohort"),
            needs_r=True,
            speed="very_slow",
            default_active=False,
            description=(
                "Per-gene variance attributed to each annotation column by a mixed "
                "model, via R's variancePartition."
            ),
            polarity="lower for batch, higher for biology",
            citation="Hoffman & Schadt 2016, BMC Bioinformatics (variancePartition).",
        ),
        "G": _spec(
            "G",
            "Graph connectivity",
            g.compute_group_g,
            "graph_connectivity_{bio}",
            column_roles=("bio",),
            needs_embedding=("pca",),
            speed="moderate",
            default_active=False,
            description=(
                "Whether each biology group forms a single connected kNN component "
                "rather than one island per batch."
            ),
            polarity="higher is better",
            citation="Luecken et al. 2022, Nat Methods (scIB graph connectivity).",
        ),
        "H": _spec(
            "H",
            "Pairwise distances",
            g.compute_group_h,
            "dist_ratio_{batch}",
            column_roles=("batch", "bio", "cohort"),
            speed="moderate",
            default_active=False,
            description=(
                "Average within-group against between-group Euclidean distance in full "
                "gene space."
            ),
            polarity="lower ratio for batch, higher for biology",
            citation="Standard cluster-separation statistic; see Luecken et al. 2022.",
        ),
        "I": _spec(
            "I",
            "WaterMelon score",
            g.compute_group_i,
            "wm_{batch}",
            column_roles=("batch", "bio"),
            needs_embedding=("pca",),
            speed="very_slow",
            default_active=False,
            description=(
                "Information a Ward dendrogram recovers about each label, against a "
                "permutation null."
            ),
            polarity="higher is better for biology",
            citation="Ward 1963, JASA; permutation null as in Luecken et al. 2022.",
        ),
        "J": _spec(
            "J",
            "Per-PC variance",
            g.compute_group_j,
            "pct_var_pc1",
            needs_annotation=False,
            needs_embedding=("pca",),
            speed="fast",
            description="Share of total variance carried by each leading component.",
            polarity="descriptive; no direction",
            citation="Standard PCA scree diagnostics.",
        ),
        "K": _spec(
            "K",
            "NA retention",
            g.compute_group_k,
            "n_genes_noNA",
            needs_annotation=False,
            speed="fast",
            description=(
                "Genes and samples free of missing values after harmonization — the "
                "cheapest way to catch a method that silently destroyed the matrix."
            ),
            polarity="higher is better",
            citation="ComboBatch-specific integrity check; no external source.",
        ),
        "L": _spec(
            "L",
            "Profile preservation",
            g.compute_group_l,
            "mk_rho_mean_all_genes",
            column_roles=("cohort",),
            needs_reference=True,
            needs_panel=True,
            speed="moderate",
            default_active=False,
            description=(
                "Within-cohort Spearman correlation of each gene before and after "
                "harmonization. Saturates at ~1.0 for any per-batch monotone transform, "
                "so read it as a guard rail, not a ranking."
            ),
            polarity="higher is better, but saturates",
            citation="ComboBatch-specific biology-preservation guard rail; see docs/METRICS.md.",
        ),
        "M": _spec(
            "M",
            "Cross-batch rank agreement",
            g.compute_group_m,
            "xb_rank_agree",
            column_roles=("batch", "bio"),
            needs_panel=True,
            speed="moderate",
            default_active=False,
            description=(
                "Whether same-biology samples agree across batches more than "
                "different-biology ones. The margin is the discriminative signal."
            ),
            polarity="higher margin is better",
            citation="ComboBatch-specific cross-batch discrimination score.",
        ),
        "N": _spec(
            "N",
            "Predictive validation",
            g.compute_group_n,
            "pv_lobo3_f1_macro_mean",
            column_roles=("batch", "bio"),
            speed="very_slow",
            default_active=False,
            description=(
                "Leave-one-batch-out classification of biology, against a label "
                "permutation control."
            ),
            polarity="higher is better, above the permutation null",
            citation="Leave-one-batch-out validation as in Luecken et al. 2022, Nat Methods.",
        ),
    }


METRIC_GROUP_REGISTRY: dict[str, MetricGroupSpec] = _build_registry()

# The default set: everything fast plus the two that catch outright failure. Deliberately
# cheap, so `--metrics` on a large cross-product does not quietly dominate the runtime.
DEFAULT_GROUPS: frozenset[str] = frozenset(
    letter for letter, spec in METRIC_GROUP_REGISTRY.items() if spec.default_active
)


def get_group(letter: str) -> MetricGroupSpec:
    """
    Look up one group by letter.

    Raises
    ------
    KeyError
        Naming the letter and listing the valid ones.
    """
    key = letter.strip().upper()
    if key not in METRIC_GROUP_REGISTRY:
        raise KeyError(
            f"unknown metric group {letter!r}; valid groups are "
            f"{''.join(sorted(METRIC_GROUP_REGISTRY))}"
        )
    return METRIC_GROUP_REGISTRY[key]


def ordered_groups(letters: frozenset[str] | set[str]) -> list[MetricGroupSpec]:
    """Return the requested groups in execution order, cheapest and most structural first."""
    return [
        METRIC_GROUP_REGISTRY[letter]
        for letter in GROUP_LETTERS
        if letter in letters and letter in METRIC_GROUP_REGISTRY
    ]


def required_embeddings(letters: frozenset[str] | set[str]) -> set[str]:
    """Return the union of embeddings the requested groups need."""
    needed: set[str] = set()
    for letter in letters:
        spec = METRIC_GROUP_REGISTRY.get(letter)
        if spec is not None:
            needed.update(spec.needs_embedding)
    return needed


__all__ = [
    "DEFAULT_GROUPS",
    "GROUP_LETTERS",
    "METRIC_GROUP_REGISTRY",
    "MetricColumns",
    "MetricGroupSpec",
    "get_group",
    "ordered_groups",
    "required_embeddings",
]
