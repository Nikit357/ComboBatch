"""Run configuration: frozen dataclasses, a YAML loader, and validation.

Everything the donor pipeline hardcoded as a module constant — batch and biology column
names, the S3 bucket and prefix, the reference batch, the metric column lists — lives
here instead, so the tool has no opinion about anyone's dataset.

``RunConfig.validate()`` is meant to be called **before the first download**. The donor
discovered a mistyped strategy name only after pulling a 1.9 GB matrix.

Note on naming: :class:`MethodSelection` here is the *user's choice* of a method plus its
parameters. It is distinct from ``combobatch.methods.base.MethodSpec``, which is the
*registry entry* describing what a method is and which parameters it accepts.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any, Iterable, Mapping

from combobatch import params as params_mod

# Metric groups are single letters A-N. The registry in combobatch.metrics is the source
# of truth once it exists; this set lets configuration be validated before it is built.
METRIC_GROUP_LETTERS = frozenset("ABCDEFGHIJKLMN")

# A deliberately cheap default: the fast groups, no R, no permutation-heavy work.
DEFAULT_METRIC_GROUPS = frozenset("ABEJK")


class ConfigError(ValueError):
    """Raised when a configuration is malformed or inconsistent with the data."""


@dataclass(frozen=True)
class ColumnSpec:
    """Which annotation columns carry batch, biology and cohort identity."""

    batch: str
    bio: str
    cohort: str | None = None
    # Metric-side column selection. Empty means "just the primary column", which keeps a
    # minimal invocation minimal; naming several reproduces the donor's multi-column
    # behaviour without hardcoding its column names.
    metric_batch_cols: tuple[str, ...] = ()
    metric_bio_cols: tuple[str, ...] = ()
    predict_class_col: str | None = None

    def __post_init__(self) -> None:
        if not self.metric_batch_cols:
            object.__setattr__(self, "metric_batch_cols", (self.batch,))
        if not self.metric_bio_cols:
            object.__setattr__(self, "metric_bio_cols", (self.bio,))
        if self.predict_class_col is None:
            object.__setattr__(self, "predict_class_col", self.bio)

    def all_columns(self) -> tuple[str, ...]:
        """Every annotation column this configuration refers to, deduplicated."""
        names = [self.batch, self.bio, *self.metric_batch_cols, *self.metric_bio_cols]
        if self.cohort:
            names.append(self.cohort)
        if self.predict_class_col:
            names.append(self.predict_class_col)
        return tuple(dict.fromkeys(names))


@dataclass(frozen=True)
class ImputationSelection:
    """One imputation strategy, with its hyperparameters."""

    name: str
    params: dict[str, Any] = field(default_factory=dict)
    # None means "whatever this imputer declares", which is 0.0 for strict and 0.20 for
    # the rest. A shared literal default here would silently give strict a 0.20 ceiling
    # and stop it being strict.
    max_na_frac: float | None = None
    param_tag: str | None = None

    def tag(self) -> str | None:
        """Return the filename tag for this selection, or None when at defaults."""
        if self.param_tag is not None:
            return params_mod.validate_param_tag(self.param_tag)
        return params_mod.encode_param_tag(self.params)


@dataclass(frozen=True)
class MethodSelection:
    """One harmonization method, with its hyperparameters."""

    name: str
    params: dict[str, Any] = field(default_factory=dict)
    param_tag: str | None = None

    def tag(self) -> str | None:
        """Return the filename tag for this selection, or None when at defaults."""
        if self.param_tag is not None:
            return params_mod.validate_param_tag(self.param_tag)
        return params_mod.encode_param_tag(self.params)


@dataclass(frozen=True)
class MetricsSpec:
    """Which quality metrics to compute, and how hard to work at them."""

    enabled: bool = False
    groups: frozenset[str] = DEFAULT_METRIC_GROUPS
    skip_slow: bool = True
    skip_wm: bool = True
    n_perm: int = 20
    marker_panel: str | None = None
    predict_classes: tuple[str, ...] | None = None
    force_groups: frozenset[str] = frozenset()

    # Tuning knobs that were unreachable module constants upstream. The defaults
    # reproduce the donor's behaviour exactly; naming them here is what makes a
    # departure from it deliberate rather than a source edit.
    n_pcs: int = 10
    min_test_n: int = 20
    min_cohort_n: int = 20
    wm_permutations: int = 200
    wm_max_samples: int = 2000
    xb_max_samples: int = 8000
    dsc_permutations: int = 999
    asw_max_samples: int = 2000
    ks_genes: int = 1000
    collect_detail: bool = False


@dataclass(frozen=True)
class PostRemovalSpec:
    """Post-harmonization removal of the most PCA-deviant batches."""

    enabled: bool = False
    n_batches: int = 1
    min_batch_size: int = 20
    n_pcs: int = 2


@dataclass(frozen=True)
class RunConfig:
    """A complete, validated description of one ComboBatch run."""

    exp_uri: str
    ann_uri: str
    out_uri: str
    columns: ColumnSpec
    imputations: tuple[ImputationSelection, ...]
    methods: tuple[MethodSelection, ...]
    metrics: MetricsSpec = field(default_factory=MetricsSpec)
    post_removal: PostRemovalSpec = field(default_factory=PostRemovalSpec)
    reference_batch: str | None = None
    subset_query: str | None = None
    random_seed: int = 42
    endpoint_url: str | None = None
    work_dir: str | None = None
    run_id: str = "local"

    # ---------------------------------------------------------------- construction --

    @classmethod
    def from_mapping(cls, data: Mapping[str, Any]) -> "RunConfig":
        """
        Build a configuration from a parsed YAML mapping.

        Parameters
        ----------
        data
            The mapping, in the layout documented in ``examples/configs/``.

        Returns
        -------
        An unvalidated :class:`RunConfig`; call :meth:`validate` before using it.
        """
        try:
            source = data.get("input", {})
            columns_raw = data.get("columns", {})
            columns = ColumnSpec(
                batch=columns_raw["batch"],
                bio=columns_raw["bio"],
                cohort=columns_raw.get("cohort"),
                metric_batch_cols=tuple(columns_raw.get("metric_batch_cols", ())),
                metric_bio_cols=tuple(columns_raw.get("metric_bio_cols", ())),
                predict_class_col=columns_raw.get("predict_class_col"),
            )
        except KeyError as exc:
            raise ConfigError(f"missing required configuration key: {exc}") from exc

        imputations = tuple(
            ImputationSelection(
                name=entry["name"],
                params=dict(entry.get("params", {})),
                max_na_frac=(
                    float(entry["max_na_frac"])
                    if entry.get("max_na_frac") is not None
                    else None
                ),
                param_tag=entry.get("param_tag"),
            )
            for entry in _as_entries(data.get("imputations", []), "imputations")
        )
        methods = tuple(
            MethodSelection(
                name=entry["name"],
                params=dict(entry.get("params", {})),
                param_tag=entry.get("param_tag"),
            )
            for entry in _as_entries(data.get("methods", []), "methods")
        )

        metrics_raw = data.get("metrics", {}) or {}
        metrics = MetricsSpec(
            enabled=bool(metrics_raw.get("enabled", False)),
            groups=parse_metric_groups(
                metrics_raw.get("groups"), DEFAULT_METRIC_GROUPS
            ),
            skip_slow=bool(metrics_raw.get("skip_slow", True)),
            skip_wm=bool(metrics_raw.get("skip_wm", True)),
            n_perm=int(metrics_raw.get("n_perm", 20)),
            marker_panel=metrics_raw.get("marker_panel"),
            predict_classes=(
                tuple(metrics_raw["predict_classes"])
                if metrics_raw.get("predict_classes")
                else None
            ),
            force_groups=parse_metric_groups(
                metrics_raw.get("force_groups"), frozenset()
            ),
            **{
                name: int(metrics_raw[name])
                for name in (
                    "n_pcs",
                    "min_test_n",
                    "min_cohort_n",
                    "wm_permutations",
                    "wm_max_samples",
                    "xb_max_samples",
                    "dsc_permutations",
                    "asw_max_samples",
                    "ks_genes",
                )
                if metrics_raw.get(name) is not None
            },
            collect_detail=bool(metrics_raw.get("collect_detail", False)),
        )

        post_raw = data.get("post_removal", {}) or {}
        if isinstance(post_raw, bool):
            post_raw = {"enabled": post_raw}
        post_removal = PostRemovalSpec(
            enabled=bool(post_raw.get("enabled", False)),
            n_batches=int(post_raw.get("n_batches", 1)),
            min_batch_size=int(post_raw.get("min_batch_size", 20)),
            n_pcs=int(post_raw.get("n_pcs", 2)),
        )

        try:
            return cls(
                exp_uri=source["expression"],
                ann_uri=source["annotation"],
                out_uri=data["output"],
                columns=columns,
                imputations=imputations,
                methods=methods,
                metrics=metrics,
                post_removal=post_removal,
                reference_batch=data.get("reference_batch"),
                subset_query=data.get("subset_query"),
                random_seed=int(data.get("random_seed", 42)),
                endpoint_url=data.get("endpoint_url"),
                work_dir=data.get("work_dir"),
                run_id=str(data.get("run_id", "local")),
            )
        except KeyError as exc:
            raise ConfigError(f"missing required configuration key: {exc}") from exc

    @classmethod
    def from_yaml(cls, path: str) -> "RunConfig":
        """Load a configuration from a YAML file."""
        import yaml

        with open(path) as handle:
            data = yaml.safe_load(handle)
        if not isinstance(data, Mapping):
            raise ConfigError(f"{path}: expected a YAML mapping at the top level")
        return cls.from_mapping(data)

    def with_overrides(self, **overrides: Any) -> "RunConfig":
        """
        Return a copy with the given fields replaced, ignoring ``None`` values.

        This is how CLI flags overlay a YAML file: an unset flag arrives as ``None`` and
        must leave the file's value alone.
        """
        applied = {key: value for key, value in overrides.items() if value is not None}
        return replace(self, **applied) if applied else self

    # ------------------------------------------------------------------ validation --

    def validate(
        self,
        ann_df: Any = None,
        *,
        known_methods: Iterable[str] | None = None,
        known_imputers: Iterable[str] | None = None,
        known_groups: Iterable[str] | None = None,
    ) -> None:
        """
        Validate the configuration, as far as the supplied information allows.

        Static checks always run. Registry checks run when the corresponding name set is
        supplied, and data checks run when ``ann_df`` is supplied — so the CLI can
        validate names at startup and columns immediately after the annotation loads,
        both before any expensive work.

        Parameters
        ----------
        ann_df
            The loaded annotation, for column and reference-batch checks.
        known_methods, known_imputers, known_groups
            Valid names from the registries, when available.

        Raises
        ------
        ConfigError
            On the first problem found, naming it specifically.
        """
        self._validate_static()
        self._validate_registries(known_methods, known_imputers, known_groups)
        if ann_df is not None:
            self._validate_against_annotation(ann_df)

    def _validate_static(self) -> None:
        """Checks that need neither the registries nor the data."""
        if not self.imputations:
            raise ConfigError("no imputation strategies selected")
        if not self.methods:
            raise ConfigError("no harmonization methods selected")

        # Keyed on (name, tag), not name alone: listing a method twice at different
        # hyperparameters is the sweep this tool exists to make possible, and the two
        # runs write to different files. Only selections that would land on the same
        # output key are a genuine duplicate.
        duplicates = _duplicates(
            [
                f"{sel.name} ({sel.tag()})" if sel.tag() else sel.name
                for sel in self.methods
            ]
        )
        if duplicates:
            raise ConfigError(
                f"method selection(s) listed more than once: {', '.join(duplicates)}. "
                f"Two entries for one method need different parameters, or a "
                f"distinguishing param_tag."
            )
        duplicates = _duplicates(
            [
                f"{sel.name} ({sel.tag()})" if sel.tag() else sel.name
                for sel in self.imputations
            ]
        )
        if duplicates:
            raise ConfigError(
                f"imputer selection(s) listed more than once: {', '.join(duplicates)}. "
                f"Two entries for one imputer need different parameters, or a "
                f"distinguishing param_tag."
            )

        for selection in self.imputations:
            if selection.max_na_frac is None:
                continue
            if not 0.0 <= selection.max_na_frac <= 1.0:
                raise ConfigError(
                    f"{selection.name}: max_na_frac must be between 0 and 1, "
                    f"got {selection.max_na_frac}"
                )

        bad_letters = sorted(set(self.metrics.groups) - METRIC_GROUP_LETTERS)
        if bad_letters:
            raise ConfigError(
                f"unknown metric group(s): {', '.join(bad_letters)}. "
                f"Valid groups are {''.join(sorted(METRIC_GROUP_LETTERS))}"
            )
        if self.metrics.n_perm < 0:
            raise ConfigError(f"n_perm must be non-negative, got {self.metrics.n_perm}")

        if self.post_removal.enabled and self.post_removal.n_batches < 1:
            raise ConfigError(
                "post_removal.n_batches must be at least 1 when post-removal is enabled"
            )
        if self.post_removal.n_pcs < 1:
            raise ConfigError("post_removal.n_pcs must be at least 1")

        # Tag construction validates any user-supplied --param-tag as a side effect.
        for selection in (*self.imputations, *self.methods):
            selection.tag()

    def _validate_registries(
        self,
        known_methods: Iterable[str] | None,
        known_imputers: Iterable[str] | None,
        known_groups: Iterable[str] | None,
    ) -> None:
        """Checks against the registries, when they are available."""
        if known_methods is not None:
            _reject_unknown(
                "method", [sel.name for sel in self.methods], set(known_methods)
            )
        if known_imputers is not None:
            _reject_unknown(
                "imputer", [sel.name for sel in self.imputations], set(known_imputers)
            )
        if known_groups is not None:
            _reject_unknown(
                "metric group", sorted(self.metrics.groups), set(known_groups)
            )

    def _validate_against_annotation(self, ann_df: Any) -> None:
        """Checks that need the annotation table."""
        available = set(ann_df.columns)
        missing = [name for name in self.columns.all_columns() if name not in available]
        if missing:
            import difflib

            details = []
            for name in missing:
                close = difflib.get_close_matches(name, sorted(available), n=1)
                hint = f" (did you mean {close[0]!r}?)" if close else ""
                details.append(f"{name!r}{hint}")
            raise ConfigError(
                f"annotation is missing column(s): {', '.join(details)}. "
                f"It has {len(available)} columns."
            )

        if self.reference_batch is not None:
            batches = set(ann_df[self.columns.batch].astype(str))
            if self.reference_batch not in batches:
                import difflib

                close = difflib.get_close_matches(
                    self.reference_batch, sorted(batches), n=1
                )
                hint = f" Did you mean {close[0]!r}?" if close else ""
                raise ConfigError(
                    f"reference_batch {self.reference_batch!r} is not a value of "
                    f"{self.columns.batch!r} ({len(batches)} distinct values).{hint}"
                )

        if self.subset_query is not None:
            from combobatch.subset import SubsetQueryError, validate_subset_query

            try:
                validate_subset_query(self.subset_query, ann_df)
            except SubsetQueryError as exc:
                raise ConfigError(str(exc)) from exc

    # ------------------------------------------------------------------- iteration --

    def combinations(self) -> list[tuple[ImputationSelection, MethodSelection]]:
        """Return the full imputation x method cross-product, in declaration order."""
        return [(imp, method) for imp in self.imputations for method in self.methods]


def _as_entries(raw: Any, field_name: str) -> list[dict[str, Any]]:
    """Normalise a list of names or mappings into a list of mappings."""
    if raw is None:
        return []
    if not isinstance(raw, list):
        raise ConfigError(f"{field_name}: expected a list, got {type(raw).__name__}")
    entries: list[dict[str, Any]] = []
    for item in raw:
        if isinstance(item, str):
            entries.append({"name": item})
        elif isinstance(item, Mapping):
            if "name" not in item:
                raise ConfigError(f"{field_name}: entry {item!r} has no 'name'")
            entries.append(dict(item))
        else:
            raise ConfigError(f"{field_name}: unexpected entry {item!r}")
    return entries


def parse_metric_groups(raw: Any, default: frozenset[str]) -> frozenset[str]:
    """Normalise a group specification into a set of upper-case single letters."""
    if raw is None:
        return default
    if isinstance(raw, str):
        letters = [part.strip() for part in raw.replace(",", "").strip()]
    else:
        letters = [str(part).strip() for part in raw]
    return frozenset(letter.upper() for letter in letters if letter)


def _duplicates(names: list[str]) -> list[str]:
    """Return names appearing more than once, in first-seen order."""
    seen: dict[str, int] = {}
    for name in names:
        seen[name] = seen.get(name, 0) + 1
    return [name for name, count in seen.items() if count > 1]


def _reject_unknown(kind: str, requested: list[str], known: set[str]) -> None:
    """Raise a ConfigError naming any requested value that is not in ``known``."""
    unknown = [name for name in requested if name not in known]
    if not unknown:
        return

    import difflib

    details = []
    for name in unknown:
        close = difflib.get_close_matches(name, sorted(known), n=1)
        hint = f" (did you mean {close[0]!r}?)" if close else ""
        details.append(f"{name!r}{hint}")
    raise ConfigError(
        f"unknown {kind}(s): {', '.join(details)}. "
        f"Available: {', '.join(sorted(known))}"
    )
