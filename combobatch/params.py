"""Hyperparameter parsing and the output-key grammar.

Two responsibilities that are deliberately kept together, because they are two halves of
the same problem: once hyperparameters are tunable, ``{imp}__{method}`` stops identifying
a result uniquely. Every non-default parameter set therefore gets a tag encoded into the
output filename, and ``run_manifest.json`` maps that tag back to the full resolved dict.

The donor pipeline split output keys on ``"__"`` and required exactly four parts, silently
dropping anything else from job enumeration. Here an unparseable key raises.
"""

from __future__ import annotations

import hashlib
import json
import re
from dataclasses import dataclass

# The tag is a filename component, so "__" must not survive into it: it is the field
# separator of the key grammar, and an unescaped one would make keys ambiguous.
KEY_SEPARATOR = "__"
MAX_PARAM_TAG_LENGTH = 64
_HASH_LENGTH = 8

_POST_RE = re.compile(r"^post([01])$")
_TRUE_LITERALS = frozenset({"true", "yes", "t", "y", "1"})
_FALSE_LITERALS = frozenset({"false", "no", "f", "n", "0"})


class ParamParseError(ValueError):
    """Raised when a ``--method-params`` / ``--impute-params`` string is malformed."""


class KeyParseError(ValueError):
    """Raised when an output key does not match the documented grammar."""


# --------------------------------------------------------------------------------------
# Value coercion
# --------------------------------------------------------------------------------------


def coerce_value(raw: str) -> object:
    """
    Convert a command-line parameter value to int, float, bool, list or str.

    A ``+``-joined value becomes a list, mirroring how lists are encoded into tags.
    Coercion is by inspection; when a declared type is available the caller should use
    :func:`coerce_with_spec` instead.

    Parameters
    ----------
    raw
        The literal text to the right of ``=``.

    Returns
    -------
    The coerced value.
    """
    text = raw.strip()
    if "+" in text:
        return [coerce_value(part) for part in text.split("+")]

    lowered = text.lower()
    if lowered in _TRUE_LITERALS and lowered not in {"1"}:
        return True
    if lowered in _FALSE_LITERALS and lowered not in {"0"}:
        return False
    if lowered in {"none", "null"}:
        return None

    for caster in (int, float):
        try:
            return caster(text)
        except ValueError:
            continue
    return text


def coerce_with_spec(name: str, raw: str, declared_type: type | None) -> object:
    """
    Coerce a value to a declared type, falling back to inspection when none is given.

    Parameters
    ----------
    name
        Parameter name, used only in the error message.
    raw
        The literal text to the right of ``=``.
    declared_type
        The type from the method's or imputer's ``HyperParam`` declaration.

    Returns
    -------
    The coerced value.

    Raises
    ------
    ParamParseError
        If the value cannot be represented as ``declared_type``.
    """
    if declared_type is None:
        return coerce_value(raw)
    if declared_type is bool:
        lowered = raw.strip().lower()
        if lowered in _TRUE_LITERALS:
            return True
        if lowered in _FALSE_LITERALS:
            return False
        raise ParamParseError(f"{name}: expected a boolean, got {raw!r}")
    try:
        return declared_type(raw)
    except (TypeError, ValueError) as exc:
        type_name = getattr(declared_type, "__name__", str(declared_type))
        raise ParamParseError(f"{name}: expected {type_name}, got {raw!r}") from exc


# --------------------------------------------------------------------------------------
# CLI parsing
# --------------------------------------------------------------------------------------


def parse_params_arg(raw: str) -> dict[str, dict[str, object]]:
    """
    Parse a ``--method-params`` / ``--impute-params`` string into per-target dicts.

    Grammar::

        '<target>:<name>=<value>[,<name>=<value>...][;<target>:...]'

    For example ``'10_mnn:k=50;38_harman:limit=0.05'``.

    Parameters
    ----------
    raw
        The raw flag value.

    Returns
    -------
    Mapping of target name to its parameter dict. An empty or whitespace-only input
    yields an empty mapping.

    Raises
    ------
    ParamParseError
        On a missing ``:``, a missing ``=``, an empty name, or a repeated parameter.
    """
    result: dict[str, dict[str, object]] = {}
    if not raw or not raw.strip():
        return result

    for block in raw.split(";"):
        block = block.strip()
        if not block:
            continue
        if ":" not in block:
            raise ParamParseError(
                f"expected '<target>:<name>=<value>' but found {block!r} "
                "(no ':' separating the target from its parameters)"
            )
        target, _, assignments = block.partition(":")
        target = target.strip()
        if not target:
            raise ParamParseError(f"empty target name in {block!r}")

        params = result.setdefault(target, {})
        for assignment in assignments.split(","):
            assignment = assignment.strip()
            if not assignment:
                continue
            if "=" not in assignment:
                raise ParamParseError(
                    f"{target}: expected '<name>=<value>' but found {assignment!r}"
                )
            name, _, value = assignment.partition("=")
            name = name.strip()
            if not name:
                raise ParamParseError(
                    f"{target}: empty parameter name in {assignment!r}"
                )
            if name in params:
                raise ParamParseError(f"{target}: parameter {name!r} given twice")
            params[name] = coerce_value(value)

    return result


def validate_param_names(
    target: str, params: dict[str, object], declared: set[str]
) -> None:
    """
    Reject parameters the target does not declare, suggesting the closest match.

    Parameters
    ----------
    target
        Method or imputer key, used in the error message.
    params
        The parsed parameter dict.
    declared
        Names the target declares.

    Raises
    ------
    ParamParseError
        If any supplied name is undeclared.
    """
    unknown = sorted(set(params) - declared)
    if not unknown:
        return

    import difflib

    details = []
    for name in unknown:
        close = difflib.get_close_matches(name, sorted(declared), n=1)
        hint = f" (did you mean {close[0]!r}?)" if close else ""
        details.append(f"{name!r}{hint}")
    known = ", ".join(sorted(declared)) or "none"
    raise ParamParseError(
        f"{target}: unknown parameter(s): {', '.join(details)}. Declared: {known}"
    )


# --------------------------------------------------------------------------------------
# Parameter tags
# --------------------------------------------------------------------------------------


def _escape(text: str) -> str:
    """Percent-escape characters that would break a filename or the key grammar."""
    out = text.replace("%", "%25").replace("/", "%2F").replace(" ", "%20")
    return out.replace(KEY_SEPARATOR, "%5F%5F")


def _encode_value(value: object) -> str:
    """Render one parameter value in the compact, human-readable tag form."""
    if isinstance(value, bool):
        return "T" if value else "F"
    if value is None:
        return "None"
    if isinstance(value, (list, tuple)):
        return "+".join(_encode_value(item) for item in value)
    if isinstance(value, float):
        return repr(value)
    return _escape(str(value))


def params_hash(params: dict[str, object]) -> str:
    """
    Return a short, stable hash of a parameter dict.

    Stable across dict ordering and across processes, so the same parameters always
    produce the same tag.
    """
    canonical = json.dumps(params, sort_keys=True, default=str, separators=(",", ":"))
    return hashlib.sha256(canonical.encode()).hexdigest()[:_HASH_LENGTH]


def encode_param_tag(
    params: dict[str, object], *, max_length: int = MAX_PARAM_TAG_LENGTH
) -> str | None:
    """
    Encode non-default parameters into a readable filename component.

    ``{"k": 50, "limit": 0.05}`` becomes ``"k50-limit0.05"``. When the readable form
    would exceed ``max_length`` it is truncated and a stable hash of the *full* dict is
    appended, so long parameter sets stay unique without producing an unusable filename.

    Parameters
    ----------
    params
        The non-default parameters only. Defaults are recorded in ``run_manifest.json``,
        not in the filename.
    max_length
        Maximum length of the returned tag.

    Returns
    -------
    The tag, or ``None`` when ``params`` is empty — which is what keeps a default run's
    filename free of any tag segment at all.
    """
    if not params:
        return None

    readable = "-".join(
        f"{_escape(name)}{_encode_value(params[name])}" for name in sorted(params)
    )
    if len(readable) <= max_length:
        return readable

    keep = max_length - _HASH_LENGTH - 1
    return f"{readable[:keep].rstrip('-')}-{params_hash(params)}"


def validate_param_tag(tag: str) -> str:
    """
    Validate a user-supplied ``--param-tag``.

    Parameters
    ----------
    tag
        The requested tag.

    Returns
    -------
    The tag unchanged.

    Raises
    ------
    ParamParseError
        If the tag is empty, contains the key separator, or contains a path separator.
    """
    if not tag or not tag.strip():
        raise ParamParseError("--param-tag must not be empty")
    if KEY_SEPARATOR in tag:
        raise ParamParseError(
            f"--param-tag must not contain {KEY_SEPARATOR!r} — it separates key fields"
        )
    if "/" in tag or " " in tag:
        raise ParamParseError("--param-tag must not contain '/' or spaces")
    if _POST_RE.match(tag):
        raise ParamParseError(
            f"--param-tag {tag!r} is reserved — it would be read as a post-removal flag"
        )
    return tag


# --------------------------------------------------------------------------------------
# Output key grammar
# --------------------------------------------------------------------------------------


@dataclass(frozen=True)
class OutputKey:
    """One parsed output key: ``{imp}__{method}[__{ptag}][__post{0|1}]``."""

    imputation: str
    method: str
    param_tag: str | None = None
    post_rm: bool | None = None

    def __str__(self) -> str:
        return build_output_key(
            self.imputation, self.method, self.param_tag, self.post_rm
        )


def _reject_separator(field: str, value: str) -> None:
    """Raise if a key field is empty or contains the field separator."""
    if not value:
        raise KeyParseError(f"{field} must not be empty")
    if KEY_SEPARATOR in value:
        raise KeyParseError(
            f"{field} {value!r} must not contain {KEY_SEPARATOR!r} — it separates fields"
        )


def build_output_key(
    imputation: str,
    method: str,
    param_tag: str | None = None,
    post_rm: bool | None = None,
) -> str:
    """
    Build an output key from its parts.

    Parameters
    ----------
    imputation
        Imputer key, e.g. ``"strict"``.
    method
        Method key, e.g. ``"10_mnn"``.
    param_tag
        Encoded non-default parameters, or ``None`` to omit the segment.
    post_rm
        ``None`` omits the post-removal segment entirely, which is what a run with
        post-removal disabled produces.

    Returns
    -------
    The key, without any file extension.
    """
    _reject_separator("imputation", imputation)
    _reject_separator("method", method)

    parts = [imputation, method]
    if param_tag is not None:
        parts.append(validate_param_tag(param_tag))
    if post_rm is not None:
        parts.append("post1" if post_rm else "post0")
    return KEY_SEPARATOR.join(parts)


def parse_output_key(stem: str) -> OutputKey:
    """
    Parse an output key, raising on anything that does not match the grammar.

    Accepts two to four ``__``-separated segments. A third segment is read as a
    post-removal flag when it looks like ``post0``/``post1`` and as a parameter tag
    otherwise; this is unambiguous because :func:`validate_param_tag` refuses tags of
    that shape.

    Parameters
    ----------
    stem
        The key with any file extension already stripped.

    Returns
    -------
    The parsed :class:`OutputKey`.

    Raises
    ------
    KeyParseError
        On the wrong number of segments, an empty segment, or a malformed post flag.
    """
    if not stem or not stem.strip():
        raise KeyParseError("output key must not be empty")

    parts = stem.split(KEY_SEPARATOR)
    if any(not part for part in parts):
        raise KeyParseError(f"output key {stem!r} has an empty segment")
    if not 2 <= len(parts) <= 4:
        raise KeyParseError(
            f"output key {stem!r} has {len(parts)} segment(s); expected 2-4 matching "
            "'{imp}__{method}[__{ptag}][__post0|post1]'"
        )

    imputation, method = parts[0], parts[1]
    rest = parts[2:]

    post_rm: bool | None = None
    if rest:
        match = _POST_RE.match(rest[-1])
        if match:
            post_rm = match.group(1) == "1"
            rest = rest[:-1]

    if len(rest) > 1:
        raise KeyParseError(
            f"output key {stem!r} has more than one parameter-tag segment"
        )
    param_tag = rest[0] if rest else None

    return OutputKey(
        imputation=imputation, method=method, param_tag=param_tag, post_rm=post_rm
    )


def build_prepared_key(imputation: str, param_tag: str | None, kind: str) -> str:
    """
    Build a key for a cached post-imputation matrix: ``{imp}[__{ptag}]__{kind}``.

    Parameters
    ----------
    imputation
        Imputer key.
    param_tag
        Encoded imputer parameters, or ``None``.
    kind
        Either ``"exp"`` or ``"ann"``.

    Returns
    -------
    The key, without any file extension.
    """
    if kind not in {"exp", "ann"}:
        raise KeyParseError(f"prepared key kind must be 'exp' or 'ann', got {kind!r}")
    _reject_separator("imputation", imputation)

    parts = [imputation]
    if param_tag is not None:
        parts.append(validate_param_tag(param_tag))
    parts.append(kind)
    return KEY_SEPARATOR.join(parts)
