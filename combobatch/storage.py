"""Storage backends: a local filesystem and S3 behind one interface.

Chosen by URI scheme, so ``--out ./run/`` and ``--out s3://bucket/run/`` produce an
identical layout and every caller stays storage-agnostic.

Two deliberate departures from the pipeline this was ported from:

* One boto3 client per backend instance, created once. The donor code constructed a
  fresh client at fifteen separate call sites.
* :meth:`StorageBackend.exists` distinguishes 404 from 403. The donor swallowed every
  ``ClientError``, so a permissions failure was indistinguishable from "not computed
  yet" — and a run would silently recompute, or silently skip, depending on the caller.
"""

from __future__ import annotations

import shutil
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Iterator

S3_SCHEME = "s3://"

# Error codes that genuinely mean "this object is not there". Anything else — most
# importantly 403/AccessDenied — is a real failure and must surface.
_MISSING_CODES = frozenset({"404", "NoSuchKey", "NotFound", "NoSuchBucket"})


class StorageBackend(ABC):
    """A key-value store rooted at some prefix, addressed by relative keys."""

    @property
    @abstractmethod
    def root(self) -> str:
        """Return the backend's root as a URI or path, for logging and error messages."""

    @abstractmethod
    def exists(self, key: str) -> bool:
        """Report whether ``key`` exists, without transferring its contents."""

    @abstractmethod
    def read_bytes(self, key: str) -> bytes:
        """Return the full contents of ``key``."""

    @abstractmethod
    def write_bytes(self, key: str, data: bytes) -> None:
        """Write ``data`` to ``key``, creating any intermediate directories."""

    @abstractmethod
    def list_keys(self, prefix: str = "") -> Iterator[str]:
        """Yield every key under ``prefix``, relative to the backend root."""

    @abstractmethod
    def delete(self, key: str) -> None:
        """Delete ``key``. Deleting a missing key is not an error."""

    def uri(self, key: str) -> str:
        """Return the fully qualified location of ``key``, for logs and manifests."""
        return f"{self.root.rstrip('/')}/{key}"


class LocalBackend(StorageBackend):
    """Filesystem-backed storage rooted at a directory."""

    def __init__(self, root: str | Path) -> None:
        self._root = Path(root).expanduser()

    @property
    def root(self) -> str:
        return str(self._root)

    def _path(self, key: str) -> Path:
        return self._root / key

    def exists(self, key: str) -> bool:
        return self._path(key).is_file()

    def read_bytes(self, key: str) -> bytes:
        return self._path(key).read_bytes()

    def write_bytes(self, key: str, data: bytes) -> None:
        path = self._path(key)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)

    def list_keys(self, prefix: str = "") -> Iterator[str]:
        base = self._root / prefix if prefix else self._root
        if not base.exists():
            return
        for path in sorted(base.rglob("*")):
            if path.is_file():
                yield str(path.relative_to(self._root))

    def delete(self, key: str) -> None:
        self._path(key).unlink(missing_ok=True)

    def rmtree(self) -> None:
        """Remove the whole root. Only used by tests and by explicit cleanup."""
        shutil.rmtree(self._root, ignore_errors=True)


class S3Backend(StorageBackend):
    """S3-backed storage rooted at ``s3://{bucket}/{prefix}``."""

    def __init__(
        self, bucket: str, prefix: str = "", *, endpoint_url: str | None = None
    ) -> None:
        import boto3

        self._bucket = bucket
        self._prefix = prefix.strip("/")
        self._client = boto3.client("s3", endpoint_url=endpoint_url)

    @property
    def root(self) -> str:
        if self._prefix:
            return f"{S3_SCHEME}{self._bucket}/{self._prefix}"
        return f"{S3_SCHEME}{self._bucket}"

    def _key(self, key: str) -> str:
        return f"{self._prefix}/{key}" if self._prefix else key

    def exists(self, key: str) -> bool:
        import botocore.exceptions

        try:
            self._client.head_object(Bucket=self._bucket, Key=self._key(key))
            return True
        except botocore.exceptions.ClientError as exc:
            code = str(exc.response.get("Error", {}).get("Code", ""))
            status = exc.response.get("ResponseMetadata", {}).get("HTTPStatusCode")
            if code in _MISSING_CODES or status == 404:
                return False
            # 403 and friends reach here on purpose: a permissions problem must never
            # be mistaken for a missing output.
            raise

    def read_bytes(self, key: str) -> bytes:
        response = self._client.get_object(Bucket=self._bucket, Key=self._key(key))
        return response["Body"].read()

    def write_bytes(self, key: str, data: bytes) -> None:
        self._client.put_object(Bucket=self._bucket, Key=self._key(key), Body=data)

    def list_keys(self, prefix: str = "") -> Iterator[str]:
        full_prefix = self._key(prefix) if prefix else self._prefix
        paginator = self._client.get_paginator("list_objects_v2")
        strip = len(self._prefix) + 1 if self._prefix else 0
        for page in paginator.paginate(Bucket=self._bucket, Prefix=full_prefix):
            for item in page.get("Contents", []):
                yield item["Key"][strip:]

    def delete(self, key: str) -> None:
        self._client.delete_object(Bucket=self._bucket, Key=self._key(key))


def is_s3_uri(uri: str) -> bool:
    """Report whether ``uri`` addresses S3."""
    return uri.startswith(S3_SCHEME)


def _split_s3_uri(uri: str) -> tuple[str, str]:
    """Split ``s3://bucket/key`` into its bucket and key."""
    remainder = uri[len(S3_SCHEME) :]
    bucket, _, key = remainder.partition("/")
    if not bucket:
        raise ValueError(f"malformed S3 URI (no bucket): {uri!r}")
    return bucket, key


def backend_for_root(uri: str, *, endpoint_url: str | None = None) -> StorageBackend:
    """
    Return a backend rooted at ``uri``, for use as an output root.

    Parameters
    ----------
    uri
        A local directory path or an ``s3://bucket/prefix`` URI.
    endpoint_url
        Alternative S3 endpoint, for MinIO or a mocked service.

    Returns
    -------
    A backend whose keys are interpreted relative to ``uri``.
    """
    if is_s3_uri(uri):
        bucket, prefix = _split_s3_uri(uri)
        return S3Backend(bucket, prefix, endpoint_url=endpoint_url)
    return LocalBackend(uri)


def resolve(uri: str, *, endpoint_url: str | None = None) -> tuple[StorageBackend, str]:
    """
    Split a single-file URI into a backend rooted at its parent, plus the filename.

    Parameters
    ----------
    uri
        A local file path or an ``s3://bucket/key`` URI.
    endpoint_url
        Alternative S3 endpoint.

    Returns
    -------
    ``(backend, key)`` such that ``backend.read_bytes(key)`` reads ``uri``.
    """
    if is_s3_uri(uri):
        bucket, key = _split_s3_uri(uri)
        if not key:
            raise ValueError(f"S3 URI addresses a bucket, not an object: {uri!r}")
        prefix, _, name = key.rpartition("/")
        return S3Backend(bucket, prefix, endpoint_url=endpoint_url), name

    path = Path(uri).expanduser()
    return LocalBackend(path.parent), path.name
