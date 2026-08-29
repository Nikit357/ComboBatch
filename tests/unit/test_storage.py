"""Storage backends: local round-trips, S3 via moto, and the 403/404 distinction."""

from __future__ import annotations

import os

import pytest

from combobatch import storage

BUCKET = "combobatch-test"


@pytest.fixture
def local(tmp_path):
    return storage.LocalBackend(tmp_path)


class TestLocalBackend:
    def test_write_then_read(self, local):
        local.write_bytes("a/b.txt", b"payload")
        assert local.read_bytes("a/b.txt") == b"payload"

    def test_write_creates_intermediate_directories(self, local, tmp_path):
        local.write_bytes("deep/er/still/x.txt", b"1")
        assert (tmp_path / "deep/er/still/x.txt").is_file()

    def test_exists(self, local):
        assert local.exists("nope.txt") is False
        local.write_bytes("yes.txt", b"1")
        assert local.exists("yes.txt") is True

    def test_list_keys_is_relative_and_recursive(self, local):
        local.write_bytes("exp/a.tsv", b"1")
        local.write_bytes("exp/b.tsv", b"1")
        local.write_bytes("metrics/c.json", b"1")
        assert sorted(local.list_keys("exp")) == ["exp/a.tsv", "exp/b.tsv"]
        assert len(list(local.list_keys())) == 3

    def test_list_keys_on_missing_prefix_is_empty(self, local):
        assert list(local.list_keys("nothing-here")) == []

    def test_delete_is_idempotent(self, local):
        local.write_bytes("x.txt", b"1")
        local.delete("x.txt")
        local.delete("x.txt")
        assert local.exists("x.txt") is False


class TestUriDispatch:
    def test_local_root(self, tmp_path):
        backend = storage.backend_for_root(str(tmp_path))
        assert isinstance(backend, storage.LocalBackend)

    def test_local_file_resolves_to_parent_and_name(self, tmp_path):
        backend, key = storage.resolve(str(tmp_path / "sub" / "file.tsv"))
        assert isinstance(backend, storage.LocalBackend)
        assert key == "file.tsv"
        assert backend.root.endswith("sub")

    def test_is_s3_uri(self):
        assert storage.is_s3_uri("s3://bucket/key") is True
        assert storage.is_s3_uri("./local/path") is False

    def test_malformed_s3_uri_raises(self):
        with pytest.raises(ValueError, match="no bucket"):
            storage.backend_for_root("s3://")

    def test_s3_uri_without_object_raises_on_resolve(self):
        with pytest.raises(ValueError, match="bucket, not an object"):
            storage.resolve("s3://bucket")

    def test_round_trip_through_a_local_backend_preserves_layout(self, tmp_path):
        backend = storage.backend_for_root(str(tmp_path / "run"))
        backend.write_bytes("exp/strict__01_raw__post0.tsv.gz", b"data")
        assert backend.exists("exp/strict__01_raw__post0.tsv.gz")
        assert backend.uri("exp/x").endswith("/exp/x")


@pytest.fixture
def aws_credentials(monkeypatch):
    """moto must never fall through to real credentials."""
    for var in (
        "AWS_ACCESS_KEY_ID",
        "AWS_SECRET_ACCESS_KEY",
        "AWS_SECURITY_TOKEN",
        "AWS_SESSION_TOKEN",
    ):
        monkeypatch.setenv(var, "testing")
    monkeypatch.setenv("AWS_DEFAULT_REGION", "us-east-1")


@pytest.fixture
def s3(aws_credentials):
    from moto import mock_aws

    with mock_aws():
        import boto3

        boto3.client("s3").create_bucket(Bucket=BUCKET)
        yield storage.S3Backend(BUCKET, "prefix")


class TestS3Backend:
    def test_write_then_read(self, s3):
        s3.write_bytes("a/b.txt", b"payload")
        assert s3.read_bytes("a/b.txt") == b"payload"

    def test_exists(self, s3):
        assert s3.exists("missing.txt") is False
        s3.write_bytes("there.txt", b"1")
        assert s3.exists("there.txt") is True

    def test_list_keys_strips_the_root_prefix(self, s3):
        s3.write_bytes("exp/a.tsv", b"1")
        s3.write_bytes("exp/b.tsv", b"1")
        assert sorted(s3.list_keys("exp")) == ["exp/a.tsv", "exp/b.tsv"]

    def test_delete(self, s3):
        s3.write_bytes("x.txt", b"1")
        s3.delete("x.txt")
        assert s3.exists("x.txt") is False

    def test_root_and_uri(self, s3):
        assert s3.root == f"s3://{BUCKET}/prefix"
        assert s3.uri("exp/x") == f"s3://{BUCKET}/prefix/exp/x"

    def test_backend_without_prefix(self, aws_credentials):
        from moto import mock_aws

        with mock_aws():
            import boto3

            boto3.client("s3").create_bucket(Bucket=BUCKET)
            backend = storage.S3Backend(BUCKET)
            backend.write_bytes("top.txt", b"1")
            assert backend.root == f"s3://{BUCKET}"
            assert list(backend.list_keys()) == ["top.txt"]


class TestExistsDistinguishes403From404:
    """The donor swallowed every ClientError, so a permissions failure looked exactly
    like "output not computed yet" — and the run silently did the wrong thing."""

    @staticmethod
    def _client_error(code: str, status: int):
        import botocore.exceptions

        return botocore.exceptions.ClientError(
            {
                "Error": {"Code": code, "Message": code},
                "ResponseMetadata": {"HTTPStatusCode": status},
            },
            "HeadObject",
        )

    @pytest.mark.parametrize("code,status", [("404", 404), ("NoSuchKey", 404)])
    def test_missing_returns_false(self, s3, code, status):
        def raise_missing(**_kwargs):
            raise self._client_error(code, status)

        s3._client.head_object = raise_missing
        assert s3.exists("whatever") is False

    @pytest.mark.parametrize(
        "code,status", [("403", 403), ("AccessDenied", 403), ("500", 500)]
    )
    def test_other_errors_propagate(self, s3, code, status):
        import botocore.exceptions

        def raise_other(**_kwargs):
            raise self._client_error(code, status)

        s3._client.head_object = raise_other
        with pytest.raises(botocore.exceptions.ClientError):
            s3.exists("whatever")
