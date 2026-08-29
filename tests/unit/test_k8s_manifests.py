"""The manifests encode constraints that fail silently when they are wrong.

A `nodeSelector` on a label no node carries leaves the pod `Pending` with no error; a
missing `capacity-type` lands a six-hour job on a Spot node that can be reclaimed with two
minutes' notice; a re-introduced requirements ConfigMap drifts from `requirements.txt`
without anything noticing. None of that is caught by `kubectl apply`, so it is caught here.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

REPO_ROOT = Path(__file__).resolve().parents[2]
K8S_DIR = REPO_ROOT / "k8s"
POD = K8S_DIR / "pod-combobatch.yaml"
JOB = K8S_DIR / "job-combobatch.yaml"
PVC = K8S_DIR / "create_pvc.yaml"

IMAGE = "ghcr.io/nikit357/combobatch:latest"
THREAD_VARS = {
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
}


def load(path: Path) -> dict:
    """The single resource in a manifest."""
    docs = [doc for doc in yaml.safe_load_all(path.read_text()) if doc]
    assert len(docs) == 1, f"{path.name} should hold exactly one resource"
    return docs[0]


def pod_spec(path: Path) -> dict:
    """The PodSpec, whether the file is a Pod or a Job wrapping one."""
    resource = load(path)
    if resource["kind"] == "Job":
        return resource["spec"]["template"]["spec"]
    return resource["spec"]


WORKLOADS = pytest.mark.parametrize("path", [POD, JOB], ids=["pod", "job"])


class TestScheduling:
    @WORKLOADS
    def test_on_demand_capacity_is_requested(self, path):
        # Spot's two-minute SIGTERM cannot save a six-hour missForest job.
        selector = pod_spec(path)["nodeSelector"]
        assert selector["karpenter.sh/capacity-type"] == "on-demand"

    @WORKLOADS
    def test_the_capacity_type_label_has_the_core_karpenter_prefix(self, path):
        # `karpenter.k8s.aws/capacity-type` (what the plan wrote) matches no node in this
        # cluster, and the only symptom is a pod that stays Pending forever.
        selector = pod_spec(path)["nodeSelector"]
        assert "karpenter.k8s.aws/capacity-type" not in selector

    @WORKLOADS
    def test_the_c6a_family_is_requested(self, path):
        assert (
            pod_spec(path)["nodeSelector"]["karpenter.k8s.aws/instance-family"] == "c6a"
        )

    @WORKLOADS
    def test_karpenter_may_not_consolidate_the_node_away(self, path):
        resource = load(path)
        metadata = (
            resource["spec"]["template"]["metadata"]
            if resource["kind"] == "Job"
            else resource["metadata"]
        )
        assert metadata["annotations"]["karpenter.sh/do-not-disrupt"] == "true"

    @WORKLOADS
    def test_the_node_group_toleration_and_affinity_are_paired(self, path):
        spec = pod_spec(path)
        assert any(t["key"] == "node-group" for t in spec["tolerations"])
        terms = spec["affinity"]["nodeAffinity"][
            "requiredDuringSchedulingIgnoredDuringExecution"
        ]["nodeSelectorTerms"]
        keys = {expr["key"] for term in terms for expr in term["matchExpressions"]}
        assert "node-group" in keys


class TestImageAndSecrets:
    @WORKLOADS
    def test_the_public_ghcr_image_is_used(self, path):
        container = pod_spec(path)["containers"][0]
        assert container["image"] == IMAGE
        assert container["imagePullPolicy"] == "Always"

    @WORKLOADS
    def test_no_image_pull_secret_is_needed(self, path):
        # The package is public; a pull secret here would hide a private package behind a
        # working pull and reintroduce the donors' node-IAM dependency.
        assert "imagePullSecrets" not in pod_spec(path)

    @WORKLOADS
    def test_nothing_is_installed_at_startup(self, path):
        # The donors' pods apt-installed R and ~25 packages at boot: 30-40 minutes, and
        # drifting from requirements.txt. The image is the environment now.
        text = path.read_text()
        for forbidden in ("apt-get install", "pip install", "BiocManager::install"):
            assert forbidden not in text, f"{path.name} installs software at startup"

    @WORKLOADS
    def test_no_configmap_is_mounted(self, path):
        # A requirements ConfigMap is the thing that used to drift out of sync.
        volumes = pod_spec(path).get("volumes", [])
        assert not [v for v in volumes if "configMap" in v]

    def test_the_ssh_key_comes_from_a_secret_and_is_not_committed(self):
        spec = pod_spec(POD)
        key_volume = next(v for v in spec["volumes"] if v["name"] == "ssh-pubkey")
        assert key_volume["secret"]["secretName"] == "danya-nikitin-ssh-pubkey"
        # The donor's pod-ssh.yaml pasted the public key inline in a ConfigMap.
        assert "ssh-ed25519 " not in POD.read_text()

    @WORKLOADS
    def test_aws_credentials_are_mounted_read_only(self, path):
        mount = next(
            m
            for m in pod_spec(path)["containers"][0]["volumeMounts"]
            if m["name"] == "aws-credentials"
        )
        assert mount["readOnly"] is True


class TestRuntime:
    @WORKLOADS
    def test_every_thread_variable_is_pinned(self, path):
        env = {
            e["name"]: e.get("value") for e in pod_spec(path)["containers"][0]["env"]
        }
        assert THREAD_VARS <= set(env)
        assert all(env[var] == "1" for var in THREAD_VARS)

    @WORKLOADS
    def test_requests_equal_limits(self, path):
        # Guaranteed QoS: a memory spike in one worker must not get the pod evicted.
        resources = pod_spec(path)["containers"][0]["resources"]
        assert resources["requests"] == resources["limits"]

    @WORKLOADS
    def test_the_workspace_pvc_is_mounted(self, path):
        spec = pod_spec(path)
        volume = next(v for v in spec["volumes"] if v["name"] == "work")
        assert volume["persistentVolumeClaim"]["claimName"] == "danya-nikitin-fl"
        mount = next(
            m for m in spec["containers"][0]["volumeMounts"] if m["name"] == "work"
        )
        assert mount["mountPath"] == "/workspace"

    def test_the_job_invokes_the_cli_explicitly(self):
        container = pod_spec(JOB)["containers"][0]
        assert container["command"] == ["combobatch"]
        assert container["args"][0] == "dispatch"

    def test_the_job_does_not_silently_retry(self):
        # A failed six-hour run should be inspected, not repeated.
        assert load(JOB)["spec"]["backoffLimit"] == 0
        assert pod_spec(JOB)["restartPolicy"] == "Never"

    def test_the_job_leaves_room_for_the_sigterm_handler(self):
        # The dispatcher flushes its failed-jobs log on SIGTERM so --retry-failed works.
        assert pod_spec(JOB)["terminationGracePeriodSeconds"] >= 30


class TestPvc:
    def test_the_claim_matches_what_the_workloads_mount(self):
        claim = load(PVC)
        assert claim["metadata"]["name"] == "danya-nikitin-fl"
        assert claim["spec"]["accessModes"] == ["ReadWriteOnce"]

    def test_every_manifest_names_the_namespace(self):
        for path in (POD, JOB, PVC):
            assert load(path)["metadata"]["namespace"] == "rnd-sandbox"
