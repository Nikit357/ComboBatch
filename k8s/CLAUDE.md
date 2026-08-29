# CLAUDE.md — `k8s/`

## Invariants

- **The image is the environment.** No requirements ConfigMap, no startup install script.
  The donors' pods installed R and ~25 packages at boot, taking 30–40 minutes and drifting
  from `requirements.txt`; that entire class of problem is designed out. Do not reintroduce
  it. The only startup work in `pod-combobatch.yaml` is enabling `sshd`, which is why
  `openssh-server` is in the image rather than apt-installed at boot.
- **`capacity-type: on-demand` is required, not optional.** Spot reclamation gives a
  two-minute SIGTERM, which cannot save a multi-hour missForest job. Pair it with the
  `karpenter.sh/do-not-disrupt` annotation.
- **The capacity-type label lives under `karpenter.sh/`, not `karpenter.k8s.aws/`.** The
  plan and an old donor note both write `karpenter.k8s.aws/capacity-type`, which matches no
  node in this cluster — a `nodeSelector` on it leaves the pod `Pending` forever with no
  error anywhere. Only *instance* attributes (`instance-family`, `instance-size`, …) carry
  the AWS-provider prefix. Verified against live nodes 2026-08-29:

  ```bash
  kubectl get nodes -L karpenter.sh/capacity-type,karpenter.k8s.aws/instance-family
  # CAPACITY-TYPE is populated (on-demand / spot); the karpenter.k8s.aws/ spelling is empty
  ```

- **The GHCR package is public, so no `imagePullSecrets`.** A new GHCR package is private by
  default — if the pod fails to pull with a 403 that looks like a missing image, that is the
  cause.
- **The AWS credentials Secret is for data only** — reading `s3://` inputs and writing
  outputs. It is never involved in pulling the image.
- **Thread-pinning env vars in the manifest are not sufficient.** They do not reach processes
  launched over SSH. The real pinning happens inside the Python entry point; the manifest
  values are belt and braces.
- **The PVC is `ReadWriteOnce`** — one pod at a time can mount it, so the interactive pod and
  the dispatch Job are mutually exclusive. Its reclaim policy is `Delete`: deleting the claim
  destroys the data.
- **Secrets are never committed.** Both the SSH public key and the AWS credentials are
  created out of band with `kubectl create secret`; the donor's `pod-ssh.yaml` embedded the
  key inline, which is what this replaces.
- Redirect long-run logs to the PVC so they survive a pod restart — `kubectl logs` is
  otherwise the only copy, and it disappears with the pod.

## Shape

`pod-combobatch.yaml` overrides `command`, so the image entry point is bypassed and the
container simply stays alive. `job-combobatch.yaml` sets `command: ["combobatch"]`
explicitly rather than relying on the entry point to infer the prefix.

Both request 32 vCPU / 64 GiB with `requests == limits` (Guaranteed QoS, so a memory spike
in one worker cannot get the pod evicted). That is exactly a `c6a.8xlarge`; raise both
together to keep the family's 2 GiB-per-vCPU ratio.

## Verifying a change

`kubectl apply --dry-run=client` validates against the live cluster's schema and is safe:

```bash
kubectl apply --dry-run=client -f k8s/ -n rnd-sandbox
```
