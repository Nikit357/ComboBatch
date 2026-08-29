# Kubernetes

Namespace `rnd-sandbox`. The image carries the whole environment, so a pod starts in the
time it takes to pull it — there is no install script and no requirements ConfigMap.

| File | Role |
|---|---|
| `pod-combobatch.yaml` | Interactive pod with SSH, for exploratory runs. |
| `job-combobatch.yaml` | `batch/v1` Job for unattended dispatcher runs. |
| `create_pvc.yaml` | The `/workspace` PVC. Already exists in this cluster. |

## Before the first launch

Two Secrets must exist in the namespace. Neither is committed — the donor's `pod-ssh.yaml`
had the SSH key inline in a ConfigMap, which put a personal identifier in the repository
and made the manifest usable by exactly one person.

```bash
kubectl create secret generic danya-nikitin-ssh-pubkey \
    --from-file=authorized_keys=$HOME/.ssh/id_ed25519.pub -n rnd-sandbox

kubectl create secret generic aws-credentials \
    --from-file=credentials=$HOME/.aws/credentials \
    --from-file=config=$HOME/.aws/config -n rnd-sandbox
```

Both already exist in `rnd-sandbox` (verified 2026-08-29), as does the PVC
`danya-nikitin-fl` — Bound, 100 GiB, `ebs-gp3-delete`.

The image must also be **public** on GHCR. There is no `imagePullSecrets` anywhere in
these manifests; a private package fails the pull with a 403 that reads like a missing
image. See [`../docker/BUILD_AND_PUSH.md`](../docker/BUILD_AND_PUSH.md) step 7.

## Interactive pod

```bash
kubectl apply -f k8s/pod-combobatch.yaml -n rnd-sandbox
kubectl get pod danya-nikitin-combobatch -n rnd-sandbox -w      # Pending → Running
kubectl logs danya-nikitin-combobatch -n rnd-sandbox            # versions, then "ready"
kubectl exec -it danya-nikitin-combobatch -n rnd-sandbox -- bash
```

Karpenter has to provision a `c6a` on-demand node first, so the first `Pending` phase is a
couple of minutes; the pull is a few more. If it is still `Pending` after that,
`kubectl describe pod` names the unsatisfied constraint.

The startup log prints the R version. **It must say 4.5.x** — 4.6 means the image pin
failed and rpy2 will segfault around `33_amdbnorm`.

For SSH (VS Code Remote, `rsync`, long sessions in `tmux`):

```bash
kubectl port-forward pod/danya-nikitin-combobatch 2222:22 -n rnd-sandbox
ssh -p 2222 root@localhost
```

Then verify the environment and run something:

```bash
combobatch selftest                       # Phase 9
combobatch list-methods                   # 34 rows
combobatch dispatch --config /workspace/configs/run.yaml --dry-run
```

## Unattended run

Put a config on the PVC, point the Job's `args` at it, then:

```bash
kubectl apply -f k8s/job-combobatch.yaml -n rnd-sandbox
kubectl logs -f job/danya-nikitin-combobatch-dispatch -n rnd-sandbox
kubectl delete job danya-nikitin-combobatch-dispatch -n rnd-sandbox   # when done
```

`backoffLimit: 0` — a failed six-hour run is inspected, not silently repeated. Re-running
after a fix is cheap because `--skip-if-exists` keeps finished combinations.

**The PVC is ReadWriteOnce**, so the Job and the interactive pod cannot run at the same
time. Delete one before starting the other.

## Cleaning up

Karpenter deprovisions the node once nothing is scheduled on it, so deleting the pod is
what stops the bill.

```bash
kubectl delete pod danya-nikitin-combobatch -n rnd-sandbox
```

Do **not** delete the PVC to tidy up: its storage class reclaim policy is `Delete`, so the
volume and everything on it goes with it.
