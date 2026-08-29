# Scripts

Maintainer utilities. None of these are part of the installed package.

| Script | Purpose |
|---|---|
| `build_and_push_image.sh` | `docker buildx` for `linux/amd64`, tagged `latest` and by date, pushed to GHCR. |
| `smoke_test.sh` | The image's acceptance gate: build, R 4.5.x, Octave, the `matlab` wrapper, zero backend skips. |
| `build_example_dataset.py` | Regenerates `examples/example_dataset/` from the source annotation and expression matrix. Run once; the output is committed. *(Phase 8.)* |

```bash
COMBOBATCH_PUSH=0 bash scripts/build_and_push_image.sh   # build locally, do not push
bash scripts/build_and_push_image.sh v0.1.0              # tag :v0.1.0 as well as :latest

bash scripts/smoke_test.sh                  # builds combobatch:test, then checks it
bash scripts/smoke_test.sh my-image:tag     # checks an image that already exists
```

`smoke_test.sh` reports `combobatch selftest` as skipped until that subcommand lands in
Phase 9; every other check is live.
