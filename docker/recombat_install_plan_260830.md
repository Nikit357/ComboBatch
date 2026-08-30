# Plan — reComBat cannot be pip-installed as declared

**Status:** awaiting approval. No code has been changed.
**Date:** 2026-08-30
**Trigger:** `bash scripts/smoke_test.sh`, log at `docker/docker_build_logs_260830_1.txt`.
Failed at step 14/20, `pip install -r requirements.txt`, after 177 s.

---

## Overview

**The VPN was off and the staging worked.** Every R stage came back `CACHED`, `github` and
`verify` completed, Octave installed, and the build reached the Python layer for the first
time — four steps further than any previous attempt. Nothing in the last three plans
regressed.

The Python layer fails on one line of `requirements.txt`:

```
git+https://github.com/BorgwardtLab/reComBat   # 30_recombat
```

reComBat 0.1.4 declares `sklearn = "^0.0"` — the **deprecated PyPI shim**, a stub package
whose only remaining behaviour is to raise:

```
The 'sklearn' PyPI package is deprecated, use 'scikit-learn' rather than 'sklearn'
error: metadata-generation-failed
```

This is drift of the same kind as `fs` 2.x: the shim began hard-erroring after the donor's
image was built, and the donor's `requirements.txt` carries the identical line, so it would
fail the same way today.

**Reading reComBat's own metadata turned up a second problem that would have bitten next:**

```toml
[tool.poetry.dependencies]
python = ">=3.8,<3.11"     ← excludes the version this project pins exactly
sklearn = "^0.0"           ← the immediate failure
```

reComBat declares that it **does not support Python 3.11**, and Python 3.11 exactly is the
first non-negotiable constraint in this project's `CLAUDE.md`. Past the `sklearn` error, pip
would have refused it on `Requires-Python` instead.

Both are stale metadata, not real constraints — verified, not assumed. Installed past its
own declarations on Python 3.11.16, reComBat runs and does its job:

```
python 3.11.16   out DataFrame (24, 40)   finite: True
batch gap before: 2.975  →  after: -0.011
```

So the fix is to install reComBat **alone**, past both stale declarations, and to supply its
real dependencies ourselves:

```
pip install --no-deps --ignore-requires-python \
    "reComBat @ git+https://github.com/BorgwardtLab/reComBat@b02bd025c967…"
```

`--no-deps` is a global pip flag with no per-requirement form, which is why this must become
its own install step rather than a line in `requirements.txt`.

---

## Background — the evidence

### What `--no-deps` discards, audited one by one

reComBat's declared `requires_dist`, from PyPI:

| Declared | Real? | Handled by |
|---|---|---|
| `sklearn (>=0.0,<0.1)` | **No** — a stub that installs no module and raises on build | Real `scikit-learn>=1.3,<1.6` is already pinned; `reComBat.py` imports `sklearn.linear_model`, which resolves to it |
| `fire (>=0.4.0,<0.5.0)` | **No** — CLI only. `reComBat/__init__.py` imports `.reComBat` and never `.cli` | Nothing. Verified: with `fire` absent, `import fire` fails and the library still runs |
| `numpy (>=1.20.0,<2.0.0)` | Yes | Already `numpy>=1.24,<2.0` — strictly tighter |
| `pandas (>=1.3.4,<2.0.0)` | Yes | Already `pandas>=1.5,<2.0` — strictly tighter |
| `tqdm (>=4.62.3,<5.0.0)` | **Yes** — `reComBat.py` line 13 imports it | **Not declared anywhere in this project.** It arrives only as a transitive dependency of `umap-learn` |

So `--no-deps` loses exactly one real requirement, `tqdm`, and it must be declared.

### Why not the PyPI wheel — a trap worth documenting

reComBat **is** on PyPI, as `reComBat 0.1.4`, a pure-Python wheel. Our own comment says
otherwise (`# 30_recombat; not on PyPI`), and that comment is wrong. The wheel installs in
**1.6 seconds** against **53 seconds** for the git clone plus build, needs no `git`, and is a
fixed artefact rather than a moving branch. It looks strictly better.

It is not. Diffing the published wheel against git `master`:

```
diff -r /wheel/reComBat /git/reComBat
59c59
<                  model='linear',          ← published wheel 0.1.4
---
>                  model='elastic_net',     ← git master, what every build so far installed
```

**One line, and it is the default model.** Both are versioned `0.1.4`. Switching source
would silently change every `30_recombat` result — visible in the test above as a batch gap
of −0.035 instead of −0.011 — and would break comparability with the dissertation's
published `metrics_comprehensive.csv`, which is the Phase 9 full-scale validation target.

So: **stay on git, and pin the commit.** The current line tracks `master` with no pin, which
is how a silent behavioural change would arrive unannounced. The build log names the commit
that has been installed all along: `b02bd025c9671e75f37dc513d7f466b05448364d`.

### What was measured

| Check | Result |
|---|---|
| `pip install git+…/reComBat` (today's line) | fails, `metadata-generation-failed` on `sklearn` |
| `pip install --no-deps --ignore-requires-python git+…` | installs |
| `reComBat().fit_transform` on Python 3.11.16 | runs; batch gap 2.975 → −0.011 |
| `import fire` after `--no-deps` | `ModuleNotFoundError` — and the library still works |
| published wheel vs git `master` | one line differs: the default `model` |
| VPN state during the build | off — `en0` holds the default route |

---

## Files to change

### 1. `requirements.txt`

#### 1a. Replace the git line with a pointer

`--no-deps` cannot be scoped to one requirement, so the line has to leave this file. It must
leave a trail, because a reader who finds no reComBat here will reasonably conclude the
image lacks it.

**Before** (lines 33–34):

```
git+https://github.com/BorgwardtLab/reComBat   # 30_recombat; not on PyPI. Needs `git`
                          # in the image and network access at build time.
```

**After:**

```
# 30_recombat's reComBat is NOT installed from this file. It declares `sklearn` (the
# deprecated PyPI stub, which now only raises) and `python >=3.8,<3.11`, so pip refuses
# it twice over. Both are stale metadata - it runs correctly on 3.11 - so the Dockerfile
# installs it alone with --no-deps --ignore-requires-python, pinned to a commit. See
# docker/recombat_install_plan_260830.md.
tqdm>=4.62.3              # reComBat imports it, and --no-deps means we must say so
```

#### 1b. Correct the header's false claim

The header says the constraints here are kept identical to `pyproject.toml` and that
`tests/unit/test_docs_in_sync.py` enforces it. **No test reads `pyproject.toml`** — verified
by grep across `tests/`. Either the claim goes or the test arrives; section 5 makes it true.

### 2. `docker/Dockerfile` — split the pip layer

**Before:**

```dockerfile
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm -f /tmp/requirements.txt
```

**After:**

```dockerfile
COPY requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt && rm -f /tmp/requirements.txt

# reComBat, alone, past two stale declarations in its own metadata: `sklearn = "^0.0"`
# (the deprecated stub, which now only raises) and `python = ">=3.8,<3.11"` (it runs fine
# on 3.11 - measured). --no-deps is a global pip flag with no per-requirement form, which
# is why this cannot live in requirements.txt. Its real imports - numpy, pandas,
# scikit-learn, tqdm - are all pinned in the layer above.
#
# Pinned to a commit, and deliberately NOT to the PyPI wheel: the published 0.1.4 defaults
# to model='linear' where this commit defaults to model='elastic_net', so the wheel would
# silently change every 30_recombat result and break comparison with the published
# dissertation metrics.
ARG RECOMBAT_COMMIT=b02bd025c9671e75f37dc513d7f466b05448364d
RUN pip install --no-cache-dir --no-deps --ignore-requires-python \
        "reComBat @ git+https://github.com/BorgwardtLab/reComBat@${RECOMBAT_COMMIT}" \
    && python -c "from reComBat import reComBat; print('reComBat import OK')"
```

Its own layer also means a network failure here costs seconds, not the whole Python layer —
the same reasoning as the staged R install.

### 3. `pyproject.toml` — remove reComBat from the `methods` extra

```toml
    "reComBat @ git+https://github.com/BorgwardtLab/reComBat",
```

`pip install -e ".[methods]"` and `".[all]"` fail today for exactly the reason above, and a
dependency needing `--no-deps` cannot be expressed in a dependency list at all. It is
removed, `tqdm>=4.62.3` is added to `methods`, and a comment records the manual command.

This is why the laptop suite never noticed: `dev` does not include `methods`, so `.venv` has
no reComBat and `30_recombat` skips locally.

### 4. `combobatch/methods/python_methods.py` — the install hint becomes wrong

`_require("reComBat", "30_recombat", "methods")` produces:

> `30_recombat needs the Python package 'reComBat' … Install it with pip install 'combobatch[methods]'.`

After section 3 that command will no longer install it — a message confidently naming a fix
that does not work, which is the defect class this project keeps correcting. `_require`
gains an optional `hint` parameter for the one package that needs a non-standard command:

```python
def _require(
    package: str, method_name: str, extra: str, hint: str | None = None
) -> None:
```

with `normalize_recombat` passing the real command.

### 5. `tests/unit/`

| Test | Guards against |
|---|---|
| `test_docker_assets.py:125` — currently `assert "reComBat" in requirements.txt` | **Must be rewritten**: it will fail once the line moves. Re-point it at the Dockerfile's dedicated install step, keeping its intent — reComBat arrives through pip, never through the R script |
| the reComBat install pins a commit and uses `--no-deps --ignore-requires-python` | a rebuild silently picking up a moved `master`, or the flags being dropped and the build failing again |
| `tqdm` is declared in `requirements.txt` | reComBat's only real dependency reverting to a transitive accident of `umap-learn` |
| **new** — every constraint shared by `requirements.txt` and `pyproject.toml` matches | makes the header's claim in 1b true instead of aspirational |

### 6. `docker/BUILD_AND_PUSH.md`

Troubleshooting rows for `The 'sklearn' PyPI package is deprecated` and for
`requires a different Python: 3.11 … not in '>=3.8,<3.11'`, both pointing at the dedicated
install step and warning against "fixing" it by adding `sklearn` to `requirements.txt`.

### 7. `docker/CLAUDE.md`

One invariant: **reComBat is installed alone, with `--no-deps --ignore-requires-python`, at
a pinned commit.** Never move it into `requirements.txt`, never switch it to the PyPI wheel
(different default model), never add `sklearn` or `fire` to satisfy its metadata.

### 8. `combobatch_tool_plan_260828.md`

One clause under Phase 6, naming this as the fourth diagnosis.

---

## Adjacent finding — `30_recombat` has no hyperparameters

**Not part of the build fix. Flagged rather than silently fixed, because it is a real defect
and its scope is yours to decide.**

`normalize_recombat` calls `reComBat()` bare and ignores `**kw`, and its `MethodSpec`
declares no `hyperparams`. The constructor accepts:

```python
def __init__(self, parametric=True, model='elastic_net', config=None,
             conv_criterion=1e-4, max_iter=1000, n_jobs=None,
             mean_only=False, optimize_params=True, reference_batch=None, verbose=True):
```

Three consequences:

- **`model` is unreachable**, so the `elastic_net`/`linear` difference this plan turns on
  cannot be pinned from a config — the only way to change it is to change which source the
  image installs.
- **`reference_batch` is unreachable**, although this project has `resolve_reference_batch()`
  and `uses_reference_batch` metadata precisely for it.
- **`n_jobs=None` becomes `cpu_count()`** inside every dispatcher worker. That is the
  "Shambhala parallelism multiplies" problem again — outer workers × inner threads — against
  a `CLAUDE.md` invariant that thread counts are pinned in the Python entry point.

`combobatch/methods/CLAUDE.md` states the rule this breaks: *"A parameter the function
accepts but does not declare is invisible and will never be reachable — that is precisely
the bug this package exists to fix."*

Fixing it is roughly: declare `model`, `parametric`, `max_iter`, `n_jobs`, `reference_batch`
as `HyperParam`s, forward them, default `n_jobs=1`, and let `docs/HYPERPARAMETERS.md`
regenerate. `tests/unit/test_hyperparams.py` would gain a case asserting that changing
`model` changes the output. Say the word and it becomes its own plan.

---

## Files that do NOT need to change

| File | Why not |
|---|---|
| `docker/install_r_packages.R`, the R stages | Fully cached and correct in this build; the failure is three layers later |
| the apt layer | `git` is still needed — `remotes::install_github` uses it for five R packages |
| `combobatch/methods/__init__.py` | Only if the adjacent finding above is taken up |
| `scripts/smoke_test.sh`, `scripts/build_and_push_image.sh` | Not reached, and unaffected by where pip gets one package |
| the three registries | `30_recombat` is correctly registered; the package was unbuildable, not misdeclared |

---

## Side effects and caveats

- **The Python layers rebuild; the R layers do not.** Everything up to and including Octave
  stays cached, so the next attempt starts at `pip install` — a few minutes, not an hour.
- **Pinning the commit changes nothing today and everything later.** It is the same commit
  the last three builds resolved to; it stops the fourth from quietly getting a different one.
- **`--ignore-requires-python` is a real override, and it is load-bearing here.** It is
  justified by the measurement above, not by optimism. If reComBat ever gains a genuine
  3.11 incompatibility, this flag will hide it — which is why the Dockerfile step ends with
  an import check, and why the acceptance suite exercises `30_recombat` inside the image.
- **`pip install -e ".[all]"` starts working again**, having been broken for as long as the
  `sklearn` stub has been raising.
- **`tqdm` becomes explicit.** It was already installed via `umap-learn`; nothing new lands
  in the image.
- **This does not touch the numeric behaviour of `30_recombat`** — that is the entire reason
  for staying on the git commit rather than taking the faster wheel.

---

## Verification commands

```bash
# 1. Reproduce the failure as it stands (expect metadata-generation-failed on sklearn)
docker run --rm python:3.11-slim-bookworm sh -c \
    'apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq --no-install-recommends git >/dev/null 2>&1; \
     pip install --no-cache-dir -q git+https://github.com/BorgwardtLab/reComBat'

# 2. The fix installs and imports
docker run --rm python:3.11-slim-bookworm sh -c \
    'apt-get update -qq >/dev/null 2>&1; apt-get install -y -qq --no-install-recommends git >/dev/null 2>&1; \
     pip install -q "numpy<2.0" "pandas<2.0" "scikit-learn<1.6" tqdm; \
     pip install -q --no-deps --ignore-requires-python \
       "reComBat @ git+https://github.com/BorgwardtLab/reComBat@b02bd025c9671e75f37dc513d7f466b05448364d"; \
     python -c "from reComBat import reComBat; print(\"import OK\")"'

# 3. The default model is still elastic_net — the reason for not taking the PyPI wheel
docker run --rm python:3.11-slim-bookworm sh -c \
    '…install as above…; python -c "from reComBat import reComBat; print(reComBat().model)"'
# expect: elastic_net

# 4. Laptop suite and lint
python -m pytest tests/unit -q
black --check .
docker buildx build --check -f docker/Dockerfile .

# 5. `[all]` installs again
python3.11 -m venv /tmp/vcheck && /tmp/vcheck/bin/pip install -q -e ".[all]" && echo "all extra OK"

# 6. The build, VPN off
bash scripts/smoke_test.sh

# 7. In the built image, the method actually runs
docker run --rm --platform linux/amd64 combobatch:test \
    pytest tests/integration/test_methods_all_backends.py -q -k recombat
```

---

## Two deviations found while implementing

**`awscli` is in the image and in no pyproject extra.** The new drift test found it
immediately. It is not a defect: `awscli` is a command-line tool for working inside the
pod, not a library this package imports, so it has no business in an extra a laptop install
would resolve. Recorded as a named exception with that reason, rather than silenced.

**The first version of the commit-pin test could not see the flags.** It collected
Dockerfile lines containing `reComBat`, but `--no-deps --ignore-requires-python` sits on the
`RUN pip install` line while the package name is on its continuation — so the assertion
failed against a correct Dockerfile. Rewritten to capture the whole block from
`ARG RECOMBAT_COMMIT` to the blank line.

---

## TODO

- [x] `requirements.txt` — replace the git line with the pointer comment; add `tqdm>=4.62.3`
- [x] `requirements.txt` — correct the header's claim about the pyproject drift test
      (now names `test_requirements_match_pyproject.py`, which exists)
- [x] `docker/Dockerfile` — dedicated reComBat layer, commit-pinned, with an import check
- [x] `pyproject.toml` — drop reComBat from the `methods` extra; add `tqdm`; comment the
      manual command
- [x] `combobatch/methods/python_methods.py` — `_require` gains `hint`; `normalize_recombat`
      passes the real install command
- [x] `tests/unit/test_docker_assets.py` — rewrite the `reComBat in requirements.txt`
      assertion around the new install step
- [x] `tests/unit/test_docker_assets.py` — assert the commit pin and both pip flags
      (six assertions in `TestReComBatInstall`, including "not from PyPI" and "no sklearn
      stub", the two ways someone would plausibly "fix" this and break it)
- [x] `tests/unit/test_docker_assets.py` — assert `tqdm` is declared
- [x] `tests/unit/test_requirements_match_pyproject.py` — new: shared constraints agree
- [x] `docker/BUILD_AND_PUSH.md` — two troubleshooting rows
- [x] `docker/CLAUDE.md` — the reComBat install invariant
- [x] `combobatch_tool_plan_260828.md` — name this plan under Phase 6
- [x] Run verification steps 1–5 — the pinned commit installs, imports, and reports
      `default model: elastic_net`; `.[all]` resolves again; 860 unit tests pass; Black
      clean; `buildx --check` clean
- [ ] **VPN off**, run `bash scripts/smoke_test.sh`
- [ ] Run verification step 7 against the built image
- [ ] *(separate decision)* the `30_recombat` hyperparameters finding above
