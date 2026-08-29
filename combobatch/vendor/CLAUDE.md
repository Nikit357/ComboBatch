# CLAUDE.md — `combobatch/vendor/`

## What this subpackage owns

Copies of third-party code, kept as close to their upstream form as possible.

## Invariants

- **Do not refactor vendored code to match this project's style.** Renaming, reformatting
  or "improving" it makes the next diff against upstream unreadable, which is the entire
  cost this directory pays for the benefit of being able to fix things at all.
- **Every deliberate change carries a `VENDOR FIX` comment at the site**, naming what was
  wrong. A change without one is indistinguishable from an upstream quirk.
- **Every vendored package's `README.md` lists its provenance and its complete set of
  changes.** If you fix something here, add it there in the same commit.
- **Adapt at the boundary, not inside.** Orientation flips, parameter resolution and
  configuration belong in the glue that calls this code — for Shambhala, that is
  `combobatch/methods/shambhala_method.py`.
- **Licence compatibility is a precondition for vendoring at all.** `39_procrustes` is not
  here for exactly this reason: BostonGene proprietary, non-commercial use only, and its
  licence forbids altering copyright notices. Nothing whose licence is incompatible with
  MIT may be added to this directory.
