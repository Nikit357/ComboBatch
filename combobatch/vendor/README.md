# Vendored code

Third-party implementations copied into this repository rather than depended on, because
the upstream is not installable, not versioned, or needed a fix that could not wait for a
release.

| Package | Upstream | Why it is vendored |
|---|---|---|
| `shambhala/` | the Shambhala2 pure-Python + Octave implementation | Not published to PyPI; needed three defect fixes; ships `.m` files that a normal install would drop. |

Each vendored package documents its own provenance and every change made to it in its
`README.md`, and marks each change at the site with a `VENDOR FIX` comment.
