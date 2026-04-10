## Vendoring notes (RepOrtR)

This directory contains a vendored copy of upstream **satMiner**.

- Upstream: `https://github.com/fjruizruano/satminer`
- Revision pinned: `fd4f40cbd3b4926ff230f516bfc8e9c827bc146d`
- Vendored date: 2026-03-30

We do not run upstream scripts directly in the RepOrtR pipeline because many are Python2-era.
Instead, RepOrtR provides Python3 adapter scripts that replicate the required behavior while
keeping satMiner semantics.

