# platt

R package: lineage graph construction from time-series + perturbation scRNA-seq. Also carries the
DE / empirical-FDR gate and phenotype-assignment code that mcclintock and the portal call.
`Depends: hooke` — so platt sits **downstream** of hooke.

## Commands (verified in Makefile)
```bash
make fast     # PLATT_FAST=1 Rscript -e 'testthat::test_local(reporter="summary", stop_on_failure=TRUE)'
make test     # same invocation without PLATT_FAST
make check    # R CMD check . --no-manual
```
11 test files under `tests/testthat/` (efdr gates, empirical FDR, background counts, gene filtering,
gsea enrichment, phenotype-assignment utils, …).

## Gotchas
- **`make fast` and `make test` are behaviourally identical.** `PLATT_FAST` is exported by the
  Makefile and read by **zero** test files. Don't rely on it to skip slow paths; if you want that
  tier, you have to add the guards first.
- **A green `make fast` can mean "ran nothing".** Both targets `quit(status=0)` if `PLNmodels`
  *or* `hooke` is unavailable, and again if `tests/testthat` has no `test*.R`. Read the output.
- `Remotes: cole-trapnell-lab/speedglm, bioc::fgsea` — neither resolves from CRAN.
- No CI in this repo. `make check` before a PR is the only gate.
- No `contract` target (root `make fast` does not require one for platt).
- `empirical_p` is the gate, not `empirical_fdr`.

## Release notes
No `NEWS.md` yet. If you bump `Version:` in `DESCRIPTION`, create one in the same commit —
see the root `CLAUDE.md` convention and `repos/monocle3/NEWS.md` for the format.

## Architecture
`R/` only. `Depends: Biobase, monocle3, PLNmodels, hooke`; roxygen 7.3.3, testthat edition 3.
