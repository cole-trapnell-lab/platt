## Testing quickstart

Scope is intentionally minimal: only small unit tests guarding the gene-filtering helper used by DEG profiling.

- Run tests: `R -q -e "devtools::test()"` (or `Rscript -e "testthat::test_dir('tests/testthat')"`).
- `R CMD check` will also pick up the tests via `tests/testthat.R`.
- Tests use tiny synthetic `dgCMatrix` inputs and finish in under a couple seconds; no large fixtures or end-to-end DEG runs are included.
