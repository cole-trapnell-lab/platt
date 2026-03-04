# Codex instructions (Platt)

## What this repo is
Platt is an R package used in the Trapnell Lab stack for statistical modeling and analysis of single-cell perturbation experiments.

## Cost & safety guardrails
- Safe to run R package unit tests and lightweight checks.
- Do NOT run long simulations, benchmarks, or analyses on full datasets unless explicitly requested.
- Do NOT download large datasets or install system-level dependencies automatically.

## Change expectations
- Prefer small, reviewable PRs.
- If changing exported functions or model behavior:
  - update roxygen documentation
  - add or update unit tests
  - update NEWS.md if present
- Avoid breaking changes unless explicitly requested.
- Use `make fast` for quick validation during development and `make test` before opening a PR.

## Fast-test conventions
- `make fast` runs the same `testthat::test_local()` invocation as `make test` but exports `PLATT_FAST=1` so tests can detect the fast tier.
- Guard expensive test blocks with `if (identical(Sys.getenv("PLATT_FAST"), "1")) testthat::skip("slow")` (or similar) so `make fast` remains cheap and deterministic.
- If `tests/testthat` is missing or contains no `test*.R` files, `make fast`/`make test` quietly skip so empty packages do not fail.

## Validation tiers
### FAST (default)
- `make fast` (sets `PLATT_FAST=1`)
- Runs cheap, deterministic unit tests only
- Should finish in minutes or less

### SMOKE (opt-in)
- Small toy examples or vignettes using synthetic data

### FULL (manual only)
- Full simulations, benchmarks, or real dataset runs
