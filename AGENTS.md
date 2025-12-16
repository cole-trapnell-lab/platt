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

## Validation tiers
### FAST (default)
- `make fast`
- Runs cheap, deterministic unit tests only
- Should finish in minutes or less

### SMOKE (opt-in)
- Small toy examples or vignettes using synthetic data

### FULL (manual only)
- Full simulations, benchmarks, or real dataset runs
