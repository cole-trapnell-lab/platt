SHELL := /bin/bash
RSCRIPT ?= Rscript

.PHONY: help fast test check

help:
	@echo "Targets:"
	@echo "  fast   Run cheap unit tests"
	@echo "  test   Run testthat::test_local"
	@echo "  check  Run R CMD check (heavier)"

fast: test

test:
	@$(RSCRIPT) -e 'if (!requireNamespace("testthat", quietly=TRUE)) {message("testthat not installed; skipping"); quit(status=0)}; testthat::test_local(reporter="summary", stop_on_failure=TRUE)'

check:
	@R CMD check . --no-manual
