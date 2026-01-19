SHELL := /bin/bash
RSCRIPT ?= Rscript
FAST_ENV_VAR ?= PLATT_FAST
TEST_DIR ?= tests/testthat
TEST_FILE_CHECK ?= find "$(TEST_DIR)" -maxdepth 1 -name 'test*.R' -print -quit | grep -q .

# Shared testthat invocation used by both fast and test targets
TESTTHAT_CMD = ok_pln <- tryCatch(requireNamespace("PLNmodels", quietly=TRUE), error=function(e) FALSE); ok_hooke <- tryCatch(requireNamespace("hooke", quietly=TRUE), error=function(e) FALSE); if (!isTRUE(ok_pln) || !isTRUE(ok_hooke)) { message("Skipping platt tests: PLNmodels/hooke not available"); quit(status=0) }; if (!requireNamespace("testthat", quietly=TRUE)) {message("testthat not installed; skipping"); quit(status=0)}; testthat::test_local(reporter="summary", stop_on_failure=TRUE)

.PHONY: help fast test check

help:
	@echo "Targets:"
	@echo "  fast   Run cheap unit tests (sets $$PLATT_FAST=1 so tests can skip heavy paths)"
	@echo "  test   Run testthat::test_local"
	@echo "  check  Run R CMD check (heavier)"

fast:
	@if [ ! -d "$(TEST_DIR)" ]; then \
		echo "No tests in $(TEST_DIR); skipping fast target."; \
	elif ! ($(TEST_FILE_CHECK)); then \
		echo "No testthat files found under $(TEST_DIR); skipping fast target."; \
	else \
		$(FAST_ENV_VAR)=1 $(RSCRIPT) -e '$(TESTTHAT_CMD)'; \
	fi

test:
	@if [ ! -d "$(TEST_DIR)" ]; then \
		echo "No tests in $(TEST_DIR); skipping test target."; \
	elif ! ($(TEST_FILE_CHECK)); then \
		echo "No testthat files found under $(TEST_DIR); skipping test target."; \
	else \
		$(RSCRIPT) -e '$(TESTTHAT_CMD)'; \
	fi

check:
	@R CMD check . --no-manual
