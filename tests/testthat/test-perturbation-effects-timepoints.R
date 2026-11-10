# get_perturbation_effects(): which stages get a contrast (mcclintock#63) and
# what the stage column is called (platt#18). make_contrast() is mocked to echo
# the newdata row it is asked for, so no model is fitted.

setClass("fake_pe_ccm", representation(ccs = "ANY"))

.fake_ccm <- function(stages, interval_col = "timepoint") {
  cd <- S4Vectors::DataFrame(stats::setNames(list(stages), interval_col))
  new("fake_pe_ccm", ccs = SummarizedExperiment::SummarizedExperiment(colData = cd))
}

.echo_contrast <- function(ccm, newdata = tibble::tibble(), adjust_q_values = FALSE) {
  tibble::tibble(cell_group = c("a", "b"), seen = as.numeric(newdata[[1]][1]))
}

test_that("by default every stage in the cell count set gets a contrast", {
  local_mocked_bindings(make_contrast = .echo_contrast)
  # 36 hpf sampled only for controls: still contrasted, as before
  out <- get_perturbation_effects(.fake_ccm(c(24, 24, 36, 48, 72)))
  expect_setequal(unique(out$timepoint), c(24, 36, 48, 72))
  expect_equal(out$seen, out$timepoint)
})

test_that("timepoints restricts contrasts to the given stages", {
  local_mocked_bindings(make_contrast = .echo_contrast)
  out <- get_perturbation_effects(.fake_ccm(c(24, 36, 48, 72)), timepoints = c(24, 48, 72))
  expect_setequal(unique(out$timepoint), c(24, 48, 72))
  expect_equal(nrow(out), 3L * 2L)
})

test_that("requested stages missing from the cell count set are skipped with a warning", {
  local_mocked_bindings(make_contrast = .echo_contrast)
  expect_warning(
    out <- get_perturbation_effects(.fake_ccm(c(24, 48)), timepoints = c(24, 48, 96)),
    "96"
  )
  expect_setequal(unique(out$timepoint), c(24, 48))
})

test_that("the stage column follows interval_col (platt#18)", {
  local_mocked_bindings(make_contrast = .echo_contrast)
  out <- get_perturbation_effects(.fake_ccm(c(11, 48), "timepointNumeric"), interval_col = "timepointNumeric")
  expect_true("timepointNumeric" %in% names(out))
  expect_false("timepoint" %in% names(out))
  expect_setequal(out$seen, c(11, 48))
})
