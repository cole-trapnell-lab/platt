# Stages where neither control nor perturbed embryos recovered a cell
# (mcclintock#134): the contrast there is interpolation on both sides.

test_that("a row with no cells observed on either side is AN, never a call", {
  # sox10-like: a -6.9 loss at q = 1e-28 with 0 cells in 24 control + 8 perturbed embryos
  expect_equal(assign_abundance_code(-6.9, 1e-28, se = 0.3, df = 40, n_observed = 0), "AN Not assessed")
  expect_equal(assign_abundance_code(-6.9, 1e-28, se = 0.3, df = 40, n_observed = 5), "A3 Near-loss")
  # absent counts (older tables) leave the code as before
  expect_equal(assign_abundance_code(-6.9, 1e-28, se = 0.3, df = 40), "A3 Near-loss")
  expect_equal(assign_abundance_code(-6.9, 1e-28, se = 0.3, df = 40, n_observed = NA), "A3 Near-loss")
  # vectorised
  expect_equal(
    assign_abundance_code(c(-6.9, 1.2), c(1e-28, 1e-6), se = 0.3, df = 40, n_observed = c(0, 12)),
    c("AN Not assessed", "A1 Expansion")
  )
})

test_that("a no-cells row is never graded", {
  expect_equal(assign_abundance_severity(-6.9, 1e-28, n_observed = 0), "none")
  expect_equal(assign_abundance_severity(-6.9, 1e-28, n_observed = 3), "severe")
  expect_equal(assign_abundance_severity(-6.9, 1e-28), "severe")
})

.counted_dact <- function() {
  tibble::tibble(
    cell_group           = c("ct", "ct", "gone"),
    timepoint_x          = c(18, 24, 24),
    present_above_thresh = TRUE,
    percent_max_abund    = c(0.5, 0.9, 0.8),
    delta_log_abund      = c(-6.9, -0.8, -4.0),
    delta_log_abund_se   = 0.3,
    delta_q_value        = c(1e-28, 0.01, 1e-10),
    df_resid             = 40,
    n_cells_x            = c(0, 40, 0),
    n_cells_y            = c(0, 15, 0)
  )
}

test_that("dacts_when_abundant picks the code row from an observed stage", {
  got <- dacts_when_abundant(.counted_dact(), percent_max_thresh = 0)
  # the 18 hpf row (min q, no cells) is skipped; 24 hpf, which has cells, is kept
  expect_equal(got$timepoint_x[got$cell_group == "ct"], 24)
  # a cell type observed at no stage is dropped, i.e. AN via the cell-type universe
  expect_false("gone" %in% got$cell_group)
})

test_that("without count columns dacts_when_abundant is unchanged", {
  d <- .counted_dact()[, setdiff(names(.counted_dact()), c("n_cells_x", "n_cells_y"))]
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x[got$cell_group == "ct"], 18)
  expect_true("gone" %in% got$cell_group)
})

test_that("the alias carries the observed count to the code assignment", {
  a <- alias_change_when_present(.counted_dact())
  expect_equal(a$change_when_present_n_observed, c(0, 55, 0))
  expect_false("change_when_present_n_observed" %in%
    names(alias_change_when_present(.counted_dact()[, 1:8])))
})
