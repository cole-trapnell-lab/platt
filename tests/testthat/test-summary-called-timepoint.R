# dacts_when_abundant(): a called timepoint beats an uncalled one; otherwise the
# minimum q decides, as before.

# Two timepoints: 24 hpf is very significant but small (codes A0), 48 hpf is
# less significant but passes the call rule (codes A2).
.hidden_call_dact <- function() {
  tibble::tibble(
    cell_group           = "ct",
    timepoint_x          = c(24, 48),
    present_above_thresh = TRUE,
    percent_max_abund    = c(0.9, 0.5),
    delta_log_abund      = c(0.2, -1.2),
    delta_log_abund_se   = c(0.03, 0.30),
    delta_q_value        = c(1e-6, 0.01),
    df_resid             = 40
  )
}

.code_at <- function(row) {
  assign_abundance_code(row$delta_log_abund, row$delta_q_value,
    se = row$delta_log_abund_se, df = row$df_resid
  )
}

test_that("the fixture's minimum-q row is A0 and the other row is called", {
  expect_equal(.code_at(.hidden_call_dact()), c("A0 No change", "A2 Depletion"))
})

test_that("a called timepoint wins over a more significant uncalled one", {
  got <- dacts_when_abundant(.hidden_call_dact(), percent_max_thresh = 0)
  expect_equal(nrow(got), 1L)
  expect_equal(got$timepoint_x, 48)
  expect_equal(.code_at(got), "A2 Depletion")
  expect_identical(got, dacts_when_abundant(.hidden_call_dact()[2:1, ], percent_max_thresh = 0))
})

test_that("among called timepoints the minimum q still wins, not the larger change", {
  d <- .hidden_call_dact()
  d$delta_log_abund <- c(-0.8, -2.5)
  d$delta_q_value <- c(1e-6, 0.02)
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x, 24)
  expect_equal(.code_at(got), "A2 Depletion")
})

test_that("an expansion and a depletion: the more significant call wins", {
  d <- .hidden_call_dact()
  d$delta_log_abund <- c(0.9, -0.6)
  d$delta_q_value <- c(0.03, 0.001)
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x, 48)
  expect_equal(.code_at(got), "A2 Depletion")
})

test_that("called rows tied on q go to peak wild-type abundance", {
  d <- .hidden_call_dact()
  d$delta_log_abund <- c(-1.2, -0.9)
  d$delta_q_value <- 0.01
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x, 24)
  ties <- dacts_when_abundant(d, percent_max_thresh = 0, with_ties = TRUE)
  expect_equal(sort(ties$timepoint_x), c(24, 48))
})

test_that("with_ties = TRUE only returns called rows when there is a call", {
  ties <- dacts_when_abundant(.hidden_call_dact(), percent_max_thresh = 0, with_ties = TRUE)
  expect_equal(ties$timepoint_x, 48)
})

test_that("with no called timepoint the minimum q is still kept", {
  d <- .hidden_call_dact()
  d$delta_log_abund <- c(0.2, -0.3)
  d$delta_q_value <- c(1e-6, 0.2)
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x, 24)
  expect_equal(.code_at(got), "A0 No change")
})

test_that("a call where the cell type is rare in wild type is not used", {
  d <- .hidden_call_dact()
  d$percent_max_abund <- c(0.9, 0.05)
  got <- dacts_when_abundant(d) # default percent_max_thresh = 0.10 neutralises 48 hpf
  expect_equal(got$timepoint_x, 24)
})

test_that("a call at a stage with no cells on either side is not used", {
  d <- .hidden_call_dact()
  d$n_cells_x <- c(50, 0)
  d$n_cells_y <- c(40, 0)
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x, 24)
})

test_that("each cell type follows its own rule, one row each, in cell_group order", {
  other <- .hidden_call_dact()
  other$cell_group <- "a_ct"
  other$delta_log_abund <- c(0.2, -0.3)
  other$delta_q_value <- c(1e-6, 0.2)
  d <- dplyr::bind_rows(.hidden_call_dact(), other)
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$cell_group, c("a_ct", "ct"))
  expect_equal(got$timepoint_x, c(24, 48))
  expect_identical(got, dacts_when_abundant(d[c(4, 1, 3, 2), ], percent_max_thresh = 0))
})

test_that("tables without se/df still find the called timepoint", {
  d <- .hidden_call_dact()[, setdiff(names(.hidden_call_dact()), c("delta_log_abund_se", "df_resid"))]
  got <- dacts_when_abundant(d, percent_max_thresh = 0)
  expect_equal(got$timepoint_x, 48)
})
