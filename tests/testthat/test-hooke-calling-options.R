# Opt-in hooke calling options (hooke >= 0.0.6): get_perturbation_effects() and
# make_contrast() pass them through to hooke::compare_abundances() untouched,
# and an empty q-value from a contrast that was not tested does not break the
# narrative. No model is fitted: the hooke calls are mocked.

setClass("fake_co_ccm", representation(ccs = "ANY"))

.fake_co_ccm <- function(stages) {
  cd <- S4Vectors::DataFrame(timepoint = stages)
  new("fake_co_ccm", ccs = SummarizedExperiment::SummarizedExperiment(colData = cd))
}

test_that("make_contrast() forwards extra arguments to compare_abundances()", {
  seen <- NULL
  local_mocked_bindings(
    estimate_abundances = function(ccm, newdata, ...) tibble::tibble(cell_group = "a", knockout = newdata$knockout[1]),
    compare_abundances = function(ccm, cond_x, cond_y, adjust_q_values = FALSE, ...) {
      seen <<- c(list(adjust_q_values = adjust_q_values), list(...))
      tibble::tibble(cell_group = "a")
    },
    .package = "hooke"
  )
  make_contrast(.fake_co_ccm(24))
  expect_equal(seen, list(adjust_q_values = FALSE))
  make_contrast(.fake_co_ccm(24), test = "lrt_diagonal", contrast_se = "covariance", support_by = c("knockout", "timepoint"))
  expect_equal(seen, list(adjust_q_values = FALSE, test = "lrt_diagonal", contrast_se = "covariance", support_by = c("knockout", "timepoint")))
})

test_that("get_perturbation_effects() forwards extra arguments at every stage", {
  calls <- list()
  local_mocked_bindings(make_contrast = function(ccm, newdata = tibble::tibble(), adjust_q_values = FALSE, ...) {
    calls[[length(calls) + 1]] <<- list(...)
    tibble::tibble(cell_group = "a")
  })
  get_perturbation_effects(.fake_co_ccm(c(24, 48)))
  expect_equal(calls, list(list(), list()))
  calls <- list()
  get_perturbation_effects(.fake_co_ccm(c(24, 48)), timepoints = c(24, 48), test = "lrt_diagonal")
  expect_equal(calls, list(list(test = "lrt_diagonal"), list(test = "lrt_diagonal")))
})

test_that("the defaults add no arguments", {
  expect_true("..." %in% names(formals(get_perturbation_effects)))
  expect_true("..." %in% names(formals(make_contrast)))
})

test_that("the narrative reports an untested contrast instead of failing on its empty q-value", {
  dact <- function(q) tibble::tibble(cell_group = "a", delta_q_value = q, delta_log_abund = -1, power_at_margin = 0.9)
  degs <- tibble::tibble(
    cell_group = character(), perturb_to_ctrl_p_value = numeric(), perturb_to_ctrl_effect = numeric(),
    gene_short_name = character(), id = character(), empirical_p = numeric(), power = numeric()
  )
  abundance <- function(q) {
    suppressMessages(summarize_cell_type_impact("a",
      perturbation_description = "x", target_gene_expression = NULL, dact_results = dact(q), degs = degs,
      ref_expression = NULL, ontologies = NULL, ai_notes_path = tempdir()
    ))$abundance
  }
  expect_equal(abundance(NA_real_), "unknown (not tested)")
  # tested rows read as before
  expect_equal(abundance(0.01), "depleted")
  expect_equal(abundance(0.5), "no_change")
})
