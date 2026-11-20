# Knockout-vs-control fit (compare_gene_expression_within_node): closed form vs monocle3::fit_models.

make_knockout_pb_cds <- function(ng = 150, n_ctrl = 24, n_ko = 8, seed = 2) {
  set.seed(seed)
  pert <- c(rep("Control", n_ctrl), rep("KO", n_ko))
  ns <- length(pert)
  cells <- round(runif(ns, 10, 80))
  mu <- outer(exp(runif(ng, -6, 0)), cells) * exp(rnorm(ng * ns, 0, 0.4))
  mu[1:10, pert == "KO"] <- mu[1:10, pert == "KO"] * 4                 # real effects
  Y <- matrix(rpois(ng * ns, mu), ng, ns, dimnames = list(paste0("g", 1:ng), paste0("pb", 1:ns)))
  Y[11:16, pert == "KO"] <- 0                                          # lost in the knockout
  Y[17:20, pert == "Control"] <- 0                                     # gained from zero
  suppressWarnings(monocle3::new_cell_data_set(methods::as(Y, "CsparseMatrix"),
    cell_metadata = data.frame(cell_group = "ct1", perturbation = pert, num_cells_in_group = cells,
                               timepoint = rep(c(24, 36, 48, 72), length.out = ns), row.names = colnames(Y)),
    gene_metadata = data.frame(id = rownames(Y), gene_short_name = rownames(Y), row.names = rownames(Y))))
}

run_knockout_fit <- function(pb_cds, ...) {
  out <- suppressWarnings(suppressMessages(utils::capture.output(res <- compare_gene_expression_within_node(
    cell_group = "ct1", ccs = NULL, pb_cds = pb_cds, control_ids = "Control", perturbation_ids = "KO",
    cv_threshold = 100, cores = 1, ...))))
  tidyr::unnest(res, perturb_effects) |> dplyr::arrange(id)
}

test_that("knockout closed-form fit matches the monocle fit", {
  skip_if_not_installed("monocle3")
  pb <- make_knockout_pb_cds()
  ref <- run_knockout_fit(pb, fit_method = "monocle")
  cf <- run_knockout_fit(pb, fit_method = "closed_form")
  expect_equal(cf$id, ref$id)
  expect_equal(cf$perturb_to_ctrl_raw_lfc, ref$perturb_to_ctrl_raw_lfc, tolerance = 1e-6)
  # standard errors differ only through the dispersion trend (speedglm stops at acc = 0.001)
  expect_equal(cf$perturb_to_ctrl_raw_lfc_se, ref$perturb_to_ctrl_raw_lfc_se, tolerance = 1e-2)
  expect_equal(cf$log_mean_expression, ref$log_mean_expression, tolerance = 1e-6)
  # genes with no counts in one arm take the same thresholded branch in contrast_helper
  zero_arm <- cf$id %in% paste0("g", 11:20)
  expect_true(all(is.finite(cf$perturb_to_ctrl_raw_lfc_se[zero_arm])))
  expect_equal(sign(cf$perturb_to_ctrl_raw_lfc[zero_arm]), sign(ref$perturb_to_ctrl_raw_lfc[zero_arm]))
})

test_that("closed form ignores max_simultaneous_genes and fits the trend once", {
  skip_if_not_installed("monocle3")
  pb <- make_knockout_pb_cds()
  one <- run_knockout_fit(pb, fit_method = "closed_form")
  blocks <- run_knockout_fit(pb, fit_method = "closed_form", max_simultaneous_genes = 40)
  expect_identical(blocks$perturb_to_ctrl_raw_lfc_se, one$perturb_to_ctrl_raw_lfc_se)
})

test_that("a nuisance term keeps the monocle fit", {
  skip_if_not_installed("monocle3")
  pb <- make_knockout_pb_cds(ng = 60)
  ref <- run_knockout_fit(pb, fit_method = "monocle", nuisance_model_formula_str = "timepoint")
  cf <- run_knockout_fit(pb, fit_method = "closed_form", nuisance_model_formula_str = "timepoint")
  expect_identical(cf$perturb_to_ctrl_raw_lfc, ref$perturb_to_ctrl_raw_lfc)
  expect_identical(cf$perturb_to_ctrl_raw_lfc_se, ref$perturb_to_ctrl_raw_lfc_se)
})
