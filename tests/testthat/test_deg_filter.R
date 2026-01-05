deg_env <- new.env(parent = baseenv())
sys.source(testthat::test_path("..", "..", "R", "deg.R"), envir = deg_env)
select_genes_for_deg <- deg_env$select_genes_for_deg

test_that("select_genes_for_deg supports global filtering and rejects unsupported modes", {
  expr <- Matrix::Matrix(c(
    1, 1, 0, 0, # gene1: detected in both control pseudobulks
    0, 0, 1, 0, # gene2: detected only in one perturb pseudobulk
    0, 0, 0, 0  # gene3: undetected everywhere
  ), nrow = 3, byrow = TRUE, sparse = TRUE)

  perturbation_labels <- c("Control", "Control", "Perturb", "Perturb")

  global_res <- select_genes_for_deg(
    expr_over_thresh = expr,
    detection_mat = expr > 0,
    perturbation_labels = perturbation_labels,
    detection_min_samples = 2,
    filter_mode = "global"
  )
  expect_equal(global_res$genes_to_test, 1)
  expect_equal(global_res$stats$n_genes_all_zero_both, 1)

  expect_error(
    select_genes_for_deg(
      expr_over_thresh = expr,
      detection_mat = expr > 0,
      perturbation_labels = perturbation_labels,
      detection_min_samples = 1,
      filter_mode = "by_condition"
    ),
    "deg_filtering_mode must be one of: global, by_background_count"
  )
})
