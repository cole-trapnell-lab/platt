test_that("compute_background_thresholds returns thresholds without densifying", {
  expr <- Matrix::Matrix(c(
    10, 0, 0, 0,
    1,  1, 1, 1,
    0.1,0.1,0.1,0.1
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  colnames(expr) <- paste0("pb", 1:4)
  rownames(expr) <- c("marker", "house", "soup")
  cell_types <- c("ct1", "ct1", "ct2", "ct2")
  perturb <- c("Control", "Treat", "Control", "Treat")
  res <- platt:::compute_background_thresholds(
    expr_mat = expr,
    cell_types = cell_types,
    perturbation_labels = perturb,
    bottom_frac = 0.5,
    threshold_type = "add",
    delta = 0.25,
    mult = 2
  )
  expect_length(res$thresholds, 3)
  expect_true(res$thresholds["marker"] < res$thresholds["house"])
  expect_true(res$thresholds["soup"] < res$thresholds["house"])
})

test_that("count_above_threshold_by_mask excludes genes below background", {
  expr <- Matrix::Matrix(c(
    2, 0,
    0, 2,
    0.1, 0.1
  ), nrow = 3, byrow = TRUE, sparse = TRUE)
  colnames(expr) <- c("c1", "p1")
  thresholds <- c(0.5, 0.5, 0.5)
  ctrl_mask <- c(TRUE, FALSE)
  pert_mask <- c(FALSE, TRUE)
  ctrl_counts <- platt:::count_above_threshold_by_mask(expr, thresholds, ctrl_mask)
  pert_counts <- platt:::count_above_threshold_by_mask(expr, thresholds, pert_mask)
  expect_equal(ctrl_counts, c(1, 0, 0))
  expect_equal(pert_counts, c(0, 1, 0))
  keep_k1 <- (ctrl_counts >= 1) | (pert_counts >= 1)
  keep_k2 <- (ctrl_counts >= 2) | (pert_counts >= 2)
  expect_identical(which(keep_k1), c(1L, 2L))
  expect_identical(length(which(keep_k2)), 0L)
})
