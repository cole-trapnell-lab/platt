test_that("background counts filter excludes soup singletons and keeps markers", {
  skip_on_cran()
  library(Matrix)

  counts <- Matrix::Matrix(
    rbind(
      c(0, 0, 1, 0),   # soup singles
      c(10, 0, 0, 0),  # ctrl-only marker
      c(0, 12, 0, 0),  # perturb-only marker
      c(5, 5, 4, 4)    # housekeeping
    ),
    sparse = TRUE
  )
  rownames(counts) <- paste0("g", seq_len(nrow(counts)))
  colnames(counts) <- c("ct1_ctrl", "ct1_pert", "ct2_ctrl", "ct2_pert")

  perturbation_labels <- c("Control", "myod1", "Control", "myod1")
  cell_types <- c("ct1", "ct1", "ct2", "ct2")
  lib_sizes <- Matrix::colSums(counts)

  res <- select_genes_for_deg(
    expr_over_thresh = counts, # not used for counts-based filtering
    detection_mat = counts > 0,
    perturbation_labels = perturbation_labels,
    min_samples_detected = 1,
    condition_min_samples_detected = 1,
    filter_mode = "by_background_counts",
    cell_types = cell_types,
    counts_mat = counts,
    library_sizes = lib_sizes,
    background_bottom_frac = 0.5,
    background_quantile_p = 0.99,
    background_count_floor = 2,
    background_min_samples_over_threshold = 1
  )

  kept_genes <- rownames(counts)[res$genes_to_test]
  expect_setequal(kept_genes, c("g2", "g3", "g4"))
  expect_equal(res$stats$n_genes_excluded_by_background_counts, 1)
})

test_that("background count floor guards against scattered single UMIs", {
  skip_on_cran()
  library(Matrix)

  counts <- Matrix::Matrix(
    rbind(
      c(1, 1),  # singletons
      c(5, 6)   # robust
    ),
    sparse = TRUE
  )
  rownames(counts) <- c("soup", "marker")
  perturbation_labels <- c("Control", "myod1")
  cell_types <- c("ct1", "ct1")
  lib_sizes <- Matrix::colSums(counts)

  res <- select_genes_for_deg(
    expr_over_thresh = counts,
    detection_mat = counts > 0,
    perturbation_labels = perturbation_labels,
    min_samples_detected = 1,
    condition_min_samples_detected = 1,
    filter_mode = "by_background_counts",
    cell_types = cell_types,
    counts_mat = counts,
    library_sizes = lib_sizes,
    background_bottom_frac = 1.0,
    background_quantile_p = 0.99,
    background_count_floor = 2,
    background_min_samples_over_threshold = 1
  )

  kept_genes <- rownames(counts)[res$genes_to_test]
  expect_setequal(kept_genes, "marker")
})
