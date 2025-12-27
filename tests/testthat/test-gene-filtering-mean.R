deg_env <- new.env(parent = baseenv())
sys.source(testthat::test_path("..", "..", "R", "deg.R"), envir = deg_env)
genes_to_test_by_mean <- deg_env$genes_to_test_by_mean

conditions <- c("control", "control", "perturb", "perturb")

make_mean_matrix <- function() {
  Matrix::sparseMatrix(
    i = c(
      2, 2,            # g2: control only
      3, 4             # g3/g4: perturb only vs both
    ),
    j = c(
      1, 2,            # control samples
      3, 4             # perturb samples
    ),
    x = c(5, 5, 10, 1),  # g3 high in perturb, g4 low
    dims = c(4, 4),
    dimnames = list(paste0("g", 1:4), NULL)
  )
}

test_that("by_mean_expression keeps genes with high mean in one arm", {
  expr <- make_mean_matrix()
  genes <- genes_to_test_by_mean(expr, conditions, min_mean = 3, condition_min_samples_detected = 0)
  expect_setequal(rownames(expr)[genes], c("g2", "g3"))
})

test_that("by_mean_expression drops genes below mean in both arms", {
  expr <- make_mean_matrix()
  genes <- genes_to_test_by_mean(expr, conditions, min_mean = 6, condition_min_samples_detected = 0)
  expect_false("g3" %in% rownames(expr)[genes]) # perturb mean 2.75 < 6
  expect_false("g2" %in% rownames(expr)[genes]) # control mean is 5 < 6
})

test_that("by_mean_expression honors optional detection K", {
  expr <- make_mean_matrix()
  # make g3 detected in only one perturb sample by zeroing one value
  expr[3, 4] <- 0
  genes <- genes_to_test_by_mean(expr, conditions, min_mean = 3, condition_min_samples_detected = 2)
  expect_false("g3" %in% rownames(expr)[genes]) # fails detection guard
  expect_true("g2" %in% rownames(expr)[genes])  # passes detection in control
})
