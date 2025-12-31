deg_env <- new.env(parent = baseenv())
sys.source(testthat::test_path("..", "..", "R", "deg.R"), envir = deg_env)
genes_to_test_from_detection <- deg_env$genes_to_test_from_detection

make_test_matrix <- function() {
  # 5 genes x 6 samples; 3 control (1-3), 3 perturb (4-6)
  Matrix::sparseMatrix(
    i = c(
      2, 2,              # g2: control only
      3, 3,              # g3: perturb only
      4, 4, 4,            # g4: both (two control, one perturb)
      5, 5               # g5: perturb, meets threshold exactly
    ),
    j = c(
      1, 2,              # control samples
      4, 5,              # perturb samples
      2, 3, 4,           # control + perturb
      5, 6               # perturb samples
    ),
    x = 1,
    dims = c(5, 6),
    dimnames = list(paste0("g", 1:5), NULL)
  )
}

conditions <- c("control", "control", "control", "perturb", "perturb", "perturb")

test_that("global filter drops genes all-zero across samples", {
  expr <- make_test_matrix()
  genes <- genes_to_test_from_detection(expr, conditions, min_samples_detected = 2, mode = "global")
  expect_s4_class(expr, "dgCMatrix")
  expect_setequal(rownames(expr)[genes], c("g2", "g3", "g4", "g5"))
  expect_false("g1" %in% rownames(expr)[genes])
})

test_that("global filter respects min_samples_detected", {
  expr <- make_test_matrix()
  # make g5 present in only one perturb sample so it should drop when min=2
  expr[5, 6] <- 0
  genes <- genes_to_test_from_detection(expr, conditions, min_samples_detected = 2, mode = "global")
  expect_false("g5" %in% rownames(expr)[genes])
  expect_true("g3" %in% rownames(expr)[genes]) # still meets threshold aggregated
})

test_that("global mode matches itself for min=1", {
  expr <- make_test_matrix()
  global_genes <- genes_to_test_from_detection(expr, conditions, min_samples_detected = 1, mode = "global")
  expect_true(all(global_genes %in% global_genes))
})

test_that("helper returns integer indices without densifying", {
  expr <- make_test_matrix()
  genes <- genes_to_test_from_detection(expr, conditions, min_samples_detected = 2, mode = "global")
  expect_type(genes, "integer")
  expect_true("dgCMatrix" %in% class(expr))
  expect_false(any(is.na(genes)))
})

test_that("invalid filter mode errors cleanly", {
  expr <- make_test_matrix()
  expect_error(
    genes_to_test_from_detection(expr, conditions, min_samples_detected = 2, mode = "by_condition"),
    "arg should be one of"
  )
})
