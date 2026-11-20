# Over-graph fast paths: mirrored reverse contrasts, count gates, closed-form state fit.

make_pem <- function(ng = 400, seed = 1) {
  set.seed(seed)
  base <- log(runif(ng, 1e-4, 0.2))
  PEM <- cbind(a = base + rnorm(ng, 0, 0.3), b = base + rnorm(ng, 0, 0.3), c = base)
  PSEM <- matrix(runif(3 * ng, 0.05, 0.6), ng, 3, dimnames = list(NULL, colnames(PEM)))
  PEM[1:20, "a"] <- log(1e-3) - 20; PSEM[1:20, "a"] <- Inf          # a side with no counts
  rownames(PEM) <- rownames(PSEM) <- paste0("g", seq_len(ng))
  list(PEM = PEM, PSEM = PSEM)
}

test_that("mirror_contrast reproduces the reverse contrast", {
  m <- make_pem()
  for (pr in list(c("a", "b"), list("a", c("b", "c")))) {
    fwd <- contrast_helper(pr[[1]], pr[[2]], m$PEM, m$PSEM, prefix = "x_to_y", n = 500, cv_threshold = 100)
    rev <- contrast_helper(pr[[2]], pr[[1]], m$PEM, m$PSEM, prefix = "y_to_x", n = 500, cv_threshold = 100)
    mir <- mirror_contrast(fwd, "x_to_y", "y_to_x")
    expect_equal(names(mir), names(rev))
    for (col in c("y_to_x_raw_lfc", "y_to_x_raw_lfc_se", "y_to_x_raw_p_value", "log_mean_expression"))
      expect_equal(mir[[col]], rev[[col]], tolerance = 1e-12)
    for (col in c("y_to_x_shrunken_lfc", "y_to_x_shrunken_lfc_se", "y_to_x_p_value"))
      expect_equal(mir[[col]], rev[[col]], tolerance = 1e-5)              # ashr optimizer tolerance
  }
})

test_that("contrast_helper survives a side with no reliable effect", {
  m <- make_pem()
  m$PSEM[, "a"] <- Inf                                                   # every effect on side a unreliable
  res <- contrast_helper("a", "b", m$PEM, m$PSEM, n = 500, cv_threshold = 100)
  expect_false(any(is.nan(res$raw_lfc)))
})

test_that("state_contrast_measurable applies the knockout count gates", {
  g <- c("deep", "thin_zero", "one_embryo", "absent")
  sup <- list(
    detections = cbind(s = c(500, 0, 30, 0), p = c(800, 40, 30, 0)),
    embryos    = cbind(s = c(20, 0, 1, 0),   p = c(30, 10, 0, 0)),
    size = c(s = 100, p = 200), pseudobulks = c(s = 20, p = 30))
  rownames(sup$detections) <- rownames(sup$embryos) <- g
  ok <- state_contrast_measurable(sup, "s", "p", g)
  # deep: expected 100 * 5 = 500; thin_zero: expected 100 * 0.2 = 20 and 10 embryos detect it on one side;
  # one_embryo: one embryo on either side; absent: nothing anywhere
  expect_equal(unname(ok), c(TRUE, TRUE, FALSE, FALSE))
  sup1 <- sup; sup1$pseudobulks["s"] <- 1                                # a one-embryo state: no contrast
  expect_false(any(state_contrast_measurable(sup1, "s", "p", g)))
})

test_that("closed_form_state_coefficients matches the monocle fit", {
  skip_if_not_installed("monocle3")
  set.seed(1)
  ng <- 60; ns <- 120; grp <- factor(rep(paste0("s", 1:4), each = 30)); sf <- runif(ns, 5, 50)
  mu <- outer(exp(runif(ng, -6, 0)), sf) * exp(rnorm(ng * ns, 0, 0.5))
  Y <- matrix(rpois(ng * ns, mu), ng, ns, dimnames = list(paste0("g", 1:ng), paste0("p", 1:ns)))
  cds <- suppressWarnings(monocle3::new_cell_data_set(methods::as(Y, "CsparseMatrix"),
    cell_metadata = data.frame(cell_group = grp, row.names = colnames(Y)),
    gene_metadata = data.frame(id = rownames(Y), gene_short_name = rownames(Y), row.names = rownames(Y))))
  SummarizedExperiment::colData(cds)$Size_Factor <- sf
  fits <- monocle3::fit_models(cds, "~ 0 + cell_group", cores = 1) |>
    dplyr::select(gene_short_name, id, model, model_summary, status)
  ref <- suppressMessages(collect_coefficients_for_shrinkage(cds, fits, 1e-3, "cell_group"))
  cf <- closed_form_state_coefficients(cds, "cell_group", 1e-3)
  ok <- is.finite(ref$stdev.unscaled) & is.finite(cf$stdev.unscaled[rownames(ref$coefficients), ])
  expect_equal(cf$coefficients[rownames(ref$coefficients), ][ok], ref$coefficients[ok], tolerance = 1e-4)
  expect_equal(unname(cf$est_dispersion[rownames(ref$coefficients)]), unname(ref$est_dispersion), tolerance = 1e-4)   # speedglm stops at acc = 0.001
  expect_equal(cf$stdev.unscaled[rownames(ref$coefficients), ][ok], ref$stdev.unscaled[ok], tolerance = 1e-2)
})
