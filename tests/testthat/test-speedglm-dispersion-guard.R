# update_summary() recomputes model summaries with a chosen dispersion. Only the lab speedglm fork
# accepts one; with CRAN speedglm it is ignored without a warning, so Platt must stop instead.

make_speedglm_model_tbl <- function(n_genes = 3, seed = 7) {
  set.seed(seed)
  df <- data.frame(g = factor(rep(c("a", "b"), each = 15)))
  models <- lapply(seq_len(n_genes), function(i) {
    df$y <- rpois(nrow(df), lambda = ifelse(df$g == "a", 6, 12) * exp(rnorm(nrow(df), 0, 0.4)))
    speedglm::speedglm(y ~ 0 + g, data = df, family = stats::quasipoisson(), acc = 1e-3, model = FALSE, y = FALSE)
  })
  tibble::tibble(
    id = paste0("gene", seq_len(n_genes)),
    model = models,
    dispersion = vapply(models, function(m) summary(m)$dispersion, numeric(1)),
    disp_fit = c(0.5, 2, 8)[seq_len(n_genes)]
  )
}

summary_se <- function(model_summary) model_summary$coefficients[, "Std. Error"]

test_that("update_summary stops when speedglm cannot take a dispersion", {
  skip_if_not_installed("speedglm")
  tbl <- make_speedglm_model_tbl()
  testthat::local_mocked_bindings(.speedglm_takes_dispersion = function() FALSE)
  expect_error(update_summary(tbl, dispersion_type = "fitted"), "no `dispersion` argument")
  expect_error(update_summary(tbl, dispersion_type = "fitted"), "cole-trapnell-lab/speedglm")
  expect_error(update_summary(tbl, dispersion_type = "max"), "no `dispersion` argument")
  # each model's own dispersion is what any speedglm uses by default, so this is not checked
  expect_no_error(update_summary(tbl, dispersion_type = "estimated"))
})

test_that("update_summary does not check when there is no speedglm model", {
  tbl <- tibble::tibble(id = c("gene1", "gene2"), model = list(NA, NA), dispersion = c(1, 1), disp_fit = c(2, 2))
  testthat::local_mocked_bindings(.speedglm_takes_dispersion = function() FALSE)
  expect_no_error(update_summary(tbl, dispersion_type = "fitted"))
})

test_that("with the speedglm fork, standard errors use the dispersion that was asked for", {
  skip_if_not_installed("speedglm")
  skip_if_not(.speedglm_takes_dispersion(), "installed speedglm has no dispersion argument (not the lab fork)")
  tbl <- make_speedglm_model_tbl()
  raw <- lapply(tbl$model, function(m) summary_se(summary(m)))
  fitted <- lapply(update_summary(tbl, dispersion_type = "fitted")$model_summary, summary_se)
  estimated <- lapply(update_summary(tbl, dispersion_type = "estimated")$model_summary, summary_se)
  for (i in seq_len(nrow(tbl))) {
    # SE scales with the square root of the dispersion
    expect_equal(fitted[[i]], raw[[i]] * sqrt(tbl$disp_fit[i] / tbl$dispersion[i]), tolerance = 1e-10)
    expect_equal(estimated[[i]], raw[[i]], tolerance = 1e-10)
  }
  # "max" never goes below the raw dispersion or min_dispersion
  maxed <- lapply(update_summary(tbl, dispersion_type = "max")$model_summary, summary_se)
  for (i in seq_len(nrow(tbl))) {
    d <- max(tbl$disp_fit[i], tbl$dispersion[i], 1)
    expect_equal(maxed[[i]], raw[[i]] * sqrt(d / tbl$dispersion[i]), tolerance = 1e-10)
  }
})
