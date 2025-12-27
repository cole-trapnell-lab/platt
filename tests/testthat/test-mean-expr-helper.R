test_that("compute_model_mean_expr uses fitted values when available", {
  skip_if_not_installed("speedglm")
  set.seed(123)
  df <- data.frame(
    x = rnorm(20),
    y = rpois(20, lambda = exp(0.5 + rnorm(20)))
  )
  m <- speedglm::speedglm(y ~ x, data = df, family = poisson(), clean = FALSE)
  if (is.null(m$fitted.values)) {
    expect_equal(
      platt:::compute_model_mean_expr(m, new_data = df),
      mean(speedglm:::predict.speedglm(m, newdata = df, type = "response")),
      tolerance = 1e-12
    )
  } else {
    expect_equal(
      platt:::compute_model_mean_expr(m, new_data = df),
      mean(m$fitted.values),
      tolerance = 1e-12
    )
  }
})

test_that("compute_model_mean_expr falls back to predict", {
  df <- data.frame(
    x = rnorm(10),
    y = rnorm(10)
  )
  m <- lm(y ~ x, data = df)
  expect_equal(
    platt:::compute_model_mean_expr(m, new_data = df),
    mean(predict(m, newdata = df)),
    tolerance = 1e-12
  )
})
