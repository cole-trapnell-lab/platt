test_that("assert_power_columns_present validates power columns", {
  expect_no_error(
    platt:::assert_power_columns_present(tibble::tibble(cell_group = "a", power = 0.8))
  )
  expect_no_error(
    platt:::assert_power_columns_present(tibble::tibble(cell_group = "a", loss_when_present_power = 0.8))
  )
  expect_error(
    platt:::assert_power_columns_present(tibble::tibble(cell_group = "a")),
    "Expected at least one power column"
  )
})
