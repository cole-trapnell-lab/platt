test_that("discordant config surface defaults and mapping are stable", {
  if (!exists("normalize_discordant_pruning_config", mode = "function")) {
    source_candidates <- c("R/graph_assembly.R", "../../R/graph_assembly.R")
    source_path <- source_candidates[file.exists(source_candidates)][1]
    expect_true(!is.na(source_path))
    source(source_path, local = TRUE)
  }

  cfg_default <- normalize_discordant_pruning_config()
  expect_identical(cfg_default$prune_mode, "existing")
  expect_equal(cfg_default$K_paths, 1)
  expect_equal(cfg_default$power_threshold, 0)
  expect_equal(cfg_default$lambda_edge, 0)

  cfg_custom <- normalize_discordant_pruning_config(
    discordant_config = list(
      power_threshold = 0.7,
      prune_mode = "greedy",
      K_paths = 3,
      lambda_edge = 0.5
    )
  )
  expect_identical(cfg_custom$prune_mode, "greedy")
  expect_equal(cfg_custom$K_paths, 3)
  expect_equal(cfg_custom$power_threshold, 0.7)
  expect_equal(cfg_custom$lambda_edge, 0.5)
})
