test_that("discordant loss pairs decrease as power threshold increases", {
  suppressPackageStartupMessages(library(dplyr))

  if (!exists("compute_discordant_loss_pairs", mode = "function")) {
    source_candidates <- c("R/graph_assembly.R", "../../R/graph_assembly.R")
    source_path <- source_candidates[file.exists(source_candidates)][1]
    expect_true(!is.na(source_path))
    source(source_path, local = TRUE)
  }

  earliest_loss_tbl <- tibble::tibble(
    cell_group = c("lostA", "u_hi", "u_mid", "u_low"),
    is_lost_at_peak = c(TRUE, FALSE, FALSE, FALSE),
    peak_time_in_ctrl_within_perturb_time_range = c(TRUE, TRUE, TRUE, TRUE),
    loss_when_present_power = c(0.95, 0.90, 0.55, 0.20)
  )

  discordant_pairs_no_power_gate <- compute_discordant_loss_pairs(
    earliest_loss_tbl = earliest_loss_tbl,
    power_threshold = 0
  )
  discordant_pairs_power_gate <- compute_discordant_loss_pairs(
    earliest_loss_tbl = earliest_loss_tbl,
    power_threshold = 0.6
  )

  expect_gt(nrow(discordant_pairs_no_power_gate), nrow(discordant_pairs_power_gate))
  expect_setequal(discordant_pairs_power_gate$unaffected_cell_groups, c("u_hi"))
})
