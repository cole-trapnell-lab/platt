test_that("edge discordance penalty is power-weighted and decreases with threshold", {
  suppressPackageStartupMessages(library(dplyr))
  suppressPackageStartupMessages(library(igraph))

  if (!exists("compute_edge_discordance_penalty", mode = "function")) {
    source_candidates <- c("R/graph_assembly.R", "../../R/graph_assembly.R")
    source_path <- source_candidates[file.exists(source_candidates)][1]
    expect_true(!is.na(source_path))
    source(source_path, local = TRUE)
  }

  state_transition_graph <- igraph::graph_from_data_frame(
    data.frame(
      from = c("A", "A"),
      to = c("B", "C")
    ),
    directed = TRUE
  )

  perturbation_ccm_tbl <- tibble::tibble(
    perturb_name = c("p1", "p2"),
    perturb_summary_tbl = list(
      tibble::tibble(
        cell_group = c("A", "B", "C"),
        is_lost_when_present = c(TRUE, FALSE, FALSE),
        loss_when_present = c(-2.0, 0.0, 0.0),
        loss_when_present_power = c(0.9, 0.9, 0.2)
      ),
      tibble::tibble(
        cell_group = c("A", "B", "C"),
        is_lost_when_present = c(TRUE, FALSE, TRUE),
        loss_when_present = c(-1.0, 0.0, -0.5),
        loss_when_present_power = c(0.9, 0.5, 0.9)
      )
    )
  )

  penalty_no_gate <- compute_edge_discordance_penalty(
    perturbation_ccm_tbl = perturbation_ccm_tbl,
    state_transition_graph = state_transition_graph,
    power_threshold = 0
  )
  penalty_power_gate <- compute_edge_discordance_penalty(
    perturbation_ccm_tbl = perturbation_ccm_tbl,
    state_transition_graph = state_transition_graph,
    power_threshold = 0.6
  )

  expect_true(all(c("discordance_penalty", "num_discordant_perturbs", "mean_discordant_power") %in% colnames(penalty_no_gate)))
  expect_gt(sum(penalty_no_gate$discordance_penalty), sum(penalty_power_gate$discordance_penalty))
  expect_gt(
    penalty_no_gate$discordance_penalty[penalty_no_gate$from == "A" & penalty_no_gate$to == "C"],
    penalty_power_gate$discordance_penalty[penalty_power_gate$from == "A" & penalty_power_gate$to == "C"]
  )
})
