test_that("prune mode reporting helpers work on toy graphs", {
  suppressPackageStartupMessages(library(igraph))

  if (!exists("count_discordant_reachability_violations", mode = "function")) {
    source_candidates <- c("R/graph_assembly.R", "../../R/graph_assembly.R")
    source_path <- source_candidates[file.exists(source_candidates)][1]
    expect_true(!is.na(source_path))
    source(source_path, local = TRUE)
  }

  g <- igraph::graph_from_data_frame(
    data.frame(from = c("A", "B"), to = c("B", "C")),
    directed = TRUE
  )
  igraph::E(g)$total_path_score_supporting <- c(2, 3)
  igraph::E(g)$discordance_penalty <- c(0.1, 0.2)

  pairs <- tibble::tibble(from = c("A", "A", "C"), to = c("C", "D", "A"))
  expect_equal(count_discordant_reachability_violations(g, pairs), 1L)

  summary_tbl <- summarize_edge_support_distributions(g, mode_label = "toy")
  expect_true(all(c("mode", "metric", "n_edges", "mean", "median", "p90", "max") %in% colnames(summary_tbl)))
  expect_true(any(summary_tbl$metric == "total_path_score_supporting"))
  expect_true(any(summary_tbl$metric == "discordance_penalty"))
})
