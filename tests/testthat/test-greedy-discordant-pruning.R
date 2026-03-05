test_that("greedy discordant pruning breaks all forbidden pairs on toy graph", {
  suppressPackageStartupMessages(library(igraph))
  suppressPackageStartupMessages(library(dplyr))

  if (!exists("prune_discordant_paths_greedily", mode = "function")) {
    source_candidates <- c("R/graph_assembly.R", "../../R/graph_assembly.R")
    source_path <- source_candidates[file.exists(source_candidates)][1]
    expect_true(!is.na(source_path))
    source(source_path, local = TRUE)
  }

  g <- igraph::graph_from_data_frame(
    data.frame(
      from = c("A", "A", "B", "C", "B"),
      to = c("B", "C", "D", "D", "C"),
      weight = c(1, 1, 1, 1, 0.5)
    ),
    directed = TRUE
  )
  igraph::E(g)$total_path_score_supporting <- c(2, 3, 2, 2, 1)

  forbidden_pairs <- tibble::tibble(
    from = c("A", "B"),
    to = c("D", "D"),
    pair_weight = c(1, 1)
  )

  pruned <- prune_discordant_paths_greedily(
    state_graph = g,
    discordant_pairs = forbidden_pairs,
    k_paths = 2,
    deletion_cost_attr = "total_path_score_supporting",
    traversal_weight_attr = "weight"
  )

  expect_gt(nrow(pruned$removed_edges), 0)

  for (i in seq_len(nrow(forbidden_pairs))) {
    from_node <- forbidden_pairs$from[[i]]
    to_node <- forbidden_pairs$to[[i]]
    paths <- suppressWarnings(igraph::all_simple_paths(pruned$graph, from = from_node, to = to_node, mode = "out"))
    expect_equal(length(paths), 0)
  }
})
