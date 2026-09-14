testthat::test_that("get_all_parents returns recursive ancestors", {
  graph <- igraph::graph_from_data_frame(
    data.frame(
      from = c("A", "B", "A"),
      to = c("B", "C", "D"),
      stringsAsFactors = FALSE
    ),
    directed = TRUE
  )

  parents_c <- platt:::get_all_parents(graph, "C")
  testthat::expect_setequal(parents_c, c("B", "A"))

  parents_d <- platt:::get_all_parents(graph, "D")
  testthat::expect_setequal(parents_d, "A")
})

testthat::test_that("get_all_parents handles cycles and unknown nodes", {
  cyclic_graph <- igraph::graph_from_data_frame(
    data.frame(
      from = c("A", "B", "B"),
      to = c("B", "A", "C"),
      stringsAsFactors = FALSE
    ),
    directed = TRUE
  )

  parents_c <- platt:::get_all_parents(cyclic_graph, "C")
  testthat::expect_setequal(parents_c, c("B", "A"))

  testthat::expect_identical(
    platt:::get_all_parents(cyclic_graph, "MISSING"),
    character(0)
  )
})
