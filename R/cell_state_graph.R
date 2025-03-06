setOldClass(c("igraph"), prototype = structure(list(), class = "igraph"))
setOldClass(c("PLNnetworkfamily"), prototype = structure(list(), class = "PLNnetworkfamily"))
setOldClass(c("PLNnetworkfit"), prototype = structure(list(), class = "PLNnetworkfit"))
setOldClass(c("PLNfamily"), prototype = structure(list(), class = "PLNfamily"))
setOldClass(c("PLNfit"), prototype = structure(list(), class = "PLNfit"))

#' The main class used by Platt to store state graphs
#'
#' @field ccs cell_count_set, the Hooke cell count set
#' @field graph igraph, the graph structure representing cell states
#' @field layout_info list, information about the layout of the graph
#' @field g data.frame, additional data associated with the graph
#' @field metadata list, metadata associated with the cell state graph
#' @field genetic_requirements data.frame, genetic requirements for the cell states
#' @name cell_state_graph
#' @rdname cell_state_graph
#' @aliases cell_state_graph-class
#' @import igraph
#' @exportClass cell_state_graph
setClass("cell_state_graph",
  # contains = "cell_count_set",
  slots = c(
    ccs = "cell_count_set",
    graph = "igraph",
    layout_info = "list",
    g = "data.frame",
    metadata = "list",
    genetic_requirements = "data.frame"
  )
)

# setMethod("is.na", "cell_state_graph", function(x) FALSE)
# currently causing problems -- ask cole why need this again

#' Create a new cell state graph object.
#'
#' This function creates a new cell state graph object from an input graph and cell count set.
#'
#' @param state_graph An input graph of class 'igraph'.
#' @param ccs An input cell count set of class 'cell_count_set'.
#' @param color_nodes_by Optional parameter to color nodes by a specific attribute.
#' @param label_nodes_by Optional parameter to label nodes by a specific attribute.
#' @param group_nodes_by Optional parameter to group nodes by a specific attribute.
#' @param edge_weights Optional parameter to specify edge weights.
#' @param min_edge_size Minimum edge size for visualization. Default is 0.1.
#' @param max_edge_size Maximum edge size for visualization. Default is 2.
#' @param hide_unlinked_nodes Logical parameter to hide unlinked nodes. Default is FALSE.
#' @param genetic_requirements A data frame specifying genetic requirements. Default is an empty data frame.
#' @param num_layers Number of layers for the graph layout. Default is 1.
#' @return A new cell_state_graph object.
#' @export
new_cell_state_graph <- function(state_graph,
                                 ccs,
                                 color_nodes_by = NULL,
                                 label_nodes_by = NULL,
                                 group_nodes_by = NULL,
                                 edge_weights = NULL,
                                 min_edge_size = 0.1,
                                 max_edge_size = 2,
                                 hide_unlinked_nodes = F,
                                 genetic_requirements = data.frame(),
                                 num_layers = 1) {
  assertthat::assert_that(is(state_graph, "igraph"))
  assertthat::assert_that(is(ccs, "cell_count_set"))

  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  layout_res <- get_graph_layout(ccs,
    state_graph,
    color_nodes_by = color_nodes_by,
    label_nodes_by = label_nodes_by,
    group_nodes_by = group_nodes_by,
    hide_unlinked_nodes = hide_unlinked_nodes,
    num_layers = num_layers
  )

  layout_info <- layout_res$layout_info
  bezier_df <- layout_info$bezier_df

  if (is.null(edge_weights) == FALSE) {
    bezier_df <- dplyr::left_join(bezier_df, edges)
    bezier_df <- bezier_df %>% mutate(
      edge_score = (weight - min(weight, na.rm = TRUE)) / max(weight, na.rm = TRUE),
      edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
      unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
      edge_thickness = tidyr::replace_na(edge_thickness, min_edge_size)
    )
  } else {
    bezier_df$edge_thickness <- (max_edge_size + min_edge_size) / 2
    bezier_df$unsupported_edge <- FALSE
  }
  layout_info$bezier_df <- bezier_df


  state_graph <- methods::new("cell_state_graph",
    graph = state_graph,
    ccs = ccs,
    layout_info = layout_info,
    g = layout_res$g,
    metadata = list(
      color_nodes_by = color_nodes_by,
      label_nodes_by = label_nodes_by,
      group_nodes_by = group_nodes_by,
      bezier_df = bezier_df
    ),
    genetic_requirements = genetic_requirements
  )

  return(state_graph)
}

#' Get Graph Layout
#'
#' This function generates a layout for a cell state graph.
#'
#' @param ccs An object containing cell state information.
#' @param state_graph An igraph object or a data frame representing the state graph.
#' @param color_nodes_by A vector specifying how to color the nodes. Default is NULL.
#' @param label_nodes_by A vector specifying how to label the nodes. Default is NULL.
#' @param group_nodes_by A vector specifying how to group the nodes. Default is NULL.
#' @param hide_unlinked_nodes A logical value indicating whether to hide unlinked nodes. Default is TRUE.
#' @param arrow_gap A numeric value specifying the gap for arrows. Default is 0.03.
#' @param num_layers An integer specifying the number of layers for the layout. Default is 1.
#'
#' @return A list containing the layout information and the ggnetwork object.
#'
#' @import dplyr
#' @import igraph
#' @import ggnetwork
#'
#' @examples
#' \dontrun{
#' ccs <- load_cell_state_data()
#' state_graph <- create_state_graph(ccs)
#' layout <- get_graph_layout(ccs, state_graph)
#' print(layout$g)
#' }
#' @export
get_graph_layout <- function(ccs,
                             state_graph,
                             color_nodes_by = NULL,
                             label_nodes_by = NULL,
                             group_nodes_by = NULL,
                             hide_unlinked_nodes = TRUE,
                             arrow_gap = 0.03,
                             num_layers = 1) {
  if (is(state_graph, "igraph")) {
    edges <- state_graph %>% igraph::as_data_frame()
  } else {
    edges <- state_graph
  }

  if (is.null(color_nodes_by)) {
    color_nodes_by <- ccs@info$cell_group
  }

  if (is.null(label_nodes_by)) {
    label_nodes_by <- ccs@info$cell_group
  }

  if (is.null(group_nodes_by)) {
    group_nodes_by <- ccs@info$cell_group
  }

  node_metadata <- collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)
  if (hide_unlinked_nodes) {
    node_metadata <- node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }

  edges <- edges %>% dplyr::ungroup()
  edges <- edges %>% dplyr::select(from, to)
  edges$weight <- 1
  edges <- edges %>% dplyr::filter(from %in% node_metadata$id & to %in% node_metadata$id)

  G <- edges %>%
    dplyr::distinct() %>%
    igraph::graph_from_data_frame(directed = TRUE, vertices = node_metadata)

  layout_info <- layout_state_graph(
    G,
    node_metadata,
    edge_labels = NULL,
    weighted = FALSE,
    num_layers = num_layers
  )

  gvizl_coords <- layout_info$gvizl_coords

  g <- ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow_gap = arrow_gap, scale = F)

  return(list(layout_info = layout_info, g = g, arrow_gap = arrow_gap))
}


#' Update the layout information of a cell state graph
#'
#' This function updates the layout information of a given cell state graph
#' based on various parameters such as node color, node label, and node grouping.
#'
#' @param cell_state_graph An object representing the cell state graph.
#' @param color_nodes_by A parameter to specify how to color the nodes (default is NULL).
#' @param label_nodes_by A parameter to specify how to label the nodes (default is NULL).
#' @param group_nodes_by A parameter to specify how to group the nodes (default is NULL).
#' @param arrow_gap A numeric value specifying the gap for arrows (default is 0.03).
#' @param arrow_unit A numeric value specifying the unit for arrows (default is 7).
#' @param bar_unit A numeric value specifying the unit for bars (default is 0.075).
#' @param node_size A numeric value specifying the size of the nodes (default is 2).
#' @param min_edge_size A numeric value specifying the minimum size of the edges (default is 0.1).
#' @param max_edge_size A numeric value specifying the maximum size of the edges (default is 2).
#' @param fract_expr A numeric value specifying the fraction expression (default is 0.0).
#' @param mean_expr A numeric value specifying the mean expression (default is 0.0).
#' @param con_colour A character string specifying the color for connections (default is "lightgray").
#'
#' @return The updated cell state graph with new layout information.
#'
#' @examples
#' # Example usage:
#' updated_graph <- update_graph_layout(cell_state_graph, color_nodes_by = "type", label_nodes_by = "name")
update_graph_layout <- function(cell_state_graph,
                                color_nodes_by = NULL,
                                label_nodes_by = NULL,
                                group_nodes_by = NULL,
                                arrow_gap = 0.03,
                                arrow_unit = 7,
                                bar_unit = .075,
                                node_size = 2,
                                min_edge_size = 0.1,
                                max_edge_size = 2,
                                fract_expr = 0.0,
                                mean_expr = 0.0,
                                con_colour = "lightgray") {
  state_graph <- cell_state_graph@graph
  ccs <- cell_state_graph@ccs

  graph_layout_res <- get_graph_layout(ccs, state_graph,
    color_nodes_by = color_nodes_by,
    label_nodes_by = label_nodes_by,
    group_nodes_by = group_nodes_by,
    arrow_gap = arrow_gap
  )

  cell_state_graph@layout_info <- graph_layout_res$layout_info
  cell_state_graph@g <- graph_layout_res$g

  return(cell_state_graph)
}
