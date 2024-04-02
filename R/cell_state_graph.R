setOldClass(c("igraph"), prototype=structure(list(), class="igraph"))
setOldClass(c("PLNnetworkfamily"), prototype=structure(list(), class="PLNnetworkfamily"))
setOldClass(c("PLNnetworkfit"), prototype=structure(list(), class="PLNnetworkfit"))
setOldClass(c("PLNfamily"), prototype=structure(list(), class="PLNfamily"))
setOldClass(c("PLNfit"), prototype=structure(list(), class="PLNfit"))

#' The cell_state_graph class
#'
#' The main class used by Platt to store state graphs
#'
#' @field ccs cell_count_set, the Hooke cell count set
#' @field
#' @name cell_state_graph
#' @rdname cell_state_graph
#' @aliases cell_state_graph-class
#' @import igraph
#' @exportClass cell_state_graph
setClass("cell_state_graph",
         slots = c(ccs = "cell_count_set", 
                   graph = "igraph", 
                   layout_info = "list", 
                   g = "data.frame")
)

# setMethod("is.na", "cell_state_graph", function(x) FALSE)
# currently causing problems -- ask cole why need this again 


#' Create a new cell state graph object.
#'
#' @param graph input graph
#' @param ccs input cell count set
#' @name cell_state_graph
#' @return a new cell_state_graph object
#' @export cell_state_graph
new_cell_state_graph <- function(state_graph, ccs) {
    assertthat::assert_that(is(state_graph, 'igraph'))
    assertthat::assert_that(is(ccs, 'cell_count_set'))
    
    if (is(state_graph, "igraph")){
      edges = state_graph %>% igraph::as_data_frame()
    } else{
      edges = state_graph
    }
    
    layout_res = get_graph_layout(ccs, state_graph)
    layout_info = layout_res$layout_info
    gvizl_coords = layout_info$gvizl_coords
    bezier_df = layout_info$bezier_df
    
    state_graph <- methods::new("cell_state_graph",
                                graph = state_graph,
                                ccs = ccs, 
                                layout_info = layout_info, 
                                g = layout_res$g)
    
    return(state_graph)
}


get_graph_layout = function(ccs, 
                            state_graph,
                            color_nodes_by = NULL,
                            label_nodes_by= NULL,
                            group_nodes_by= NULL,
                            arrow.gap=0.03) {
  
  
  if (is(state_graph, "igraph")){
    edges = state_graph %>% igraph::as_data_frame()
  } else{
    edges = state_graph
  }
  
  if (is.null(color_nodes_by)) { 
    color_nodes_by = ccs@info$cell_group
  }
  
  if (is.null(label_nodes_by)) { 
    label_nodes_by = ccs@info$cell_group
  }
  
  if (is.null(group_nodes_by)) { 
    group_nodes_by = ccs@info$cell_group
  }
  
  node_metadata = platt:::collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)
  node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  
  edges = edges %>% dplyr::ungroup()
  edges = edges %>% filter(from %in% node_metadata$id & to %in% node_metadata$id)
  
  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)
  
  layout_info = platt:::layout_state_graph(G, node_metadata, NULL, weighted=FALSE)
  
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  
  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)
  
  return(list(layout_info = layout_info, g = g, arrow.gap = arrow.gap))  
  
}



# update anything about the layout information about a graph 
update_graph_layout = function(cell_state_graph, 
                               color_nodes_by = NULL,
                               label_nodes_by= NULL,
                               group_nodes_by= NULL,
                               arrow.gap=0.03,
                               arrow_unit = 7,
                               bar_unit = .075,
                               node_size = 2,
                               min_edge_size=0.1,
                               max_edge_size=2,
                               fract_expr=0.0,
                               mean_expr=0.0, 
                               con_colour = "darkgrey") {
  
  state_graph = cell_state_graph@graph
  ccs = cell_state_graph@ccs
  
  graph_layout_res = get_graph_layout(ccs, state_graph, 
                                      color_nodes_by = color_nodes_by,
                                      label_nodes_by= label_nodes_by,
                                      group_nodes_by= group_nodes_by,
                                      arrow.gap=arrow.gap)

  cell_state_graph@layout_info <- graph_layout_res$layout_info
  cell_state_graph@g <- graph_layout_res$g
  
  return(cell_state_graph)

}
  

# bar_unit = .075,
# 
# min_edge_size=0.1,
# max_edge_size=2,


