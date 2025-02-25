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
#' @name cell_state_graph
#' @rdname cell_state_graph
#' @aliases cell_state_graph-class
#' @import igraph
#' @exportClass cell_state_graph
setClass("cell_state_graph",
         # contains = "cell_count_set",
         slots = c(ccs = "cell_count_set", 
                   graph = "igraph", 
                   layout_info = "list", 
                   g = "data.frame", 
                   metadata = "list", 
                   genetic_requirements = "data.frame")
)

# setMethod("is.na", "cell_state_graph", function(x) FALSE)
# currently causing problems -- ask cole why need this again 


#' Create a new cell state graph object.
#'
#' @param graph input graph
#' @param ccs input cell count set
#' @return a new cell_state_graph object
#' @export 
new_cell_state_graph <- function(state_graph, 
                                 ccs, 
                                 color_nodes_by = NULL,
                                 label_nodes_by= NULL,
                                 group_nodes_by= NULL, 
                                 edge_weights=NULL, 
                                 min_edge_size=0.1,
                                 max_edge_size=2, 
                                 hide_unlinked_nodes = F, 
                                 genetic_requirements = data.frame(),
                                 num_layers = 1) {
    assertthat::assert_that(is(state_graph, 'igraph'))
    assertthat::assert_that(is(ccs, 'cell_count_set'))
    
    if (is(state_graph, "igraph")){
      edges = state_graph %>% igraph::as_data_frame()
    } else{
      edges = state_graph
    }
    
    layout_res = get_graph_layout(ccs, 
                                  state_graph, 
                                  color_nodes_by = color_nodes_by, 
                                  label_nodes_by = label_nodes_by, 
                                  group_nodes_by = group_nodes_by, 
                                  hide_unlinked_nodes = hide_unlinked_nodes,
                                  num_layers = num_layers)
    
    layout_info = layout_res$layout_info
    gvizl_coords = layout_info$gvizl_coords
    bezier_df = layout_info$bezier_df
    
    
    if (is.null(edge_weights) == FALSE){
      bezier_df = left_join(bezier_df, edges)
      bezier_df = bezier_df %>% mutate(edge_score =  (weight - min(weight, na.rm=TRUE)) / max(weight, na.rm=TRUE),
                                       edge_thickness = ((max_edge_size - min_edge_size) * edge_score) + min_edge_size,
                                       unsupported_edge = ifelse(is.na(weight), TRUE, FALSE),
                                       edge_thickness = replace_na(edge_thickness, min_edge_size))
    }else{
      bezier_df$edge_thickness = (max_edge_size + min_edge_size) / 2
      bezier_df$unsupported_edge = FALSE
    }
    layout_info$bezier_df = bezier_df
    
    
    state_graph <- methods::new("cell_state_graph",
                                graph = state_graph,
                                ccs = ccs, 
                                layout_info = layout_info, 
                                g = layout_res$g,
                                metadata = list(color_nodes_by = color_nodes_by, 
                                                label_nodes_by= label_nodes_by,
                                                group_nodes_by= group_nodes_by, 
                                                bezier_df = bezier_df), 
                                genetic_requirements = genetic_requirements)
    
    return(state_graph)
}


#'
#' @param ccs
#' @param state_graph
#' @param color_nodes_by
#' @param label_nodes_by
#' @param group_nodes_by
get_graph_layout = function(ccs, 
                            state_graph,
                            color_nodes_by = NULL,
                            label_nodes_by= NULL,
                            group_nodes_by= NULL,
                            hide_unlinked_nodes = TRUE,
                            arrow.gap=0.03,
                            num_layers=1) {
  
  
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
  
  node_metadata = collect_psg_node_metadata(ccs, color_nodes_by, label_nodes_by, group_nodes_by)
  if (hide_unlinked_nodes){
    node_metadata = node_metadata %>% filter(id %in% edges$from | id %in% edges$to)
  }
  
  edges = edges %>% dplyr::ungroup()
  edges = edges %>% select(from, to)
  edges$weight = 1
  edges = edges %>% filter(from %in% node_metadata$id & to %in% node_metadata$id)
  
  G = edges %>% distinct() %>% igraph::graph_from_data_frame(directed = T, vertices=node_metadata)
  
  layout_info = layout_state_graph(G, node_metadata, edge_labels=NULL, weighted=FALSE, num_layers=num_layers)
  
  gvizl_coords = layout_info$gvizl_coords
  bezier_df = layout_info$bezier_df
  
  g = ggnetwork::ggnetwork(G, layout = gvizl_coords, arrow.gap = arrow.gap, scale=F)
  
  return(list(layout_info = layout_info, g = g, arrow.gap = arrow.gap))  
  
}



#' update anything about the layout information about a graph 
#' @param cell_state_graph
#' @param color_nodes_by description
#' @param label_nodes_by description
#' @param group_nodes_by 
update_graph_layout = function(cell_state_graph, 
                               color_nodes_by = NULL,
                               label_nodes_by = NULL,
                               group_nodes_by = NULL,
                               arrow.gap = 0.03,
                               arrow_unit = 7,
                               bar_unit = .075,
                               node_size = 2,
                               min_edge_size = 0.1,
                               max_edge_size = 2,
                               fract_expr = 0.0, 
                               mean_expr = 0.0, 
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
  
