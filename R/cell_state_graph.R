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
                   info = "SimpleList")
)
setMethod("is.na", "cell_state_graph", function(x) FALSE)


#' Create a new cell state graph object.
#'
#' @param graph input graph
#' @param ccs input cell count set
#' @name cell_state_graph
#' @return a new cell_state_graph object
#' @export cell_state_graph
new_cell_state_graph <- function(graph, ccs) {
    assertthat::assert_that(is(graph, 'igraph'))
    assertthat::assert_that(is(ccs, 'cell_count_set'))
    state_graph <- methods::new("cell_state_graph",
                                graph = graph,
                                ccs = ccs, 
                                info = SimpleList(sample_group = ccs@info$sample_group,
                                                  cell_group = ccs@info$cell_group))
    return(state_graph)
}








