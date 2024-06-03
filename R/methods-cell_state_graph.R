


# # plot a basic version of the graph
# # with a bezier layout



# # get edge info
# edge_info = function(state_graph) {
#     stopifnot( methods::is( state_graph, "cell_state_graph" ) )
# }

# # get node ifo

# node_info = function(state_graph) {
#     stopifnot( methods::is( state_graph, "cell_state_graph" ) )
# }


# # 
# supported_edges = function() {

# }

# nonsupported_edges = function() {

# }


# # 

# add_edge = function()
# remove_edge = function()
# add_node = function()
# remove_node = function()


# #' get the parent(s) of a state in a state transition graph
# #' @noRd
# get_parents = function(state_graph, cell_state){
#   parents = igraph::neighbors(state_graph, cell_state, mode="in")
#   if (length(parents) > 0)
#     return (parents$name)
#   else
#     return (c())
# }

# #' get the children of a state in a state transition graph
# #' @noRd
# get_children = function(state_graph, cell_state){
#   children = igraph::neighbors(state_graph, cell_state, mode="out")
#   if (length(children) > 0)
#     return (children$name)
#   else
#     return (c())
# }

# #' get the siblings of a state in a state transition graph
# #' @noRd
# get_siblings = function(state_graph, cell_state){
#   parents = get_parents(state_graph, cell_state)
#   if (length(parents) > 0){
#     siblings = igraph::neighbors(state_graph, parents, mode="out")
#     siblings = setdiff(siblings$name, cell_state) #exclude self
#     return(siblings)
#   } else{
#     return (c())
#   }

# }
