#'
#' @param res
build_whitelist = function(res, graph_type = c("wt", "mt")) {
  # build a whitelist
  
  graph_type = match.arg(graph_type)
  
  if (graph_type == "wt"){
    edge_whitelist = do.call(igraph::union, res %>% filter(is.na(wt_graph) == FALSE) %>% pull(wt_graph))
  } else {
    edge_whitelist = do.call(igraph::union, res %>% filter(is.na(mt_graph) == FALSE) %>% pull(mt_graph))
  }
  
  edge_whitelist = igraph::as_data_frame(edge_whitelist)
  edge_whitelist = edge_whitelist %>% select(from, to) %>% distinct()
  
  return(edge_whitelist)
  
}


#'
#' @param cds
#' @param res
transfer_cell_states = function(cds, res, colname="cell_state") {
  colData(cds)$cds_row_id = colnames(cds)
  colData(cds)[[colname]] = left_join(colData(cds) %>% as.data.frame(),
                                      res %>% tidyr::unnest(data) %>% select(cds_row_id, subassembly_group),
                                      by = "cds_row_id") %>% pull(subassembly_group)
  return(cds)
}


#' after fitting on sub partitions, fit a model on the entire cds
#' 
#' @export
fit_global_wt_model = function(cds, 
                               res, 
                               sample_group = "embryo",
                               cell_group = "cell_state",
                               main_model_formula_str = NULL,
                               start_time = 18,
                               stop_time = 72,
                               interval_col = "timepoint",
                               ctrl_ids = c("Control"),
                               perturbation_col = "gene_target",
                               component_col = "partition", 
                               nuisance_model_formula_str = "~1", 
                               num_threads=1, 
                               assembly_backend = "nlopt", 
                               num_time_breaks = 4,
                               vhat_method = "bootstrap", 
                               sparsity_factor = 0.01) {
  
  
  global_wt_graph_edge_whitelist = build_whitelist(res)
  
  # transfer the cell states back to global cds
  
  colData(cds)$cds_row_id = colnames(cds)
  
  colData(cds)$cell_state = left_join(colData(cds) %>% as.data.frame(),
                                      res %>% tidyr::unnest(data) %>% select(cds_row_id, subassembly_group),
                                      by = "cds_row_id") %>% pull(subassembly_group)
  
  
  # Fit a single wild-type cell count timeseries model to all the cell states at once
  global_wt_ccm = fit_wt_model(cds,
                               sample_group = sample_group,
                               cell_group = cell_group,
                               start_time = start_time,
                               stop_time = stop_time,
                               interval_col = interval_col,
                               vhat_method = vhat_method,
                               num_time_breaks = num_time_breaks,
                               nuisance_model_formula_str = nuisance_model_formula_str,
                               ctrl_ids = ctrl_ids,
                               sparsity_factor = sparsity_factor,
                               perturbation_col = perturbation_col,
                               edge_whitelist = global_wt_graph_edge_whitelist,
                               keep_cds = TRUE,
                               num_threads = num_threads,
                               backend = assembly_backend,
                               verbose = TRUE,
                               penalize_by_distance = TRUE,
                               pln_num_penalties = 30)
  
  return(global_wt_ccm)
  
  # global_wt_graph = assemble_wt_graph(cds,
  #                                     global_wt_ccm,
  #                                     sample_group = sample_group,
  #                                     cell_group = cell_group,
  #                                     main_model_formula_str = NULL,
  #                                     start_time = assembly_start_time,
  #                                     stop_time = assembly_stop_time,
  #                                     interval_col = interval_col,
  #                                     ctrl_ids = control_genotypes,
  #                                     sparsity_factor = 0.01,
  #                                     perturbation_col = perturbation_col,
  #                                     component_col = component_col,
  #                                     edge_whitelist = global_wt_graph_edge_whitelist,
  #                                     verbose=TRUE)
  
  # return(list(ccm = global_wt_ccm,
  #             graph = global_wt_graph))
  
}


#' fit a global perturbation model 
#' @param cds input cds
#' @param global_wt_ccm ccm fit on all the umaps
#' @export
fit_global_mt_model = function(cds, 
                               global_wt_ccm, 
                               res, 
                               sample_group = "embryo",
                               cell_group = "cell_state",
                               main_model_formula_str = NULL,
                               start_time = assembly_start_time,
                               stop_time = assembly_stop_time,
                               interval_col = "timepoint",
                               ctrl_ids = control_genotypes,
                               perturbation_col = "gene_target",
                               component_col="partition", 
                               nuisance_model_formula_str = "~1", 
                               num_time_breaks=3,
                               sparsity_factor = 0.01, 
                               vhat_method="bootstrap") {
  
  
  # Learn a single graph on all states at once (using the subgraphs as a whitelist/prior)
  
  
  # Fit cell count models to each mutant vs control
  global_perturb_models_tbl = fit_mt_models(cds,
                                            sample_group = sample_group,
                                            cell_group = cell_group,
                                            main_model_formula_str = NULL,
                                            start_time = assembly_start_time,
                                            stop_time = assembly_stop_time,
                                            interval_col= interval_col,
                                            num_time_breaks = num_time_breaks,
                                            ctrl_ids = control_genotypes,
                                            mt_ids = mt_genotypes,
                                            sparsity_factor = sparsity_factor,
                                            perturbation_col = perturbation_col,
                                            keep_cds=FALSE,
                                            num_threads=num_threads,
                                            backend=assembly_backend,
                                            vhat_method=vhat_method,
                                            penalize_by_distance=TRUE)
  
  # # Build a whitelist of edges by collecting the edges in the subassemblies from the mutants
  # global_mt_graph_edge_whitelist = build_whitelist(res)
  # 
  # # Build a global assembly from all the mutant models, using the subassembly as a whitelist
  # # NOTE: this graph should really only have edges that are directly supported by
  # # genetic perturbations, and may therefore be somewhat sparse.
  # global_mt_graph = assemble_mt_graph(global_wt_ccm,
  #                                             global_perturb_models_tbl,
  #                                             start_time = assembly_start_time,
  #                                             stop_time = assembly_stop_time,
  #                                             interval_col = interval_col,
  #                                             perturbation_col = perturbation_col,
  #                                             component_col = component_col, 
  #                                             edge_whitelist = global_mt_graph_edge_whitelist,
  #                                             q_val=0.1,
  #                                             verbose=TRUE)
  
  return(global_perturb_models)
  
  
}


#' to be filled in
#' @export
fit_global_models = function(res,
                             cds = cds,
                             sample_group = "embryo",
                             cell_group = "cell_state",
                             main_model_formula_str = NULL,
                             start_time = assembly_start_time,
                             stop_time = assembly_stop_time,
                             interval_col = "timepoint",
                             ctrl_ids = control_genotypes,
                             perturbation_col = "gene_target") {
  
  
  
  
  # 
  # # make the annotated graph
  # global_wt_graph_edges = igraph::as_data_frame(global_wt_graph)
  # global_mt_graph_edges = igraph::as_data_frame(global_mt_graph)
  # mt_only = setdiff(global_mt_graph_edges %>% select(from, to), global_wt_graph_edges %>% select(from, to))
  # 
  # global_graph_annotated = left_join(global_wt_graph_edges, global_mt_graph_edges)
  # global_graph_annotated = global_graph_annotated %>% select(-support)
  # global_graph_annotated = rbind(global_graph_annotated,
  #                                global_mt_graph_edges %>% inner_join(mt_only))
  # global_graph_annotated = igraph::graph_from_data_frame(global_graph_annotated)
  # 
  # return(global_graph_annotated)
  
  
}


#' annotate graph by wt and mt graphs
#' @param global_wt_graph wt graph
#' @param global_mt_graph mt graph 
#' @export
get_annotated_graph = function(global_wt_graph,
                               global_mt_graph) {
  
  global_wt_graph_edges = igraph::as_data_frame(global_wt_graph)
  global_mt_graph_edges = igraph::as_data_frame(global_mt_graph)
  
  global_wt_graph_edges = global_wt_graph_edges %>% filter(to %in% global_mt_graph_edges$to == FALSE) #exclude WT edges to nodes that have at least one mutant-supported parent
  
  global_wt_graph_nodes = igraph::as_data_frame(global_wt_graph, what="vertices")
  global_mt_graph_nodes = igraph::as_data_frame(global_mt_graph, what="vertices")
  global_annotated_graph_nodes = data.frame(name=row.names(global_wt_ccm@ccs))
  global_annotated_graph_nodes = left_join(global_annotated_graph_nodes, global_mt_graph_nodes)
  
  mt_only = setdiff(global_mt_graph_edges %>% select(from, to), global_wt_graph_edges %>% select(from, to)) %>% as.data.frame
  
  global_graph_annotated = left_join(global_wt_graph_edges, global_mt_graph_edges)
  global_graph_annotated = global_graph_annotated %>% select(-support)
  global_graph_annotated = rbind(global_graph_annotated,
                                 global_mt_graph_edges %>% inner_join(mt_only))
  global_graph_annotated = igraph::graph_from_data_frame(global_graph_annotated, vertices = global_annotated_graph_nodes)
  global_graph_annotated = platt:::break_cycles_in_state_transition_graph(global_graph_annotated, "total_perturb_path_score_supporting")
  return(global_graph_annotated)
  
}
