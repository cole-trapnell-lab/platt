library(argparse)
library(hooke)
library(tidyr)
library(dplyr)
library(splines)
library(monocle3)
#library(data.table)
# -----------------------------------------------------------------------------

# FIXME: set this to number of CPUs on machine (or in cluster session)
RhpcBLASctl::omp_set_num_threads(8)

outdir = "/Users/maddyduran/UW/Trapnell/zebrafish-atlas-assembly/"
setwd("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/")
# currently using the
# remotes::install_github("pln-team/PLNmodels", "88b18dd875c2037cf592012ca4731474ee2ee567")

# this is failing because of PLN model version PLNmodels_0.11.7-9720
# PLNmodels_0.11.7-9721 on cluster
# remotes::install_github("pln-team/PLNmodels", "03aed9519e219cbbafa26dda01cbd365f6550bc6")

# devtools::load_all("~/OneDrive/UW/Trapnell/hooke_sandwich/")
# remotes::install_github("pln-team/PLNmodels", "28f753a2")

# this doesnt have pln params?
# remotes::install_github("pln-team/PLNmodels", "7c102365792a4d85cb01ef523ec7be1dea7e4f30")

# load data -------------------------------------------------------------------

# wt_cds = readRDS("~/OneDrive/UW/Trapnell/zf-atlas-portal/data/full_gap_hf_ctrl_ref_mito-filt_50k_sub_anno_cds_updated.RDS")
# mt_cds = readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/R_objects/gap16_no-ctrls_projected_major-group-anno_100k_cds.RDS")


#wt_cds = readRDS("/net/trapnell/vol1/home/sanjays/projects/GAP/COMB_GAP/maddy_updates/R_objects/full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")
# wt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")

# copy it down

# wt_cds = readRDS("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/filter_doublets/ref_single_reprocessed.rds")
# wt_cds = readRDS("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/filter_doublets/")

#FIXME: NEED TO COMPUTE STAGING ADJUSTMENTS AND USE THEM INSTEAD OF RAW TIME

# colData(wt_cds)$embryo = colData(wt_cds)$Oligo

source("./R/assembly_utils.R")

# do it on combined data ------------------------------------------------------

# mt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/gap16_no-ctrls_projected_major-group-anno_clean_cds.RDS")
# mt_cds = readRDS("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/filter_doublets/")

# comb_cds = combine_cds(list(wt_cds, mt_cds), keep_reduced_dims = T)

comb_cds = readRDS("filter_doublets/cns_comb_cds_reprocessed.rds")

colData(comb_cds)$timepoint = as.numeric(colData(comb_cds)$timepoint)

control_genotypes = unique(colData(comb_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]

# FIXME: DUMP THIS TO RUN ON ALL MUTANTS
selected_mt_ids = c("tbx16", "tbxta", "noto")

# FIXME: set comb_sample_cds = comb_cds to run on all cells
comb_sample_cds = comb_cds[,colData(comb_cds)$gene_target %in% c(control_genotypes, selected_mt_ids)]
comb_sample_cds = comb_sample_cds[,sample(ncol(comb_sample_cds), 250000)]

rm (comb_cds)

#comb_res = sub_partition_cds(comb_sample_cds)
#saveRDS(comb_res, "comb_res.rds")

colData(comb_sample_cds)$partition = NULL
colData(comb_sample_cds)$cluster = NULL
colData(comb_sample_cds)$cell_state = NULL

# Make sub-umaps/sub-CDS objects automatically:

debug(sub_partition_assembly)
# Constructing sub-UMAP for partition 6_1_1
# No preprocess_method specified, and aligned coordinates have been computed previously. Using preprocess_method = 'Aligned'
# Clustering all cells
# Clustering 6355 cells at resolution = 7.34094486767569e-05
# 64 NAs found in sample group. Dropping NAs.



comb_res = sub_partition_assembly(comb_sample_cds,
                                  ctrl_ids = control_genotypes,
                                  mt_ids = selected_mt_ids,
                                  verbose=TRUE)

comb_res = bind_rows(comb_res)
cell_state_assignments = comb_res %>% select(data) %>% tidyr::unnest(data)
cell_state_assignments$cluster = as.character(cell_state_assignments$cluster)
cell_state_assignments = cell_state_assignments %>% as.data.frame(stringsAsFactors=FALSE)
row.names(cell_state_assignments) = cell_state_assignments$cell
#comb_res = sub_partition_cds(comb_sample_cds)

# Write the cell states back to the main CDS
colData(comb_sample_cds)$cell_state = cell_state_assignments[colData(comb_sample_cds)$cell,]$cell_state

wt_sample_cds = comb_sample_cds[, colData(comb_sample_cds)[["gene_target"]] %in% control_genotypes]

# Need to convert the cluster IDs in each graph to cell_state IDs (which are what we'll use in the final model)
convert_graph_ids = function(state_graph, partition, cluster_to_state_id_tbl, support_col="support"){
  if (igraph::is.igraph(state_graph) == FALSE)
    return (NA)

  cluster_to_state_id_tbl = cluster_to_state_id_tbl %>% filter(grepl(paste("^",partition,"-", sep=""), cell_state))

  #cluster_to_state_id_tbl = colData(cds)[,c("cluster", "cell_state")] %>% as_tibble() %>% distinct()
  state_graph = igraph::as_data_frame(state_graph)
  state_graph$to = as.character(state_graph$to)
  state_graph$from = as.character(state_graph$from)
  state_graph = left_join(state_graph, cluster_to_state_id_tbl, by=c("from"="cluster")) %>%
    dplyr::select(-from) %>% dplyr::rename(from=cell_state)
  state_graph = left_join(state_graph, cluster_to_state_id_tbl, by=c("to"="cluster")) %>%
    dplyr::select(-to) %>% dplyr::rename(to=cell_state)
  state_graph = state_graph[,c("from", "to", support_col)]
  state_graph = igraph::graph_from_data_frame(state_graph, directed=TRUE)
}
#debug(convert_graph_ids)

c_to_s_id_tbl = cell_state_assignments %>% select(cluster, cell_state) %>% distinct()

comb_res = comb_res %>%
  mutate(mt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = mt_graph,
                                          .y=partition,
                                          cluster_to_state_id_tbl=c_to_s_id_tbl,
                                          support_col=NULL))

comb_res = comb_res %>%
  mutate(wt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = wt_graph,
                                          .y=partition,
                                          c_to_s_id_tbl))




# Save everything again so we capture all the cell count models and state graphs
#saveRDS(res_graphs, paste(outdir, "wt_res_graphs_cds.rds", sep="/"))

# Build a whitelist of graph edges from all the sub-CDS state graphs
global_wt_graph_edge_whitelist = do.call(igraph::union, comb_res %>% filter(is.na(wt_graph_converted) == FALSE) %>% pull(wt_graph_converted))
global_wt_graph_edge_whitelist = igraph::as_data_frame(global_wt_graph_edge_whitelist)
global_wt_graph_edge_whitelist = global_wt_graph_edge_whitelist %>% select(from, to) %>% distinct()

gc()

# Fit a single cell count model to all the cell states at once
global_wt_ccm = fit_wt_model(wt_sample_cds,
                             sample_group = "embryo",
                             cell_group = "cell_state",
                             main_model_formula_str = NULL,
                             start_time = 18,
                             stop_time = 48,
                             #nuisance_model_formula_str = "~expt",
                             ctrl_ids = control_genotypes,
                             sparsity_factor = 0.01,
                             perturbation_col = "gene_target",
                             keep_cds = FALSE)

# Learn a single graph on all states at once (using the subgraphs as a whitelist/prior)
global_wt_graph = assemble_wt_graph(wt_sample_cds,
                                    global_wt_ccm,
                                    sample_group = "embryo",
                                    cell_group = "cell_state",
                                    main_model_formula_str = NULL,
                                    start_time = 18,
                                    stop_time = 48,
                                    #nuisance_model_formula_str = "~expt",
                                    ctrl_ids = control_genotypes,
                                    sparsity_factor = 0.01,
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_wt_graph_edge_whitelist,
                                    verbose=TRUE)


# Plot the giant graph
plot_state_graph_annotations(global_wt_ccm, global_wt_graph, color_nodes_by = "timepoint", group_nodes_by="cell_type_sub",
                             edge_weights = "support",
                             hide_unlinked_nodes = FALSE)
ggplot2::ggsave(paste(outdir, "global_wt_graph.png", sep="/"), width=20, height=20)


global_mt_graph_edge_whitelist = do.call(igraph::union, comb_res %>% filter(is.na(mt_graph_converted) == FALSE) %>% pull(mt_graph_converted))
global_mt_graph_edge_whitelist = igraph::as_data_frame(global_mt_graph_edge_whitelist)
global_mt_graph_edge_whitelist = global_mt_graph_edge_whitelist %>% select(from, to) %>% distinct()

global_perturb_models_tbl = fit_mt_models(comb_sample_cds,
                                          sample_group = "embryo",
                                          cell_group = "cell_state",
                                          main_model_formula_str = NULL,
                                          start_time = 18,
                                          stop_time = 48,
                                          #nuisance_model_formula_str = "~expt",
                                          ctrl_ids = control_genotypes,
                                          mt_ids = selected_mt_ids,
                                          sparsity_factor = 0.01,
                                          perturbation_col = "gene_target",
                                          keep_cds=FALSE)

global_mt_graph = assemble_mt_graph(global_wt_ccm,
                                    global_perturb_models_tbl,
                                    start_time = 18,
                                    stop_time = 48,
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_mt_graph_edge_whitelist,
                                    verbose=TRUE)

plot_state_graph_annotations(global_wt_ccm, global_mt_graph,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "max_path_score_supporting",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_mt_graph.png", sep="/"), width=30, height=15)

saveRDS(global_wt_ccm, paste(outdir, "global_wt_ccm.rds", sep="/"))
saveRDS(global_perturb_models_tbl, paste(outdir, "global_perturb_models_tbl.rds", sep="/"))
saveRDS(global_mt_graph, paste(outdir, "global_mt_graph.rds", sep="/"))


#####################################

get_subset_by_genotype = function(cds_list, index, perturbation_col = "gene_target", genotypes=NULL) {

  cds = cds_list[[index]]

  # if genotypes is not specified, just return wild-type and controls
  if (is.null(genotypes)) {
    genotypes = unique(colData(cds)[[perturbation_col]])
    genotypes = genotypes[grepl("wt|ctrl", genotypes)]
  }

  cds = cds[, colData(cds)[[perturbation_col]] %in% genotypes]


  return(cds)
}




# Let's process the wt data first:
res_graphs = tibble("index" = 1:length(comb_res)) %>%
  mutate(cds = purrr::map(.f = get_subset_by_genotype,
                          .x = index,
                          cds_list = comb_res))

# res_graphs = readRDS("res_graphs.rds")
# fit_wt_model(res_graphs$cds[[1]],
#              sample_group = "embryo",
#              cell_group = "cluster",
#              main_model_formula_str = NULL,
#              start_time = 18,
#              stop_time = 48,
#              #nuisance_model_formula_str = "~expt",
#              ctrl_ids = control_genotypes,
#              sparsity_factor = 0.01,
#              perturbation_col = "gene_target")



# Fit cell count models to each subset of cells
res_graphs = res_graphs %>%
  mutate(wt_ccm = purrr::map(.f = fit_wt_model,
                             .x = cds,
                             sample_group = "embryo",
                             cell_group = "cluster",
                             main_model_formula_str = NULL,
                             start_time = 18,
                             stop_time = 48,
                             #nuisance_model_formula_str = "~expt",
                             ctrl_ids = control_genotypes,
                             sparsity_factor = 0.01,
                             perturbation_col = "gene_target"))

# Learn a state graph for each subset of cells
res_graphs = res_graphs %>%
  mutate(wt_graph = purrr::map2(.f = assemble_wt_graph,
                                .x = cds,
                                .y = wt_ccm,
                                sample_group = "embryo",
                                cell_group = "cluster",
                                main_model_formula_str = NULL,
                                start_time = 18,
                                stop_time = 48,
                                #nuisance_model_formula_str = "~expt",
                                ctrl_ids = control_genotypes,
                                sparsity_factor = 0.01,
                                perturbation_col = "gene_target"))


# Now let's plot all the state graphs
res_graphs = res_graphs %>%
  mutate(state_graph_plot = purrr::map2(.f = purrr::possibly(plot_state_graph_annotations_wrapper, NA_real_),
                                        .x = wt_ccm,
                                        .y = wt_graph))
plotted_graphs = res_graphs %>% filter(is.na(state_graph_plot) == FALSE)
cowplot::plot_grid(plotlist = plotted_graphs$state_graph_plot, labels=plotted_graphs$index)
ggplot2::ggsave(paste(outdir, "wt_all_sub_state_graphs_by_cell_type_sub.png", sep="/"), width=20, height=20)

# Need to convert the cluster IDs in each graph to cell_state IDs (which are what we'll use in the final model)
convert_graph_ids = function(cds, state_graph, support_col="support"){
  if (igraph::is.igraph(state_graph) == FALSE)
    return (NA)

  cluster_to_state_id_tbl = tibble(cell_state = as.character(colData(cds)$cell_state),
                                   cluster = as.character(clusters(cds))) %>%
    distinct()

  #cluster_to_state_id_tbl = colData(cds)[,c("cluster", "cell_state")] %>% as_tibble() %>% distinct()
  state_graph = igraph::as_data_frame(state_graph)
  state_graph$to = as.character(state_graph$to)
  state_graph$from = as.character(state_graph$from)
  state_graph = left_join(state_graph, cluster_to_state_id_tbl, by=c("from"="cluster")) %>%
    dplyr::select(-from) %>% dplyr::rename(from=cell_state)
  state_graph = left_join(state_graph, cluster_to_state_id_tbl, by=c("to"="cluster")) %>%
    dplyr::select(-to) %>% dplyr::rename(to=cell_state)
  state_graph = state_graph[,c("from", "to", support_col)]
  state_graph = igraph::graph_from_data_frame(state_graph, directed=TRUE)
}
#debug(convert_graph_ids)

res_graphs = res_graphs %>%
  mutate(wt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = cds,
                                          .y = wt_graph))

# Save everything again so we capture all the cell count models and state graphs
#saveRDS(res_graphs, paste(outdir, "wt_res_graphs_cds.rds", sep="/"))

# Build a whitelist of graph edges from all the sub-CDS state graphs
global_wt_graph_edge_whitelist = do.call(igraph::union, res_graphs %>% filter(is.na(wt_graph_converted) == FALSE) %>% pull(wt_graph_converted))
global_wt_graph_edge_whitelist = igraph::as_data_frame(global_wt_graph_edge_whitelist)
global_wt_graph_edge_whitelist = global_wt_graph_edge_whitelist %>% select(from, to) %>% distinct()

gc()

# Fit a single cell count model to all the cell states at once
global_wt_ccm = fit_wt_model(wt_sample_cds,
                             sample_group = "embryo",
                             cell_group = "cell_state",
                             main_model_formula_str = NULL,
                             start_time = 18,
                             stop_time = 48,
                             #nuisance_model_formula_str = "~expt",
                             ctrl_ids = control_genotypes,
                             sparsity_factor = 0.01,
                             perturbation_col = "gene_target")

# Learn a single graph on all states at once (using the subgraphs as a whitelist/prior)
global_wt_graph = assemble_wt_graph(wt_sample_cds,
                                    global_wt_ccm,
                                    sample_group = "embryo",
                                    cell_group = "cell_state",
                                    main_model_formula_str = NULL,
                                    start_time = 18,
                                    stop_time = 48,
                                    #nuisance_model_formula_str = "~expt",
                                    ctrl_ids = control_genotypes,
                                    sparsity_factor = 0.01,
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_wt_graph_edge_whitelist)

# Plot the giant graph
plot_state_graph_annotations(global_wt_ccm, global_wt_graph, color_nodes_by = "timepoint", group_nodes_by="cell_type_sub",
                             edge_weights = "support",
                             hide_unlinked_nodes = FALSE)
ggplot2::ggsave(paste(outdir, "global_wt_graph.png", sep="/"), width=20, height=20)

# Save the global model and the state graph
saveRDS(global_wt_ccm, paste(outdir, "global_wt_ccm.rds", sep="/"))
saveRDS(global_wt_graph, paste(outdir, "global_wt_graph.rds", sep="/"))

# Now let's start analyzing the sub-CDS objects
res_graphs = res_graphs %>%
  mutate(cds = purrr::map(.f = subset_index,
                          .x = index,
                          cds_list = comb_res))

# make sub-CDS UMAP plots:
res_graphs = res_graphs %>%
  mutate(umap_cell_states = purrr::map(plot_cells, .x = cds, color_cells_by="cell_state"),
         umap_time = purrr::map(plot_cells, .x = cds, color_cells_by="timepoint"),
         umap_cell_type_sub = purrr::map(plot_cells, .x = cds, color_cells_by="cell_type_sub")
  )

plotted_umaps = res_graphs %>% filter(is.na(umap_cell_states) == FALSE)
cowplot::plot_grid(plotlist = plotted_umaps$umap_cell_states)
ggplot2::ggsave(paste(outdir, "comb_cds_all_sub_umaps_by_cell_state.png", sep="/"), width=12, height=12)

# This one needs some cleanup to look good:
plotted_umaps = res_graphs %>% filter(is.na(umap_time) == FALSE)
cowplot::plot_grid(plotlist = plotted_umaps$umap_time)
ggplot2::ggsave(paste(outdir, "comb_cds_all_sub_umaps_by_time.png", sep="/"), width=12, height=12)

plotted_umaps = res_graphs %>% filter(is.na(umap_cell_type_sub) == FALSE)
cowplot::plot_grid(plotlist = plotted_umaps$umap_cell_type_sub)
ggplot2::ggsave(paste(outdir, "comb_cds_all_sub_umaps_by_cell_type_sub.png", sep="/"), width=20, height=20)

#saveRDS(res_graphs, "comb_res_graphs_cds.rds")

#res_graphs = res_graphs %>% head(3)


# res_graphs %>% head(3) %>%
#   mutate(perturb_models_tbl = purrr::map(.f = fit_mt_models,
#                                          .x = cds,
#                                          sample_group = "embryo",
#                                          cell_group = "cluster",
#                                          main_model_formula_str = NULL,
#                                          start_time = 18,
#                                          stop_time = 48,
#                                          #nuisance_model_formula_str = "~expt",
#                                          ctrl_ids = control_genotypes,
#                                          #mt_ids = selected_mt_ids,
#                                          sparsity_factor = 0.01,
#                                          perturbation_col = "gene_target"))

res_graphs = res_graphs %>%
  mutate(perturb_models_tbl = purrr::map(.f = fit_mt_models,
                                         .x = cds,
                                         sample_group = "embryo",
                                         cell_group = "cluster",
                                         main_model_formula_str = NULL,
                                         start_time = 18,
                                         stop_time = 48,
                                         #nuisance_model_formula_str = "~expt",
                                         ctrl_ids = control_genotypes,
                                         #mt_ids = selected_mt_ids,
                                         sparsity_factor = 0.01,
                                         perturbation_col = "gene_target"))

res_graphs = res_graphs %>%
  mutate(mt_graph = purrr::map2(.f = assemble_mt_graph,
                                .x = wt_ccm,
                                .y = perturb_models_tbl,
                                start_time = 18,
                                stop_time = 48,
                                perturbation_col = "gene_target"))

res_graphs = res_graphs %>%
  mutate(state_graph_plot = purrr::map2(.f = purrr::possibly(plot_state_graph_annotations_wrapper, NA_real_),
                                        .x = wt_ccm,
                                        .y = mt_graph,
                                        edge_weights="total_path_score_supporting"))

plotted_graphs = res_graphs %>% filter(is.na(state_graph_plot) == FALSE)
cowplot::plot_grid(plotlist = plotted_graphs$state_graph_plot, labels=plotted_graphs$index)
ggplot2::ggsave(paste(outdir, "mt_all_sub_state_graphs_by_cell_type_sub.png", sep="/"), width=20, height=20)


res_graphs = res_graphs %>%
  mutate(mt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = cds,
                                          .y = mt_graph,
                                          support_col="total_path_score_supporting"))

saveRDS(res_graphs, "comb_res_graphs_cds.rds")

global_mt_graph_edge_whitelist = do.call(igraph::union, res_graphs %>% filter(is.na(mt_graph_converted) == FALSE) %>% pull(mt_graph_converted))
global_mt_graph_edge_whitelist = igraph::as_data_frame(global_mt_graph_edge_whitelist)
global_mt_graph_edge_whitelist = global_mt_graph_edge_whitelist %>% select(from, to) %>% distinct()

global_perturb_models_tbl = fit_mt_models(comb_sample_cds,
                                          sample_group = "embryo",
                                          cell_group = "cell_state",
                                          main_model_formula_str = NULL,
                                          start_time = 18,
                                          stop_time = 48,
                                          #nuisance_model_formula_str = "~expt",
                                          ctrl_ids = control_genotypes,
                                          mt_ids = selected_mt_ids,
                                          sparsity_factor = 0.01,
                                          perturbation_col = "gene_target")

global_mt_graph = assemble_mt_graph(global_wt_ccm,
                                    global_perturb_models_tbl,
                                    start_time = 18,
                                    stop_time = 48,
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_mt_graph_edge_whitelist)

plot_state_graph_annotations(global_wt_ccm, global_mt_graph, color_nodes_by = "timepoint", group_nodes_by="cell_type_sub",
                             edge_weights = NULL,
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_mt_graph.png", sep="/"), width=20, height=20)

saveRDS(global_wt_ccm, paste(outdir, "global_wt_ccm.rds", sep="/"))
saveRDS(global_perturb_models_tbl, paste(outdir, "global_perturb_models_tbl.rds", sep="/"))
saveRDS(global_mt_graph, paste(outdir, "global_mt_graph.rds", sep="/"))

###### Not needed, for development, curiosity


loss_timing_wrapper = function(perturbation_ccm, perturb_time_window){
  return(hooke:::estimate_loss_timing(perturbation_ccm,
                                      start_time=min(perturb_time_window$start_time),
                                      stop_time=min(perturb_time_window$stop_time),
                                      interval_step = 1,
                                      interval_col="timepoint",
                                      log_abund_detection_thresh=-5,
                                      q_val = 0.01,
                                      delta_log_abund_loss_thresh=0))
}

#debug(loss_timing_wrapper)
# FIXME: shouldn't use hardcoded knockout here

loss_table = global_perturb_models_tbl %>%
  filter (is.na(perturb_ccm) == FALSE) %>%
  mutate(perturbation_losses = purrr::map2(.f = loss_timing_wrapper,
                                           .x=perturb_ccm,
                                           .y=perturb_time_window))

loss_table = loss_table %>% tidyr::unnest(perturbation_losses)

