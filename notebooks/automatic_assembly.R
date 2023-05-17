library(argparse)
library(hooke)
library(tidyr)
library(dplyr)
library(splines)
library(monocle3)
#library(data.table)
library(future)
library(purrr)
library(ggplot2)
library(ggtext)

plan(multicore)
options(future.fork.multithreading.enable=FALSE)
options(future.globals.maxSize = 8192 * 1024^2) # 8GB

assembly_start_time = 18
assembly_stop_time = 72
assembly_backend = "nlopt"


source("R/assembly_utils.R")

detect_outlier_cells = function(cds, prefix , k=10) {

  # build annoy index
  cds = make_cds_nn_index(cds, reduction_method = "UMAP")

  # save it
  # save_transform_models(cds, paste0(prefix, "mt_umap_nn_models"))

  # throw them away or just NA everything?

  query_ann = cds@reduce_dim_aux$UMAP$nn_index$annoy$nn_index$annoy_index
  query_dims = reducedDims(cds)[["UMAP"]]
  query_res = uwot:::annoy_search(query_dims, k = k + 1, ann = query_ann)

  df = lapply(1:nrow(query_dims), function(i) {
    # remember to ignore first index
    neighbor_indices = query_res$idx[i,2:(k+1)]
    neighbor_dists = query_res$dist[i,2:(k+1)]
    data.frame("idx"=i,
               "mean_nn_dist" = mean(neighbor_dists),
               "median_nn_dist" = median(neighbor_dists),
               "min_nn_dist" = min(neighbor_dists)
               )
  }) %>% bind_rows()

  colData(cds)$mean_nn_dist = df$mean_nn_dist
  colData(cds)$median_nn_dist = df$median_nn_dist
  colData(cds)$min_nn_dist = df$min_nn_dist

  # fig = plot_cells_3d(cds, color_cells_by = "mean")
  # saveWidget(fig, paste0("plots/", prefix,"_mean_nn_dist.html"))

  # p = colData(cds) %>% as.data.frame %>% ggplot(aes(mean)) + geom_density() + geom_vline(xintercept = 0.1)
  # ggsave(p, paste0("plots", prefix, "_mean_nn_dist.png"))

  return(cds)
}

# -----------------------------------------------------------------------------

# FIXME: set this to number of CPUs on machine (or in cluster session)
#RhpcBLASctl::omp_set_num_threads(10)

num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
print (paste("running assembly with", num_threads, "threads"))

Sys.setenv("OMP_NUM_THREADS" = 1)
Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

outdir = "//Users/coletrap/tmp_output"

#setwd("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/")

setwd("/Users/coletrap/Google Drive/My Drive/develop/zebrafish-atlas-assembly")


# load data -------------------------------------------------------------------

comb_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/full_cds_combined.rds")
colData(comb_cds)$timepoint = as.numeric(colData(comb_cds)$timepoint)


control_genotypes = unique(colData(comb_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]
selected_mt_ids = NULL

colData(comb_sample_cds)$sample = NULL
colData(comb_sample_cds)$mean = NULL


# Old loading section:
# -------------------------------------------------------------------

# wt_cds = readRDS("~/OneDrive/UW/Trapnell/zf-atlas-portal/data/full_gap_hf_ctrl_ref_mito-filt_50k_sub_anno_cds_updated.RDS")
# mt_cds = readRDS("~/OneDrive/UW/Trapnell/gap-notebook-ct-1/R_objects/gap16_no-ctrls_projected_major-group-anno_100k_cds.RDS")

#setwd("/net/trapnell/vol1/home/duran/zebrafish-atlas-assembly")


#wt_cds = readRDS("/net/trapnell/vol1/home/sanjays/projects/GAP/COMB_GAP/maddy_updates/R_objects/full_gap_hf_ctrl_ref_mito-filt_1.25M_model-update_anno_cds.RDS")
wt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/ref_doublet_filtered_reprocessed.rds")

#FIXME: NEED TO COMPUTE STAGING ADJUSTMENTS AND USE THEM INSTEAD OF RAW TIME

#colData(wt_cds)$embryo = colData(wt_cds)$Oligo
colData(wt_cds)$timepoint = as.numeric(colData(wt_cds)$timepoint)

mt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/gap16_doublet_filtered_projected_cds.rds")

colData(mt_cds)$timepoint = as.numeric(colData(mt_cds)$timepoint)

control_genotypes = unique(colData(wt_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]

selected_mt_ids = NULL

# FIXME: DUMP THIS TO RUN ON ALL MUTANTS
#selected_mt_ids = c("tbx16", "tbx16-msgn1", "tbxta", "noto", "mafba", "egr2b")
#selected_mt_ids = c("tbx16", "foxi1", "phox2a")
#mt_cds = mt_cds[, colData(mt_cds)[["gene_target"]] %in% selected_mt_ids]

gc()

comb_cds = combine_cds(list(wt_cds, mt_cds), keep_reduced_dims = T)

# FIXME: we should really remove mean as a column
comb_cds = detect_outlier_cells(comb_cds)

min_dist_quantiles = quantile(colData(comb_cds)$min_nn_dist, probs = seq(0, 1, 0.01), na.rm=TRUE)
min_dist_thresh = min_dist_quantiles["99%"]
colData(comb_cds)$outlier = colData(comb_cds)$min_nn_dist > min_dist_thresh
comb_cds = comb_cds[,colData(comb_cds)$outlier == FALSE]

#qplot(colData(comb_cds)$mean_nn_dist, log="x")

#colData(comb_cds)$timepoint = as.numeric(colData(comb_cds)$timepoint)

rm(mt_cds)
rm(wt_cds)
gc()

##### End old loading section


# staging_df = colData(comb_cds) %>% as.data.frame %>% select(cell, embryo, expt, timepoint, mean_nn_time) %>% as_tibble()
# staging_df = staging_df %>% mutate(timepoint = as.numeric(timepoint))
#
# staging_model = lm(mean_nn_time ~ as.numeric(timepoint) * expt, data=staging_df)
# staging_df$predicted_timepoint = predict(staging_model, newdata=staging_df)
#
# colData(comb_cds)$adjusted_timepoint = staging_df$predicted_timepoint

#
# # FIXME: set comb_sample_cds = comb_cds to run on all cells

comb_sample_cds = comb_cds

#
# # Comment the following block when doing the full dataset
# ####
# comb_sample_cds = comb_cds[,colData(comb_cds)$gene_target %in% c(control_genotypes, selected_mt_ids)]
# comb_sample_cds = comb_sample_cds[,sample(ncol(comb_sample_cds), 250000)]
# comb_sample_cds = suppressMessages(suppressWarnings(preprocess_cds(comb_sample_cds))) %>%
#   align_cds(residual_model_formula_str = "~log.n.umi") %>%
#   suppressMessages(suppressWarnings(reduce_dimension(max_components = max_components,
#                                                      preprocess_method="Aligned",
#                                                      umap.fast_sgd=TRUE,
#                                                      cores=num_threads)))
# ####

#comb_sample_cds = comb_cds

rm (comb_cds)
gc()

# TODO: possibly we should explicitly model this in formulae as well?
bycatch_cell_types = c("periderm")
comb_sample_cds = comb_sample_cds[,colData(comb_sample_cds)$cell_type_broad %in% bycatch_cell_types == FALSE]

#comb_sample_ccs = new_cell_count_set(comb_sample_cds, sample_group = "embryo", cell_group="cell_type_broad")
#comb_sample_ccs = preprocess_cds(comb_sample_ccs)
#comb_sample_ccs = align_cds(comb_sample_ccs, alignment_group="expt")

#comb_sample_ccs = reduce_dimension(comb_sample_ccs)
#plot_cells(comb_sample_ccs, color_cells_by="timepoint")

#comb_res = sub_partition_cds(comb_sample_cds)
#saveRDS(comb_res, "comb_res.rds")

colData(comb_sample_cds)$partition = NULL
colData(comb_sample_cds)$cluster = NULL
colData(comb_sample_cds)$cell_state = NULL

#comb_sample_cds = cluster_cells(comb_sample_cds, resolution=1e-4)



# Calculate per-embryo size factors so we can use those in sub-assemblies:
global_ccs = new_cell_count_set(comb_sample_cds,
                                sample_group = "embryo",
                                cell_group = "cell_type_broad")
global_embryo_size_factors = size_factors(global_ccs)

rm(global_ccs)
gc()

# Make sub-umaps/sub-CDS objects automatically:
system.time({comb_res = sub_partition_assembly(comb_sample_cds,
                                  ctrl_ids = control_genotypes,
                                  mt_ids = selected_mt_ids,
                                  verbose=FALSE,
                                  start_time=assembly_start_time,
                                  stop_time=assembly_stop_time,
                                  interval_col="timepoint",
                                  nuisance_model_formula_str="~expt",
                                  min_res=5e-7,
                                  max_res=5e-4,
                                  num_threads=num_threads,
                                  backend=assembly_backend,
                                  expts_excluded_from_assembly=c("GAP13"),
                                  embryo_size_factors=global_embryo_size_factors)})

comb_res = bind_rows(comb_res)
cell_state_assignments = comb_res %>% select(data) %>% tidyr::unnest(data)
cell_state_assignments$cluster = as.character(cell_state_assignments$cluster)
cell_state_assignments = cell_state_assignments %>% distinct() %>% as.data.frame(stringsAsFactors=FALSE)
row.names(cell_state_assignments) = cell_state_assignments$cds_row_id
#comb_res = sub_partition_cds(comb_sample_cds)

cell_state_assignments = cell_state_assignments %>% mutate(partition = stringr::str_split_fixed(cell_state, "-", n=2)[,1])
cell_state_assignmen_summary = cell_state_assignments %>% group_by(partition) %>%
  summarize(total_cells = n(),
            num_states = length(unique(cell_state)))

plotted_graphs = comb_res %>% filter(is.na(cell_plot_state) == FALSE)
png(paste(outdir, "all_sub_cds_umaps_by_cell_state.png", sep="/"), width=20, height=20, units="in", res=300)
cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_state, labels=plotted_graphs$partition)
dev.off()
#ggplot2::ggsave(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20)

plotted_graphs = comb_res %>% filter(is.na(cell_plot_type) == FALSE)
png(paste(outdir, "all_sub_cds_umaps_by_cell_type.png", sep="/"), width=20, height=20, units="in", res=300)
cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_type, labels=plotted_graphs$partition)
dev.off()
#ggplot2::ggsave(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20)

plotted_graphs = comb_res %>% filter(is.na(cell_plot_time) == FALSE)
png(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20, units="in", res=300)
cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_time, labels=plotted_graphs$partition)
dev.off()
#ggplot2::ggsave(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20)


plotted_graphs = comb_res %>% filter(is.null(wt_state_graph_plot) == FALSE && class(wt_graph) == "list")
pdf(paste(outdir, "sub_cds_wt_state_graphs.pdf", sep="/"), width=20, height=20)
cowplot::plot_grid(plotlist = plotted_graphs$wt_state_graph_plot, labels=plotted_graphs$partition)
dev.off()
#ggplot2::ggsave(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20)


plotted_graphs = comb_res %>% filter(is.null(mt_state_graph_plot) == FALSE && class(mt_graph) == "list")
pdf(paste(outdir, "sub_cds_mt_state_graphs.pdf", sep="/"), width=30, height=30)
cowplot::plot_grid(plotlist = plotted_graphs$mt_state_graph_plot, labels=plotted_graphs$partition)
dev.off()
#ggplot2::ggsave(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20)



plotted_graphs = comb_res %>% filter(is.null(wt_graph_blacklist_plot) == FALSE && class(mt_graph_blacklist) == "list")
pdf(paste(outdir, "sub_cds_wt_blacklist_state_graphs.pdf", sep="/"), width=20, height=20)
cowplot::plot_grid(plotlist = plotted_graphs$wt_graph_blacklist_plot, labels=plotted_graphs$partition)
dev.off()
#ggplot2::ggsave(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20)


plotted_graphs = comb_res %>% filter(is.null(mt_state_graph_plot) == FALSE && class(mt_graph_blacklist) == "list")
pdf(paste(outdir, "sub_cds_mt_blacklist_state_graphs.pdf", sep="/"), width=30, height=30)
cowplot::plot_grid(plotlist = plotted_graphs$mt_state_graph_plot, labels=plotted_graphs$partition)
dev.off()
#ggplot2::ggsave(paste(outdir, "all_sub_cds_umaps_by_cell_time.png", sep="/"), width=20, height=20)


rm (plotted_graphs)

comb_res$cell_plot_time = NULL
comb_res$cell_plot_type = NULL
comb_res$cell_plot_state = NULL
#comb_res$wt_state_graph_plot = NULL
#comb_res$mt_state_graph_plot = NULL


gc()


# Write the cell states back to the main CDS

colData(comb_sample_cds)$cell_state = cell_state_assignments[row.names(colData(comb_sample_cds) %>% as.data.frame()),]$cell_state

#wt_sample_cds = comb_sample_cds[, colData(comb_sample_cds)[["gene_target"]] %in% control_genotypes]

# admitted_wt_cell_states = colData(wt_sample_cds) %>% as.data.frame () %>%
#   select(cell_state) %>%
#   group_by(cell_state) %>%
#   summarize(total_cells = n()) %>%
#   filter (total_cells > 10000) %>% pull(cell_state) %>% unique
# wt_sample_cds = wt_sample_cds[,colData(wt_sample_cds)$cell_state %in% admitted_wt_cell_states]

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

comb_res = comb_res %>%
  mutate(wt_graph_blacklist_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = wt_graph_blacklist,
                                          .y=partition,
                                          c_to_s_id_tbl))

# comb_res = comb_res %>%
#   mutate(mt_graph_blacklist_converted = purrr::map2(.f = convert_graph_ids,
#                                                     .x = mt_graph_blacklist,
#                                                     .y=partition,
#                                                     c_to_s_id_tbl))




# Save everything again so we capture all the cell count models and state graphs
#saveRDS(res_graphs, paste(outdir, "wt_res_graphs_cds.rds", sep="/"))

# Build a whitelist of graph edges from all the sub-CDS state graphs
global_wt_graph_edge_whitelist = do.call(igraph::union, comb_res %>% filter(is.na(wt_graph_converted) == FALSE) %>% pull(wt_graph_converted))
global_wt_graph_edge_whitelist = igraph::as_data_frame(global_wt_graph_edge_whitelist)
global_wt_graph_edge_whitelist = global_wt_graph_edge_whitelist %>% select(from, to) %>% distinct()


#rm(comb_res)

gc()

# Fit a single cell count model to all the cell states at once
# global_wt_ccm = fit_wt_model(comb_sample_cds,
#                                           sample_group = "embryo",
#                                           cell_group = "cell_state",
#                                           #main_model_formula_str = NULL,
#                                           #main_model_formula_str = "~1",
#                                           main_model_formula_str = "~ ns( timepoint , knots= c(36,54) ) + expt + ns( timepoint , knots= c(36,54) ):expt",
#                                           start_time = assembly_start_time,
#                                           stop_time = assembly_stop_time,
#                                           interval_col="timepoint",
#                                           vhat_method="bootstrap",
#                                           num_time_breaks=4,
#                                           #nuisance_model_formula_str = "~expt",
#                                           ctrl_ids = control_genotypes,
#                                           sparsity_factor = 0.01,
#                                           perturbation_col = "gene_target",
#                                           edge_whitelist = global_wt_graph_edge_whitelist,
#                                           #base_penalty=1000,
#                                           keep_cds = FALSE,
#                                           num_threads=num_threads,
#                                           backend=assembly_backend,
#                                           verbose=TRUE,
#                                           penalize_by_distance=TRUE,
#                                           pln_num_penalties=30)

# Fit a single wild-type cell count timeseries model to all the cell states at once
system.time({global_wt_ccm = fit_wt_model(comb_sample_cds,
                             sample_group = "embryo",
                             cell_group = "cell_state",
                             #main_model_formula_str = NULL,
                             #main_model_formula_str = "~1",
                             #main_model_formula_str = "~ns(timepoint, df=5)",
                             start_time = assembly_start_time,
                             stop_time = assembly_stop_time,
                             interval_col="timepoint",
                             vhat_method="bootstrap",
                             num_time_breaks=4,
                             nuisance_model_formula_str = "~expt",
                             ctrl_ids = control_genotypes,
                             sparsity_factor = 0.01,
                             perturbation_col = "gene_target",
                             edge_whitelist = global_wt_graph_edge_whitelist,
                             #base_penalty=1000,
                             keep_cds = FALSE,
                             num_threads=num_threads,
                             backend=assembly_backend,
                             verbose=TRUE,
                             penalize_by_distance=TRUE,
                             pln_num_penalties=30)})



# Learn a single graph on all states at once (using the subgraphs as a whitelist/prior)
global_wt_graph = assemble_wt_graph(comb_sample_cds,
                                    global_wt_ccm,
                                    sample_group = "embryo",
                                    cell_group = "cell_state",
                                    main_model_formula_str = NULL,
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="timepoint",
                                    #nuisance_model_formula_str = "~expt",
                                    ctrl_ids = control_genotypes,
                                    sparsity_factor = 0.01,
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_wt_graph_edge_whitelist,
                                    verbose=TRUE)


# Plot the giant graph
plot_state_graph_annotations(global_wt_ccm, global_wt_graph, color_nodes_by = "timepoint", group_nodes_by="cell_type_sub",
                             edge_weights = "support",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph.png", sep="/"), width=30, height=10)

# Plot it by cell type broad just because the labels can get crazy
plot_state_graph_annotations(global_wt_ccm, global_wt_graph, color_nodes_by = "timepoint", group_nodes_by="cell_type_broad",
                             edge_weights = "support",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph_broad.png", sep="/"), width=30, height=10)

# What, if any edges, were in the timeseries subassemblies that did NOT wind up
# in the global assembly:
wt_edges_not_included_in_global_graph = setdiff (global_wt_graph_edge_whitelist %>% select(from, to),
         global_wt_graph %>% igraph::as_data_frame() %>% select(from, to))
print (paste("Of ", nrow(global_wt_graph_edge_whitelist), "edges from WT subassemblies",
             nrow(wt_edges_not_included_in_global_graph), "not included in global WT graph"))

saveRDS(global_wt_ccm, paste(outdir, "global_wt_ccm.rds", sep="/"))
saveRDS(global_wt_graph, paste(outdir, "global_wt_graph.rds", sep="/"))


# Fit cell count models to each mutant vs control
global_perturb_models_tbl = fit_mt_models(comb_sample_cds,
                                          sample_group = "embryo",
                                          cell_group = "cell_state",
                                          main_model_formula_str = NULL,
                                          start_time = assembly_start_time,
                                          stop_time = assembly_stop_time,
                                          interval_col="timepoint",
                                          num_time_breaks=3,
                                          #nuisance_model_formula_str = "~expt",
                                          ctrl_ids = control_genotypes,
                                          #mt_ids = c("tbx16", "mafba"), #selected_mt_ids,
                                          sparsity_factor = 0.01,
                                          perturbation_col = "gene_target",
                                          keep_cds=FALSE,
                                          num_threads=num_threads,
                                          backend=assembly_backend,
                                          vhat_method="bootstrap",
                                          penalize_by_distance=TRUE)


## Debugging stuff:
# selected_perturb_models_tbl = fit_mt_models(comb_sample_cds,
#                                           sample_group = "embryo",
#                                           cell_group = "cell_state",
#                                           main_model_formula_str = NULL,
#                                           start_time = assembly_start_time,
#                                           stop_time = assembly_stop_time,
#                                           interval_col="timepoint",
#                                           num_time_breaks=3,
#                                           #nuisance_model_formula_str = "~expt",
#                                           ctrl_ids = control_genotypes,
#                                           mt_ids = c("mafba", "hoxb1b", "tbx16-msgn1"), #selected_mt_ids,
#                                           sparsity_factor = 0.01,
#                                           perturbation_col = "gene_target",
#                                           keep_cds=FALSE,
#                                           num_threads=num_threads,
#                                           backend=assembly_backend,
#                                           vhat_method="bootstrap",
#                                           penalize_by_distance=TRUE,
#                                           independent_spline_for_ko=TRUE,
#                                           num_bootstraps=10)
#
# selected_perturb = "phox2a"
# selected_perturb_results = global_perturb_models_tbl %>% filter(selected_perturb == perturb_name)
# selected_ccm = selected_perturb_results$perturb_ccm[[1]]
#
#
# #debug(plot_cell_type_perturb_kinetics)
# source("R/assembly_utils.R")
# #debug(plot_cell_type_perturb_kinetics)
# kinetics_plot = plot_cell_type_perturb_kinetics(selected_ccm,
#                                                 cell_groups = selected_cell_groups,
#                                                 start_time=assembly_start_time,
#                                                 stop_time=assembly_stop_time,
#                                                 interval_col="timepoint",
#                                                 expt="GAP16") +
#   monocle3:::monocle_theme_opts()
#
#
# pdf(paste(outdir, "muscle_selected_ko_cell_type_distribution_over_time.pdf", sep="/"), height=10, width=6)
# print(kinetics_plot)
# #ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)
# dev.off()

# Build a whitelist of edges by collecting the edges in the subassemblies from the mutants
sub_mt_graph_union = do.call(igraph::union, comb_res %>% filter(is.na(mt_graph_converted) == FALSE) %>% pull(mt_graph_converted))
global_mt_graph_edge_whitelist = igraph::as_data_frame(sub_mt_graph_union)
global_mt_graph_edge_whitelist = global_mt_graph_edge_whitelist %>% select(from, to) %>% distinct()

# Build a global assembly from all the mutant models, using the subassembly as a whitelist
# NOTE: this graph should really only have edges that are directly supported by
# genetic perturbations, and may therefore be somewhat sparse.
global_mt_graph = assemble_mt_graph(global_wt_ccm,
                                    global_perturb_models_tbl,# %>% filter(perturb_name %in% c("tbx16-msgn1")),
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="timepoint",
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_mt_graph_edge_whitelist,
                                    q_val=0.1,
                                    verbose=TRUE)
# RGRAPHVIZ CRASHES R
# plot_state_graph_perturb_effects(global_wt_ccm, global_mt_graph,
#                              label_nodes_by="cell_state",
#                              #color_nodes_by = "timepoint",
#                              group_nodes_by="cell_type_sub",
#                              edge_weights = "num_perturbs_supporting",
#                              label_edges_by="support_label",
#                              hide_unlinked_nodes = TRUE)
# ggplot2::ggsave(paste(outdir, "global_mt_graph_perturb_effects.pdf", sep="/"), width=50, height=25, limitsize = FALSE)


plot_state_graph_annotations(global_wt_ccm, global_mt_graph,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_mt_graph.pdf", sep="/"), width=50, height=25, limitsize = FALSE)

plot_state_graph_annotations(global_wt_ccm, global_mt_graph,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_mt_graph_broad.pdf", sep="/"), width=50, height=25,limitsize = FALSE)

mt_edges_not_included_in_global_graph = setdiff (global_mt_graph_edge_whitelist %>% select(from, to),
                                                 global_mt_graph %>% igraph::as_data_frame() %>% select(from, to))
print (paste("Of ", nrow(global_mt_graph_edge_whitelist), "edges from mutant subassemblies",
             nrow(mt_edges_not_included_in_global_graph), "not included in global mutant graph"))

# plot_state_graph_annotations(global_wt_ccm, sub_mt_graph_union,
#                              label_nodes_by="cell_state",
#                              #color_nodes_by = "timepoint",
#                              group_nodes_by="cell_type_sub",
#                              #edge_weights = "num_perturbs_supporting",
#                              #label_edges_by="support_label",
#                              hide_unlinked_nodes = TRUE)
# ggplot2::ggsave(paste(outdir, "sub_mt_graph_union.pdf", sep="/"), width=50, height=25, limitsize = FALSE)
#
# plot_state_graph_annotations(global_wt_ccm, sub_mt_graph_union,
#                              label_nodes_by="cell_state",
#                              #color_nodes_by = "timepoint",
#                              group_nodes_by="cell_type_broad",
#                              edge_weights = "num_perturbs_supporting",
#                              label_edges_by="support_label",
#                              hide_unlinked_nodes = TRUE)
# ggplot2::ggsave(paste(outdir, "sub_mt_graph_union_broad.pdf", sep="/"), width=50, height=25,limitsize = FALSE)


# OK, let's build the final, global graph that we'll use for downstream computations
# We'll do this by annotating the wild-type graph with the graph built from perturbations,
# and adding any perturbation-supported edges that aren't in the WT graph.
# You might think all the mutant edges should exist in the WT graph, bu because we break cycles in both graphs using different heuristics, it's not quite that straightforward

# TODO: we could (and probably should) use the blacklisted wild-type graph instead of
# the WT only graph. But for some reason, that graph is missing some stuff that
# ought to be included somewhere, like the liver. Not clear why that is, and likely
# A bug somewhere in the blacklist procedure. Needs more investigating.
global_wt_graph_edges = igraph::as_data_frame(global_wt_graph)
global_mt_graph_edges = igraph::as_data_frame(global_mt_graph)
global_wt_graph_nodes = igraph::as_data_frame(global_wt_graph, what="vertices")
global_mt_graph_nodes = igraph::as_data_frame(global_mt_graph, what="vertices")
global_wt_graph_nodes = left_join(global_wt_graph_nodes,  global_mt_graph_nodes) %>% rename(id=name)

mt_only = setdiff(global_mt_graph_edges %>% select(from, to), global_wt_graph_edges %>% select(from, to))

global_graph_annotated = left_join(global_wt_graph_edges, global_mt_graph_edges)
global_graph_annotated = global_graph_annotated %>% select(-support)
global_graph_annotated = rbind(global_graph_annotated,
                                  global_mt_graph_edges %>% inner_join(mt_only))
global_graph_annotated = igraph::graph_from_data_frame(global_graph_annotated, vertices=global_wt_graph_nodes)


pdf(paste(outdir, "global_graph_annotated_perturb_effects.pdf", sep="/"), width=50, height=25)
plot_state_graph_perturb_effects(global_wt_ccm, global_graph_annotated,
                                 label_nodes_by="cell_state",
                                 #color_nodes_by = "timepoint",
                                 group_nodes_by="cell_type_sub",
                                 edge_weights = "num_perturbs_supporting",
                                 label_edges_by="support_label",
                                 hide_unlinked_nodes = TRUE)
dev.off()

#ggplot2::ggsave(paste(outdir, "global_graph_annotated_perturb_effects.pdf", sep="/"), width=50, height=25, limitsize = FALSE)

pdf(paste(outdir, "global_graph_annotated.pdf", sep="/"), width=50, height=25)
plot_state_graph_annotations(global_wt_ccm, global_graph_annotated,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
dev.off()
#ggplot2::ggsave(paste(outdir, "global_graph_annotated.pdf", sep="/"), width=50, height=25, limitsize = FALSE)

plot_state_graph_annotations(global_wt_ccm, global_graph_annotated,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_graph_annotated_broad.pdf", sep="/"), width=50, height=25,limitsize = FALSE)

plot_state_graph_annotations(global_wt_ccm, global_graph_annotated,
                             #label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             edge_weights = "num_perturbs_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             con_colour = "black",
                             label_groups=FALSE,
                             min_edge_size=0.5,
                             max_edge_size=3)
ggplot2::ggsave(paste(outdir, "global_graph_annotated_broad_small_version.pdf", sep="/"), width=16, height=9,limitsize = FALSE)

plot_state_graph_annotations(global_wt_ccm, global_graph_annotated,
                             #label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="tissue",
                             edge_weights = "num_perturbs_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             con_colour = "black",
                             label_groups=TRUE,
                             min_edge_size=0.5,
                             max_edge_size=3)
ggplot2::ggsave(paste(outdir, "global_graph_annotated_tissue_small_version.pdf", sep="/"), width=16, height=9,limitsize = FALSE)


#############
# Various additional plots over the global graph:

plot_loss_kinetics_wrapper = function(perturb_ccm,
                                      time_window,
                                      state_graph,
                                      loss_time,
                                      control_ccm=perturb_ccm,
                                      control_start_time=NULL,
                                      control_stop_time=NULL,
                                      interval_step=1,
                                      label_nodes_by="cell_state",
                                      label_groups=FALSE,
                                      group_outline=TRUE,
                                      interval_col="timepoint",
                                      log_abund_detection_thresh=-2,
                                      q_val = 1,
                                      group_nodes_by="cell_type_broad",
                                      expt="GAP16",
                                      legend_position="none",
                                      hide_unlinked_nodes = TRUE){

  perturb_start_time = min(time_window$start_time)
  perturb_stop_time = min(time_window$stop_time)

  if (is.null(control_start_time))
    control_start_time = perturb_start_time
  if (is.null(control_stop_time))
    control_stop_time = perturb_stop_time

  if (perturb_start_time == perturb_stop_time)
    return (NULL)

  #print (time_window)
  # "largest_loss", "largest_loss_time", "earliest_time", "latest_time", "peak_loss_time"
  loss_graph = hooke:::plot_state_graph_losses(perturb_ccm,
                                               state_graph,
                                               start_time=perturb_start_time,
                                               stop_time=perturb_stop_time,
                                               loss_time=loss_time,
                                               interval_step=interval_step,
                                               label_nodes_by=label_nodes_by,
                                               label_groups=label_groups,
                                               label_edges_by="support_label",
                                               group_outline=group_outline,
                                               interval_col=interval_col,
                                               log_abund_detection_thresh=log_abund_detection_thresh,
                                               q_val = q_val,
                                               group_nodes_by=group_nodes_by,
                                               expt=expt,
                                               legend_position=legend_position,
                                               hide_unlinked_nodes = hide_unlinked_nodes,
                                               control_ccm=perturb_ccm,
                                               control_start_time=control_start_time,
                                               control_stop_time=control_stop_time
  )
  #loss_graph + ggplot2::ggtitle(perturb_name)
  return (loss_graph)
}

add_title_to_loss_plot =  function(state_graph_plot, perturb_name) {
  state_graph_plot = state_graph_plot + ggplot2::ggtitle(perturb_name)
  return(state_graph_plot)
}

render_plot =  function(state_graph_plot, perturb_name, plot_prefix, outdir, width=40, height=15) {
  plot_name = paste(outdir, paste(plot_prefix, "_", perturb_name, ".pdf", sep=""), sep="/")
  pdf(plot_name, width=width, height=height)
  print(state_graph_plot)
  dev.off()
  return(plot_name)
}

perturbation_model_plots = global_perturb_models_tbl %>%
  filter(is.na(perturb_ccm) == FALSE) %>%
  mutate(largest_loss_graph = purrr::map2(.f = plot_loss_kinetics_wrapper,
                                          .x=perturb_ccm,
                                          .y=perturb_time_window,
                                          global_graph_annotated,
                                          control_ccm=global_wt_ccm,
                                          control_start_time=assembly_start_time,
                                          control_stop_time=assembly_stop_time,
                                          loss_time="largest_loss_time",
                                          q_val=0.1)) %>%
  tidyr::unnest(perturb_time_window) %>%
  filter(start_time < stop_time) %>%
  mutate(largest_loss_graph = purrr::map2(.f = add_title_to_loss_plot,
                                          .x=largest_loss_graph,
                                          .y=perturb_name),
         largest_loss_graph_plotfile = purrr::map2(.f = purrr:::possibly(render_plot, NA_real_),
                                                   .x=largest_loss_graph,
                                                   .y=perturb_name,
                                                   plot_prefix="global_graph_largest_loss_time",
                                                   outdir=outdir))



perturbation_model_plots = global_perturb_models_tbl %>%
  filter(is.na(perturb_ccm) == FALSE) %>%
  mutate(largest_loss_graph = purrr::map2(.f = plot_loss_kinetics_wrapper,
                                          .x=perturb_ccm,
                                          .y=perturb_time_window,
                                          global_graph_annotated,
                                          control_ccm=global_wt_ccm,
                                          control_start_time=assembly_start_time,
                                          control_stop_time=assembly_stop_time,
                                          loss_time="largest_loss",
                                          q_val=0.1)) %>%
  tidyr::unnest(perturb_time_window) %>%
  filter(start_time < stop_time) %>%
  mutate(largest_loss_graph = purrr::map2(.f = add_title_to_loss_plot,
                                          .x=largest_loss_graph,
                                          .y=perturb_name),
         largest_loss_graph_plotfile = purrr::map2(.f =  purrr::possibly(render_plot, NA_real_),
                                                   .x=largest_loss_graph,
                                                   .y=perturb_name,
                                                   plot_prefix="global_graph_largest_loss",
                                                   outdir=outdir))


perturbation_model_plots = global_perturb_models_tbl %>%
  filter(is.na(perturb_ccm) == FALSE) %>%
  mutate(largest_loss_graph = purrr::map2(.f = plot_loss_kinetics_wrapper,
                                         .x=perturb_ccm,
                                         .y=perturb_time_window,
                                         global_graph_annotated,
                                         control_ccm=global_wt_ccm,
                                         control_start_time=assembly_start_time,
                                         control_stop_time=assembly_stop_time,
                                         loss_time="peak_loss_time",
                                         q_val=0.1)) %>%
  tidyr::unnest(perturb_time_window) %>%
  filter(start_time < stop_time) %>%
  mutate(largest_loss_graph = purrr::map2(.f = add_title_to_loss_plot,
                                          .x=largest_loss_graph,
                                          .y=perturb_name),
         largest_loss_graph_plotfile = purrr::map2(.f = purrr::possibly(render_plot, NA_real_),
                                                   .x=largest_loss_graph,
                                                   .y=perturb_name,
                                                   plot_prefix="global_graph_peak_loss_time",
                                                   outdir=outdir))


perturbation_model_plots = global_perturb_models_tbl %>%
  filter(is.na(perturb_ccm) == FALSE) %>%
  mutate(largest_loss_graph = purrr::map2(.f = plot_loss_kinetics_wrapper,
                                          .x=perturb_ccm,
                                          .y=perturb_time_window,#tibble(start_time=assembly_start_time, stop_time=assembly_stop_time),
                                          global_graph_annotated,
                                          control_ccm=global_wt_ccm,
                                          control_start_time=assembly_start_time,
                                          control_stop_time=assembly_stop_time,
                                          loss_time="delta_log_abund_at_peak",
                                          q_val=1)) %>%
  tidyr::unnest(perturb_time_window) %>%
  filter(start_time < stop_time) %>%
  mutate(largest_loss_graph = purrr::map2(.f = add_title_to_loss_plot,
                                          .x=largest_loss_graph,
                                          .y=perturb_name),
         largest_loss_graph_plotfile = purrr::map2(.f = purrr::possibly(render_plot, NA_real_),
                                                   .x=largest_loss_graph,
                                                   .y=perturb_name,
                                                   plot_prefix="global_graph_peak_loss_no_sig_filter",
                                                   outdir=outdir))

perturbation_model_plots = global_perturb_models_tbl %>%
  filter(is.na(perturb_ccm) == FALSE) %>%
  mutate(largest_loss_graph = purrr::map2(.f = plot_loss_kinetics_wrapper,
                                         .x=perturb_ccm,
                                         .y=perturb_time_window,#tibble(start_time=assembly_start_time, stop_time=assembly_stop_time),
                                         global_graph_annotated,
                                         control_ccm=global_wt_ccm,
                                         control_start_time=assembly_start_time,
                                         control_stop_time=assembly_stop_time,
                                         loss_time="delta_log_abund_at_peak",
                                         q_val=0.1)) %>%
  tidyr::unnest(perturb_time_window) %>%
  filter(start_time < stop_time) %>%
  mutate(largest_loss_graph = purrr::map2(.f = add_title_to_loss_plot,
                                          .x=largest_loss_graph,
                                          .y=perturb_name),
         largest_loss_graph_plotfile = purrr::map2(.f = purrr::possibly(render_plot, NA_real_),
                                                   .x=largest_loss_graph,
                                                   .y=perturb_name,
                                                   plot_prefix="global_graph_peak_loss",
                                                   outdir=outdir))



perturbation_model_plots = selected_perturb_models_tbl %>%
  filter(is.na(perturb_ccm) == FALSE) %>%
  mutate(earliest_loss_graph = purrr::map2(.f = plot_loss_kinetics_wrapper,
                                           .x=perturb_ccm,
                                           .y=perturb_time_window,
                                           global_wt_graph,
                                           loss_time="earliest_time",
                                           q_val=0.1)) %>%
  tidyr::unnest(perturb_time_window) %>%
  filter(start_time < stop_time) %>%
  mutate(earliest_loss_graph = purrr::map2(.f = add_title_to_loss_plot,
                                           .x=earliest_loss_graph,
                                           .y=perturb_name),
         largest_loss_graph_plotfile = purrr::map2(.f = render_plot,
                                                   .x=earliest_loss_graph,
                                                   .y=perturb_name,
                                                   plot_prefix="wt_graph_earliest_loss",
                                                   outdir=outdir))





#########################################
# END OF IMPORTANT STUFF.
# EVERYTHING BELOW IS CODE FOR DEBUGGING,
# TESTING OUT ALTERNATIVE STRATEGIES, ETC.
##########################################


#################################
# Some other strategies for graph construction:

### Another possible strategy: build a wild-type graph from wild-type embryos,
# But only use the controls and constrain the assembly to no include edges that
# "violate" genetic perturbations.
sub_wt_blacklist_graph_union = do.call(igraph::union, comb_res %>% filter(is.na(wt_graph_blacklist_converted) == FALSE) %>% pull(mt_graph_converted))
global_wt_graph_blacklist_edge_whitelist = igraph::as_data_frame(sub_wt_blacklist_graph_union)
global_wt_graph_blacklist_edge_whitelist = global_wt_graph_blacklist_edge_whitelist %>% select(from, to) %>% distinct()

global_wt_graph_blacklist = assemble_wt_graph(comb_sample_cds,
                                    global_wt_ccm,
                                    sample_group = "embryo",
                                    cell_group = "cell_state",
                                    main_model_formula_str = NULL,
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="timepoint",
                                    #nuisance_model_formula_str = "~expt",
                                    ctrl_ids = control_genotypes,
                                    sparsity_factor = 0.01,
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_wt_graph_blacklist_edge_whitelist,
                                    verbose=TRUE)

# Plot the giant graph
plot_state_graph_annotations(global_wt_ccm, global_wt_graph_blacklist, color_nodes_by = "timepoint", group_nodes_by="cell_type_sub",
                             edge_weights = "support",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph_blacklist.png", sep="/"), width=30, height=10)

# Plot it by cell type broad just because the labels can get crazy
plot_state_graph_annotations(global_wt_ccm, global_wt_graph_blacklist, color_nodes_by = "timepoint", group_nodes_by="cell_type_broad",
                             edge_weights = "support",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph_blacklist_broad.png", sep="/"), width=30, height=10)


wt_edges_not_included_in_global_graph = setdiff (global_wt_graph_blacklist_edge_whitelist %>% select(from, to),
                                                 global_wt_graph_blacklist %>% igraph::as_data_frame() %>% select(from, to))
print (paste("Of ", nrow(global_wt_graph_blacklist_edge_whitelist), "edges from WT subassemblies",
             nrow(wt_edges_not_included_in_global_graph), "not included in global WT graph"))

plot_state_graph_annotations(global_wt_ccm, sub_wt_blacklist_graph_union,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             #edge_weights = "num_perturbs_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "sub_wt_blacklist_graph_union.pdf", sep="/"), width=50, height=25, limitsize = FALSE)

plot_state_graph_annotations(global_wt_ccm, sub_wt_blacklist_graph_union,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "sub_wt_blacklist_graph_union.pdf", sep="/"), width=50, height=25,limitsize = FALSE)


#saveRDS(global_perturb_models_tbl, paste(outdir, "global_perturb_models_tbl.rds", sep="/"))
saveRDS(global_mt_graph, paste(outdir, "global_mt_graph.rds", sep="/"))

#####################################
global_perturb_models_tbl = global_perturb_models_tbl %>%
  filter(is.na(perturb_ccm) == FALSE) %>%
  mutate(discordant_state_pairs = purrr::map2(.f = get_discordant_loss_pairs,
                                              .x = perturb_ccm,
                                              .y = perturb_time_window,
                                              control_timeseries_ccm = global_wt_ccm,
                                              control_time_window=tibble(start_time=assembly_start_time, stop_time=assembly_stop_time),
                                              interval_step=2,
                                              interval_col="timepoint",
                                              log_abund_detection_thresh=-5,
                                              q_val=0.1,
                                              expt="GAP16"))

discordant_pair_summary = global_perturb_models_tbl %>%
  tidyr::unnest(discordant_state_pairs) %>%
  select(perturb_name, from=lost_cell_groups, to=unaffected_cell_groups) %>%
  group_by(from, to) %>%
  summarize(num_perturbs_discordant = n())

mt_graph_edge_whitelist = igraph::as_data_frame(global_mt_graph)
mt_graph_edge_whitelist = mt_graph_edge_whitelist %>% select(from, to, num_perturbs_supporting) %>% distinct()

mt_graph_edge_blacklist = discordant_pair_summary %>% left_join(mt_graph_edge_whitelist)
mt_graph_edge_blacklist = mt_graph_edge_blacklist %>% mutate(num_perturbs_supporting = replace_na(num_perturbs_supporting, 0))

mt_graph_edge_blacklist = mt_graph_edge_blacklist %>%
  filter(num_perturbs_supporting == 0 & num_perturbs_supporting < num_perturbs_discordant)

print (as.data.frame(mt_graph_edge_blacklist))

mt_graph_wl_minus_bl = dplyr::setdiff(mt_graph_edge_whitelist %>% select(from, to),
                                      mt_graph_edge_blacklist %>% select(from, to))

wt_graph_edge_whitelist = igraph::as_data_frame(global_wt_graph)
wt_graph_edge_whitelist = wt_graph_edge_whitelist %>% select(from, to) %>% distinct()

wt_graph_wl_minus_bl = dplyr::setdiff(wt_graph_edge_whitelist %>% select(from, to),
                                      mt_graph_edge_blacklist %>% select(from, to))

global_wt_graph_blacklist = assemble_wt_graph(comb_sample_cds,
                                              global_wt_ccm,
                                              start_time = assembly_start_time,
                                              stop_time = assembly_stop_time,
                                              interval_col="timepoint",
                                              perturbation_col = "gene_target",
                                              edge_whitelist = wt_graph_wl_minus_bl,
                                              edge_blacklist = mt_graph_edge_blacklist,
                                              verbose=TRUE,
                                              expt="GAP16")


plot_state_graph_annotations(global_wt_ccm, global_wt_graph_blacklist, color_nodes_by = "timepoint", group_nodes_by="cell_type_sub",
                             edge_weights = "support",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph_blacklist.pdf", sep="/"), width=50, height=25, limitsize = FALSE)

plot_state_graph_annotations(global_wt_ccm, global_wt_graph_blacklist,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             #edge_weights = "num_perturbs_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph_blacklist_broad.pdf", sep="/"), width=50, height=25, limitsize = FALSE)



global_mt_graph_blacklist = assemble_mt_graph(global_wt_ccm,
                                    global_perturb_models_tbl,
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="timepoint",
                                    perturbation_col = "gene_target",
                                    edge_whitelist = mt_graph_wl_minus_bl,
                                    edge_blacklist = mt_graph_edge_blacklist,
                                    verbose=TRUE,
                                    expt="GAP16")

plot_state_graph_annotations(global_wt_ccm, global_mt_graph_blacklist,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_mt_graph_blacklist.pdf", sep="/"), width=50, height=25, limitsize = FALSE)

plot_state_graph_annotations(global_wt_ccm, global_mt_graph_blacklist,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_mt_graph_blacklist_broad.pdf", sep="/"), width=50, height=25, limitsize = FALSE)


mt_graph_not_in_bl = global_mt_graph - global_mt_graph_blacklist

plot_state_graph_annotations(global_wt_ccm, mt_graph_not_in_bl,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "mt_graph_not_in_bl.pdf", sep="/"), width=40, height=15)




# Learn a single graph on all states at once (using the subgraphs as a whitelist/prior)
global_wt_graph_blacklist = assemble_wt_graph(comb_sample_cds,
                                    global_wt_ccm,
                                    sample_group = "embryo",
                                    cell_group = "cell_state",
                                    main_model_formula_str = NULL,
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="timepoint",
                                    #nuisance_model_formula_str = "~expt",
                                    ctrl_ids = control_genotypes,
                                    sparsity_factor = 0.01,
                                    perturbation_col = "gene_target",
                                    edge_whitelist = global_mt_graph_wl_minus_bl,
                                    edge_blacklist = global_mt_graph_edge_blacklist,
                                    verbose=TRUE)


# Plot the giant graph
plot_state_graph_annotations(global_wt_ccm, global_wt_graph_blacklist, color_nodes_by = "timepoint", group_nodes_by="cell_type_sub",
                             edge_weights = "support",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph_blacklist.png", sep="/"), width=30, height=10)

# Plot it by cell type broad just because the labels can get crazy
plot_state_graph_annotations(global_wt_ccm, global_wt_graph_blacklist, color_nodes_by = "timepoint", group_nodes_by="cell_type_broad",
                             edge_weights = "support",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_wt_graph_broad_blacklist.png", sep="/"), width=30, height=10)


###### Plots for loss heatmaps:

loss_timing_wrapper = function(perturbation_ccm, perturb_time_window){
  return(hooke:::estimate_loss_timing(perturbation_ccm,
         start_time=min(perturb_time_window$start_time),
         stop_time=min(perturb_time_window$stop_time),
         interval_step = 1,
         interval_col="timepoint",
         log_abund_detection_thresh=-5,
         q_val = 0.1,
         delta_log_abund_loss_thresh=0,
         expt="GAP16"))
}
#debug(loss_timing_wrapper)
# FIXME: shouldn't use hardcoded knockout here
loss_table = global_perturb_models_tbl %>%
  filter (is.na(perturb_ccm) == FALSE) %>%
  mutate(perturbation_losses = purrr::map2(.f = loss_timing_wrapper,
                                          .x=perturb_ccm,
                                          .y=perturb_time_window))

loss_table = loss_table %>% tidyr::unnest(perturbation_losses)
loss_matrix = loss_table %>% select(perturb_name, cell_group, delta_log_abund_at_peak) %>%
  pivot_wider(names_from = perturb_name, values_from = delta_log_abund_at_peak, values_fill=0)
loss_matrix_names = loss_matrix %>% pull(cell_group)
loss_matrix = loss_matrix %>% select(-cell_group)
loss_matrix = as.matrix(loss_matrix)
row.names(loss_matrix) = loss_matrix_names
loss_matrix[loss_matrix < -2] = -2
loss_matrix[loss_matrix > 2] = 2

pheatmap::pheatmap(t(loss_matrix),
                   clustering_method="ward.D2",
                   color=colorRampPalette(c("steelblue", "white", "orangered"))(50),
                   fontsize=6,
                   show_rownames=TRUE,
                   show_colnames=FALSE,
                   treeheight_col = 0,
                   treeheight_row = 0,
                   filename=paste(outdir, "loss_matrix.png", sep="/"),
                   width=3, height=2)

######
# Plots for global assemblies:

selected_perturb = "mafba"
selected_perturb_results = global_perturb_models_tbl %>% filter(selected_perturb == perturb_name)
selected_ccm = selected_perturb_results$perturb_ccm[[1]]

hooke:::plot_state_graph_losses(selected_ccm,
                                global_mt_graph,
                                min(selected_perturb_results$perturb_time_window[[1]]$start_time),
                                min(selected_perturb_results$perturb_time_window[[1]]$stop_time),
                                interval_step=1,
                                label_nodes_by="cell_state",
                                label_groups=FALSE,
                                group_outline=TRUE,
                                interval_col="timepoint",
                                log_abund_detection_thresh=-2,
                                q_val = 1,
                                loss_time="largest_loss_time",
                                group_nodes_by="cell_type_broad",
                                experiment="GAP16",
                                legend_position="right",
                                hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "mt_graph_largest_loss_times_linked.pdf", sep="/"), width=40, height=15)



hooke:::plot_state_graph_losses(selected_ccm,
                                global_mt_graph,
                                min(selected_perturb_results$perturb_time_window[[1]]$start_time),
                                min(selected_perturb_results$perturb_time_window[[1]]$stop_time),
                                interval_step=1,
                                label_nodes_by="cell_state",
                                label_groups=FALSE,
                                group_outline=TRUE,
                                interval_col="adjusted_timepoint",
                                log_abund_detection_thresh=-2,
                                q_val = 1,
                                loss_time="largest_loss_time",
                                group_nodes_by="cell_type_broad",
                                experiment="GAP16",
                                legend_position="right",
                                hide_unlinked_nodes = FALSE)
ggplot2::ggsave(paste(outdir, "mt_graph_largest_loss_times_all.pdf", sep="/"), width=40, height=15)


hooke:::plot_state_graph_losses(selected_ccm,
                                global_mt_graph,
                                min(selected_perturb_results$perturb_time_window[[1]]$start_time),
                                min(selected_perturb_results$perturb_time_window[[1]]$stop_time),
                                interval_step=1,
                                label_nodes_by="cell_state",
                                label_groups=FALSE,
                                group_outline=TRUE,
                                interval_col="adjusted_timepoint",
                                log_abund_detection_thresh=-2,
                                q_val = 1,
                                loss_time="earliest_time",
                                group_nodes_by="cell_type_broad",
                                experiment="GAP16",
                                legend_position="right",
                                hide_unlinked_nodes = FALSE)
ggplot2::ggsave(paste(outdir, "mt_graph_earliest_loss_times_all.pdf", sep="/"), width=40, height=15)


# hooke:::plot_state_graph_kinetics(global_wt_ccm,
#                                 global_mt_graph,
#                                 assembly_start_time,
#                                 assembly_stop_time,
#                                 interval_step=1,
#                                 label_nodes_by="cell_state",
#                                 label_groups=FALSE,
#                                 group_outline=TRUE,
#                                 interval_col="adjusted_timepoint",
#                                 log_abund_detection_thresh=-2,
#                                 q_val = 1,
#                                 time_to_plot="emergence_time",
#                                 group_nodes_by="cell_type_broad",
#                                 experiment="GAP16",
#                                 legend_position="right",
#                                 hide_unlinked_nodes = FALSE)
# ggplot2::ggsave(paste(outdir, "mt_graph_emergence_time_times_all.pdf", sep="/"), width=40, height=15)
#
#
# hooke:::plot_state_graph_kinetics(global_wt_ccm,
#                                   global_mt_graph,
#                                   assembly_start_time,
#                                   assembly_stop_time,
#                                   interval_step=1,
#                                   label_nodes_by="cell_state",
#                                   label_groups=FALSE,
#                                   group_outline=TRUE,
#                                   interval_col="adjusted_timepoint",
#                                   log_abund_detection_thresh=-2,
#                                   q_val = 1,
#                                   time_to_plot="extinction_time",
#                                   group_nodes_by="cell_type_broad",
#                                   experiment="GAP16",
#                                   legend_position="right",
#                                   hide_unlinked_nodes = FALSE)
# ggplot2::ggsave(paste(outdir, "mt_graph_extinction_times_all.pdf", sep="/"), width=40, height=15)
#
# hooke:::plot_state_graph_kinetics(global_wt_ccm,
#                                   global_mt_graph,
#                                   assembly_start_time,
#                                   assembly_stop_time,
#                                   interval_step=1,
#                                   label_nodes_by="cell_state",
#                                   label_groups=FALSE,
#                                   group_outline=TRUE,
#                                   interval_col="adjusted_timepoint",
#                                   log_abund_detection_thresh=-2,
#                                   q_val = 1,
#                                   time_to_plot="peak_abund_time",
#                                   group_nodes_by="cell_type_broad",
#                                   experiment="GAP16",
#                                   legend_position="right",
#                                   hide_unlinked_nodes = FALSE)
# ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)


### Plotting cell type kinetics:


#debug(plot_cell_type_perturb_kinetics)
source("R/assembly_utils.R")
#debug(plot_cell_type_control_kinetics)

pdf(paste(outdir, "selected_ko_cell_type_distribution_over_time.pdf", sep="/"), height=49, width=49)
print(kinetics_plot)
#ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)
dev.off()


kinetics_plot = plot_cell_type_perturb_kinetics(selected_ccm, interval_col="adjusted_timepoint") +
  monocle3:::monocle_theme_opts()


pdf(paste(outdir, "selected_ko_cell_type_distribution_over_time.pdf", sep="/"), height=49, width=49)
print(kinetics_plot)
#ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)
dev.off()

selected_broad_groups = c("neurons (gabaergic, glutamatergic)",
                          "neural progenitor (hindbrain R7/8)",
                          "neural progenitor (telencephalon/diencephalon)",
                          "neural progenitor (MHB)",
                          "differentiating neuron 2",
                          "differentiating neuron (hindbrain)",
                          "differentiating neuron 1"
                          )

#selected_broad_groups = c("intestine")

selected_broad_groups = c("cranial neural crest",
                          "neuron (cranial ganglia sensory, Rohon-Beard)")

selected_cell_groups = colData(comb_sample_cds) %>%
  as.data.frame %>%
  group_by(cell_type_broad, cell_state) %>%
  summarize(total_cells = n()) %>% ungroup() %>%
  group_by(cell_state) %>%
  slice_max(total_cells, n=1) %>%
  filter(cell_type_broad %in% selected_broad_groups) %>%
  pull(cell_state) %>% unique

#debug(plot_cell_type_perturb_kinetics)
source("R/assembly_utils.R")
#debug(plot_cell_type_perturb_kinetics)
kinetics_plot = plot_cell_type_perturb_kinetics(selected_ccm,
                                                cell_groups = selected_cell_groups,
                                                interval_col="timepoint",
                                                start_time=assembly_start_time,
                                                stop_time=assembly_stop_time,
                                                expt="GAP16",
                                                q_val = 0.1) +
  monocle3:::monocle_theme_opts() #+ scale_y_continuous()


pdf(paste(outdir, "csg_selected_ko_cell_type_distribution_over_time.pdf", sep="/"), height=30, width=6)
print(kinetics_plot)
#ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)
dev.off()


selected_broad_groups = c("mature fast muscle")
selected_cell_groups = colData(comb_sample_cds) %>%
  as.data.frame %>%
  group_by(cell_type_broad, cell_state) %>%
  summarize(total_cells = n()) %>% ungroup() %>%
  group_by(cell_state) %>%
  slice_max(total_cells, n=1) %>%
  filter(cell_type_broad %in% selected_broad_groups) %>%
  pull(cell_state) %>% unique

source("R/assembly_utils.R")
#debug(plot_cell_type_control_kinetics)
kinetics_plot = plot_cell_type_control_kinetics(global_wt_ccm,
                                                cell_groups = selected_cell_groups,
                                                interval_col="timepoint") +
  monocle3:::monocle_theme_opts()
kinetics_plot


#debug(plot_cell_type_perturb_kinetics)
source("R/assembly_utils.R")
#debug(plot_cell_type_perturb_kinetics)
kinetics_plot = plot_cell_type_perturb_kinetics(selected_ccm,
                                                cell_groups = selected_cell_groups,
                                                interval_col="adjusted_timepoint") +
  monocle3:::monocle_theme_opts()


#pdf(paste(outdir, "muscle_selected_ko_cell_type_distribution_over_time.pdf", sep="/"), height=10, width=6)
print(kinetics_plot)
#ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)
#dev.off()





#
# sel_ccs_counts = normalized_counts(selected_ccm@ccs, method="size_only", pseudocount=1e-5)
# sel_ccs_counts_long = summary(sel_ccs_counts)
# sel_ccs_counts_long = data.frame(cell_group = rownames(sel_ccs_counts)[sel_ccs_counts_long$i],
#            embryo = colnames(sel_ccs_counts)[sel_ccs_counts_long$j],
#            num_cells      = sel_ccs_counts_long$x)
#
# cell_group_metadata = hooke:::collect_psg_node_metadata(selected_ccm,
#                                                         group_nodes_by="cell_type_broad",
#                                                         color_nodes_by="cell_type_broad",
#                                                         label_nodes_by="cell_type_broad") %>%
#   select(id, cell_type_broad = group_nodes_by)
#
# sel_ccs_counts_long = left_join(sel_ccs_counts_long,
#                                 cell_group_metadata,
#                                 by=c("cell_group"="id"))
#
# sel_ccs_counts_long = left_join(sel_ccs_counts_long,
#                                colData(selected_ccm@ccs) %>% as.data.frame %>% select(sample, adjusted_timepoint, knockout),
#                                 by=c("embryo"="sample"))
#
# #sel_ccs_counts_long = sel_ccs_counts_long %>% filter(expt == "GAP16")
#
# #sel_ccs_counts_long = sel_ccs_counts_long %>% filter(grepl("uscle", cell_type_broad))
# sel_ccs_counts_long = sel_ccs_counts_long %>% filter(is.na(cell_group) == FALSE & is.na(adjusted_timepoint) == FALSE)
#
# # cell_densities <- sel_ccs_counts_long %>%
# #   group_by(cell_group, knockout) %>%
# #   group_modify(~ ggplot2:::compute_density(.x$adjusted_timepoint, .x$num_cells)) %>%
# #   rename(adjusted_timepoint = x)
# # cell_densities = cell_densities %>% filter(adjusted_timepoint >= min(sel_ccs_counts_long$adjusted_timepoint) &
# #                                            adjusted_timepoint <= max(sel_ccs_counts_long$adjusted_timepoint))
# # #cell_densities = cell_densities %>% filter(adjusted_timepoint >= 18 &
# # #                                           adjusted_timepoint <= 72)
# # cell_densities = cell_densities %>% filter(is.na(cell_group) == FALSE)
#
#
# timing_table = hooke:::estimate_loss_timing(selected_ccm,
#                              start_time=min(sel_ccs_counts_long$adjusted_timepoint),
#                              stop_time=max(sel_ccs_counts_long$adjusted_timepoint),
#                              interval_step=1,
#                              interval_col="adjusted_timepoint")
#
# timing_table = timing_table %>%
#   left_join(sel_ccs_counts_long %>% select(cell_group, cell_type_broad) %>% distinct()) %>%
#   arrange(cell_type_broad, peak_wt_time) %>%
#   mutate(cell_group = forcats::fct_reorder(cell_group, peak_wt_time))
#
# plot_counts = ggplot(sel_ccs_counts_long, aes(x = adjusted_timepoint, y = num_cells, height = density, fill=knockout)) +
#   facet_wrap(cell_group, ncol=1) + monocle3:::monocle_theme_opts()
#
# cell_densities$cell_group = factor(cell_densities$cell_group, levels=levels(timing_table$cell_group))
#
# ridge_plot = ggplot(cell_densities, aes(x = adjusted_timepoint, y = cell_group, height = density, fill=knockout)) +
#   ggridges::geom_density_ridges(stat = "identity") + monocle3:::monocle_theme_opts()
#
# pdf(paste(outdir, "selected_ko_cell_type_distribution_over_time.pdf", sep="/"), height=49, width=4)
# print(ridge_plot)
# #ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)
# dev.off()
#
#
# ggplot2::ggplot(aes(x=adjusted_timepoint, y=num_cells, fill=knockout), data=sel_ccs_counts_long) + geom_boxplot() + facet_wrap(~cell_group)
#

# Debugging:

periderm_cells = comb_res %>% filter(partition == "1") %>% select(data) %>% tidyr::unnest(data)
periderm_cds = comb_sample_cds[,periderm_cells %>% pull(cds_row_id)]
periderm_states = colData(periderm_cds)$cell_state %>% unique

periderm_wt_whitelist = global_wt_graph_edge_whitelist %>% filter(from %in% periderm_states | to %in% periderm_states)
periderm_mt_whitelist = global_mt_graph_edge_whitelist %>% filter(from %in% periderm_states | to %in% periderm_states)

# Fit a single cell count model to all the cell states at once
periderm_wt_ccm = fit_wt_model(periderm_cds,
                             sample_group = "embryo",
                             cell_group = "cell_state",
                             main_model_formula_str = NULL,
                             #num_breaks=4,
                             #main_model_formula_str = "~1",
                             #main_model_formula_str = wt_main_model_formula_str,
                             start_time = assembly_start_time,
                             stop_time = assembly_stop_time,
                             interval_col="timepoint",
                             #nuisance_model_formula_str = "~expt",
                             ctrl_ids = control_genotypes,
                             #ctrl_ids = colData(periderm_cds)$gene_target %>% unique,
                             sparsity_factor = 0.01,
                             perturbation_col = "gene_target",
                             edge_whitelist = periderm_wt_whitelist,
                             #base_penalty=1000,
                             keep_cds = FALSE,
                             verbose=TRUE,
                             num_threads=num_threads,
                             backend=assembly_backend)


noto_periderm_cds = periderm_cds[,colData(periderm_cds)$gene_target %in% c("noto", control_genotypes)]
peridem_perturb_models_tbl = fit_mt_models(noto_periderm_cds,
                                          sample_group = "embryo",
                                          cell_group = "cell_state",
                                          main_model_formula_str = NULL,
                                          start_time = assembly_start_time,
                                          stop_time = assembly_stop_time,
                                          interval_col="timepoint",
                                          #nuisance_model_formula_str = "~expt",
                                          ctrl_ids = control_genotypes,
                                          mt_ids = selected_mt_ids,
                                          sparsity_factor = 0.01,
                                          perturbation_col = "gene_target",
                                          keep_cds=FALSE,
                                          num_threads=num_threads,
                                          backend=assembly_backend)

periderm_mt_graph = assemble_mt_graph(periderm_wt_ccm,
                                      peridem_perturb_models_tbl,
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="adjusted_timepoint",
                                    perturbation_col = "gene_target",
                                    edge_whitelist = periderm_mt_whitelist,
                                    verbose=TRUE)

