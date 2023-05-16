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

plan(multicore)
options(future.fork.multithreading.enable=FALSE)
options(future.globals.maxSize = 8192 * 1024^2) # 8GB

assembly_start_time = 18
assembly_stop_time = 72
assembly_backend = "nlopt"

devtools::load_all("~/OneDrive/UW/Trapnell/hooke_pln_20230104/")

# -----------------------------------------------------------------------------

# FIXME: set this to number of CPUs on machine (or in cluster session)
#RhpcBLASctl::omp_set_num_threads(10)
RhpcBLASctl::omp_set_num_threads(8)

num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
print (paste("running assembly with", num_threads, "threads"))

Sys.setenv("OMP_NUM_THREADS" = 1)
Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

outdir = "/Users/maddyduran/UW/Trapnell/zebrafish-atlas-assembly/"
setwd("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/")
source("R/assembly_utils.R")

# starting with the full cds --------------------------------------------------

wt_cds = readRDS("~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/filtered_full_cds5e-06.rds")
# mt_cds = readRDS("~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/gap16_proj_cds_doublet_filtered.rds")
# 
# # drop periderm and then combine
# 
# 
# wt_ccs = new_cell_count_set(wt_cds, 
#                             sample_group = "embryo", 
#                             cell_group = "cell_type_sub", keep_cds = F)
# mt_ccs = new_cell_count_set(mt_cds, 
#                             sample_group = "embryo", 
#                             cell_group = "cell_type_sub", keep_cds = F)
# 
# combine_ccs(ccs_list = list(wt_ccs, mt_ccs))


# ------------------------------------------------------------------------------

# the cell type sub 
# time is a factor 
# cluster resolution so that clusters have different times 
# global model have spline of time in nuisance model formula 
# subtracted out of the pcor 
# lumping things together in time  

# like slow and fast muscle 







# mt_cds = detect_outlier_cells(mt_cds)
# 
# 
# min_dist_quantiles = quantile(colData(mt_cds)$min_nn_dist, probs = seq(0, 1, 0.001), na.rm=TRUE)
# min_dist_thresh = min_dist_quantiles["99.9%"]
# colData(mt_cds)$outlier = colData(mt_cds)$min_nn_dist > min_dist_thresh
# mt_cds = mt_cds[,colData(mt_cds)$outlier == FALSE]
# 
# plot_cells_3d(mt_cds[, sample(ncol(mt_cds), 100000)], color_cells_by = "cell_type_sub")
# 
# saveRDS(mt_cds, 
#         "~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/gap16_proj_cds_doublet_filtered.rds")


# transfer partition labels 

# mt_cds = load_transform_models(mt_cds,"../hooke_manuscript/supplement/filter_doublets/filtered_full_transform_models/")
# 
# mt_cds = transfer_cell_labels(mt_cds, 
#                               reduction_method = "UMAP", 
#                               as.data.frame(colData(wt_cds)), 
#                               ref_column_name = "partition")
# 
# colData(wt_cds)$new_partition = colData(wt_cds)$partition
# 
# colData(mt_cds)$new_partition = colData(mt_cds)$partition
# 
# plot_cells_3d(wt_cds[, sample(ncol(wt_cds), 100000)], color_cells_by = "new_partition")
# 
# plot_cells_3d(mt_cds[, sample(ncol(mt_cds), 100000)], color_cells_by = "new_partition")
# 
# saveRDS(mt_cds,
#         "~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/gap16_proj_cds_doublet_filtered.rds")






# -----------------------------------------------------------------------------


# cds = cds[,colData(cds)$outlier == FALSE]


# for initial debugging let's subsample these
# wt_cds = wt_cds[, sample(ncol(wt_cds), 50000)]
# mt_cds = mt_cds[, sample(ncol(mt_cds), 50000)]
comb_cds = combine_cds(list(wt_cds, mt_cds), keep_reduced_dims = T)

# can i reconstruct my sub cds from my big cds --------------------------------
# aka pull partition umap coords from coldata
part_2_cds = get_partition_cds(comb_cds, 2)
plot_cells(part_2_cds)

part_3_cds = get_partition_cds(comb_cds, 3)
plot_cells(part_3_cds)

part_4_cds = get_partition_cds(comb_cds, 4)
plot_cells(part_4_cds)

part_7_cds = get_partition_cds(wt_cds, 7)


# test recursive clustering ----------------------------------------------------

# debug(subcluster_cds)
test_cds = subcluster_cds(part_3_cds, partition_name = 3, recursive_subcluster = TRUE)
# should be able to separate
unique(colData(test_cds)$partition)
# "3_1_1" "3_2_1"
unique(colData(test_cds)$cluster)
# 1 2
unique(colData(test_cds)$cell_state)
# "3_1-1" "3_1-2" "3_2-1

test_cds2 = subcluster_cds(part_3_cds, partition_name = 3, recursive_subcluster = FALSE)

unique(colData(test_cds2)$partition)
# "3_1" "3_2"
unique(colData(test_cds2)$cluster)
# 1 2 3
unique(colData(test_cds2)$cell_state)
# "3-1" "3-2" "3-3

# just recluster once
test_cds3 = subcluster_cds(part_3_cds, partition_name = NULL, recursive_subcluster = FALSE)
unique(colData(test_cds3)$partition)
# 1 2 
unique(colData(test_cds3)$cluster)
# 1 2 3
unique(colData(test_cds3)$cell_state)
# 1 2 3 

test_comb_cds = subcluster_cds(comb_cds, partition_name = NULL, recursive_subcluster = FALSE)
unique(colData(test_comb_cds)$partition)
# 1 - 24
unique(colData(test_comb_cds)$cluster)
# 1 - 27 
unique(colData(test_comb_cds)$cell_state)
# 1 - 27


# undebug(subcluster_cds)
test_comb_cds_2 = subcluster_cds(comb_cds, partition_name = NULL, recursive_subcluster = TRUE)
unique(colData(test_comb_cds_2)$partition)
# 1 2 
unique(colData(test_comb_cds_2)$cluster)
# 1 2 3
unique(colData(test_comb_cds_2)$cell_state)
# 1 2 3 


# should be able to plot cells and facet also 
colData(test_comb_cds_2) %>% as.data.frame %>% 
  ggplot(aes(partition_umap1, partition_umap2)) + geom_point(size=0.5) + 
  monocle3:::monocle_theme_opts() + facet_wrap(~partition, scales = "free")

plot_cells(test_comb_cds_2, color_cells_by = "partition") + 
  facet_wrap(~partition, scales = "free")


# debug partition naming ------------------------------------------------------

part_cells = colData(wt_cds) %>% as.data.frame %>% filter(partition == "5_1") %>% pull(cell)
part_cds = wt_cds[,colData(wt_cds)$cell %in% part_cells] 

reducedDims(part_cds)[["UMAP"]] = colData(part_cds) %>% as.data.frame %>% 
  select(partition_umap1, partition_umap2, partition_umap3) %>% 
  as.matrix()

plot_cells(part_cds)  

undebug(subcluster_cds)
part_cds = subcluster_cds(part_cds, recursive_subcluster = T)
unique(colData(part_cds)$partition)

part_2_cds = subcluster_cds(part_cds, recursive_subcluster = T)
unique(colData(part_2_cds)$partition)
unique(colData(part_2_cds)$cell_state)

# can you rescue clusters from making them with cell state ---------------------

# clusters(test_cds)
# 
# test_cds_part1 = test_cds[, colData(test_cds)$partition == "3_1_1"]
# unique(colData(test_cds_part1)$cluster)
# 
# partition_list = as.factor(colData(test_cds_part1)$partition)
# names(partition_list) = colnames(test_cds_part1)
# cluster_list = as.factor(colData(test_cds_part1)$cell_state)
# names(cluster_list) = colnames(test_cds_part1)

# ------------------------------------------------------------------------------



# this version starts with directly loading single cds object partition --------

partition = 3
cds_dir = "~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/R_objects/"
part_wt_cds = readRDS(paste0(cds_dir, "filtered_cds_and_transform_models/filtered_",partition,"_cds.rds"))
part_mt_cds = readRDS(paste0(cds_dir, "gap16_filtered_cds/mt_",partition,"_doublet_filtered_cds.rds"))

part_comb_cds = combine_cds(list(part_wt_cds, part_mt_cds), keep_reduced_dims = T)


plot_cells(part_comb_cds, color_cells_by = "cell_type_sub")


rm(mt_cds)
rm(wt_cds)
gc()


comb_cds = part_3_cds
# drop outliers
comb_cds = drop_outlier_cells(comb_cds)
# adjust timepoints 
comb_cds = adjust_time_stage(comb_cds)


# remove periderm -------------------------------------------------------------

remove_periderm = TRUE

if (remove_periderm) {
  
  bycatch_cell_types = c("periderm")
  comb_cds = comb_cds[,colData(comb_cds)$cell_type_broad %in% bycatch_cell_types == FALSE]
    
}


# other processing steps ------------------------------------------------------

colData(comb_cds)$sample = NULL
colData(comb_cds)$mean = NULL
# colData(comb_cds)$partition = NULL
colData(comb_cds)$cluster = NULL
colData(comb_cds)$cell_state = NULL


# what if we sub cluster beforehand --------------------------------------------

# this just runs one round of clustering
comb_cds = subcluster_cds(comb_cds,
                          # partition_name = 3,
                          num_dim = NULL,
                          max_components = 3,
                          resolution_fun = NULL,
                          max_num_cells = NULL,
                          min_res = 5e-6,
                          max_res = 1e-5)

plot_cells(comb_cds, color_cells_by = "cluster")
plot_cells(comb_cds, color_cells_by = "partition")

# save this for ease
# saveRDS(comb_cds, "R/tmp_files/comb_cds_5e06res_3_cds.rds")
comb_cds = readRDS("R/tmp_files/comb_cds_5e06res_3_cds.rds")


# run assembly directly on the cds ---------------------------------------------

unique(colData(comb_cds)$partition) #24 partitions
unique(colData(comb_cds)$cluster) #30 clusters

# run assembly on the previously run clustering
partition_results = assemble_partition_from_cds(comb_cds, 
                                                recluster = FALSE, 
                                                interval_col = "timepoint",
                                                partition_name = 3)

# if comb_cds wasn't clustered, recluster
partition_results_recluster = assemble_partition_from_cds(comb_cds, 
                                                          recluster=TRUE, 
                                                          recursive_subcluster = FALSE,
                                                          interval_col = "timepoint",
                                                          partition_name = 3)


# run assembly but subcluster
# this currently fails
partition_results_recursive = assemble_partition_from_cds(comb_cds, 
                                                          recluster = TRUE, 
                                                          recursive_subcluster = TRUE,
                                                          interval_col = "timepoint",
                                                          partition_name = 3)


# Go from coldata to cds to assembly -------------------------------------------

# this cds was clustered/has partitions, but i only want to run assembly on 1
# partition
# wrapper function for purrr use cases
run_assembly = function(cds, 
                        partition_group, 
                        interval_col = "timepoint",
                        recluster = FALSE,
                        recursive_subcluster = FALSE) {
  
  part_cds = get_partition_cds(cds, partition_group)
  partition_results = assemble_partition_from_cds(part_cds, 
                                                  recluster = recluster, 
                                                  recursive_subcluster = recursive_subcluster,
                                                  interval_col = interval_col,
                                                  partition_name = partition_group)
  return(partition_results)
}

# try running assembly on partition 1

run_assembly(comb_cds, 1)


# test different scenarios ----------------------------------------------------

# case 1: i have already done 1 round of partitioning of the cds
# i want to run the assembler on each other those WITHOUT recursively 
# partitioning, but i will need to recluster

test_df_1 = data.frame(partition=1:2) %>% 
  mutate(partition_results = purrr::map(.f = run_assembly,
                                        .x = partition, 
                                        recluster = T, 
                                        recursive_subcluster = F, 
                                        cds = comb_cds))

# can also apply this 
test_df_1_2 = data.frame(partition=1:2) %>% 
  mutate(partition_results = purrr::map(.f = run_assembly,
                                        .x = partition, 
                                        recluster = T, 
                                        recursive_subcluster = F, 
                                        cds = comb_cds))



# case 2: i have already done round 1 of partitioning of the cds, 
# but i want to further partition each of those partiions (if possible)
# i will use run assembly, recluster = T, recursive_subcluster
# but pass in the partition cds into the assembler 


test_df_2 = data.frame(partition=1:2) %>% 
          mutate(partition_results = purrr::map(.f = run_assembly,
                                                .x = partition, 
                                                recluster = T, 
                                                recursive_subcluster = T, 
                                                cds = comb_cds))



test_df_3 = data.frame("partitions" = unique(test_comb_cds_2$partition)) %>% 
  head() %>% 
  mutate(partition_results = purrr::map(.f = run_assembly,
                                        .x = partitions, 
                                        recluster = T, 
                                        recursive_subcluster = F, 
                                        cds = test_comb_cds_2))


# need to figure out how to deal with the mismatch in colnames
test_df_3$partition_results

# ------------------------------ end of testing ------------------------------ #

# now try on full data --------------------------------------------------------


rm(mt_cds)
rm(wt_cds)

# currently can't hold comb_cds in memory
comb_cds_subclustered = subcluster_cds(comb_cds,
                                       recursive_subcluster = recursive_subcluster,
                                       num_dim = NULL,
                                       max_components = 3,
                                       resolution_fun = NULL,
                                       max_num_cells = NULL,
                                       min_res = 5e-6,
                                       max_res = 1e-5)



# or run by individual partition -----------------------------------------------

# wrapper function to split it up 
# because i can't hold comb_cds in memory
run_partition_assembly = function(wt_cds, 
                                  mt_cds,
                                  partition, 
                                  coemb = FALSE){
  
  
  wt_i_cds = wt_cds[, colData(wt_cds)$partition == partition]
  mt_i_cds = mt_cds[, replace_na(colData(mt_cds)$partition == partition, F)]
  
  # my mt cds is wrong 
  
  comb_i_cds = combine_cds(list(wt_i_cds, mt_i_cds), keep_reduced_dims = T)
  
  # switch to partition umap space
  if (coemb) {
    
    comb_i_cds = comb_i_cds %>% 
      preprocess_cds(num_dim = 50) %>% 
      align_cds(residual_model_formula_str = "~log.n.umi") %>% 
      reduce_dimension(max_components = 3)
    
  } else {
    comb_i_cds = get_partition_cds(comb_i_cds, partition)
  }
  
  # remove outliers
  # comb_i_cds = drop_outlier_cells(comb_i_cds)
  
  # recursively sub cluster 
  comb_i_cds = subcluster_cds(comb_i_cds,
                              recursive_subcluster = T, 
                              partition_name = partition,
                              num_dim = NULL,
                              max_components = 3,
                              resolution_fun = NULL,
                              max_num_cells = NULL,
                              min_res = 5e-6,
                              max_res = 1e-5)
  
  partitions = unique(colData(comb_i_cds)$partition)
  unique(colData(comb_i_cds)$partition)
  
  comb_res = lapply(partitions, function(p){
    
    comb_p_cds = get_partition_cds(comb_i_cds, p)
    
    partition_results = assemble_partition_from_cds(comb_p_cds, 
                                                    recluster = T, 
                                                    recursive_subcluster = F,
                                                    interval_col = "timepoint",
                                                    partition_name = p)
    partition_results
    
  })
  
  comb_res = bind_rows(comb_res)
  
  return(comb_res)
  
}



# 



cloaca_cells = colData(kidney_cds) %>% as.data.frame %>% filter(cluster_ext_type == "Cloaca") %>% pull(cell)
length(cloaca_cells)
colData(coembed_cds)$is_cloaca = ifelse(colData(coembed_cds)$cell %in% cloaca_cells, T, F)
sum(colData(wt_cds)$cell %in% cloaca_cells) + sum(colData(mt_cds)$cell %in% cloaca_cells)
colData(wt_cds) %>% as.data.frame %>% 
  filter(cell %in% cloaca_cells) %>% 
  group_by(partition) %>% tally() %>% arrange(-n)


colData(wt_cds) %>% as.data.frame %>% 
  filter(cell %in% colData(kidney_cds)$cell) %>% 
  group_by(partition) %>% tally() %>% arrange(-n)
  


# dir.create("partition_results")
# problem w partition 2
all_partition_results = lapply(1:30, function(partition){
  
  comb_res = run_partition_assembly(wt_cds, 
                                    mt_cds, 
                                    partition)
  
  # save the partition results 
  saveRDS(comb_res, paste0("partition_results/partition_", partition, ".rds"))
  comb_res
})


all_partition_results = bind_rows(all_partition_results$partition_results)


# debug cranial sensory ganglion graph ----------------------------------------

# if identify the partition
# switch to subspace
csg_cds = get_partition_cds(comb_cds, partition_id = "")

csg_partition_results = assemble_partition_from_cds(csg_cds, 
                                                    recluster = F, 
                                                    recursive_subcluster = F,
                                                    interval_col = "timepoint",
                                                    partition_name = partition_group)




# can i reconstruct a pseudobulked cds ----------------------------------------

# given partition results
# from a full model 
# can i construct a pb_cds

cell_to_cell_state = bind_rows(partition_results$data) %>% select(cell, cell_state)

cell_to_cell_state

wt_ccs = new_cell_count_set(wt_cds, 
                            sample_group = "embryo", 
                            cell_group = "cell_state")

mt_ccs = new_cell_count_set(mt_cds, 
                            sample_group = "embryo", 
                            cell_group = "cell_state")


pb_cds = pseudobulk_ccs_for_states(wt_ccs)



# to still do: -----------------------------------------------------------------

# comb_res = bind_rows(comb_res)


cell_state_assignments = assign_cell_states(partition_results)
cell_state_assignment_summary = cell_state_assignments %>% group_by(partition) %>%
  summarize(total_cells = n(),
            num_states = length(unique(cell_state)))


plotted_graphs = partition_results %>% filter(is.na(cell_plot_state) == FALSE)


cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_state, labels=plotted_graphs$partition)

cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_type, labels=plotted_graphs$partition)

cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_time, labels=plotted_graphs$partition)

cowplot::plot_grid(plotlist = plotted_graphs$wt_state_graph_plot, labels=plotted_graphs$partition)

cowplot::plot_grid(plotlist = plotted_graphs$mt_state_graph_plot, labels=plotted_graphs$partition)


# LEFT OFF HERE

# -----------------------------------------------------------------------------

comb_res$cell_plot_time = NULL
comb_res$cell_plot_type = NULL
comb_res$cell_plot_state = NULL
comb_res$wt_state_graph_plot = NULL
comb_res$mt_state_graph_plot = NULL

# Write the cell states back to the main CDS -----------------------------------

colData(comb_cds)$cell_state = cell_state_assignments[row.names(colData(comb_cds) %>% as.data.frame()),]$cell_state



c_to_s_id_tbl = cell_state_assignments %>% select(cluster, cell_state) %>% distinct()


# this is once you have a bunch of them ----------------------------------------
# to make the global graph ----------------------------------------------------

comb_res = comb_res %>%
  mutate(mt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = mt_graph,
                                          .y=partition,
                                          cluster_to_state_id_tbl=c_to_s_id_tbl,
                                          support_col=NULL))

omb_res = comb_res %>%
  mutate(wt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = wt_graph,
                                          .y=partition,
                                          c_to_s_id_tbl))

comb_res = comb_res %>%
  mutate(wt_graph_blacklist_converted = purrr::map2(.f = convert_graph_ids,
                                                    .x = wt_graph_blacklist,
                                                    .y=partition,
                                                    c_to_s_id_tbl))

comb_res = comb_res %>%
  mutate(mt_graph_blacklist_converted = purrr::map2(.f = convert_graph_ids,
                                                    .x = mt_graph_blacklist,
                                                    .y=partition,
                                                    c_to_s_id_tbl))


# Save everything again so we capture all the cell count models and state graphs
#saveRDS(res_graphs, paste(outdir, "wt_res_graphs_cds.rds", sep="/"))

# Build a whitelist of graph edges from all the sub-CDS state graphs
global_wt_graph_edge_whitelist = do.call(igraph::union, comb_res %>% filter(is.na(wt_graph_converted) == FALSE) %>% pull(wt_graph_converted))
global_wt_graph_edge_whitelist = igraph::as_data_frame(global_wt_graph_edge_whitelist)
global_wt_graph_edge_whitelist = global_wt_graph_edge_whitelist %>% select(from, to) %>% distinct()

rm(comb_res())



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





