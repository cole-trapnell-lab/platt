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
options(future.globals.maxSize = 2* 8192 * 1024^2) # 16GB

assembly_start_time = 18
assembly_stop_time = 72
assembly_backend = "nlopt"

#devtools::load_all("~/OneDrive/UW/Trapnell/hooke_pln_20230104/")

# -----------------------------------------------------------------------------

# FIXME: set this to number of CPUs on machine (or in cluster session)
#RhpcBLASctl::omp_set_num_threads(10)
RhpcBLASctl::omp_set_num_threads(10)

num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
print (paste("running assembly with", num_threads, "threads"))

Sys.setenv("OMP_NUM_THREADS" = 1)
Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

#outdir = "/Users/maddyduran/UW/Trapnell/zebrafish-atlas-assembly/"
#setwd("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/")
outdir = "//Users/coletrap/tmp_output"
setwd("/Users/coletrap/Google Drive/My Drive/develop/zebrafish-atlas-assembly")
source("R/assembly_utils.R")

# starting with the full cds --------------------------------------------------

#wt_cds = readRDS("~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/filtered_full_cds5e-06.rds")
#mt_cds = readRDS("~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/gap16_proj_cds_doublet_filtered.rds")
wt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/filtered_full_cds5e-06.rds")
colData(wt_cds)$timepoint = as.numeric(colData(wt_cds)$timepoint)
mt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/gap16_proj_cds_doublet_filtered.rds")
colData(mt_cds)$timepoint = as.numeric(colData(mt_cds)$timepoint)

wt_cds = drop_outlier_cells(wt_cds)
mt_cds = drop_outlier_cells(mt_cds)

control_genotypes = unique(colData(wt_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]


remove_periderm = TRUE

if (remove_periderm) {

  bycatch_cell_types = c("periderm")
  wt_cds = wt_cds[,colData(wt_cds)$cell_type_broad %in% bycatch_cell_types == FALSE]
  mt_cds = mt_cds[,colData(mt_cds)$cell_type_broad %in% bycatch_cell_types == FALSE]
  gc()
}


# dir.create("partition_results")
# problem w partition 2

num_partitions = length(unique(colData(wt_cds)$partition))
all_partition_results = lapply(1:num_partitions, function(partition){

  comb_res = run_partition_assembly(wt_cds,
                                    mt_cds,
                                    partition,
                                    #coemb=TRUE,
                                    coemb=FALSE,
                                    #max_res=5e-4,
                                    max_res=1e-4,
                                    max_num_cells=(ncol(wt_cds) + ncol(mt_cds)),
                                    ctrl_ids = control_genotypes,
                                    #mt_ids = mt_genotypes,
                                    #mt_ids = c("mafba"),
                                    num_threads=num_threads)

  # save the partition results
  #saveRDS(comb_res, paste0("partition_results/partition_", partition, ".rds"))
  comb_res
})


# FIXME: HACK: THIS IS NEEDED BECAUSE OF THOSE EXTRA PARTITION NAMES:
#all_partition_results$partition = stringr::str_replace(all_partition_results$partition, "_1-", "")

all_partition_results = do.call(bind_rows, all_partition_results)

all_partition_results_fixed = lapply(all_partition_results, function(res){
  res$partition = as.character(res$partition)
  return(res)
})

all_partition_results_fixed = do.call(bind_rows, all_partition_results_fixed)
all_partition_results = all_partition_results_fixed
rm(all_partition_results_fixed)

png(paste(outdir, "subassembly_cell_state_plots.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$cell_plot_state, labels=all_partition_results$partition)
dev.off()

png(paste(outdir, "subassembly_cell_type_plots.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$cell_plot_type, labels=all_partition_results$partition)
dev.off()

png(paste(outdir, "subassembly_time_plots.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$cell_plot_time, labels=all_partition_results$partition)
dev.off()

png(paste(outdir, "subassembly_wt_state_graph_plot.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$wt_state_graph_plot, labels=all_partition_results$partition)
dev.off()

png(paste(outdir, "subassembly_mt_state_graph_plot.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$mt_state_graph_plot, labels=all_partition_results$partition)
dev.off()

all_partition_results$cell_plot_state = NULL
all_partition_results$cell_plot_time = NULL
all_partition_results$cell_plot_type = NULL
all_partition_results$wt_state_graph_plot = NULL
all_partition_results$mt_state_graph_plot = NULL

saveRDS(all_partition_results, paste(outdir, "all_partition_results.rds", sep="/"))

gc()

comb_cds = combine_cds(list(wt_cds, mt_cds), keep_reduced_dims = T)

cell_state_assignments = assign_cell_states(all_partition_results)
cell_state_assignment_summary = cell_state_assignments %>% group_by(partition) %>%
  summarize(total_cells = n(),
            num_states = length(unique(cell_state)))
colData(comb_cds)$cell_state = cell_state_assignments[row.names(colData(comb_cds) %>% as.data.frame()),]$cell_state


c_to_s_id_tbl = cell_state_assignments %>% select(cluster, cell_state) %>% distinct()

rm(wt_cds)
rm(mt_cds)





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

remove_periderm = TRUE

if (remove_periderm) {

  bycatch_cell_types = c("periderm")
  wt_cds = wt_cds[,colData(wt_cds)$cell_type_broad %in% bycatch_cell_types == FALSE]
  mt_cds = mt_cds[,colData(mt_cds)$cell_type_broad %in% bycatch_cell_types == FALSE]
  gc()
}


# dir.create("partition_results")
# problem w partition 2
num_partitions = length(unique(colData(wt_cds)$partition))
all_partition_results = lapply(1:30, function(partition){

  comb_res = run_partition_assembly(wt_cds,
                                    mt_cds,
                                    partition)

  # save the partition results
  #saveRDS(comb_res, paste0("partition_results/partition_", partition, ".rds"))
  comb_res
})




comb_cds = combine_cds(list(wt_cds, mt_cds), keep_reduced_dims = T)
control_genotypes = unique(colData(comb_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]



#rm(wt_cds)
#rm(mt_cds)
#gc()

# # can i reconstruct my sub cds from my big cds --------------------------------
# # aka pull partition umap coords from coldata
# part_2_cds = get_partition_cds(comb_cds, 2)
# plot_cells(part_2_cds)
#
# part_3_cds = get_partition_cds(comb_cds, 3)
# plot_cells(part_3_cds)
#
# part_4_cds = get_partition_cds(comb_cds, 4)
# plot_cells(part_4_cds)
#
#
# # test recursive clustering ----------------------------------------------------
#
# # debug(subcluster_cds)
# test_cds = subcluster_cds(part_3_cds, partition_name = 3, recursive_subcluster = TRUE)
# # should be able to separate
# unique(colData(test_cds)$partition)
# # "3_1_1" "3_2_1"
# unique(colData(test_cds)$cluster)
# # 1 2
# unique(colData(test_cds)$cell_state)
# # "3_1-1" "3_1-2" "3_2-1
#
# test_cds2 = subcluster_cds(part_3_cds, partition_name = 3, recursive_subcluster = FALSE)
#
# unique(colData(test_cds2)$partition)
# # "3_1" "3_2"
# unique(colData(test_cds2)$cluster)
# # 1 2 3
# unique(colData(test_cds2)$cell_state)
# # "3-1" "3-2" "3-3
#
# # just recluster once
# test_cds3 = subcluster_cds(part_3_cds, partition_name = NULL, recursive_subcluster = FALSE)
# unique(colData(test_cds3)$partition)
# # 1 2
# unique(colData(test_cds3)$cluster)
# # 1 2 3
# unique(colData(test_cds3)$cell_state)
# # 1 2 3
#
# test_comb_cds = subcluster_cds(comb_cds, partition_name = NULL, recursive_subcluster = FALSE)
# unique(colData(test_comb_cds)$partition)
# # 1 - 24
# unique(colData(test_comb_cds)$cluster)
# # 1 - 27
# unique(colData(test_comb_cds)$cell_state)
# # 1 - 27
#
#
# # undebug(subcluster_cds)
# test_comb_cds_2 = subcluster_cds(comb_cds, partition_name = NULL, recursive_subcluster = TRUE)
# unique(colData(test_comb_cds_2)$partition)
# # 1 2
# unique(colData(test_comb_cds_2)$cluster)
# # 1 2 3
# unique(colData(test_comb_cds_2)$cell_state)
# # 1 2 3
#
#
# # should be able to plot cells and facet also
# colData(test_comb_cds_2) %>% as.data.frame %>%
#   ggplot(aes(partition_umap1, partition_umap2)) + geom_point(size=0.5) +
#   monocle3:::monocle_theme_opts() + facet_wrap(~partition, scales = "free")
#
# plot_cells(test_comb_cds_2, color_cells_by = "partition") +
#   facet_wrap(~partition, scales = "free")
#
#
# # can you rescue clusters from making them with cell state ---------------------
#
# # clusters(test_cds)
# #
# # test_cds_part1 = test_cds[, colData(test_cds)$partition == "3_1_1"]
# # unique(colData(test_cds_part1)$cluster)
# #
# # partition_list = as.factor(colData(test_cds_part1)$partition)
# # names(partition_list) = colnames(test_cds_part1)
# # cluster_list = as.factor(colData(test_cds_part1)$cell_state)
# # names(cluster_list) = colnames(test_cds_part1)
#
# # ------------------------------------------------------------------------------
#
#
#
# # this version starts with directly loading single cds object partition --------
#
# partition = 3
# cds_dir = "~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/R_objects/"
# part_wt_cds = readRDS(paste0(cds_dir, "filtered_cds_and_transform_models/filtered_",partition,"_cds.rds"))
# part_mt_cds = readRDS(paste0(cds_dir, "gap16_filtered_cds/mt_",partition,"_doublet_filtered_cds.rds"))
#
# part_comb_cds = combine_cds(list(part_wt_cds, part_mt_cds), keep_reduced_dims = T)
#
#
# plot_cells(part_comb_cds, color_cells_by = "cell_type_sub")
#
#
# rm(mt_cds)
# rm(wt_cds)
# gc()
#
#
# comb_cds = part_3_cds
# # drop outliers
# comb_cds = drop_outlier_cells(comb_cds)
# # adjust timepoints
# comb_cds = adjust_time_stage(comb_cds)


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
saveRDS(comb_cds, paste(outdir, "comb_cds_5e06res_3_cds.rds", sep="/"))
comb_cds = readRDS(paste(outdir, "comb_cds_5e06res_3_cds.rds", sep="/"))

# dist_from_other_cell_type_sub_cells = function(ref_coldata, centroids_per_cts, i, cell_type_label) {
#   curr_loc = ref_coldata[i,c("umap3d_1","umap3d_2", "umap3d_3")]
#   centroid_loc = centroids_per_cts[cell_type_label,]
#   eucl.dist = dist(rbind(curr_loc,centroid_loc) %>% as.matrix)
#   return(eucl.dist[[1]])
# }
#
# # this currently needs to be run on the full cds
# clean_up_labels = function(cds,
#                            old_colname,
#                            new_colname = old_colname,
#                            k = 10,
#                            max_dist_fct = 10) {
#
#   # build a knn
#   cds = make_cds_nn_index(cds, reduction_method = "UMAP")
#   query_ann = cds@reduce_dim_aux$UMAP$nn_index$annoy$nn_index$annoy_index
#   query_dims = reducedDims(cds)[["UMAP"]]
#   query_res = uwot:::annoy_search(query_dims, k = k + 1, ann = query_ann)
#   ref_coldata = colData(cds) %>% as.data.frame()
#   ref_coldata = ref_coldata %>% mutate(rn = row_number())
#
#   ref_coldata$umap3d_1 = reducedDim(x = cds, type = "UMAP")[,1]
#   ref_coldata$umap3d_2 = reducedDim(x = cds, type = "UMAP")[,2]
#   ref_coldata$umap3d_3 = reducedDim(x = cds,type = "UMAP")[,3]
#
#
#   # where are the rest of those labels
#   centroids_per_cts = reducedDims(cds)[["UMAP"]] %>% as.data.frame
#   centroids_per_cts$cell_type_sub = colData(cds)$cell_type_sub
#   centroids_per_cts = aggregate(.~cell_type_sub, data = centroids_per_cts, FUN=mean)
#   centroids_per_cts = centroids_per_cts %>% tibble::column_to_rownames("cell_type_sub")
#   names(centroids_per_cts) = c("umap3d_1","umap3d_2", "umap3d_3")
#
#   #
#   new_labels = sapply(ref_coldata$rn, function(i) {
#     self_label = ref_coldata[query_res$idx[i,1], "cell_type_sub"]
#     neighbor_indices = query_res$idx[i,2:(k+1)]
#     neighbor_labels = ref_coldata[neighbor_indices, "cell_type_sub"]
#     neighbor_dists = query_res$dist[i,2:(k+1)]
#     # mean_dist = mean(neighbor_dists)
#     max_dist = max(neighbor_dists)
#     min_dist = min(neighbor_dists)
#     majority_label = monocle3:::which_mode(neighbor_labels)
#     eucl.dist = dist_from_other_cell_type_sub_cells(ref_coldata, centroids_per_cts, i, self_label)
#
#     if (self_label != "") {
#       eucl.dist = dist_from_other_cell_type_sub_cells(ref_coldata, i, self_label)
#     } else {
#       # we can keep these as ""
#       eucl.dist = 0
#     }
#
#
#     # only change it if far from like-things and clear majority
#     if (eucl.dist > max_dist_fct*max_dist & !is.na(majority_label)) {
#       majority_label
#     } else {
#       self_label
#     }
#
#   })
#
#   colData(cds)[[new_colname]] = new_labels
#   return(cds)
# }

# comb_cds_subclustered = comb_cds_subclustered[,is.na(colData(comb_cds_subclustered)$cell_type_sub) == FALSE]
# comb_cds_subclustered = comb_cds_subclustered[,is.na(colData(comb_cds_subclustered)$cell_type_broad) == FALSE]
#
# colData(comb_cds_subclustered)$orig_cell_type_sub = colData(comb_cds_subclustered)$cell_type_sub
# colData(comb_cds_subclustered)$orig_cell_type_broad = colData(comb_cds_subclustered)$cell_type_broad
# comb_cds_subclustered = clean_up_labels(comb_cds_subclustered, "orig_cell_type_sub", "cell_type_sub")
# comb_cds_subclustered = clean_up_labels(comb_cds_subclustered, "orig_cell_type_broad", "cell_type_broad")

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
# run_assembly = function(cds,
#                         partition_group,
#                         interval_col = "timepoint",
#                         recluster = FALSE,
#                         recursive_subcluster = FALSE,
#                         ...) {
#
#   part_cds = get_partition_cds(cds, partition_group)
#   partition_results = assemble_partition_from_cds(part_cds,
#                                                   recluster = recluster,
#                                                   recursive_subcluster = recursive_subcluster,
#                                                   interval_col = interval_col,
#                                                   partition_name = partition_group)
#   return(partition_results)
# }

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


# rm(mt_cds)
# rm(wt_cds)
#
# comb_cds = cluster_cells(comb_cds, resolution= 5e-7, k = 15)
#
# # currently can't hold comb_cds in memory
# comb_cds = subcluster_cds(comb_cds,
#                                        recursive_subcluster = T,
#                                        num_dim = NULL,
#                                        max_components = 3,
#                                        #resolution_fun = NULL,
#                                        max_num_cells = ncol(comb_cds),
#                                        min_res = 5e-7,
#                                        max_res = 5e-4,
#                           cluster_k=15)
#
# partition_assembly_df = data.frame("partitions" = unique(colData(comb_cds)$partition)) %>%
#   #head() %>%
#   mutate(partition_results = purrr::map(.f = run_assembly,
#                                         .x = partitions,
#                                         recluster = T,
#                                         recursive_subcluster = T,
#                                         cds = comb_cds,
#                                         ctrl_ids = control_genotypes,
#                                         min_res = 5e-7,
#                                         max_res = 5e-4,
#                                         max_num_cells=ncol(comb_cds)))
# partition_assembly_df = partition_assembly_df %>% tidyr::unnest(partition_results)
#
#

cell_state_assignments = assign_cell_states(partition_assembly_df)
cell_state_assignment_summary = cell_state_assignments %>% group_by(partition) %>%
  summarize(total_cells = n(),
            num_states = length(unique(cell_state)))


plotted_graphs = partition_assembly_df %>% filter(is.na(cell_plot_state) == FALSE)

png(paste(outdir, "all_sub_cds_umaps_by_cell_state.png", sep="/"), width=20, height=20, units="in", res=300)
cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_state, labels=plotted_graphs$partition)
dev.off()



# # or run by individual partition -----------------------------------------------
#
# # wrapper function to split it up
# # because i can't hold comb_cds in memory
# run_partition_assembly = function(wt_cds,
#                                   mt_cds,
#                                   partition,
#                                   coemb = FALSE,
#                                   resolution_fun = NULL,
#                                   max_num_cells = NULL,
#                                   min_res = 5e-6,
#                                   max_res = 1e-5,
#                                   ...){
#
#
#   wt_i_cds = wt_cds[, colData(wt_cds)$partition == partition]
#   mt_i_cds = mt_cds[, replace_na(colData(mt_cds)$partition == partition, F)]
#
#   # my mt cds is wrong
#
#   comb_i_cds = combine_cds(list(wt_i_cds, mt_i_cds), keep_reduced_dims = T)
#
#   # switch to partition umap space
#
#   if (coemb) {
#
#     comb_i_cds = comb_i_cds %>%
#       preprocess_cds(num_dim = 50) %>%
#       align_cds(residual_model_formula_str = "~log.n.umi") %>%
#       reduce_dimension(max_components = 3)
#
#   } else {
#     comb_i_cds = get_partition_cds(comb_i_cds, partition)
#   }
#
#   # remove outliers
#   # comb_i_cds = drop_outlier_cells(comb_i_cds)
#
#   # recursively sub cluster
#   comb_i_cds = subcluster_cds(comb_i_cds,
#                               recursive_subcluster = T,
#                               partition_name = partition,
#                               num_dim = NULL,
#                               max_components = 3,
#                               resolution_fun = resolution_fun,
#                               max_num_cells = max_num_cells,
#                               min_res = min_res,
#                               max_res = max_res)
#
#   partitions = unique(colData(comb_i_cds)$partition)
#   unique(colData(comb_i_cds)$partition)
#
#   comb_res = lapply(partitions, function(p){
#
#     comb_p_cds = get_partition_cds(comb_i_cds, p)
#
#     partition_results = assemble_partition_from_cds(comb_p_cds,
#                                                     recluster = T,
#                                                     recursive_subcluster = F,
#                                                     interval_col = "timepoint",
#                                                     partition_name = p,
#                                                     max_num_cells = max_num_cells,
#                                                     min_res = min_res,
#                                                     max_res = max_res,
#                                                     ...)
#     partition_results
#
#   })
#
#   comb_res = bind_rows(comb_res)
#
#   return(comb_res)
#
# }



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
num_partitions = length(unique(colData(wt_cds)$partition))
all_partition_results = lapply(1:30, function(partition){

  comb_res = run_partition_assembly(wt_cds,
                                    mt_cds,
                                    partition)

  # save the partition results
  #saveRDS(comb_res, paste0("partition_results/partition_", partition, ".rds"))
  comb_res
})


all_partition_results = bind_rows(all_partition_results$partition_results)


# debug cranial sensory ganglion graph ----------------------------------------

csg_assembly = run_partition_assembly(wt_cds,
                       mt_cds,
                       "7",
                       coemb=TRUE,
                       max_res=1e-2,
                       max_num_cells=(ncol(wt_cds) + ncol(mt_cds)),
                       ctrl_ids = control_genotypes)

#png(paste(outdir, "all_sub_cds_umaps_by_cell_state.png", sep="/"), width=20, height=20, units="in", res=300)
cowplot::plot_grid(plotlist = csg_assembly$cell_plot_state, labels=csg_assembly$partition)
#dev.off()

cowplot::plot_grid(plotlist = csg_assembly$cell_plot_type, labels=csg_assembly$partition)

cowplot::plot_grid(plotlist = csg_assembly$cell_plot_time, labels=csg_assembly$partition)

cowplot::plot_grid(plotlist = csg_assembly$wt_state_graph_plot, labels=csg_assembly$partition)

cowplot::plot_grid(plotlist = csg_assembly$mt_state_graph_plot, labels=csg_assembly$partition)



# if identify the partition
# switch to subspace
csg_cds = get_partition_cds(comb_cds_subclustered, partition_id = "7_1")

bad_cell_types = colData(csg_cds) %>% as.data.frame %>% group_by(cell_type_sub) %>% tally() %>% filter (n < 50) %>% pull(cell_type_sub)
colData(csg_cds)$cell_type_sub[colData(csg_cds)$cell_type_sub %in% bad_cell_types] = NA_character_

csg_cds = fix_missing_cell_labels(csg_cds, reduction_method='UMAP', from_column_name='cell_type_sub')
plot_cells(csg_cds, color_cells_by="cell_type_sub")


csg_cds = cluster_cells(csg_cds, resolution=1e-3)
plot_cells(csg_cds, color_cells_by="cluster")


#csg_cds = cluster_cells(csg_cds)

csg_partition_results = assemble_partition_from_cds(csg_cds,
                                                    recluster = F,
                                                    recursive_subcluster = F,
                                                    interval_col = "timepoint",
                                                    partition_name = "partition",
                                                    ctrl_ids = control_genotypes)


# pronephros ----------------------------------------



# # if identify the partition
# # switch to subspace
# pronephros_cds = get_partition_cds(comb_cds_subclustered, partition_id = "10_1")
# pronephros_cds = cluster_cells(pronephros_cds)
# pronephros_classifier = readRDS("./notebooks/supporting_data/pronephros_classifier.RDS")
# # pronephros_classifier <- train_cell_classifier(cds = kidney_cds,
# #                                          marker_file = "./examples/pronephros_cell_types.txt",
# #                                          db="none",#org.Dr.eg.db::org.Dr.eg.db,
# #                                          cds_gene_id_type = "SYMBOL",
# #                                          num_unknown = 50,
# #                                          marker_file_gene_id_type = "SYMBOL")
#
# # Commenting this out so I don't actually overwrite the classifier:
# #saveRDS(pronephros_classifier, "pronephros_classifier.RDS")
#
# colData(pronephros_cds)$garnett_cluster = clusters(pronephros_cds)
# pronephros_cds = classify_cells(pronephros_cds, pronephros_classifier, db="none", cluster_extend=TRUE)
# plot_cells(pronephros_cds, color_cells_by="cluster_ext_type")
#
# plot_cells(pronephros_cds, color_cells_by="cell_type_sub")
#
# proneprhos_partition_results = assemble_partition_from_cds(pronephros_cds,
#                                                     recluster = T,
#                                                     recursive_subcluster = F,
#                                                     interval_col = "timepoint",
#                                                     partition_name = "partition",
#                                                     ctrl_ids = control_genotypes)




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


cell_state_assignments = assign_cell_states(partition_assembly_df)
cell_state_assignment_summary = cell_state_assignments %>% group_by(partition) %>%
  summarize(total_cells = n(),
            num_states = length(unique(cell_state)))


plotted_graphs = partition_assembly_df %>% filter(is.na(cell_plot_state) == FALSE)

png(paste(outdir, "all_sub_cds_umaps_by_cell_state.png", sep="/"), width=20, height=20, units="in", res=300)
cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_state, labels=plotted_graphs$partition)
dev.off()

cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_type, labels=plotted_graphs$partition)

cowplot::plot_grid(plotlist = plotted_graphs$cell_plot_time, labels=plotted_graphs$partition)

pdf(paste(outdir, "sub_cds_wt_state_graphs.pdf", sep="/"), width=20, height=20)
cowplot::plot_grid(plotlist = plotted_graphs$wt_state_graph_plot, labels=plotted_graphs$partition)
dev.off()

pdf(paste(outdir, "sub_cds_mt_state_graphs.pdf", sep="/"), width=30, height=30)
cowplot::plot_grid(plotlist = plotted_graphs$mt_state_graph_plot, labels=plotted_graphs$partition)
dev.off()

# LEFT OFF HERE

# -----------------------------------------------------------------------------

partition_assembly_df$cell_plot_time = NULL
partition_assembly_df$cell_plot_type = NULL
partition_assembly_df$cell_plot_state = NULL
partition_assembly_df$wt_state_graph_plot = NULL
partition_assembly_df$mt_state_graph_plot = NULL

# Write the cell states back to the main CDS -----------------------------------

colData(comb_cds)$cell_state = cell_state_assignments[row.names(colData(comb_cds) %>% as.data.frame()),]$cell_state



c_to_s_id_tbl = cell_state_assignments %>% select(cluster, cell_state) %>% distinct()


# this is once you have a bunch of them ----------------------------------------
# to make the global graph ----------------------------------------------------

partition_assembly_df = partition_assembly_df %>%
  mutate(mt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = mt_graph,
                                          .y=partitions,
                                          cluster_to_state_id_tbl=c_to_s_id_tbl,
                                          support_col=NULL))

partition_assembly_df = partition_assembly_df %>%
  mutate(wt_graph_converted = purrr::map2(.f = convert_graph_ids,
                                          .x = wt_graph,
                                          .y=partitions,
                                          c_to_s_id_tbl))

partition_assembly_df = partition_assembly_df %>%
  mutate(wt_graph_blacklist_converted = purrr::map2(.f = convert_graph_ids,
                                                    .x = wt_graph_blacklist,
                                                    .y=partitions,
                                                    c_to_s_id_tbl))

partition_assembly_df = partition_assembly_df %>%
  mutate(mt_graph_blacklist_converted = purrr::map2(.f = convert_graph_ids,
                                                    .x = mt_graph_blacklist,
                                                    .y=partitions,
                                                    c_to_s_id_tbl))


# Save everything again so we capture all the cell count models and state graphs
#saveRDS(res_graphs, paste(outdir, "wt_res_graphs_cds.rds", sep="/"))

# Build a whitelist of graph edges from all the sub-CDS state graphs
global_wt_graph_edge_whitelist = do.call(igraph::union, partition_assembly_df %>% filter(is.na(wt_graph_converted) == FALSE) %>% pull(wt_graph_converted))
global_wt_graph_edge_whitelist = igraph::as_data_frame(global_wt_graph_edge_whitelist)
global_wt_graph_edge_whitelist = global_wt_graph_edge_whitelist %>% select(from, to) %>% distinct()

#rm(comb_res())



# Fit a single wild-type cell count timeseries model to all the cell states at once
system.time({global_wt_ccm = fit_wt_model(comb_cds,
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
global_wt_graph = assemble_wt_graph(comb_cds,
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
global_perturb_models_tbl = fit_mt_models(comb_cds,
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


# Build a whitelist of edges by collecting the edges in the subassemblies from the mutants
sub_mt_graph_union = do.call(igraph::union, partition_assembly_df %>% filter(is.na(mt_graph_converted) == FALSE) %>% pull(mt_graph_converted))
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

plot_state_graph_perturb_effects(global_wt_ccm, global_mt_graph,
                                 label_nodes_by="cell_state",
                                 #color_nodes_by = "timepoint",
                                 group_nodes_by="cell_type_sub",
                                 edge_weights = "num_perturbs_supporting",
                                 label_edges_by="support_label",
                                 hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_mt_graph_perturb_effects.pdf", sep="/"), width=50, height=25, limitsize = FALSE)


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


# OK, let's build the final, global graph that we'll use for downstream computations
# We'll do this by annotating the wild-type graph with the graph built from perturbations,
# and adding any perturbation-supported edges that aren't in the WT graph.
# You might think all the mutant edges should exist in the WT graph, bu because we break cycles in both graphs using different heuristics, it's not quite that straightforward
global_wt_graph_edges = igraph::as_data_frame(global_wt_graph)
global_mt_graph_edges = igraph::as_data_frame(global_mt_graph)
mt_only = setdiff(global_mt_graph_edges %>% select(from, to), global_wt_graph_edges %>% select(from, to))

global_graph_annotated = left_join(global_wt_graph_edges, global_mt_graph_edges)
global_graph_annotated = global_graph_annotated %>% select(-support)
global_graph_annotated = rbind(global_graph_annotated,
                               global_mt_graph_edges %>% inner_join(mt_only))
global_graph_annotated = igraph::graph_from_data_frame(global_graph_annotated)

plot_state_graph_annotations(global_wt_ccm, global_graph_annotated,
                             label_nodes_by="cell_state",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "num_perturbs_supporting",
                             label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE)
ggplot2::ggsave(paste(outdir, "global_graph_annotated.pdf", sep="/"), width=50, height=25, limitsize = FALSE)

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





