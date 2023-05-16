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
library(tidyverse)

### NOTES:
# Should be able to use weighted leiden when graph is intersect(corr_graph, paga_graph)
# really should be aiming to break up the partitions to avoid conflating things that aren't from
# the same lineage


#setwd("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/")
source("R/assembly_utils.R")

RhpcBLASctl::omp_set_num_threads(8)

num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
print (paste("running assembly with", num_threads, "threads"))

Sys.setenv("OMP_NUM_THREADS" = 1)
Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

outdir = "/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/"
setwd("/Users/coletrap/Google Drive/My Drive/develop/zebrafish-atlas-assembly")


combined_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/hooke-analysis/combined_ref_gap16_doubletfiltered_labelupdate202304_with_clusters_cds.rds")
colData(combined_cds)$timepoint = as.numeric(colData(combined_cds)$timepoint)

#combined_cds = combined_cds[,colData(combined_cds)$expt != "GAP13"]

#wt_cds = readRDS("~/OneDrive/UW/Trapnell/hooke_manuscript/supplement/filter_doublets/filtered_full_cds5e-06.rds")
# wt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/filtered_full_cds5e-06.rds")
# colData(wt_cds)$timepoint = as.numeric(colData(wt_cds)$timepoint)
# mt_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/gap16_proj_cds_doublet_filtered.rds")
# colData(mt_cds)$timepoint = as.numeric(colData(mt_cds)$timepoint)

assembly_distance_space = "UMAP"
#assembly_distance_space = "Aligned"

## Use these settings for clustering in PCA space
# combined_cds = cluster_cells(combined_cds, resolution=1e-4, reduction_method = assembly_distance_space)
# colData(combined_cds)$cell_state = clusters(combined_cds, reduction_method = assembly_distance_space)
# colData(combined_cds)$PAGA_partition = partitions(combined_cds, reduction_method = assembly_distance_space)

## Use these settings for clustering in UMAP space:
combined_cds = cluster_cells(combined_cds, resolution=1e-5, reduction_method = assembly_distance_space)
colData(combined_cds)$cell_state = clusters(combined_cds, reduction_method = assembly_distance_space)
colData(combined_cds)$PAGA_partition = partitions(combined_cds, reduction_method = assembly_distance_space)

saveRDS(combined_cds, "/Users/coletrap/dropbox_lab/Analysis/fish-mutants/hooke-analysis/combined_ref_gap16_doubletfiltered_labelupdate202304_with_clusters_cds.rds")

plot_cells_3d(combined_cds[, sample(ncol(combined_cds), 100000)], color_cells_by = "cell_state")

#
# saveRDS(combined_cds, "/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/full_cds_combined.rds")

# wt_cds = subcluster_cds(wt_cds,
#                         max_res = 5e-4,
#                         max_num_cells = (ncol(wt_cds)),
#                         partition_name = NULL,
#                         recursive_subcluster = TRUE)

# # save cell state column
# saveRDS(colData(wt_cds) %>% as.data.frame %>% select(cell_state, cell),
#         "wt_cds_cell_states_5e-4_cell_states.rds")


# colData(wt_cds)$cell_state =
#   readRDS("wt_cds_cell_states_5e-4_cell_states.rds")


# wt_cds[,colData(wt_cds)$cell_state == "11-1"] %>% plot_cells()


# # how many cell states are in here
# num_cell_states = unique(colData(wt_cds)$cell_state) %>% length() # 225
#
# p = plot_cells_3d(wt_cds[, sample(ncol(wt_cds), 100000)], color_cells_by = "cell_state")
# saveWidget(p, file = paste0("wt_cds_by_cell_state.html"))


# just batch
# just expt + batch
# expt * batch
# drop periderm
#

# ccs = new_cell_count_set(combined_cds,
#                          sample_group = "embryo",
#                          cell_group = "cell_state", keep_cds = T)
#
# cell_type_metadata =
#   hooke:::collect_psg_node_metadata(ccs,
#                                     color_nodes_by="germ_layer",
#                                     label_nodes_by="cell_type_broad",
#                                     group_nodes_by="tissue") %>%
#   dplyr::rename(germ_layer=color_nodes_by,
#                 cell_type_broad=label_nodes_by,
#                 tissue=group_nodes_by)
#
# periderm_groups = cell_type_metadata %>% filter(cell_type_broad == "periderm") %>% pull(id)
# periderm_counts = Matrix::colSums(normalized_counts(ccs, norm_method="size_only", pseudocount=0)[periderm_groups,])
# periderm_props = periderm_counts / Matrix::colSums(normalized_counts(ccs, norm_method="size_only", pseudocount=0))

ccs = new_cell_count_set(combined_cds,
                         sample_group = "embryo",
                         cell_group = "cell_state", keep_cds = T)

cell_type_metadata =
  hooke:::collect_psg_node_metadata(ccs,
                                    color_nodes_by="germ_layer",
                                    label_nodes_by="cell_type_broad",
                                    group_nodes_by="tissue") %>%
  dplyr::rename(germ_layer=color_nodes_by,
                cell_type_broad=label_nodes_by,
                tissue=group_nodes_by)

combined_cds = combined_cds[,!is.na(colData(combined_cds)$cell_type_broad) &
                              !colData(combined_cds)$expt == "GAP13" &
                              !colData(combined_cds)$cell_type_broad == "periderm"]

combined_cds = combined_cds[,colData(combined_cds)$expt == "GAP16"]

control_genotypes = unique(colData(combined_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]

ccm_null_wt_only = new_cell_count_model(ccs[,colData(ccs)$gene_target %in% control_genotypes],
                                        nuisance_model_formula_str = "~ns(timepoint, df = 3)",
                                        main_model_formula_str = "~1",
                                        #nuisance_model_formula_str = "~expt",
                                        #nuisance_model_formula_str = "~1",
                                        penalize_by_distance = FALSE,
                                        penalty_scale_exponent=1,
                                        num_threads=num_threads,
                                        reduction_method=assembly_distance_space,
                                        verbose=TRUE
                                        # penalty_matrix = NULL
)
ccm_null_wt_only = select_model(ccm_null_wt_only, "EBIC", 0.01)

cell_type_perturb_models_tbl = fit_mt_models(combined_cds,
                                             sample_group = "embryo",
                                             cell_group = "cell_state",
                                             main_model_formula_str = NULL,
                                             start_time = 18,
                                             stop_time = 72,
                                             interval_col="timepoint",
                                             num_time_breaks=3,
                                             #nuisance_model_formula_str = "~expt",
                                             ctrl_ids = control_genotypes,
                                             #mt_ids = mt_genotypes, #selected_mt_ids,
                                             sparsity_factor = 0.01,
                                             perturbation_col = "gene_target",
                                             keep_cds=TRUE,
                                             num_threads=num_threads,
                                             backend="nlopt",
                                             vhat_method="bootstrap",
                                             penalize_by_distance=TRUE,
                                             num_bootstraps=10)

cell_type_perturb_models_tbl = hooke:::assess_perturbation_effects(ccm_null_wt_only,
                                                                   cell_type_perturb_models_tbl,
                                                                   q_val=0.1,
                                                                   start_time = 18,
                                                                   stop_time = 72,
                                                                   perturbation_col="gene_target",
                                                                   interval_col="timepoint",
                                                                   batch_col="expt",
                                                                   interval_step = 2,
                                                                   min_interval = 4,
                                                                   max_interval = 24,
                                                                   log_abund_detection_thresh=log(1),
                                                                   min_lfc=0,
                                                                   verbose=TRUE,
                                                                   expt="GAP16")


cell_type_perturb_models_tbl = cell_type_perturb_models_tbl %>% mutate(
  cov_mat = purrr::map(.f=function(ccm) {
    cov_mat = model(ccm, "reduced") %>% sigma() %>% cov2cor()
    #cov_mat[cov_mat < 0] = 0
    diag(cov_mat) = 0
    return(cov_mat)
  }, .x = perturb_ccm),
  loss_summaries = purrr::map(.f=function(perturb_sum_tbl) {
    perturb_sum_tbl %>% pull(loss_when_present)
  }, .x = perturb_summary_tbl),
  gain_summaries = purrr::map(.f=function(perturb_sum_tbl) {
    perturb_sum_tbl %>% pull(gain_when_present)
  }, .x = perturb_summary_tbl),

)

loss_summaries = do.call(bind_cols, cell_type_perturb_models_tbl$loss_summaries) %>% as.data.frame
colnames(loss_summaries) = cell_type_perturb_models_tbl$perturb_name
loss_summaries[is.na(loss_summaries)] = 0
gain_summaries = do.call(bind_cols, cell_type_perturb_models_tbl$gain_summaries) %>% as.data.frame
colnames(gain_summaries) = cell_type_perturb_models_tbl$perturb_name
gain_summaries[is.na(gain_summaries)] = 0

effect_summaries = loss_summaries + gain_summaries
effect_summaries[loss_summaries != 0 & gain_summaries != 0] = loss_summaries[loss_summaries != 0 & gain_summaries != 0]

row.names(effect_summaries) = row.names(ccs)

breaklist = seq(-4, 4, by=0.25)
pheatmap::pheatmap(t(effect_summaries),
                   annotation_col=cell_type_metadata %>% select(germ_layer),
                   cutree_rows=4,
                   cutree_cols=16,
                   clustering_method="ward.D2",
                   #cutree_rows=num_cov_blocks,
                   #cutree_cols=num_cov_blocks,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaklist)),
                   breaks=breaklist,
                   show_rownames=TRUE,
                   show_colnames=FALSE,
                   fontsize=6,
                   filename=paste(outdir, "genotype_effect_summaries.pdf", sep="/"), width=8, height=3)





# Plot the covariance matrices
purrr::map2(.x=cell_type_perturb_models_tbl$perturb_name, .y=cell_type_perturb_models_tbl$perturb_ccm, .f = function(p_name, ccm){

  cell_cluster_cor_matrix = model(ccm, "reduced") %>% sigma() %>% cov2cor()
  diag(cell_cluster_cor_matrix) = 0

  pheatmap::pheatmap(cell_cluster_cor_matrix,
                     annotation_row=cell_type_metadata %>% select(germ_layer),
                     cluster_cols=TRUE,
                     cluster_rows=TRUE,
                     #cutree_rows=num_cov_blocks,
                     #cutree_cols=num_cov_blocks,
                     #clustering_method="ward.D2",
                     #cutree_rows=num_cov_blocks,
                     #cutree_cols=num_cov_blocks,
                     show_rownames=FALSE,
                     show_colnames=FALSE,
                     #color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaklist)),
                     #breaks = breaklist,
                     fontsize=4,
                     filename=paste(outdir, paste("cell_state_cor_matrix_heatmap_", p_name, ".pdf", sep=""), sep="/"))
})


wt_cor_matrix = model(ccm_null_wt_only, "reduced") %>% sigma() %>% cov2cor()
diag(wt_cor_matrix) = 0

cumulative_cor_matrix = accumulate(cell_type_perturb_models_tbl$cov_mat, .f=function(x, y) { return (x + y) })
cumulative_cor_matrix = cumulative_cor_matrix[[length(cell_type_perturb_models_tbl$cov_mat)]]

#cumulative_cor_matrix = cumulative_cor_matrix - length(cell_type_perturb_models_tbl$cov_mat) * wt_cor_matrix

#cumulative_cor_matrix[cumulative_cor_matrix < 0] = 0

cov_dendrogram = cumulative_cor_matrix %>% dist() %>% hclust(method="ward.D2")
num_cov_blocks = 5
dendro_cuts = cutree(cov_dendrogram, k=num_cov_blocks)

cell_type_metadata$dendro_cut = as.factor(dendro_cuts[row.names(cell_type_metadata)])

pheatmap::pheatmap(cumulative_cor_matrix,
                   annotation_row=cell_type_metadata %>% select(germ_layer, dendro_cut),
                   cluster_cols=cov_dendrogram,
                   cluster_rows=cov_dendrogram,
                   cutree_rows=num_cov_blocks,
                   cutree_cols=num_cov_blocks,
                   #clustering_method="ward.D2",
                   #cutree_rows=num_cov_blocks,
                   #cutree_cols=num_cov_blocks,
                   show_rownames=FALSE,
                   show_colnames=FALSE,
                   #color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaklist)),
                   #breaks = breaklist,
                   fontsize=4,
                   filename=paste(outdir, "cell_state_cor_matrix_heatmap.pdf", sep="/"))

#cell_cluster_cov_matrix = cell_cluster_cov_matrix[cov_dendrogram$order, cov_dendrogram$order]

# plot(cov_dendrogram)
#


colData(combined_cds)$pcor_cluster = NULL
cds_pdata = colData(combined_cds) %>% as.data.frame
cds_pdata$comb_cell_id = row.names(cds_pdata)
cds_pdata = left_join(cds_pdata, cell_type_metadata %>% select(cell_state=id, pcor_cluster=dendro_cut))
row.names(cds_pdata) = cds_pdata$comb_cell_id


colData(combined_cds)$pcor_cluster = cds_pdata[row.names(colData(combined_cds)),]$pcor_cluster

plot_cells_3d(combined_cds[,sample(ncol(combined_cds), 100000)], color_cells_by="pcor_cluster")













ccs = new_cell_count_set(combined_cds,
                         sample_group = "embryo",
                         cell_group = "cell_state", keep_cds = F)

colData(ccs)$periderm_prop = periderm_props[colnames(ccs)]

# ccs = new_cell_count_set(combined_cds,
#                          sample_group = "embryo",
#                          cell_group = "cell_type_broad", keep_cds = T)

# build

# ccm_time = new_cell_count_model(ccs,
#                            main_model_formula_str = "~1",
#                            nuisance_model_formula_str = "~ns(timepoint, df = 3)",
#                            penalize_by_distance = TRUE,
#                            num_threads=num_threads,
#                            verbose=TRUE
#                            # penalty_matrix = NULL
# )
#

wt_only = ccs[,colData(ccs)$gene_target %in% c("ctrl-inj") & colData(ccs)$expt == "GAP16"]
ccm_null_wt_only = new_cell_count_model(wt_only,
                                        nuisance_model_formula_str = "~ns(timepoint, df = 3) + periderm_prop",
                                        main_model_formula_str = "~1",
                                        #nuisance_model_formula_str = "~expt",
                                        #nuisance_model_formula_str = "~1",
                                        penalize_by_distance = FALSE,
                                        penalty_scale_exponent=1,
                                        num_threads=num_threads,
                                        reduction_method=assembly_distance_space,
                                        verbose=TRUE
                                        # penalty_matrix = NULL
)

#wt_and_muts = ccs[,colData(ccs)$gene_target %in% c("ctrl-inj") == FALSE & colData(ccs)$expt == "GAP16"]
#wt_and_muts = wt_and_muts[,sample(ncol(wt_and_muts), ncol(wt_only))]

wt_and_muts = ccs[,colData(ccs)$gene_target %in% c("tfap2a-foxd3", "ctrl-inj") & colData(ccs)$expt == "GAP16"]
wt_and_muts = wt_and_muts[,sample(ncol(wt_and_muts), ncol(wt_only))]


# time in nuisance model, so residuals mostly describe genetic perturbs
ccm_null = new_cell_count_model(wt_and_muts,
                                nuisance_model_formula_str = "~ns(timepoint, df = 3) + periderm_prop",
                                main_model_formula_str = "~1",
                                #nuisance_model_formula_str = "~expt",
                                #nuisance_model_formula_str = "~1",
                                penalize_by_distance = FALSE,
                                penalty_scale_exponent=1,
                                num_threads=num_threads,
                                reduction_method=assembly_distance_space,
                                verbose=TRUE
                                # penalty_matrix = NULL
)


cell_type_metadata =
  hooke:::collect_psg_node_metadata(ccm_null@ccs,
                                    color_nodes_by="germ_layer",
                                    label_nodes_by="cell_type_sub",
                                    group_nodes_by="tissue") %>%
  dplyr::rename(germ_layer=color_nodes_by,
                cell_type_sub=label_nodes_by,
                #timepoint=label_nodes_by,
                tissue=group_nodes_by)

ccm_null = select_model(ccm_null, "EBIC", 0.01)
cell_cluster_cov_matrix = model(ccm_null, "reduced") %>% sigma() %>% cov2cor()
diag(cell_cluster_cov_matrix) = 0

cov_dendrogram = cell_cluster_cov_matrix %>% dist() %>% hclust(method="ward.D2")
#cell_cluster_cov_matrix = cell_cluster_cov_matrix[cov_dendrogram$order, cov_dendrogram$order]

# plot(cov_dendrogram)
#
# dendro_cuts = cutree(cov_dendrogram, k=8)
# cell_type_metadata$dendro_cut = as.factor(dendro_cuts[row.names(cell_type_metadata)])
#

breaklist = seq(-0.5, 0.5, by = 0.1)

num_cov_blocks=9
pheatmap::pheatmap(cell_cluster_cov_matrix,
                   annotation_row=cell_type_metadata %>% select(germ_layer),
                   cluster_cols=cov_dendrogram,
                   cluster_rows=cov_dendrogram,
                   cutree_rows=num_cov_blocks,
                   cutree_cols=num_cov_blocks,
                   #clustering_method="ward.D2",
                   #cutree_rows=num_cov_blocks,
                   #cutree_cols=num_cov_blocks,
                   show_rownames=FALSE,
                   show_colnames=FALSE,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaklist)),
                   breaks = breaklist,
                   fontsize=4,
                   filename=paste(outdir, "cell_state_cor_matrix_heatmap.pdf", sep="/"))


ccm_null_wt_only = select_model(ccm_null_wt_only, "EBIC", 0.01)
cell_cluster_cov_matrix_wt_only = model(ccm_null_wt_only, "reduced") %>% sigma() %>% cov2cor()
diag(cell_cluster_cov_matrix_wt_only) = 0
#cell_cluster_cov_matrix_wt_only = cell_cluster_cov_matrix_wt_only[cov_dendrogram$order, cov_dendrogram$order]

cell_type_metadata =
  hooke:::collect_psg_node_metadata(ccm_null_wt_only@ccs,
                                    color_nodes_by="germ_layer",
                                    label_nodes_by="cell_type_sub",
                                    group_nodes_by="tissue") %>%
  dplyr::rename(germ_layer=color_nodes_by,
                cell_type_sub=label_nodes_by,
                #timepoint=label_nodes_by,
                tissue=group_nodes_by)


num_cov_blocks=9
pheatmap::pheatmap(cell_cluster_cov_matrix_wt_only,
                   annotation_row=cell_type_metadata %>% select(germ_layer),
                   #cluster_cols=FALSE, cluster_rows=FALSE,
                   cluster_cols=cov_dendrogram,
                   cluster_rows=cov_dendrogram,
                   cutree_rows=num_cov_blocks,
                   cutree_cols=num_cov_blocks,
                   #clustering_method="ward.D2",
                   #cutree_rows=num_cov_blocks,
                   #cutree_cols=num_cov_blocks,
                   show_rownames=FALSE,
                   show_colnames=FALSE,
                   color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name = "RdYlBu")))(length(breaklist)),
                   breaks = breaklist,
                   fontsize=4,
                   filename=paste(outdir, "cell_state_cor_matrix_wt_only_heatmap.pdf", sep="/"))

qplot(cell_cluster_cov_matrix_wt_only, cell_cluster_cov_matrix)

selected_cell_groups = cell_type_metadata %>% filter (tissue == "Kidney") %>% pull(id)
kinetics_plot = plot_cell_type_control_kinetics(ccm_null,
                                                cell_groups = selected_cell_groups,
                                                interval_col="timepoint",
                                                batch_col="expt",
                                                periderm_prop=0,
                                                log_abund_detection_thresh=log(1)) +
  monocle3:::monocle_theme_opts()
kinetics_plot


embryo_latent_pos = t(ccm_null@best_full_model$latent_pos)
embryo_fitted_vals = t(predict(ccm_null@best_full_model, newdata=cbind(colData(ccm_null@ccs), Offset=1)))# t(fitted(ccm_null@best_full_model))

embryo_latent_pos_resids = embryo_latent_pos - embryo_fitted_vals
latent_pos_cds = new_cell_data_set(embryo_latent_pos,
                                   cell_metadata=colData(ccm_null@ccs) %>% as.data.frame)
colData(latent_pos_cds)$Size_Factor = 1
latent_pos_cds = preprocess_cds(latent_pos_cds, norm_method="none", num_dim=5, scaling=FALSE)
latent_pos_cds = reduce_dimension(latent_pos_cds)
latent_pos_cds = cluster_cells(latent_pos_cds, resolution=1e-2)
plot_cells(latent_pos_cds, color_cells_by="timepoint")
plot_cells(latent_pos_cds, color_cells_by="gene_target")

# # time not in model, so residuals include time. Lots of correlation
# # structure due to time (not surprising)
# ccm_expt = new_cell_count_model(ccs,
#                                 #nuisance_model_formula_str = "~expt",
#                                 main_model_formula_str = "~1",
#                                 #nuisance_model_formula_str = "~expt",
#                                 #nuisance_model_formula_str = "~1",
#                                 penalize_by_distance = TRUE,
#                                 penalty_scale_exponent=1,
#                                 num_threads=num_threads,
#                                 reduction_method=assembly_distance_space,
#                                 verbose=TRUE
#                                 # penalty_matrix = NULL
# )
#
#
# ccm_expt = select_model(ccm_expt, "EBIC", 1)
# model(ccm_expt, "reduced") %>% sigma() %>% cov2cor() %>% #abs() %>%
#   pheatmap::pheatmap(annotation_row=cell_type_metadata %>% select(germ_layer, timepoint), fontsize=4)




# helper functions ------------------------------------------------------------

# calculate distance matrix

get_dist_df = function(ccs) {

  cell_group_centroids = centroids(ccs)
  dist_matrix = as.matrix(dist(cell_group_centroids[,-1], method = "euclidean", upper=T, diag = T))

  row.names(dist_matrix) <- cell_group_centroids$cell_group
  colnames(dist_matrix) <- cell_group_centroids$cell_group

  dist_df = dist_matrix %>% as.data.frame() %>%
    tibble::rownames_to_column("from") %>%
    pivot_longer(-from, names_to = "to")

  return(dist_df)

}

get_pcor_edges = function(pln_model, pcor_range = NULL ) {

  # pln_model = model(ccm, "reduced")
  cov_graph <- hooke:::return_igraph(pln_model)

  if (is.null(pcor_range)) {
    cov_edges <- igraph::as_data_frame(cov_graph, what="edges")
  } else {
    cov_edges <- igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(weight > max(pcor_range) | weight < min(pcor_range))

  }

  cov_graph <- cov_edges %>% igraph::graph_from_data_frame()
  return(cov_edges)
}

get_pcor_edges_wrapper = function(pln_model, threshold = NULL) {
  pln_model = pln_model$model[[1]]
  pcor_edges = get_pcor_edges(pln_model, threshold)
  pcor_edges
}

get_cor_edges = function(pln_model, cor_range = NULL ) {

  # pln_model = model(ccm, "reduced")
  cor_matrix <-  pln_model %>% sigma() %>% cov2cor()
  #row.names(cor_matrix) = colnames(cor_matrix) = row.names(counts(ccm@ccs))

  cov_graph =  igraph::graph_from_adjacency_matrix(cor_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)

  if (is.null(cor_range)) {
    cor_edges <- igraph::as_data_frame(cov_graph, what="edges")
  } else {
    cor_edges <- igraph::as_data_frame(cov_graph, what="edges") %>% dplyr::filter(weight > max(cor_range) | weight < min(cor_range))

  }

  #cov_graph <- cor_edges %>% igraph::graph_from_data_frame()
  return(cor_edges)
}

get_cor_edges_wrapper = function(pln_model, threshold = NULL) {
  pln_model = pln_model$model[[1]]
  cor_edges = get_cor_edges(pln_model, threshold)

  cor_edges
}


get_membership = function(g,  ccs, group_column="cell_type_sub", method="leiden", resolution_parameter=1e-1, weighted=TRUE) {

  if (weighted){
    #graph_weights = sqrt(abs(igraph::E(g)$weight))
    graph_weights = igraph::E(g)$weight
  }else{
    graph_weights = NULL
  }


  if (method=="leiden"){
    ldc =  leidenbase::leiden_find_partition(g, edge_weights=graph_weights, resolution_parameter=resolution_parameter, verbose=TRUE)
    #igraph::cluster_leiden(g, resolution_parameter=resolution_parameter, n_iterations=10)
    membership_df = data.frame(id=igraph::V(g)$name, "membership" = ldc$membership)
    colnames(membership_df)[1] = group_column
  }else if (method == "components"){
    membership_df = data.frame("membership" = igraph::clusters(g)$membership) %>%
      rownames_to_column(group_column)
  }

  membership_df = colData(ccs@cds) %>% as.data.frame %>%
    select(!!sym(group_column)) %>% distinct() %>%
    left_join(membership_df, by = group_column)

  num_new_states = membership_df$membership[is.na(membership_df$membership)] %>% length()
  num_curr_states = max(membership_df$membership, na.rm=TRUE)

  membership_df$membership[is.na(membership_df$membership)] = seq(num_curr_states+1, num_new_states+num_curr_states)

  #membership = as.character(membership)
  membership_df$membership = as.character(membership_df$membership)
  membership_df
}
#debug(get_membership)

assign_pcor_components <- function(ccm_sparsity, clustering_resolution, ccm, paga_graph=NULL, cor_range=NULL, reduction_method=assembly_distance_space){
  ccm = select_model(ccm, sparsity_factor=ccm_sparsity)

  cov_edges = get_cor_edges(model(ccm, "reduced"), cor_range=cor_range)
  cov_edges = cov_edges %>% filter(from != "" & to != "")
  cov_graph = igraph::graph_from_data_frame(cov_edges, directed=FALSE)

  if (is.null(paga_graph) == FALSE){
      #cov_graph = igraph::union(cov_graph, paga_graph)
    cov_graph = igraph::intersection(cov_graph, paga_graph)
  }

  pcor_cluster_df = get_membership(cov_graph, ccm@ccs, group_column="cell_state", method="leiden", resolution_parameter=clustering_resolution)
  return(pcor_cluster_df)
}

collect_pcor_component_stats <- function(ccm_sparsity, clustering_resolution, ccm, cell_type_metadata, paga_graph=TRUE, cor_range=NULL){
  ccm = select_model(ccm, sparsity_factor=ccm_sparsity)
  pcor_cluster_df = assign_pcor_components(ccm_sparsity, clustering_resolution, ccm, paga_graph=paga_graph, cor_range=cor_range)

  cds_pdata = left_join(cell_type_metadata, pcor_cluster_df, by=c("id"="cell_state"))

  pcor_stat_df = tibble(
    cell_type_sub_ARI = adj_rand_index(cds_pdata$membership, cds_pdata$cell_type_sub),
    cell_type_broad_ARI = adj_rand_index(cds_pdata$membership, cds_pdata$cell_type_broad),
    germ_layer_ARI = adj_rand_index(cds_pdata$membership, cds_pdata$germ_layer),
    tissue_ARI = adj_rand_index(cds_pdata$membership, cds_pdata$tissue),
    timepoint_ARI = adj_rand_index(cds_pdata$membership, cds_pdata$timepoint),
    PAGA_partition_ARI = adj_rand_index(cds_pdata$membership, cds_pdata$PAGA_partition)
  )


  return(pcor_stat_df)
}



collect_pcor_component_link_stats <- function(ccm_sparsity, clustering_resolution, ccm, cell_type_metadata, paga_graph=NULL){
  ccm = select_model(ccm, sparsity_factor=ccm_sparsity)

  pcor_cluster_df = assign_pcor_components(ccm_sparsity, clustering_resolution, ccm, paga_graph=paga_graph)

  cell_type_metadata = cell_type_metadata %>% left_join(pcor_cluster_df, by=c("id"="cell_state"))
  dist_df = get_dist_df(ccs)

  cov_edges = get_cor_edges(model(ccm, "reduced"))
  dist_vs_pcor_df = left_join(dist_df, cov_edges) %>%
    rename(umap_distance=value, pcor=weight) %>%
    mutate(pcor = ifelse(is.na(pcor), 0, pcor))

  dist_vs_pcor_df = dplyr::left_join(dist_vs_pcor_df, cell_type_metadata %>% setNames(paste0('to_', names(.))), by=c("to"="to_id")) #%>%
  dist_vs_pcor_df = dplyr::left_join(dist_vs_pcor_df, cell_type_metadata %>% setNames(paste0('from_', names(.))), by=c("from"="from_id")) #%>%
  dist_vs_pcor_df = dist_vs_pcor_df %>%
    mutate(same_cell_type_sub = from_cell_type_broad == to_cell_type_sub,
           same_cell_type_broad = from_cell_type_broad == to_cell_type_broad,
           same_tissue = from_tissue == to_tissue,
           same_pcor_component = from_membership == to_membership)
  dist_vs_pcor_df = dist_vs_pcor_df %>% filter(from != to)
  return(dist_vs_pcor_df)
}

# from mclust package
adj_rand_index <- function (x, y)
{
  x <- as.vector(x)
  y <- as.vector(y)
  if (length(x) != length(y))
    stop("arguments must be vectors of the same length")
  tab <- table(x, y)
  if (all(dim(tab) == c(1, 1)))
    return(1)
  a <- sum(choose(tab, 2))
  b <- sum(choose(rowSums(tab), 2)) - a
  c <- sum(choose(colSums(tab), 2)) - a
  d <- choose(sum(tab), 2) - a - b - c
  ARI <- (a - (a + b) * (a + c)/(a + b + c + d))/((a + b +
                                                     a + c)/2 - (a + b) * (a + c)/(a + b + c + d))
  return(ARI)
}
#debug(adj_rand_index)

#debug(assign_pcor_components)
cell_type_metadata =
  hooke:::collect_psg_node_metadata(ccm_null@ccs,
                                    color_nodes_by="germ_layer",
                                    label_nodes_by="cell_type_broad",
                                    group_nodes_by="tissue") %>%
  dplyr::rename(germ_layer=color_nodes_by,
                cell_type_broad=label_nodes_by,
                tissue=group_nodes_by) %>% left_join(
  hooke:::collect_psg_node_metadata(ccm_null@ccs,
                                    color_nodes_by="PAGA_partition",
                                    label_nodes_by="timepoint",
                                    group_nodes_by="cell_type_sub")) %>%
  dplyr::rename(PAGA_partition=color_nodes_by,
                timepoint=label_nodes_by,
                cell_type_sub=group_nodes_by)


# pcor_component_params = expand.grid(sparsity=c(0.1, seq(1,10, by=2)), resolution=10^seq(-5,-1)) %>% as_tibble() # including sparsity 50 basically represents PAGA only
# pcor_component_stats = pcor_component_params %>%
#   mutate(
#     pcor_component_df = purrr::map2(.f=assign_pcor_components,
#                              .x=sparsity,
#                              .y=resolution,
#                              ccm_null,
#                              pcor_range=c(0, Inf),
#                              include_paga = FALSE)) %>%
#   # mutate(
#   # pcor_link_stats = purrr::map2(.f=collect_pcor_component_link_stats,
#   #                          .x=sparsity,
#   #                          .y=resolution,
#   #                          ccm_null,
#   #                          cell_type_metadata)) %>%
#   mutate(
#     pcor_stats = purrr::map2(.f=collect_pcor_component_stats,
#                              .x=sparsity,
#                              .y=resolution,
#                              ccm_null,
#                              cell_type_metadata,
#                              pcor_range=c(-0.01, Inf),
#                              include_paga = TRUE))
# pcor_component_stats = pcor_component_stats %>% tidyr::unnest(pcor_stats)

paga_graph = hooke:::get_paga_graph(ccm_null@ccs@cds, reduction_method=assembly_distance_space)

pcor_component_params = expand.grid(sparsity=c(0.01, 0.1, seq(1,10, by=2)), resolution=10^seq(-3,0)) %>% as_tibble() # including sparsity 50 basically represents PAGA only
#pcor_component_params = expand.grid(sparsity=10^seq(-3,0), resolution=10^seq(-3,0)) %>% as_tibble() # including sparsity 50 basically represents PAGA only

cor_filter_range = c(-Inf, 0.05)
pcor_component_stats = pcor_component_params %>%
  mutate(
    pcor_component_df = purrr::map2(.f=assign_pcor_components,
                                    .x=sparsity,
                                    .y=resolution,
                                    ccm=ccm_null,
                                    cor_range=cor_filter_range,
                                    paga_graph = paga_graph)) %>%
  mutate(
    pcor_stats = purrr::map2(.f=collect_pcor_component_stats,
                             .x=sparsity,
                             .y=resolution,
                             ccm=ccm_null,
                             cell_type_metadata=cell_type_metadata,
                             cor_range=cor_filter_range,
                             paga_graph = paga_graph))
pcor_component_stats = pcor_component_stats %>% tidyr::unnest(pcor_stats)

pcor_component_stats = pcor_component_stats %>%
  tidyr::unnest(pcor_component_df) %>%
  group_by(sparsity, resolution) %>%
  mutate(num_pcor_components = length(unique(membership))) %>%
  select(-membership, -cell_state)

ARI_table = pcor_component_stats %>% select(sparsity,
                                            resolution,
                                            num_pcor_components,
                                            cell_type_sub_ARI,
                                            cell_type_broad_ARI,
                                            germ_layer_ARI,
                                            tissue_ARI,
                                            timepoint_ARI,
                                            PAGA_partition_ARI)

# ARI_table = ARI_table %>% pivot_longer(cols = ends_with("ARI")) %>% rename(ARI_type=name, ARI=value)
# qplot(sparsity, ARI, color=as.factor(resolution), geom="line", data=ARI_table) + facet_wrap(~ARI_type)

ARI_table = ARI_table %>% pivot_longer(cols = c("tissue_ARI", "germ_layer_ARI", "num_pcor_components", "cell_type_broad_ARI")) %>% rename(ARI_type=name, ARI=value)
qplot(sparsity, ARI, color=as.factor(resolution), geom="line", data=ARI_table) + facet_wrap(~ARI_type, scales="free_y")


# pcor_component_stat_summary = pcor_component_stats %>%
#   group_by(sparsity, resolution) %>%
#   summarize(num_nz_pcors = sum(pcor != 0),
#             num_components = length(unique(c(to_membership, from_membership))))
# qplot(sparsity, num_components, color=as.factor(resolution), geom="line", data=pcor_component_stat_summary)
# qplot(sparsity, num_nz_pcors, color=as.factor(resolution), geom="line", data=pcor_component_stat_summary)

#pcor_comp_df = assign_pcor_components(5, 0.1, ccm_null, pcor_range=c(-0.01, Inf), include_paga = TRUE)

# This worked well with PCA distances (i.e. components corresponded to germ layers or large lineages even across PAGA partitions)
pcor_comp_df = assign_pcor_components(0.01, 0.001, ccm_null, cor_range=cor_filter_range, paga_graph = paga_graph)

#pcor_comp_df = assign_pcor_components(1, 0.1, ccm_null, pcor_range=c(-0.01, 0.01), include_paga = FALSE)


cds_pdata = colData(combined_cds) %>% as.data.frame
cds_pdata$comb_cell_id = row.names(cds_pdata)

cds_pdata = left_join(cds_pdata, pcor_comp_df, by="cell_state")
row.names(cds_pdata) = cds_pdata$comb_cell_id

colData(combined_cds)$pcor_cluster = cds_pdata[row.names(colData(combined_cds)),]$membership

#plot_cells(combined_cds[,sample(ncol(combined_cds), 100000)], color_cells_by="pcor_cluster")

pheatmap::pheatmap(table(colData(combined_cds)$pcor_cluster, colData(combined_cds)$germ_layer), scale="column")
plot_cells_3d(combined_cds[,sample(ncol(combined_cds), 100000)], color_cells_by="pcor_cluster")

pheatmap::pheatmap(table(partitions(combined_cds), colData(combined_cds)$tissue), scale="column")



saveRDS(combined_cds, "/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/full_cds_combined.rds")

# Some extra stuff for exploring
select_pcor_group_cds = combined_cds[,colData(combined_cds)$pcor_cluster == "11"]

select_pcor_group_cds = select_pcor_group_cds %>%
  preprocess_cds(num_dim = 50) %>%
  align_cds(residual_model_formula_str = "~log.n.umi") %>%
  reduce_dimension(max_components = 3)
select_pcor_group_cds = cluster_cells(select_pcor_group_cds, resolution = 5e-4)
plot_cells_3d(select_pcor_group_cds, color_cells_by="cell_type_sub")


# Do cell types that emerge at the same time or that are most abundant at the same time have higher pcor?

# Do cell types from the same tissue or germ layer have higher pcor?

# run hooke --------------------------------------------------------------------


# for each cluster which sparsity parameter that causes that node to become isolated

# scan over regularization path 30 modles that

# once a cell state becomes isolated at a given lambda -- it is also isolated

# what kinds of cell types become isolated early vs late
  # small vs late
  # early in time vs late
  # etc

# stars method?



ccs = new_cell_count_set(wt_cds,
                         sample_group = "embryo",
                         cell_group = "cell_state", keep_cds = T)

# build

ccm = new_cell_count_model(ccs,
                           main_model_formula_str = "~1",
                           nuisance_model_formula_str = "~ns(timepoint, df = 3)*expt ",
                           penalize_by_distance = FALSE,
                           num_threads=num_threads
                           # penalty_matrix = NULL
                           )


ccm_null = new_cell_count_model(ccs,
                           main_model_formula_str = "~1",
                           nuisance_model_formula_str = "~1",
                           # penalize_by_distance = FALSE,
                           # penalty_matrix = NULL
)

ccm_inter = new_cell_count_model(ccs,
                                main_model_formula_str = "~1",
                                nuisance_model_formula_str = "~ ns(timepoint, df = 3)*expt",
                                # penalize_by_distance = FALSE,
                                # penalty_matrix = NULL
)

# ccm = select_model(ccm, criterion = "EBIC", sparsity_factor = 0.1)


get_penalty_df = function(ccm) {

  penalty_list = lapply(1:30, function(i) {
    ccm@reduced_model_family$models[[i]]$penalty
  })

  penalty_df = data.frame("penalty" = unlist(penalty_list))
  penalty_df$model = unlist(ccm@reduced_model_family$models)
  penalty_df = penalty_df %>% tidyr::nest(data=model)

  penalty_df = penalty_df %>%
    mutate(cov_edges = purrr::map(.f = get_cor_edges_wrapper,
                                  .x = data,
                                  threshold = NULL)) %>%
    mutate(graph = purrr::map(.f = igraph::graph_from_data_frame,
                              .x = cov_edges)) %>%
    mutate(graph_plot = purrr::map(.f = plot,
                                   .x = graph,
                                   edge.arrow.size = 0.25,
                                   vertex.label.dist = 1.5,
                                   vertex.size = 5,
                                   vertex.label = NA)) %>%
    mutate(membership_df = purrr::map2(.f = get_membership,
                                       .x = graph,
                                       .y = cov_edges,
                                       ccs = ccs,
                                       wt_cds = wt_cds)) %>%
    tidyr::unnest(membership_df)

  return(penalty_df)


}


# color by penalty where cluster is first isolated --------------------------------------


color_cds_by_penalty = function(cds, penalty_df, colname="penalty"){
  min_penalty = penalty_df %>%
    group_by(penalty, membership) %>%
    mutate(is_isolated = ifelse(n() == 1, T, F)) %>%
    ungroup %>%
    filter(is_isolated) %>%
    group_by(cell_state) %>%
    top_n(wt = -penalty, n = 1) %>%
    select(cell_state, penalty)

  colData(cds)[[colname]] = colData(wt_cds) %>%
    as.data.frame() %>%
    left_join(min_penalty, by = "cell_state") %>%
    pull(penalty)
  return(cds)
}

#

penalty_df = get_penalty_df(ccm,)
penalty_null_df = get_penalty_df(ccm_null)
penalty_inter_df = get_penalty_df(ccm_inter)

wt_cds = color_cds_by_penalty(wt_cds, ccm, penalty_df, colname = "penalty")
wt_cds = color_cds_by_penalty(wt_cds, ccm_null, penalty_null_df, colname = "penalty_null")
wt_cds = color_cds_by_penalty(wt_cds, ccm_inter, penalty_inter_df, colname = "penalty_inter")


p = plot_cells_3d(wt_cds[, sample(ncol(wt_cds), 100000)], color_cells_by = "penalty")
saveWidget(p, file = paste0("wt_cds_by_min_penalty.html"))

p = plot_cells_3d(wt_cds[, sample(ncol(wt_cds), 100000)], color_cells_by = "penalty_null")
saveWidget(p, file = paste0("wt_cds_by_min_penalty_null.html"))

p = plot_cells_3d(wt_cds[, sample(ncol(wt_cds), 100000)], color_cells_by = "penalty_inter")
saveWidget(p, file = paste0("wt_cds_by_min_penalty_inter.html"))



# get cell type info

total_counts = counts(ccs) %>%
  as.matrix() %>%
  Matrix::rowSums() %>%
  as.data.frame() %>%
  rownames_to_column("cell_state") %>%
  rename("total_count" = ".")




penalty_df %>%
  group_by(penalty, membership) %>%
  mutate(is_isolated = ifelse(n() == 1, T, F)) %>%
  left_join(total_counts, by = "cell_state") %>%
  ggplot(aes(x = penalty, y = total_count, color = is_isolated)) + geom_point() +
  monocle3:::monocle_theme_opts()


penalty_df %>%
  group_by(penalty, membership) %>%
  mutate(is_isolated = ifelse(n() == 1, T, F)) %>%
  left_join(total_counts, by = "cell_state") %>%
  ungroup %>%
  group_by(penalty) %>%
  summarise(num_isolated = sum(is_isolated),
            num_not_isolated =
            avg_isolated_count = mean(total_count)) %>%
  ggplot(aes(x = penalty, y = num_isolated)) + geom_point()


# plot num components

penalty_df %>%
  group_by(penalty) %>%
  summarize(num_components = max(membership)) %>%
  ggplot(aes(x = penalty, y = num_components)) + geom_point()






# ------------------------------------------------------------------------------

# kidney stuff

kid_cs = colData(wt_cds) %>% as.data.frame %>%
  filter(tissue == "Kidney")%>%
  group_by(cell_state, cell_type_sub) %>% tally() %>% filter(n > 10) %>% select(-n)

kid_cs

membership %>% filter(cell_state %in% kid_cs$cell_state)
kid_cs %>% filter(cell_state == "9-6")
kid_cs %>% filter(cell_state == "9-5")


# repartition by partial corr graph --------------------------------------------

membership_coldatas = lapply(unique(membership$membership), function(m) {

  if (!is.na(m)) {

    cell_states = membership[membership$membership == m,"cell_state"]
    wt_sub_cds = wt_cds[, colData(wt_cds)$cell_state %in% cell_states]

    # print(paste0(m, "total cells", ncol(wt_sub_cds))
    wt_sub_cds = wt_sub_cds %>%
      preprocess_cds(num_dim = 30) %>%
      align_cds(residual_model_formula_str = "~log.n.umi") %>%
      reduce_dimension(max_components = 3)

    plot_cells_3d(wt_sub_cds, color_cells_by = "cell_state")
    plot_cells_3d(wt_sub_cds, color_cells_by = "cell_type_sub")
    plot_cells_3d(wt_sub_cds, color_cells_by = "partition")

    colData(wt_sub_cds)$auto_partition_umap3d_1 = reducedDims(wt_sub_cds)[["UMAP"]][,1]
    colData(wt_sub_cds)$auto_partition_umap3d_2 = reducedDims(wt_sub_cds)[["UMAP"]][,2]
    colData(wt_sub_cds)$auto_partition_umap3d_3 = reducedDims(wt_sub_cds)[["UMAP"]][,3]
    as.data.frame(colData(wt_sub_cds))
  }


})

membership_df = do.call(rbind, membership_coldatas)


# # where are the kidney cells
#
# kidney_cds <- readRDS("~/OneDrive/UW/Trapnell/hooke/examples/R_objects/pronephros.RDS")
# kidney_cells = colData(kidney_cds)$cell
#
#
# # pheatmap::pheatmap(sigma(ccm@best_reduced_model) %>% cov2cor())
#
#
#
# # which cell types correspond to which partition or cell state
#
# cts_to_st = colData(wt_cds) %>%
#   as.data.frame() %>%
#   # filter(cell %in% kidney_cells) %>%
#   group_by(cell_type_sub, cell_state) %>%
#   tally() %>% top_n(n=1, wt=n)
#
#
#
#
# colData(wt_cds) %>%
#   as.data.frame() %>%
#   filter(cell %in% kidney_cells) %>%
#   group_by(cell_state, cell_type_sub) %>%
#   tally() %>% ungroup() %>%
#   group_by(cell_type_sub) %>%
#   top_n(n = 1, wt=n)
#
# #3_2-9
# #4-41
# # 9-1, 9-2, 9-3, 9-5, 9-5, 9-6
# kidney_partitions = c("9-1", "9-2", "9-3", "9-5", "9-6")
#
# cov_edge_df = cov_edges %>%
#   mutate(from_partition = stringr::str_split(from, "-") %>% map_chr(., 1)) %>%
#   mutate(to_partition = stringr::str_split(to, "-") %>% map_chr(., 1)) %>%
#   mutate(same_partition = (from_partition == to_partition)) %>%
#   arrange(-weight)
#
#
# cov_edge_df %>%
#   filter(to == "3_2-9" | from == "3_2-9") %>%
#   group_by(from, to, same_partition) %>%
#   top_n(1, weight)
#
#
# cov_edge_df %>%
#   filter(to == "4-41" | from == "4-41") %>%
#   group_by(from, to, same_partition) %>%
#   top_n(1, weight)
#
# cov_edge_df %>%
#   filter(to %in% kidney_partitions | to %in% kidney_partitions) %>%
#   group_by(from, to, same_partition) %>%
#   top_n(1, weight)
#
#
#
#
#
# # compare the distances
#
# # 3_2-9
# # 4-41
# # 9-1, 9-2, 9-3, 9-5, 9-5, 9-6
#
# left_join(cov_edges, dist_df, by = c("from", "to")) %>%
#   mutate(c = case_when(
#     to == "3_2-9" & from == "4-14" ~T,
#     from == "3_2-9" & to == "4-14" ~ T,
#     to %in% kidney_partitions & from %in% c("3_2-9", "4-14") ~ T,
#     from %in% kidney_partitions & to %in% c("3_2-9", "4-14") ~ T,
#     TRUE ~ F)) %>%
#   filter(weight > 0) %>%
#   ggplot(aes(weight, value, color=c)) + geom_point()
#
#
# dist_df %>%
#   mutate(p = case_when(
#   to == "3_2-9" & from == "4-14" ~T,
#   from == "3_2-9" & to == "4-14" ~ T,
#   to %in% kidney_partitions & from %in% c("3_2-9", "4-14") ~ T,
#   from %in% kidney_partitions & to %in% c("3_2-9", "4-14") ~ T,
#   TRUE ~ F)) %>%
#   filter(p)
