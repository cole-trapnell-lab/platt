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
#outdir = "//Users/coletrap/tmp_output"
outdir = "/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/"
setwd("/Users/coletrap/Google Drive/My Drive/develop/zebrafish-atlas-assembly")
source("R/assembly_utils.R")

# load data -------------------------------------------------------------------

#combined_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/full_cds_combined.rds")
combined_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/hooke-analysis/combined_ref_gap16_doubletfiltered_labelupdate202304_with_clusters_cds.rds")
colData(combined_cds)$timepoint = as.numeric(colData(combined_cds)$timepoint)

#combined_cds = clean_up_labels(combined_cds, "germ_layer", "germ_layer")

# Need to have run:
# combined_cds = cluster_cells(combined_cds, resolution=1e-5, reduction_method = "UMAP")

control_genotypes = unique(colData(combined_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]

#selected_mt_ids = c("tbxta", "cdx4", "noto", "tfap2a", "foxd3")
selected_mt_ids = NULL

combined_cds = combined_cds[,!is.na(colData(combined_cds)$germ_layer) & colData(combined_cds)$expt %in% c("GAP13") == FALSE]
#combined_cds = combined_cds[,colData(combined_cds)$gene_target %in% c(control_genotypes, selected_mt_ids)]
#combined_cds = combined_cds
#rm(combined_cds)
#combined_cds = combined_cds[,sample(ncol(combined_cds), 250000)]

colData(combined_cds)$sample = NULL
colData(combined_cds)$mean = NULL


assemble_component = function(cds,
                                  partition,
                                  partition_col="pcor_component",
                                  coemb = FALSE,
                                  resolution_fun = NULL,
                                  max_num_cells = NULL,
                                  min_res=5e-7,
                                  max_res=5e-4,
                              num_threads=1,
                                  ...){

  partition_cds = cds[, !is.na(colData(cds)[,partition_col]) & colData(cds)[,partition_col] == partition]

  # switch to partition umap space

  if (coemb) {

    partition_cds = partition_cds %>%
      preprocess_cds(num_dim = 50) %>%
      align_cds(residual_model_formula_str = "~log.n.umi") %>%
      reduce_dimension(max_components = 3,
                       cores=num_threads)

    partition_cds = drop_outlier_cells(partition_cds)
    gc()


  } else {
    partition_cds = get_partition_cds(partition_cds, partition, partition_col)
  }


  partition_results = assemble_partition_from_cds(partition_cds,
                                                  recluster = T,
                                                  recursive_subcluster = F,
                                                  interval_col = "timepoint",
                                                  partition_name = partition,
                                                  max_num_cells = max_num_cells,
                                                  min_res = min_res,
                                                  max_res = max_res,
                                                  num_threads=num_threads,
                                                  component_col=partition_col,
                                                  ...)
  partition_results


  partition_results

  return(partition_results)

}

#xxx_cds = cluster_cells(combined_cds[,sample(ncol(combined_cds), 25000)], resolution=1e-5)
#plot_cells(xxx_cds)


subassembly_group_ccs = new_cell_count_set(combined_cds, "embryo", "subassembly_group")
cell_type_metadata =
  hooke:::collect_psg_node_metadata(subassembly_group_ccs,
                                    color_nodes_by="germ_layer",
                                    label_nodes_by="cell_type_broad",
                                    group_nodes_by="tissue") %>%
  dplyr::rename(germ_layer=color_nodes_by,
                cell_type_broad=label_nodes_by,
                tissue=group_nodes_by)

periderm_groups = cell_type_metadata %>% filter(cell_type_broad == "periderm") %>% pull(id)
periderm_counts = Matrix::colSums(normalized_counts(subassembly_group_ccs[periderm_groups,], norm_method="size_only", pseudocount=0))
periderm_props = periderm_counts / Matrix::colSums(normalized_counts(subassembly_group_ccs, norm_method="size_only", pseudocount=0))
colData(combined_cds)$periderm_prop = periderm_props[colData(combined_cds)$embryo]

# Remove periderm
combined_cds = combined_cds[,colData(combined_cds)$subassembly_group %in% periderm_groups == FALSE]

# recompute size factors with periderm removed:
subassembly_group_ccs = new_cell_count_set(combined_cds, "embryo", "subassembly_group")
embryo_sz_fcts = size_factors(subassembly_group_ccs)
rm(subassembly_group_ccs)
gc()


#### Debugging:
# Create a count set, grouping cells by subassembly_group

# # Notochord
# notochord_subassembly = assemble_component(combined_cds,
#                                            "11",
#                                            partition_col = "pcor_cluster",
#                                            coemb=TRUE,
#                                            #coemb=FALSE,
#                                            min_res=1e-7,
#                                            max_res=5e-4,
#                                            start_time=assembly_start_time,
#                                            stop_time=assembly_stop_time,
#                                            expts_excluded_from_assembly=c("GAP13"),
#                                            #nuisance_model_formula_str="~expt", # FIXME: put this back if we start working on more than just GAP16!!!
#                                            max_num_cells=ncol(combined_cds),
#                                            ctrl_ids = control_genotypes,
#                                            mt_ids = selected_mt_ids,
#                                            #mt_ids = c("zc4h2"),
#                                            num_threads=num_threads,
#                                            num_bootstraps=10,
#                                            embryo_size_factors = embryo_sz_fcts
# )

# Neural crest
nc_subassembly = assemble_component(combined_cds,
                                           "mesoderm/neural crest",
                                           partition_col = "germ_layer",
                                           coemb=TRUE,
                                           #coemb=FALSE,
                                           min_res=1e-7,
                                           max_res=5e-4,
                                           start_time=assembly_start_time,
                                           stop_time=assembly_stop_time,
                                           expts_excluded_from_assembly=c("GAP13"),
                                           #nuisance_model_formula_str="~expt", # FIXME: put this back if we start working on more than just GAP16!!!
                                           max_num_cells=ncol(combined_cds),
                                           ctrl_ids = control_genotypes,
                                           mt_ids = selected_mt_ids,
                                           #mt_ids = c("zc4h2"),
                                           num_threads=num_threads,
                                           num_bootstraps=10
)


# Mesoderm
mesoderm_subassembly = assemble_component(combined_cds,
                                    "mesoderm",
                                    partition_col = "germ_layer",
                                    coemb=TRUE,
                                    #coemb=FALSE,
                                    min_res=1e-7,
                                    max_res=5e-4,
                                    start_time=assembly_start_time,
                                    stop_time=assembly_stop_time,
                                    expts_excluded_from_assembly=c("GAP13"),
                                    #nuisance_model_formula_str="~expt", # FIXME: put this back if we start working on more than just GAP16!!!
                                    max_num_cells=ncol(combined_cds),
                                    ctrl_ids = control_genotypes,
                                    mt_ids = selected_mt_ids,
                                    #mt_ids = c("zc4h2"),
                                    num_threads=num_threads,
                                    num_bootstraps=10
)





partitions_for_assembly = unique(colData(combined_cds)$germ_layer)

#partitions_for_assembly = unique(colData(combined_cds)$pcor_cluster)
# partitions_for_assembly =  c("6",
#                                                   "18", "13", # Cranial NC + pigment
#                                                   "17", # cranial muscle
#                                                   "27", # schwann cells (and retina)
#                                                   "22", # vascular endothelium/aorta
#                                                   "40", # chondrocranium
#                                                   "31", # jaw chondrocytes
#                                                   "23" # cardiomyocytes
# )
#
# partitions_for_assembly =  c("18", "13", # Cranial NC + pigment
#                              "23" # cardiomyocytes
# )


num_partitions = length(partitions_for_assembly)

all_partition_results = lapply(partitions_for_assembly, function(partition){

  comb_res = assemble_component(combined_cds,
                                    partition,
                                    partition_col = "germ_layer",
                                    coemb=TRUE,
                                    #coemb=FALSE,
                                min_res=1e-7,
                                max_res=5e-4,
                                start_time=assembly_start_time,
                                stop_time=assembly_stop_time,
                                expts_excluded_from_assembly=c("GAP13"),
                                #nuisance_model_formula_str="~expt", # FIXME: put this back once we figure out how to also add periderm fraction as a nuisance into the model!
                                    max_num_cells=ncol(combined_cds),
                                    ctrl_ids = control_genotypes,
                                    mt_ids = selected_mt_ids,
                                    #mt_ids = c("mafba"),
                                    num_threads=num_threads#,
                                #embryo_size_factors = embryo_sz_fcts
                                )

  # save the partition results
  #saveRDS(comb_res, paste0("partition_results/partition_", partition, ".rds"))
  comb_res
})


all_partition_results = do.call(bind_rows, all_partition_results)
all_partition_results$partition = partitions_for_assembly

png(paste(outdir, "subassembly_cell_state_plots.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$cell_plot_state, labels=all_partition_results$partition)
dev.off()

png(paste(outdir, "subassembly_cell_type_plots.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$cell_plot_type, labels=all_partition_results$partition)
dev.off()

png(paste(outdir, "subassembly_time_plots.png", sep="/"), width=40, height=40,  units="in", res=300)
cowplot::plot_grid(plotlist = all_partition_results$cell_plot_time, labels=all_partition_results$partition)
dev.off()

# png(paste(outdir, "subassembly_wt_state_graph_plot.png", sep="/"), width=40, height=40,  units="in", res=300)
# cowplot::plot_grid(plotlist = all_partition_results$wt_state_graph_plot)
# dev.off()
#
# png(paste(outdir, "subassembly_mt_state_graph_plot.png", sep="/"), width=40, height=40,  units="in", res=300)
# cowplot::plot_grid(plotlist = all_partition_results$mt_state_graph_plot, labels=all_partition_results$partition)
# dev.off()

all_partition_results$cell_plot_state = NULL
all_partition_results$cell_plot_time = NULL
all_partition_results$cell_plot_type = NULL
all_partition_results$wt_state_graph_plot = NULL
all_partition_results$mt_state_graph_plot = NULL

saveRDS(all_partition_results, paste(outdir, "all_partition_results.rds", sep="/"))


# For now, we need to load this up to attached the subassembly info to the CDS
all_partition_results = readRDS(paste(outdir, "all_partition_results.rds", sep="/"))
cell_state_assignments = all_partition_results %>% select(data) %>% tidyr::unnest(data)
cell_state_assignments$subassembly_group = as.character(cell_state_assignments$subassembly_group)
# Do this if you want to use all the cells:
cell_downsample = cell_state_assignments
cell_downsample = cell_downsample %>% distinct() %>% as.data.frame(stringsAsFactors=FALSE)
row.names(cell_downsample) = cell_downsample$cds_row_id
combined_cds = combined_cds[,cell_downsample$cds_row_id]

colData(combined_cds)$subassembly_group = cell_downsample[row.names(colData(combined_cds) %>% as.data.frame()),]$subassembly_group
colData(combined_cds)$pcor_cluster = stringr::str_split_fixed(colData(combined_cds)$subassembly_group, "-", 2)[,1]

#
# # collect the assignments to subassembly groups and write them back to the main CDS:
# cell_state_assignments = all_partition_results %>% select(data) %>% tidyr::unnest(data)
# cell_state_assignments$subassembly_group = as.character(cell_state_assignments$subassembly_group)
# cell_state_assignments = cell_state_assignments %>% distinct() %>% as.data.frame(stringsAsFactors=FALSE)
# row.names(cell_state_assignments) = cell_state_assignments$cds_row_id
# colData(combined_cds)$subassembly_group = cell_state_assignments[row.names(colData(combined_cds) %>% as.data.frame()),]$subassembly_group

# Combine all the subassembly graphs into a single graph object
wt_subassembly_graphs = all_partition_results %>%
  #filter(is.na(wt_graph) == FALSE ) %>%
  pull(wt_graph)

wt_subassembly_graphs = lapply(wt_subassembly_graphs, function(x) {
  if(is.null(x) == FALSE && is.na(x) == FALSE) return(igraph::as_data_frame(x))})
subassembly_graph_wt_union = do.call(bind_rows, wt_subassembly_graphs)
subassembly_graph_wt_union = igraph::graph_from_data_frame(subassembly_graph_wt_union, vertices = unique(colData(combined_cds)$subassembly_group))


# Combine all the subassembly graphs into a single graph object
mt_subassembly_graphs = all_partition_results %>%
  #filter(is.na(mt_graph) == FALSE ) %>%
  pull(mt_graph)

mt_subassembly_graphs = lapply(mt_subassembly_graphs, function(x) {
  if(is.null(x) == FALSE && is.na(x) == FALSE) return(igraph::as_data_frame(x))})
subassembly_graph_mt_union = do.call(bind_rows, mt_subassembly_graphs)
subassembly_graph_mt_union = igraph::graph_from_data_frame(subassembly_graph_mt_union, vertices = unique(colData(combined_cds)$subassembly_group))

rm(all_partition_results)
gc()

# Create a count set, grouping cells by subassembly_group
subassembly_group_ccs = new_cell_count_set(combined_cds, "embryo", "subassembly_group")

pdf(paste(outdir, "subassembly_wt_graph_cell_type_sub.pdf", sep="/"), width=40, height=20)
plot_state_graph_annotations(subassembly_group_ccs, subassembly_graph_wt_union,
                             label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             edge_weights = "support",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=TRUE)
dev.off()

pdf(paste(outdir, "subassembly_wt_graph_cell_type_sub_small.pdf", sep="/"), width=6, height=6)
plot_state_graph_annotations(subassembly_group_ccs, subassembly_graph_wt_union,
                             #label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="tissue",
                             edge_weights = "support",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=FALSE)
dev.off()

pdf(paste(outdir, "subassembly_mt_graph_cell_type_sub.pdf", sep="/"), width=40, height=20)
plot_state_graph_annotations(subassembly_group_ccs, subassembly_graph_mt_union,
                             label_nodes_by="subassembly_group",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_sub",
                             edge_weights = "total_perturb_path_score_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=TRUE)
dev.off()

pdf(paste(outdir, "subassembly_mt_graph_cell_type_sub_small.pdf", sep="/"), width=6, height=6)
plot_state_graph_annotations(subassembly_group_ccs, subassembly_graph_mt_union,
                             #label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="tissue",
                             edge_weights = "total_perturb_path_score_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=FALSE)
dev.off()

pb_cds = hooke:::pseudobulk_ccs_for_states(subassembly_group_ccs)

plan(sequential) # I don't know why by future does not play nice with monocle's threading:
#debug(classify_genes_over_graph)
gene_patterns_over_subassembly_graph = classify_genes_over_graph(subassembly_group_ccs, subassembly_graph_mt_union, abs_expr_thresh=1e-2, log_fc_thresh=2, cores=8)
marker_scores_over_subassembly_graph = top_markers(pb_cds, group_cells_by = "cell_group", cores=8)

top_subassembly_group_specific_markers = marker_scores_over_subassembly_graph %>%
  filter(fraction_expressing >= 0.10) %>%
  group_by(cell_group) %>%
  top_n(1, pseudo_R2)

top_specific_marker_ids = unique(top_subassembly_group_specific_markers %>% pull(gene_id))

pdf(paste(outdir, "subassembly_group_top_markers_dotplot.pdf", sep="/"), width=20, height=40)
plot_genes_by_group(pb_cds,
                    top_specific_marker_ids,
                    group_cells_by="cell_group",
                    ordering_type="maximal_on_diag",
                    max.size=3)
dev.off()

plan(multicore)

# Now contract that single graph object on cell type
cell_type_graph_wt = contract_state_graph(subassembly_group_ccs, subassembly_graph_wt_union, group_nodes_by="cell_type_sub")
#cell_type_graph = igraph::delete_vertices(cell_type_graph, "Unknown")
cell_type_graph_wt = hooke:::break_cycles_in_state_transition_graph(cell_type_graph_wt, "support")
saveRDS(cell_type_graph_wt, paste(outdir,"cell_type_graph_wt.rds", sep="/"))


# Now contract that single graph object on cell type
cell_type_graph_mt = contract_state_graph(subassembly_group_ccs, subassembly_graph_mt_union, group_nodes_by="cell_type_sub")
#cell_type_graph = igraph::delete_vertices(cell_type_graph, "Unknown")
cell_type_graph_mt = hooke:::break_cycles_in_state_transition_graph(cell_type_graph_mt, "total_perturb_path_score_supporting")
saveRDS(cell_type_graph_mt, paste(outdir,"cell_type_graph_mt.rds", sep="/"))

cell_type_graph_mt = readRDS(paste(outdir,"cell_type_graph_mt.rds", sep="/"))
cell_type_graph_wt = readRDS(paste(outdir,"cell_type_graph_wt.rds", sep="/"))

cell_type_ccs = new_cell_count_set(combined_cds, "embryo", "cell_type_sub", keep_cds=FALSE)

pdf(paste(outdir, "global_wt_graph_cell_type_sub.pdf", sep="/"), width=40, height=10)
plot_state_graph_annotations(cell_type_ccs, cell_type_graph_wt,
                             label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             #edge_weights = "support",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=TRUE)
dev.off()

pdf(paste(outdir, "global_mt_graph_cell_type_sub.pdf", sep="/"), width=40, height=10)
plot_state_graph_annotations(cell_type_ccs, cell_type_graph_mt,
                             label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             edge_weights = "total_perturb_path_score_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=TRUE)
dev.off()

pdf(paste(outdir, "global_mt_graph_cell_type_sub_small.pdf", sep="/"), width=40, height=10)
plot_state_graph_annotations(cell_type_ccs, cell_type_graph_mt,
                             #label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="tissue",
                             edge_weights = "total_perturb_path_score_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = FALSE,
                             #con_colour = "black",
                             label_groups=TRUE)
dev.off()


# Now fit models on cell types so we can plot the genetic requirements over the cell type graph
system.time({cell_type_ccm = fit_wt_model(combined_cds,
                                          sample_group = "embryo",
                                          cell_group = "cell_type_sub",
                                          #main_model_formula_str = NULL,
                                          #main_model_formula_str = "~1",
                                          #main_model_formula_str = "~ns(timepoint, df=5)",
                                          start_time = assembly_start_time,
                                          stop_time = assembly_stop_time,
                                          interval_col="timepoint",
                                          vhat_method="bootstrap",
                                          num_time_breaks=4,
                                          #nuisance_model_formula_str = "~expt",
                                          ctrl_ids = control_genotypes,
                                          sparsity_factor = 0.01,
                                          perturbation_col = "gene_target",
                                          edge_whitelist = igraph::as_data_frame(cell_type_graph_wt),
                                          #base_penalty=1000,
                                          keep_cds = TRUE,
                                          num_threads=num_threads,
                                          backend=assembly_backend,
                                          verbose=TRUE,
                                          penalize_by_distance=TRUE,
                                          pln_num_penalties=30)})

global_graph_edge_whitelist = igraph::as_data_frame(cell_type_graph_wt)
global_graph_edge_whitelist = global_graph_edge_whitelist %>% filter(to %in% row.names(cell_type_ccm@ccs) &
                                                                    from %in% row.names(cell_type_ccm@ccs))

# Learn a single graph on all states at once (using the subgraphs as a whitelist/prior)
global_cell_type_wt_graph = assemble_wt_graph(combined_cds,
                                              cell_type_ccm,
                                    sample_group = "embryo",
                                    cell_group = "cell_type_sub",
                                    main_model_formula_str = NULL,
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="timepoint",
                                    #nuisance_model_formula_str = "~expt",
                                    ctrl_ids = control_genotypes,
                                    sparsity_factor = 0.01,
                                    perturbation_col = "gene_target",
                                    edge_whitelist =  global_graph_edge_whitelist,
                                    component_col="pcor_cluster",
                                    verbose=TRUE)

pdf(paste(outdir, "global_wt_graph_cell_type_sub_cross_partion_assembly.pdf", sep="/"), width=40, height=10)
plot_state_graph_annotations(cell_type_ccm@ccs, global_cell_type_wt_graph,
                             label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="cell_type_broad",
                             edge_weights = "support",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=TRUE)
dev.off()

pdf(paste(outdir, "global_wt_graph_cell_type_sub_cross_partion_assembly_small.pdf", sep="/"), width=40, height=10)
plot_state_graph_annotations(cell_type_ccs, global_cell_type_wt_graph,
                             #label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="tissue",
                             edge_weights = "support",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = FALSE,
                             #con_colour = "black",
                             label_groups=TRUE)
dev.off()




# Fit cell count models to each mutant vs control
cell_type_perturb_models_tbl = fit_mt_models(combined_cds,
                                             sample_group = "embryo",
                                             cell_group = "cell_type_sub",
                                             main_model_formula_str = NULL,
                                             start_time = assembly_start_time,
                                             stop_time = assembly_stop_time,
                                             interval_col="timepoint",
                                             num_time_breaks=3,
                                             #nuisance_model_formula_str = "~expt",
                                             ctrl_ids = control_genotypes,
                                             #mt_ids = mt_genotypes, #selected_mt_ids,
                                             sparsity_factor = 0.01,
                                             perturbation_col = "gene_target",
                                             edge_whitelist = igraph::as_data_frame(cell_type_graph_mt),
                                             keep_cds=TRUE,
                                             num_threads=num_threads,
                                             backend=assembly_backend,
                                             vhat_method="bootstrap",
                                             penalize_by_distance=TRUE,
                                             num_bootstraps=10)

cell_type_perturb_models_tbl = hooke:::assess_perturbation_effects(cell_type_ccm,
                                                                   cell_type_perturb_models_tbl,
                                                                   q_val=0.1,
                                                                   start_time = assembly_start_time,
                                                                   stop_time = assembly_stop_time,
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


# Learn a single graph on all states at once (using the subgraphs as a whitelist/prior)
global_cell_type_mt_graph = assemble_mt_graph(cell_type_ccm,
                                              cell_type_perturb_models_tbl,# %>% filter(perturb_name %in% c("tbx16-msgn1")),
                                    start_time = assembly_start_time,
                                    stop_time = assembly_stop_time,
                                    interval_col="timepoint",
                                    perturbation_col = "gene_target",
                                    break_cycles=FALSE,
                                    edge_whitelist = global_graph_edge_whitelist,
                                    q_val=0.1,
                                    verbose=TRUE)
global_cell_type_mt_graph = hooke:::break_cycles_in_state_transition_graph(global_cell_type_mt_graph, "total_perturb_path_score_supporting")


pdf(paste(outdir, "global_mt_graph_cell_type_sub.pdf", sep="/"), width=40, height=10)
plot_state_graph_annotations(cell_type_ccm@ccs, global_cell_type_mt_graph,
                             label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="tissue",
                             #group_nodes_by="cell_type_broad",
                             edge_weights = "total_perturb_path_score_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=FALSE)
dev.off()

pdf(paste(outdir, "global_mt_graph_cell_type_sub_small.pdf", sep="/"), width=6, height=6)
plot_state_graph_annotations(cell_type_ccm@ccs, cell_type_graph,
                             #label_nodes_by="cell_type_sub",
                             #color_nodes_by = "timepoint",
                             group_nodes_by="tissue",
                             edge_weights = "total_perturb_path_score_supporting",
                             #label_edges_by="support_label",
                             hide_unlinked_nodes = TRUE,
                             #con_colour = "black",
                             label_groups=FALSE)
dev.off()

saveRDS(cell_type_graph, paste(outdir,  "cell_type_graph.RDS", sep="/"))
saveRDS(cell_type_ccm, paste(outdir, "cell_type_ccm.RDS", sep="/"))
saveRDS(notochord_cds, paste(outdir,  "notochord_cds.RDS", sep="/"))


cell_type_group_ccs = new_cell_count_set(combined_cds, "embryo", "cell_type_sub")
pb_cds = hooke:::pseudobulk_ccs_for_states(ccs)

gene_patterns_over_cell_type_graph = classify_genes_over_graph(cell_type_group_ccs, cell_type_graph, abs_expr_thresh=1e-2, log_fc_thresh=2, cores=8)
marker_scores_over_cell_type_graph = top_markers(pb_cds, group_cells_by = "cell_type_sub", cores=8)

