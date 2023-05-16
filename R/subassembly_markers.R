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

#plan(multicore)
options(future.fork.multithreading.enable=FALSE)
options(future.globals.maxSize = 2* 8192 * 1024^2) # 16GB

assembly_start_time = 18
assembly_stop_time = 72
assembly_backend = "nlopt"


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
outdir = "/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/"
setwd("/Users/coletrap/Google Drive/My Drive/develop/zebrafish-atlas-assembly")
source("R/assembly_utils.R")

all_partition_results = readRDS(paste(outdir, "all_partition_results.rds", sep="/"))

combined_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/hooke-analysis/combined_ref_gap16_doubletfiltered_labelupdate202304_cds.rds")
colData(combined_cds)$timepoint = as.numeric(colData(combined_cds)$timepoint)

# collect the assignments to subassembly groups and write them back to the main CDS:
cell_state_assignments = all_partition_results %>% select(data) %>% tidyr::unnest(data)
cell_state_assignments$subassembly_group = as.character(cell_state_assignments$subassembly_group)

# Downsample for now because even analyzing pseudobulks takes a lot of memory
#cell_downsample = cell_state_assignments %>% group_by(subassembly_group) %>% dplyr::slice_sample(n=2000, replace=TRUE)
#cell_downsample = cell_downsample %>% distinct() %>% as.data.frame(stringsAsFactors=FALSE)
#row.names(cell_downsample) = cell_downsample$cds_row_id

# Do this if you want to use all the cells:
cell_downsample = cell_state_assignments
cell_downsample = cell_downsample %>% distinct() %>% as.data.frame(stringsAsFactors=FALSE)
row.names(cell_downsample) = cell_downsample$cds_row_id

combined_cds = combined_cds[,cell_downsample$cds_row_id]

colData(combined_cds)$subassembly_group = cell_downsample[row.names(colData(combined_cds) %>% as.data.frame()),]$subassembly_group
#colData(combined_cds)$pcor_cluster = stringr::str_split_fixed(colData(combined_cds)$subassembly_group, "-", 2)[,1]


# Create a count set, grouping cells by subassembly_group
subassembly_group_ccs = new_cell_count_set(combined_cds, "embryo", "subassembly_group")


# Combine all the subassembly graphs into a single graph object
mt_subassembly_graphs = all_partition_results %>%
  #filter(is.na(mt_graph) == FALSE ) %>%
  pull(mt_graph)

mt_subassembly_graphs = lapply(mt_subassembly_graphs, function(x) {
  if(is.null(x) == FALSE && is.na(x) == FALSE) return(igraph::as_data_frame(x))})
subassembly_graph_mt_union = do.call(bind_rows, mt_subassembly_graphs)
subassembly_graph_mt_union = igraph::graph_from_data_frame(subassembly_graph_mt_union, vertices = unique(colData(combined_cds)$subassembly_group))

pb_cds = hooke:::pseudobulk_ccs_for_states(subassembly_group_ccs)
#marker_scores_over_subassembly_graph = top_markers(pb_cds, group_cells_by = "cell_group", cores=8)
#saveRDS(marker_scores_over_subassembly_graph, paste(outdir, "marker_scores_over_subassembly_graph.rds", sep="/"))
marker_scores_over_subassembly_graph = readRDS(paste(outdir, "marker_scores_over_subassembly_graph.rds", sep="/"))

subassembly_group_markers = marker_scores_over_subassembly_graph %>%
  filter(pseudo_R2 > 0.25) %>%
  select(cell_group, gene_short_name) %>%
  dplyr::rename(id=cell_group, marker_name=gene_short_name)
subassembly_group_markers$marker_type = "Specifically activated"

# pdf(paste(outdir, "subassembly_state_graph_top_markers_summary.pdf", sep="/"), width=40, height=20)
# plot_state_graph_marker_genes(subassembly_group_ccs,
#                               subassembly_graph_mt_union,
#                               subassembly_group_markers,
#                               label_nodes_by="subassembly_group",
#                               #color_nodes_by = "tissue",
#                               group_nodes_by="pcor_cluster",
#                               #edge_weights = "total_perturb_path_score_supporting",
#                               num_top_genes = 6,
#                               #label_edges_by="support_label",
#                               label_groups=TRUE,
#                               hide_unlinked_nodes = TRUE)
# dev.off()


subassemblies = unique(all_partition_results$partition)
plot_subassembly_markers = function(subassembly){


  subassembly_graph = all_partition_results %>% filter(partition == subassembly) %>% pull(mt_graph)
  subassembly_graph = subassembly_graph[[1]]

  if (is.null(subassembly_graph) == FALSE && is.logical(subassembly_graph) == FALSE){

    # marker_state_graph_plot_name = paste("subassembly_", subassembly, "_state_graph_marker_plot.pdf", sep="")
    # p = plot_state_graph_marker_genes(subassembly_group_ccs,
    #                               subassembly_graph,
    #                               subassembly_group_markers,
    #                               label_nodes_by="subassembly_group",
    #                               #color_nodes_by = "tissue",
    #                               group_nodes_by="cell_type_sub",
    #                               #edge_weights = "total_perturb_path_score_supporting",
    #                               num_top_genes = 6,
    #                               #label_edges_by="support_label",
    #                               label_groups=TRUE,
    #                               hide_unlinked_nodes = TRUE)
    # pdf(paste(outdir, marker_state_graph_plot_name, sep="/"), width=10, height=10)
    # print (p)
    # dev.off()

    subassembly_perturb_effects = all_partition_results %>% filter(partition == subassembly) %>% pull(perturbation_effects)
    subassembly_perturb_effects = subassembly_perturb_effects[[1]]

    observed_genetic_requirements = categorize_genetic_requirements(subassembly_perturb_effects, subassembly_graph)
    observed_genetic_requirements = observed_genetic_requirements #%>%
    #mutate(perturb_name = stringr::str_replace(perturb_name, "-mut", "")) %>% distinct()
    genetic_requirements = observed_genetic_requirements


    marker_state_graph_plot_name = paste("subassembly_", subassembly, "_state_graph_requirements_plot.pdf", sep="")
    p = plot_state_graph_perturb_effects(subassembly_group_ccs,
                                      subassembly_graph,
                                      genetic_requirements,
                                      label_nodes_by="subassembly_group",
                                      #color_nodes_by = "tissue",
                                      group_nodes_by="cell_type_sub",
                                      edge_weights = "total_perturb_path_score_supporting",
                                      num_top_genes = 6,
                                      #label_edges_by="support_label",
                                      label_groups=TRUE,
                                      hide_unlinked_nodes = TRUE)
    pdf(paste(outdir, marker_state_graph_plot_name, sep="/"), width=20, height=10)
    print (p)
    dev.off()
  } else {
    message(paste("No graph for subassembly", subassembly))
  }

  # subassembly_groups = all_partition_results %>% filter(partition == subassembly) %>% tidyr::unnest(data) %>%
  #   pull(subassembly_group) %>% unique
  # if (length(subassembly_groups) > 1){
  #   subassembly_markers = subassembly_group_markers %>% filter(id %in% subassembly_groups) %>% pull(marker_name) %>% unique
  #   marker_dot_plot_name = paste("subassembly_", subassembly, "_dot_plot.pdf", sep="")
  #   subgroup_pb_cds = pb_cds[,colData(pb_cds)$cell_group %in% subassembly_groups]
  #   p = plot_genes_by_group(subgroup_pb_cds,
  #                           group_cells_by="cell_group",
  #                           markers=subassembly_markers)
  #   pdf(paste(outdir, marker_dot_plot_name, sep="/"), width=10, height=10)
  #   print (p)
  #   dev.off()
  # }

}
#debug(plot_subassembly_markers)
lapply(subassemblies, plot_subassembly_markers)

#### Now just go through and fix up a bunch of annotations:

colData(combined_cds)$updated_cell_type = colData(combined_cds)$cell_type_sub

source("R/update_subassembly_annotations.R")
combined_cds = udpate_annotations(combined_cds)

colData(combined_cds)$cell_type_sub = colData(combined_cds)$updated_cell_type

saveRDS(colData(combined_cds), paste(outdir, "updated_coldata.rds", sep="/"))

xxx_cds = combined_cds[,sample(ncol(combined_cds), 100000)]

# Exploring specific sets of cells:

#################################
# EAR:
selected_subassembly = c("36", "44")

selected_subassembly_cds = combined_cds[,colData(combined_cds)$pcor_cluster %in% selected_subassembly]
selected_subassembly_cds = selected_subassembly_cds %>%
  preprocess_cds(num_dim = 50) %>%
  align_cds(residual_model_formula_str = "~log.n.umi") %>%
  reduce_dimension(max_components = 3)
plot_cells_3d(selected_subassembly_cds, color_cells_by="updated_cell_type")

# Bunch of important ear genes from ZFIN:
plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers = c("bmp2b", "bmp4", "bmp7a", "bmp7b",
                                "cdh1", "cdh2", "cdh11",
                                "col2a1",
                                "dla", "dlb", "dld", "jag2b",
                                "dlx2", "dlx3", "dlx4", "dlx6", "dlx7",
                                "emx2",
                                "eya1",
                                "fgf8a", "fgf8b",
                                "hsp47",
                                "anos1a", "anos1b",
                                "msx3", "msx2b", "msx1a",
                                "otx1",
                                "pax2a", "pax2b", "pax5", "pax8",
                                "prox1a",
                                "epha4a", ### NOTE!!! <- we have a crispant for this. Check losses in support cells
                                "snai1b",
                                "tbx2b",
                                "wnt4a"))

shi_raible_elife_markers = c("matn4",
                             "col2a1a",
                             "fgfr2",
                             "ifng",
                             "myo6b",
                             "dla",
                             "otogl",
                             "tectb",
                             "zpld1a",
                             "pard3bb",
                             "neo1a",
                             "fat1a",
                             "cx30.3",
                             "chrna7",
                             "slc1a2b",
                             "zlpd1b",
                             "glis1b",
                             "slc25a22b",
                             "klhl30",
                             "espnla",
                             "dld",
                             "cabp1b",
                             "tmc1",
                             "ush1c",
                             "rtn4r",
                             "rem1",
                             "dpp6b",
                             "cabp2b",
                             "pvalb9",
                             "dlx3b", "dlx7b", "eya1", "six4b", # early otic placode
                             "pax8",
                             "pax5",
                             "pax2a",
                             "pax2b",
                             "cldna",
                             "cldnb",
                             "dacha",
                             "dachb",
                             "dachc",
                             "msx3",
                             "msx2b",
                             "bmp2b",
                             "bmp4",
                             "prox1a",
                             "bmp2a",
                             "odf3l2b",
                             "adam11",
                             "glis1b",
                             "gpc3"
)

# Markers from Raible recent eLife paper:
plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers=shi_raible_elife_markers)


pdf(paste(outdir, "all_cells_otic_vesicle_progenitors.pdf", sep="/"), width=150, height=4)
p= plot_genes_by_group(pb_cds,
                    group_cells_by="cell_group",
                    markers=c("pard3bb", "igsf3", "fgfr2", "neo1a", "fat1a"))
print (p)
dev.off()


#updated_coldata = colData(combined_cds)



plot_cells_3d(selected_subassembly_cds, color_cells_by="updated_cell_type")

selected_subassembly_pb_cds = pb_cds[,colData(pb_cds)$cell_group %in% colData(selected_subassembly_cds)$subassembly_group]


pb_cds = hooke:::pseudobulk_ccs_for_states(subassembly_group_ccs)

marker_scores_over_selected_subassembly = top_markers(selected_subassembly_pb_cds, group_cells_by = "cell_group", cores=8)

selected_subassembly_group_markers = marker_scores_over_selected_subassembly %>%
  filter(pseudo_R2 > 0.25) %>%
  select(cell_group, gene_short_name) %>%
  dplyr::rename(id=cell_group, marker_name=gene_short_name)

selected_subassembly_markers = selected_subassembly_group_markers %>% pull(marker_name) %>% unique

pdf(paste(outdir, "selected_subassembly_top_markers.pdf", sep="/"), width=10, height=20)
p= plot_genes_by_group(selected_subassembly_pb_cds,
                        group_cells_by="cell_group",
                        markers=selected_subassembly_markers)
print(p)
dev.off()


###############

#################################
# Cranial muscle, pec fin muscle?
selected_subassembly = c("17")

selected_subassembly_cds = combined_cds[,colData(combined_cds)$pcor_cluster %in% selected_subassembly]
bogus_cell_types = colData(selected_subassembly_cds) %>% as.data.frame %>% pull(cell_type_sub) %>% table
bogus_cell_types = names(bogus_cell_types)[which(bogus_cell_types < 100)]
colData(selected_subassembly_cds)$cell_type_sub[colData(selected_subassembly_cds)$cell_type_sub %in% bogus_cell_types] = NA
selected_subassembly_cds = selected_subassembly_cds %>%
  preprocess_cds(num_dim = 50) %>%
  align_cds(residual_model_formula_str = "~log.n.umi") %>%
  reduce_dimension(max_components = 3)
plot_cells_3d(selected_subassembly_cds, color_cells_by="cell_type_sub")

# Genes that L&S used to annotated originally:
LS_markers = c("msc", "gfra3", "si:ch211-263k4.2", # cranial muscle early
                 "tcf21", "msc", # cranial muscle (progenitor)
                 "erbb4a", # cranial muscle mid
                 "myhb", "fhl2b", "tnni2a.2") # cranial muscle late)


# Genes that regulate dorsal, intermediate, ventral PA patterning:
DIV_markers = c("eya1", "six1a", "six1b", "nr2fs",  "jag1a", "jag1b", "hey1", "grem2b", "pou3f3a", "pou3f3b", # Dorsal
                "bapx1", "gdf5", # intermediate
                "mef2c", "edn1", "dlx3b", "dlx4a", "dlx4b", "dlx5a", "dlx6a", "bmp4", "hand2", "msx1a", "msx1b") #ventral


plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers = DIV_markers)

# Create a count set, grouping cells by subassembly_group
subassembly_group_ccs = new_cell_count_set(selected_subassembly_cds, "embryo", "subassembly_group")
selected_subassembly_pb_cds = hooke:::pseudobulk_ccs_for_states(subassembly_group_ccs)

marker_scores_over_selected_subassembly = top_markers(selected_subassembly_pb_cds, group_cells_by = "cell_group", cores=8)

selected_subassembly_group_markers = marker_scores_over_selected_subassembly %>%
  #filter(pseudo_R2 > 0.25) %>%
  arrange(cell_group, desc(marker_score)) %>%
  group_by(cell_group) %>%
  slice_head(n=5) %>%
  select(cell_group, gene_short_name) %>%
  dplyr::rename(id=cell_group, marker_name=gene_short_name)

selected_subassembly_markers = selected_subassembly_group_markers %>% pull(marker_name) %>% unique


p= plot_genes_by_group(selected_subassembly_pb_cds,
                       group_cells_by="cell_group",
                       markers=selected_subassembly_markers)
pdf(paste(outdir, "cranial_muscle_subassembly_group_top_markers.pdf", sep="/"), width=6, height=15)
print (p)
dev.off()


#################################
# Neural crest + Pharyngeal arches:

# Note possibly the NC derived neurons are mixed into the neural progenitor partition, there are some foxd3+ cells in there
selected_subassembly = c("6",
                         "18", "13", # Cranial NC + pigment
                         "17", # cranial muscle
                         "27", # schwann cells (and retina)
                         "22", # vascular endothelium/aorta
                         "40", # chondrocranium
                         "31", # jaw chondrocytes
                         "23" # cardiomyocytes
                         )

#colData(selected_subassembly_cds) %>% as.data.frame %>% filter(grepl("pharyngeal muscle (contains muscle, early", cell_type_sub)) %>% pull(subassembly_group) %>% table

selected_subassembly_cds = combined_cds[,colData(combined_cds)$pcor_cluster %in% selected_subassembly]
bogus_cell_types = colData(selected_subassembly_cds) %>% as.data.frame %>% pull(cell_type_sub) %>% table
bogus_cell_types = names(bogus_cell_types)[which(bogus_cell_types < 100)]
colData(selected_subassembly_cds)$cell_type_sub[colData(selected_subassembly_cds)$cell_type_sub %in% bogus_cell_types] = NA
selected_subassembly_cds = selected_subassembly_cds %>%
  preprocess_cds(num_dim = 50) %>%
  align_cds(residual_model_formula_str = "~log.n.umi") %>%
  reduce_dimension(max_components = 3, cores=8)
plot_cells_3d(selected_subassembly_cds, color_cells_by="cell_type_sub")

selected_subassembly_cds = cluster_cells(selected_subassembly_cds)


p = plot_cells(selected_subassembly_cds, color_cells_by="cluster")
png(paste(outdir, "pharyngeal_arch_subassembly_group_by_cluster_UMAP.png", sep="/"), width=6, height=6, units="in", res=300)
print (p)
dev.off()


p = plot_cells(selected_subassembly_cds, color_cells_by="updated_cell_type")
png(paste(outdir, "pharyngeal_arch_subassembly_group_by_type_UMAP.png", sep="/"), width=6, height=6, units="in", res=300)
print (p)
dev.off()

p = plot_cells(selected_subassembly_cds, color_cells_by="subassembly_group")
png(paste(outdir, "pharyngeal_arch_subassembly_group_by_subassembly_group_UMAP.png", sep="/"), width=6, height=6, units="in", res=300)
print (p)
dev.off()

p = plot_cells(selected_subassembly_cds, color_cells_by="timepoint")
png(paste(outdir, "pharyngeal_arch_subassembly_group_by_timepoint_UMAP.png", sep="/"), width=6, height=6, units="in", res=300)
print (p)
dev.off()


NC_markers = c("acana",  "matn1", "col10a1a", # jaw chondrocyte
               #"myl1", "myhz2", "myom2a", "obscnb", "tcf21", "msc", "myod1", # muscle
               "lmx1bb", "pitx2", "eya2", # periocular mesenchyme
               "mki67",
               "hand2",
               "prdm1a", "grem2b", "dlx2a", "fli1a", # posterior arch NC
               "tfap2a", "crestin", "sox9a", # general NC
               "snai1b", "prickle1a", "wnt11", # migrating NC
               "mitfa", "gpr143", "trpm1b", # pigment
               "sox10", "ntm", "foxd3", "cdh1", "cdh2", "gfra4a", "znf536") #cranial neural crest


# these are from https://www.science.org/doi/10.1126/science.aas9536
mouse_trunk_NC_markers = c("zic1", "zic3", "msx1a", "msx1b", "gdf7", "wnt3a", "mafba", # neural-tube and pre-migratory crest
                           "cdh11", "arhgef7a", "arhgef7b", "sox9a", "sox9b", "dlx5a", "dlx5b", # delaminating NC
                           "smtnl2", "sfrp5", "heyl", "nkain4",
                           "neurog1", "neurog3", "neurog4", "car11", # sensory neuron fated
                           "aslc1a", "aslc1b", "itga4", "phox2b", "mcama", "mcamb", # autonomic neuron fated
                           "egflam", "mstna", "mstnb", #  glia fated (both sensory and autonomic)
                           "meox1", "foxc1", "prrx1", "twist1a", "fli1", "six2a", "six2b" # mesenchyme fated
)


plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers = NC_markers)



p = plot_genes_by_group(combined_cds,
                    group_cells_by="subassembly_group",
                    markers = NC_markers) + coord_flip()

pdf(paste(outdir, "all_subassembly_ncc_markers_dotplot.pdf", sep="/"), width=10, height=75)
print (p)
dev.off()

p = plot_cells(selected_subassembly_cds, genes=NC_markers)
png(paste(outdir, "pharyngeal_arch_subassembly_group_NC_genes_UMAP.png", sep="/"), width=6, height=6, units="in", res=300)
print (p)
dev.off()

# Genes that regulate dorsal, intermediate, ventral PA patterning:
DIV_markers = c("eya1", "six1a", "six1b", "nr2fs",  "jag1a", "jag1b", "hey1", "grem2b", "pou3f3a", "pou3f3b", # Dorsal
                "bapx1", "gdf5", # intermediate
                "mef2c", "edn1", "dlx3b", "dlx4a", "dlx4b", "dlx5a", "dlx6a", "bmp4", "hand2", "msx1a", "msx1b") #ventral

p= plot_cells(selected_subassembly_cds,
               genes=DIV_markers)
png(paste(outdir, "pharyngeal_arch_subassembly_group_DIV_markers_UMAP.png", sep="/"), width=15, height=15, units="in", res=300)
print (p)
dev.off()

plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers = DIV_markers)

hox_genes = rowData(selected_subassembly_cds)$gene_short_name[grepl("^hox", rowData(selected_subassembly_cds)$gene_short_name)]
pbx_genes = rowData(selected_subassembly_cds)$gene_short_name[grepl("^pbx", rowData(selected_subassembly_cds)$gene_short_name)]
dlx_genes = rowData(selected_subassembly_cds)$gene_short_name[grepl("^dlx", rowData(selected_subassembly_cds)$gene_short_name)]
otx_genes = rowData(selected_subassembly_cds)$gene_short_name[grepl("^otx", rowData(selected_subassembly_cds)$gene_short_name)]

plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers = c(hox_genes, pbx_genes, dlx_genes, otx_genes))

plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers = c(hox_genes))


p= plot_cells(selected_subassembly_cds,
              genes=c(hox_genes, pbx_genes, dlx_genes, otx_genes))
png(paste(outdir, "pharyngeal_arch_subassembly_group_hox_pbx_otx_pbx_UMAP.png", sep="/"), width=25, height=25, units="in", res=300)
print (p)
dev.off()


# Create a count set, grouping cells by subassembly_group
subassembly_group_ccs = new_cell_count_set(selected_subassembly_cds, "embryo", "subassembly_group")
selected_subassembly_pb_cds = hooke:::pseudobulk_ccs_for_states(subassembly_group_ccs)

marker_scores_over_selected_subassembly = top_markers(selected_subassembly_pb_cds, group_cells_by = "cell_group", cores=8)

selected_subassembly_group_markers = marker_scores_over_selected_subassembly %>%
  #filter(pseudo_R2 > 0.25) %>%
  arrange(cell_group, desc(marker_score)) %>%
  group_by(cell_group) %>%
  slice_head(n=5) %>%
  select(cell_group, gene_short_name) %>%
  dplyr::rename(id=cell_group, marker_name=gene_short_name)

selected_subassembly_markers = selected_subassembly_group_markers %>% pull(marker_name) %>% unique


p= plot_genes_by_group(selected_subassembly_pb_cds,
                       group_cells_by="cell_group",
                       markers=selected_subassembly_markers)
pdf(paste(outdir, "cranial_muscle_subassembly_group_top_markers.pdf", sep="/"), width=6, height=15)
print (p)
dev.off()



#################################
# Roof plate:
selected_subassembly = c("42")

selected_subassembly_cds = combined_cds[,colData(combined_cds)$pcor_cluster %in% selected_subassembly]
selected_subassembly_cds = selected_subassembly_cds %>%
  preprocess_cds(num_dim = 50) %>%
  align_cds(residual_model_formula_str = "~log.n.umi") %>%
  reduce_dimension(max_components = 3)
plot_cells_3d(selected_subassembly_cds, color_cells_by="updated_cell_type")

# Bunch of important ear genes from ZFIN:
zfin_markers = c()

plot_genes_by_group(selected_subassembly_cds,
                    group_cells_by="subassembly_group",
                    markers = zfin_markers)


pdf(paste(outdir, "all_roof_plate_selectedmarkers.pdf", sep="/"), width=150, height=4)
p= plot_genes_by_group(pb_cds,
                       group_cells_by="cell_group",
                       markers=c("pard3bb", "igsf3", "fgfr2", "neo1a", "fat1a"))
print (p)
dev.off()


marker_scores_over_selected_subassembly = top_markers(selected_subassembly_pb_cds, group_cells_by = "cell_group", cores=8)

selected_subassembly_group_markers = marker_scores_over_selected_subassembly %>%
  filter(pseudo_R2 > 0.25) %>%
  select(cell_group, gene_short_name) %>%
  dplyr::rename(id=cell_group, marker_name=gene_short_name)

selected_subassembly_markers = selected_subassembly_group_markers %>% pull(marker_name) %>% unique

pdf(paste(outdir, "selected_subassembly_top_markers.pdf", sep="/"), width=10, height=20)
p= plot_genes_by_group(selected_subassembly_pb_cds,
                       group_cells_by="cell_group",
                       markers=selected_subassembly_markers)
print(p)
dev.off()


#updated_coldata = colData(combined_cds)



plot_cells_3d(selected_subassembly_cds, color_cells_by="updated_cell_type")

selected_subassembly_pb_cds = pb_cds[,colData(pb_cds)$cell_group %in% colData(selected_subassembly_cds)$subassembly_group]

marker_scores_over_selected_subassembly = top_markers(selected_subassembly_pb_cds, group_cells_by = "cell_group", cores=8)

selected_subassembly_group_markers = marker_scores_over_selected_subassembly %>%
  filter(pseudo_R2 > 0.25) %>%
  select(cell_group, gene_short_name) %>%
  dplyr::rename(id=cell_group, marker_name=gene_short_name)

selected_subassembly_markers = selected_subassembly_group_markers %>% pull(marker_name) %>% unique

pdf(paste(outdir, "selected_subassembly_top_markers.pdf", sep="/"), width=10, height=20)
p= plot_genes_by_group(selected_subassembly_pb_cds,
                       group_cells_by="cell_group",
                       markers=selected_subassembly_markers)
print(p)
dev.off()


###############

