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


setwd("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/")
source("R/assembly_utils.R")

# outdir = "/Users/coletrap/dropbox_lab/Analysis/fish-mutants/gap-notebook-ct-1/"
outdir = "~/OneDrive/UW/Trapnell/hooke_manuscript/main_figures/figure3_components/"

RhpcBLASctl::omp_set_num_threads(8)

num_threads = as.numeric(Sys.getenv("OMP_NUM_THREADS"))
print (paste("running assembly with", num_threads, "threads"))

Sys.setenv("OMP_NUM_THREADS" = 1)
Sys.setenv("OPENBLAS_NUM_THREADS" = 1)

RhpcBLASctl::omp_set_num_threads(1)
RhpcBLASctl::blas_set_num_threads(1)

# combined_cds = readRDS("/Users/coletrap/dropbox_lab/Analysis/fish-mutants/hooke-analysis/combined_ref_gap16_doubletfiltered_labelupdate202304_cds.rds")
combined_cds = readRDS("~/OneDrive/UW/Trapnell/zebrafish-atlas-assembly/R_objects/combined_ref_gap16_doubletfiltered_labelupdate202304_cds.rds")
colData(combined_cds)$timepoint = as.numeric(colData(combined_cds)$timepoint)


plot_cells_3d(combined_cds[, sample(ncol(combined_cds), 100000)], color_cells_by = "cell_type_sub")


ccs = new_cell_count_set(combined_cds,
                         sample_group = "embryo",
                         cell_group = "cell_type_sub", keep_cds = T)

cell_type_metadata =
  hooke:::collect_psg_node_metadata(ccs,
                                    color_nodes_by="germ_layer",
                                    label_nodes_by="cell_type_broad",
                                    group_nodes_by="tissue") %>%
  dplyr::rename(germ_layer=color_nodes_by,
                cell_type_broad=label_nodes_by,
                tissue=group_nodes_by)

periderm_groups = cell_type_metadata %>% filter(cell_type_broad == "periderm") %>% pull(id)
periderm_counts = Matrix::colSums(normalized_counts(ccs, norm_method="size_only", pseudocount=0)[periderm_groups,])
periderm_props = periderm_counts / Matrix::colSums(normalized_counts(ccs, norm_method="size_only", pseudocount=0))

combined_cds = combined_cds[,!is.na(colData(combined_cds)$cell_type_broad) &
                             !colData(combined_cds)$expt == "GAP13" &
                             !colData(combined_cds)$cell_type_broad == "periderm"]

#
# ccs = new_cell_count_set(combined_cds,
#                          sample_group = "embryo",
#                          cell_group = "cell_type_sub", keep_cds = T)
#
# colData(ccs)$periderm_prop = periderm_props[row.names(colData(ccs))]

control_genotypes = unique(colData(combined_cds)[["gene_target"]])
control_genotypes = control_genotypes[grepl("wt|ctrl", control_genotypes)]

# time in nuisance model, so residuals mostly describe genetic perturbs
# ccm_controls_all_by_time = new_cell_count_model(ccs[,colData(ccs)$gene_target %in% control_genotypes],
#                                 nuisance_model_formula_str = "~ns(timepoint, df = 3) + expt + periderm_prop",
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


system.time({cell_type_ccm = fit_wt_model(combined_cds,
                                          sample_group = "embryo",
                                          cell_group = "cell_type_sub",
                                          #main_model_formula_str = NULL,
                                          #main_model_formula_str = "~1",
                                          #main_model_formula_str = "~ns(timepoint, df=5)",
                                          start_time = 18,
                                          stop_time = 72,
                                          interval_col="timepoint",
                                          vhat_method="bootstrap",
                                          num_time_breaks=4,
                                          #nuisance_model_formula_str = "~expt",
                                          ctrl_ids = control_genotypes,
                                          sparsity_factor = 0.01,
                                          perturbation_col = "gene_target",
                                          #base_penalty=1000,
                                          keep_cds = TRUE,
                                          num_threads=8,
                                          backend="nlopt",
                                          verbose=TRUE,
                                          penalize_by_distance=TRUE,
                                          pln_num_penalties=30)})

cell_type_metadata =
  hooke:::collect_psg_node_metadata(cell_type_ccm@ccs,
                                    color_nodes_by="germ_layer",
                                    label_nodes_by="cell_type_broad",
                                    group_nodes_by="tissue") %>%
  dplyr::rename(germ_layer=color_nodes_by,
                cell_type_broad=label_nodes_by,
                tissue=group_nodes_by)


#selected_cell_types = cell_type_metadata %>% filter(tissue == "Kidney") %>% pull(id)

selected_cell_types = c(
  "myoblast",
  "fast-committed myocyte (pre-fusion)",
  "fast-committed myocyte (fusing)",
  #"cranial muscle (non-somitic, fast-twitch)",
  #"anterior migratory muscle",
  #"mature slow muscle 3",
  #"satellite cells",
  #"cranial muscle (late)",
  "mature fast muscle 1",
  "mature fast muscle 2",
  "mature fast muscle 3",
  "mature fast muscle 4",
  #"mature slow muscle 2",
  #"mature slow muscle 1",
  "mature fast muscle 5",
  "mature fast muscle 6"#,
  #"slow-committed myocyte",
  #"unknown (dcn+, col6+)"
)

colData(cell_type_ccm@ccs)$total_cells = Matrix::colSums(counts(cell_type_ccm@ccs))

ggplot(aes(expt, total_cells, fill=expt), data=colData(cell_type_ccm@ccs) %>% as.data.frame) + geom_boxplot()

wt_kinetic_plot = plot_cell_type_control_kinetics(cell_type_ccm,
                                cell_groups=selected_cell_types,
                                log_abund_detection_thresh=log(1),
                                batch_col="expt") + monocle3:::monocle_theme_opts()

pdf(paste(outdir, "muscle_wt_kinetics.pdf", sep="/"), width=4, height=12)
print(wt_kinetic_plot)
dev.off()

rm(cell_type_ccm)

# Fit cell count models to each mutant vs control
cell_type_perturb_models_tbl = fit_mt_models(combined_cds,
                                             sample_group = "embryo",
                                             cell_group = "cell_type_sub",
                                             main_model_formula_str = NULL,
                                             start_time = 18,
                                             stop_time = 72,
                                             interval_col="timepoint",
                                             num_time_breaks=3,
                                             #nuisance_model_formula_str = "~expt",
                                             ctrl_ids = "ctrl-inj",
                                             mt_ids = c("tbx16", "tbx16-msgn1"), #selected_mt_ids,
                                             sparsity_factor = 0.01,
                                             perturbation_col = "gene_target",
                                             keep_cds=TRUE,
                                             num_threads=8,
                                             backend="nlopt",
                                             vhat_method="bootstrap",
                                             penalize_by_distance=TRUE,
                                             num_bootstraps=10)


selected_ccm = cell_type_perturb_models_tbl %>% filter(perturb_name == "tbx16-msgn1")
selected_ccm = selected_ccm$perturb_ccm[[1]]
kinetics_plot = plot_cell_type_perturb_kinetics(selected_ccm,
                                                cell_groups = selected_cell_types,
                                                interval_col="timepoint",
                                                start_time=18,
                                                stop_time=36,
                                                expt="GAP16",
                                                q_val = 0.1) +
  monocle3:::monocle_theme_opts() #+ scale_y_continuous()


pdf(paste(outdir, "muscle_tbx16_msgn1_cell_type_distribution_over_time.pdf", sep="/"), height=10, width=4)
print(kinetics_plot)
#ggplot2::ggsave(paste(outdir, "mt_graph_peak_abund_times_all.pdf", sep="/"), width=40, height=15)
dev.off()







