Reference
=========

All functions
-------------------------------

### Kinetic functions

[`fit_genotype_ccm()`](fit_genotype_ccm)

Fits a cell count model (CCM) for a given genotype using a specified dataset of cell counts (CCS). 

[`fit_wt_model()`](fit_wt_model)
Fits a wild type (WT) model to a cell dataset (CDS) by estimating cell count dynamics over time and accounting for nuisance variables.

[`fit_mt_models()`](fit_mt_models)

[`plot_cell_type_control_kinetics()`](plot_cell_type_control_kinetics)
Generates a kinetic plot of cell type control data over time.

[`plot_cell_type_perturb_kinetics()`](plot_cell_type_perturb_kinetics)
Generates a plot to visualize the kinetics of cell type perturbations over time.

### Graph functions

[`new_cell_state_graph()`](new_cell_state_graph)
Creates a new cell state graph object from an input graph and a cell count set.

[`assemble_partition()`](assemble_partition)
Assembles a partition for a given cell dataset (CDS) by fitting wild-type (WT) and mutant (MT) models, constructing state transition graphs, and assessing perturbation effects.

### Plotting functions

[`plot_annotations()`](plot_annotations)
Generates a plot of cell state graphs with various customization options.

[`plot_gene_expr()`](plot_gene_expr)
Plots gene expression data on a cell state graph.

[`plot_abundance_changes()`](plot_abundance_changes)
Generates a plot to visualize changes in cell state abundances.

[`plot_degs()`](plot_degs)
Plots differentially expressed genes (DEGs) on a cell state graph, with customizable appearance and layout options.


