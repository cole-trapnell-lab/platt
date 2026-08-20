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

Fits models for analyzing multi-timepoint perturbation data using a cell dataset (`cds`).

[`plot_cell_type_control_kinetics()`](plot_cell_type_control_kinetics)

Generates a kinetic plot of cell type control data over time.

[`plot_cell_type_perturb_kinetics()`](plot_cell_type_perturb_kinetics)

Generates a plot to visualize the kinetics of cell type perturbations over time.

### Graph functions

[`new_cell_state_graph()`](new_cell_state_graph)

Creates a new cell state graph object from an input graph and a cell count set.

[`run_wildtype_assembly()`](run_wildtype_assembly)

Fits the wild-type model and assembles a state transition graph for a partition of the data.

[`contract_state_graph()`](contract_state_graph)

Simplifies a directed state transition graph by contracting nodes according to a specified grouping variable.

[`get_all_parents()`](get_all_parents)

Recursively retrieves all upstream ancestor parents for a given cell state in a state graph.

### DEG functions

[`compare_genes_over_graph()`](compare_genes_over_graph)

Compares gene expression over a given state transition graph, scoring genes according to their expression pattern across states.

[`compare_genes_within_state_graph()`](compare_genes_within_state_graph)

Compares gene expression within a state graph, contrasting a perturbation against control(s) within each cell state.

[`compare_genes_in_cell_state()`](compare_genes_in_cell_state)

Classifies a single cell state's gene expression pattern relative to its parents, children, and siblings in the state graph.

[`annotate_empirical_fdr()`](annotate_empirical_fdr)

Decorates a within-state DEG table with an empirical FDR for the low-expression, control-sampling artifact in control-heavy designs.

[`efdr_perturbation_support()`](efdr_perturbation_support)

Per-cell-type perturbation-arm replicate support, computed from a `cell_data_set`.

[`efdr_support_from_coldata()`](efdr_support_from_coldata)

Per-cell-type perturbation-arm support, from a plain coldata table.

[`build_empirical_null()`](build_empirical_null)

Builds a control-split empirical null for the DEG artifact, for training the full eFDR model.

[`train_efdr_model()`](train_efdr_model)

Trains the empirical-FDR null model (2-D tail surfaces over expression x thin-arm cell count) from control-split null draws.

[`annotate_empirical_fdr_model()`](annotate_empirical_fdr_model)

Annotates a DEG table with per-gene `empirical_p`/`empirical_fdr` using a trained `efdr_model` — the fuller, arm-size-conditioned version of `annotate_empirical_fdr()`.

### Phenotyping functions

[`categorize_genetic_requirements()`](categorize_genetic_requirements)

Classifies lost cell states as directly or indirectly required for a genotype, based on state graph topology.

[`impact_to_phenos()`](impact_to_phenos)

Translates a hand- or upstream-assembled impact table into the standardized fields the phenotype-plotting functions expect.

[`plot_phenotypes_from_impact()`](plot_phenotypes_from_impact)

Converts an impact table with `impact_to_phenos()` and renders the annotated phenotype calls on a cell state graph in one call.

[`plot_phenotypes_glyphs()`](plot_phenotypes_glyphs)

Draws phenotype calls (abundance fill + identity glyph) on a cell state graph from a table already run through `impact_to_phenos()`.

### Plotting functions

[`plot_annotations()`](plot_annotations)

Generates a plot of cell state graphs with various customization options.

[`plot_gene_expression()`](plot_gene_expression)

Plots gene expression data on a cell state graph.

[`plot_abundance_changes()`](plot_abundance_changes)

Generates a plot to visualize changes in cell state abundances.

[`plot_degs()`](plot_degs)

Plots differentially expressed genes (DEGs) on a cell state graph, with customizable appearance and layout options.


