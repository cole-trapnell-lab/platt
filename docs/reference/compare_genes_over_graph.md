# `compare_genes_over_graph`

Compares gene expression over a given state transition graph, scoring genes according to their expression pattern across the states (nodes) of the graph.

```r
compare_genes_over_graph(
  ccs,
  state_graph,
  gene_ids = NULL,
  group_nodes_by = NULL,
  assembly_group = NULL,
  label_nodes_by = "cell_state",
  states_to_assess = list(),
  nuisance_model_formula_str = "0",
  log_fc_thresh = 1,
  abs_expr_thresh = 0.001,
  sig_thresh = 0.05,
  min_samples_detected = 2,
  min_cells_per_pseudobulk = 3,
  cores = 1,
  cv_threshold = 100,
  ...
)
```

## Arguments

- **ccs**  
  *cell_count_set*  
  Cell count summary object.

- **state_graph**  
  *igraph*  
  A graph object representing the state graph.

- **gene_ids**  
  *character vector*  
  A vector of gene IDs to be assessed. Default is `NULL`.

- **group_nodes_by**  
  *character*  
  A column name to group nodes by. Default is `NULL`.

- **assembly_group**  
  *character*  
  A specific assembly group to assess. Default is `NULL`.

- **label_nodes_by**  
  *character*  
  A label for nodes. Default is `"cell_state"`.

- **states_to_assess**  
  *list*  
  A list of states to assess. Default is an empty list.

- **nuisance_model_formula_str**  
  *character*  
  A string representing the nuisance model formula. Default is `"0"`.

- **log_fc_thresh**  
  *numeric*  
  Log fold change threshold. Default is `1`.

- **abs_expr_thresh**  
  *numeric*  
  Absolute expression threshold. Default is `0.001`.

- **sig_thresh**  
  *numeric*  
  Significance threshold. Default is `0.05`.

- **min_samples_detected**  
  *numeric*  
  Minimum number of samples detected. Default is `2`.

- **min_cells_per_pseudobulk**  
  *numeric*  
  Minimum number of cells per pseudobulk. Default is `3`.

- **cores**  
  *numeric*  
  Number of cores to use for parallel processing. Default is `1`.

- **cv_threshold**  
  *numeric*  
  Coefficient of variation threshold. Default is `100`.

- **...**  
  Additional arguments.

## Value

A tibble containing gene class scores for each cell state.

## Details

This function compares gene expression over a given state graph, in contrast to `compare_genes_within_state_graph()` which compares gene expression between a perturbation and a control within each state.

## Examples

```r
pf_graph_degs <- compare_genes_over_graph(
  pf_ccs,
  state_graph,
  gene_ids = c("gene1", "gene2")
)
```
