# `compare_genes_in_cell_state`

Compares gene expression levels in a given cell state with its parent, sibling, and child states, and interprets the resulting patterns.

```r
compare_genes_in_cell_state(
  cell_state,
  state_graph,
  estimate_matrix,
  stderr_matrix,
  n,
  state_term = "cell_group",
  log_fc_thresh = 1,
  abs_expr_thresh = 0.001,
  sig_thresh = 0.05,
  cores = 1,
  expected_effect_mode_interval = c(-10, 10),
  cv_threshold = NULL
)
```

## Arguments

- **cell_state**  
  *character*  
  The cell state to be analyzed.

- **state_graph**  
  *igraph*  
  The state graph.

- **estimate_matrix**  
  *matrix*  
  Per-cell-state coefficient estimates (genes x cell states).

- **stderr_matrix**  
  *matrix*  
  Standard errors for `estimate_matrix`, same shape.

- **n**  
  *numeric*  
  Number of pseudobulk samples the models were fit on.

- **state_term**  
  *character*  
  Column identifying cell states. Default is `"cell_group"`.

- **log_fc_thresh**  
  *numeric*  
  Log fold change threshold. Default is `1`.

- **abs_expr_thresh**  
  *numeric*  
  Absolute expression threshold. Default is `1e-3`.

- **sig_thresh**  
  *numeric*  
  Significance threshold. Default is `0.05`.

- **cores**  
  *numeric*  
  Number of cores for parallel processing. Default is `1`.

- **expected_effect_mode_interval**  
  *numeric vector*  
  Expected effect mode interval. Default is `c(-10, 10)`.

- **cv_threshold**  
  *numeric*  
  Coefficient of variation threshold. Default is `NULL`.

## Value

A tibble containing gene expression comparisons and interpretations.

## Details

Identifies the parent, sibling, and child states of `cell_state` in `state_graph`, computes gene expression levels and significance for all of them, and compares the target state against each neighbor to classify each gene's expression pattern.

## Examples

```r
condensate_genes <- compare_genes_in_cell_state(
  cell_state = "pectoral fin condensate",
  state_graph = pf_cell_state_graph@graph,
  estimate_matrix = pb_coeffs$coefficients,
  stderr_matrix = pb_coeffs$stdev.unscaled,
  n = ncol(pb_cds)
)
```
