# `plot_deg_heatmap`

Builds a clustered heatmap of log fold changes for a supplied set of genes, typically a set of common DEGs shared across perturbations or cell groups.

```r
plot_deg_heatmap(
  deg_table,
  genes,
  cell_groups = NULL,
  sig_p_value_thresh = 0.05,
  min_log_mean_expression = -2,
  gene_col = "id",
  gene_label_col = "gene_short_name",
  cell_group_col = "cell_group",
  lfc_col = "perturb_to_ctrl_shrunken_lfc",
  p_value_col = "perturb_to_ctrl_p_value",
  present_col = "present_above_thresh",
  log_mean_expression_col = "log_mean_expression"
)
```

## Arguments

- **deg_table**  
  *data frame*  
  A DEG results table. By default it must include gene identifiers, gene labels, cell-group labels, shrunken log fold changes, p-values, a detection flag, and log mean expression values.

- **genes**  
  *character vector*  
  Genes or IDs to include in the heatmap. Matching is performed against `gene_col`.

- **cell_groups**  
  *character vector, optional*  
  Cell groups to include. If `NULL`, cell groups are selected automatically by requiring at least one gene from `genes` to satisfy the detection, significance, and mean expression thresholds.

- **sig_p_value_thresh**  
  Significance threshold applied to `p_value_col` when selecting cell groups and when outlining significant tiles in the plot.

- **min_log_mean_expression**  
  Minimum value in `log_mean_expression_col` required for a DEG to count toward retaining a cell group.

- **gene_col**, **gene_label_col**, **cell_group_col**, **lfc_col**, **p_value_col**, **present_col**, **log_mean_expression_col**  
  Column-name overrides for `deg_table`, for tables that don't use the default naming.

## Value

A `ggplot` object showing clustered log fold changes for the selected genes across retained cell groups.

## Details

Cell groups are retained only if they contain at least one gene from `genes` that is detected, significant, and above `min_log_mean_expression`. The retained gene-by-cell-group log fold-change matrix is then hierarchically clustered along both axes and plotted with significant tiles outlined.

## Examples

```r
plot_deg_heatmap(
  deg_table = tbx16_degs,
  genes = common_degs,
  cell_groups = c("Neurons", "Glia"),
  sig_p_value_thresh = 0.05,
  min_log_mean_expression = -2
)
```
