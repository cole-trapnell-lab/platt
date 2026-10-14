# `plot_degs`

Plots differentially expressed genes (DEGs) on a cell state graph, with customizable appearance and layout options.

```r
plot_degs(
  cell_state_graph,
  deg_table,
  perturb_table = NULL,
  facet_group = "term",
  arrow_unit = 7,
  node_size = 1,
  arrow_color = "darkgrey",
  fract_expr = 0.0,
  mean_expr = 0.0,
  legend_position = "none",
  fc_limits = c(-3, 3),
  plot_labels = TRUE,
  label_size = 3,
  group_label_size = 1
)
```

## Arguments

- **cell_state_graph**  
  *object*  
  Contains the cell state graph and layout information.

- **deg_table**  
  *data.frame*  
  Data frame with DEGs including log fold changes and p-values.

- **perturb_table**  
  *data.frame* (optional)  
  Perturbation information. Default is `NULL`.

- **facet_group**  
  *string*  
  Facet group for plotting. Default is `"term"`.

- **arrow_unit**  
  *numeric*  
  Arrow size. Default is `7`.

- **node_size**  
  *numeric*  
  Node size. Default is `1`.

- **arrow_color**  
  *string*  
  Connection color. Default is `"darkgrey"`.

- **fract_expr**  
  *numeric*  
  Minimum fraction of expression. Default is `0.0`.

- **mean_expr**  
  *numeric*  
  Minimum mean expression. Default is `0.0`.

- **legend_position**  
  *string*  
  Legend position. Default is `"none"`.

- **fc_limits**  
  *numeric vector*  
  Limits for fold change scale. Default is `c(-3, 3)`.

- **plot_labels**  
  *logical*  
  Whether to display labels. Default is `TRUE`.

- **label_size**  
  *numeric*  
  Font size of node labels. Default is `3`.

- **group_label_size**  
  *numeric*  
  Font size for group labels. Default is `1`.

## Value

A `ggplot` object representing the cell state graph with DEGs plotted.

## Examples

```r
# Example usage:
plot_degs(cell_state_graph, deg_table, perturb_table = NULL)
```
