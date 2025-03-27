# `plot_gene_expr`

Plots gene expression data on a cell state graph.

```r
plot_gene_expr(
  cell_state_graph,
  genes,
  arrow_unit = 7,
  con_colour = "lightgrey",
  fract_expr = 0.0,
  mean_expr = 0.0,
  legend_position = "right",
  plot_labels = TRUE,
  aggregate = FALSE,
  scale_to_range = FALSE,
  log_expr = FALSE,
  node_size,
  pseudocount = 1e-5,
  expr_limits = NULL,
  node_label_width = 50,
  group_label_font_size = 1
)
```

## Arguments

- **cell_state_graph**  
  *object*  
  Cell state graph data.

- **genes**  
  *vector*  
  Gene names to plot.

- **arrow_unit**  
  *numeric*  
  Size of arrows in the plot. Default is `7`.

- **con_colour**  
  *character*  
  Color for connections. Default is `"lightgrey"`.

- **fract_expr**  
  *numeric*  
  Minimum fraction of cells expressing the gene. Default is `0.0`.

- **mean_expr**  
  *numeric*  
  Minimum mean expression level. Default is `0.0`.

- **legend_position**  
  *character*  
  Legend position. Default is `"right"`.

- **plot_labels**  
  *logical*  
  Whether to show labels. Default is `TRUE`.

- **aggregate**  
  *logical*  
  Aggregate gene expression data. Default is `FALSE`.

- **scale_to_range**  
  *logical*  
  Scale expression data to range. Default is `FALSE`.

- **log_expr**  
  *logical*  
  Log-transform expression data. Default is `FALSE`.

- **node_size**  
  *numeric*  
  Offset used to prevent clipping when saving plots.

- **pseudocount**  
  *numeric*  
  Pseudocount for log transformation. Default is `1e-5`.

- **expr_limits**  
  *numeric vector*  
  Length 2 vector specifying expression limits. Default is `NULL`.

- **node_label_width**  
  *numeric*  
  Width of node labels. Default is `50`.

- **group_label_font_size**  
  *numeric*  
  Font size for group labels. Default is `1`.

## Value

A `ggplot2` object representing the gene expression on the cell state graph.

## Examples

```r
# Example usage:
# plot_gene_expr(cell_state_graph, genes = c("Gene1", "Gene2"))
```
