# `plot_gene_expression`

Plots gene expression data on a cell state graph.

```r
plot_gene_expression(
  cell_state_graph,
  genes,
  arrow_unit = 7,
  node_size = 2,
  arrow_color = "lightgrey",
  fract_expr = 0,
  mean_expr = 0,
  legend_position = "right",
  plot_labels = TRUE,
  label_size = 3,
  aggregate = FALSE,
  scale_to_range = FALSE,
  log_expr = FALSE,
  pseudocount = 1e-05,
  expr_limits = NULL,
  group_label_size = 1
)
```

## Arguments

- **cell_state_graph**  
  *cell_state_graph*  
  An object containing the cell state graph data.

- **genes**  
  *character vector*  
  A vector of gene names to plot.

- **arrow_unit**  
  *numeric*  
  Size of the arrows in the plot. Default is `7`.

- **node_size**  
  *numeric*  
  Does not actually control the size of the nodes; used to specify the offset for invisible points that help prevent the plot from getting clipped upon saving. Default is `2`.

- **arrow_color**  
  *character*  
  Colour for the connections in the plot. Default is `"lightgrey"`.

- **fract_expr**  
  *numeric*  
  Minimum fraction of cells expressing the gene to be considered. Default is `0`.

- **mean_expr**  
  *numeric*  
  Minimum mean expression level to be considered. Default is `0`.

- **legend_position**  
  *character*  
  Position of the legend in the plot. Default is `"right"`.

- **plot_labels**  
  *logical*  
  Whether to plot labels. Default is `TRUE`.

- **label_size**  
  *numeric*  
  Font size of node labels. Default is `3`.

- **aggregate**  
  *logical*  
  Whether to aggregate gene expression data. Default is `FALSE`.

- **scale_to_range**  
  *logical*  
  Whether to scale expression data to range. Default is `FALSE`.

- **log_expr**  
  *logical*  
  If `TRUE`, the expression values will be log-transformed. Default is `FALSE`.

- **pseudocount**  
  *numeric*  
  Pseudocount to add when log-transforming expression data. Default is `1e-05`.

- **expr_limits**  
  *numeric vector*  
  Vector of length 2 specifying the limits for expression values. Default is `NULL`.

- **group_label_size**  
  *numeric*  
  Font size of group labels. Default is `1`.

## Value

A ggplot2 object representing the gene expression on the cell state graph.

## Examples

```r
plot_gene_expression(cell_state_graph, genes = c("Gene1", "Gene2"))
```
