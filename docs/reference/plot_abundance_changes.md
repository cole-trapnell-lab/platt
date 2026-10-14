# `plot_abundance_changes`

Generates a plot to visualize changes in cell state abundances.

```r
plot_abundance_changes(
  cell_state_graph,
  comp_abund_table,
  facet_group = NULL,
  scale_node = FALSE,
  plot_labels = TRUE,
  label_size = 3,
  node_size = 2,
  node_scale = 1,
  arrow_unit = 7,
  arrow_color = "darkgrey",
  legend_position = "right",
  group_label_size = 1,
  fc_limits = c(-3, 3)
)
```

## Arguments

- **cell_state_graph**  
  *object*  
  A cell state graph containing the graph and layout info.

- **comp_abund_table**  
  *data.frame*  
  Comparison abundance table.

- **facet_group**  
  *character*  
  Optional column in `comp_abund_table` to facet the plot.

- **scale_node**  
  *logical*  
  Whether to scale nodes by significance.

- **plot_labels**  
  *logical*  
  Whether to label nodes in the plot.

- **label_size**  
  *numeric*  
  Font size of node labels. Default is `3`.

- **node_size**  
  *numeric*  
  Base size of nodes.

- **node_scale**  
  *numeric*  
  Scaling factor for node sizes.

- **arrow_unit**  
  *numeric*  
  Size of arrows in the plot.

- **arrow_color**  
  *character*  
  Color of the connections between nodes. Default is `"darkgrey"`.

- **legend_position**  
  *character*  
  Position of the legend.

- **group_label_size**  
  *numeric*  
  Font size for group labels. Default is `1`.

- **fc_limits**  
  *numeric vector*  
  Length-2 vector for fold change color scale limits.

## Value

A `ggplot` object representing the abundance changes.

## Examples

```r
# Example usage:
# plot_abundance_changes(cell_state_graph, comp_abund_table)
```
