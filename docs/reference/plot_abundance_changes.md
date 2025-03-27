# `plot_abundance_changes`

Generates a plot to visualize changes in cell state abundances.

```r
plot_abundance_changes(
  cell_state_graph,
  comp_abund_table,
  facet_group = NULL,
  scale_node,
  plot_labels,
  arrow_unit,
  node_size,
  node_scale,
  con_colour,
  legend_position,
  node_label_width,
  group_label_font_size,
  fc_limits
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

- **arrow_unit**  
  *numeric*  
  Size of arrows in the plot.

- **node_size**  
  *numeric*  
  Base size of nodes.

- **node_scale**  
  *numeric*  
  Scaling factor for node sizes.

- **con_colour**  
  *character*  
  Color for node connections.

- **legend_position**  
  *character*  
  Position of the legend.

- **node_label_width**  
  *numeric*  
  Width of node labels.

- **group_label_font_size**  
  *numeric*  
  Font size for group labels.

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
