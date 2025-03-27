# `plot_annotations`

Generates a plot of cell state graphs with various customization options.

```r
plot_annotations(
  cell_state_graph,
  color_nodes_by = NULL,
  label_nodes_by = NULL,
  arrow_unit = 7,
  node_size = 2,
  con_colour = "darkgrey",
  legend_position = "none",
  min_edge_size = 0.1,
  max_edge_size = 2,
  edge_weights = NULL,
  plot_labels = TRUE,
  group_label_font_size = 1,
  node_label_width = 50
)
```

## Arguments

- **cell_state_graph**  
  *object*  
  The cell state graph data.

- **color_nodes_by**  
  *character*  
  Attribute to color nodes by. Default is `NULL`.

- **label_nodes_by**  
  *character*  
  Attribute to label nodes by. Default is `NULL`.

- **arrow_unit**  
  *numeric*  
  Size of the arrows. Default is `7`.

- **node_size**  
  *numeric*  
  Size of the nodes. Default is `2`.

- **con_colour**  
  *character*  
  Color of the connections. Default is `"darkgrey"`.

- **legend_position**  
  *character*  
  Position of the legend. Default is `"none"`.

- **min_edge_size**  
  *numeric*  
  Minimum edge size. Default is `0.1`.

- **max_edge_size**  
  *numeric*  
  Maximum edge size. Default is `2`.

- **edge_weights**  
  *numeric vector*  
  Weights for edges. Default is `NULL`.

- **plot_labels**  
  *logical*  
  Whether to show labels. Default is `TRUE`.

- **group_label_font_size**  
  *numeric*  
  Font size for group labels. Default is `1`.

- **node_label_width**  
  *numeric*  
  Width of node labels. Default is `50`.

## Value

A `ggplot` object representing the cell state graph with annotations.

## Examples

```r
# Example usage:
# plot_annotations(cell_state_graph)
```
