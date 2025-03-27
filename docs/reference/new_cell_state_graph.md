# `new_cell_state_graph`

Creates a new cell state graph object from an input graph and a cell count set.

```r
new_cell_state_graph(
  state_graph,
  ccs,
  color_nodes_by = NULL,
  label_nodes_by = NULL,
  group_nodes_by = NULL,
  edge_weights = NULL,
  min_edge_size = 0.1,
  max_edge_size = 2,
  hide_unlinked_nodes = FALSE,
  genetic_requirements = data.frame(),
  num_layers = 1
)
```

## Arguments

- **state_graph**  
  *igraph object*  
  Input graph of class `'igraph'`.

- **ccs**  
  *cell_count_set*  
  Input cell count set of class `'cell_count_set'`.

- **color_nodes_by**  
  *optional*  
  Attribute to color nodes by.

- **label_nodes_by**  
  *optional*  
  Attribute to label nodes by.

- **group_nodes_by**  
  *optional*  
  Attribute to group nodes by.

- **edge_weights**  
  *optional*  
  Values to specify edge weights.

- **min_edge_size**  
  *numeric*  
  Minimum edge size for visualization. Default is `0.1`.

- **max_edge_size**  
  *numeric*  
  Maximum edge size for visualization. Default is `2`.

- **hide_unlinked_nodes**  
  *logical*  
  Whether to hide unlinked nodes. Default is `FALSE`.

- **genetic_requirements**  
  *data.frame*  
  Genetic requirements. Default is an empty data frame.

- **num_layers**  
  *integer*  
  Number of layers for the graph layout. Default is `1`.

## Value

A new `cell_state_graph` object.

