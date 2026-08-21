# `get_children`

Retrieves the children of a given cell state in a state graph.

```r
get_children(state_graph, cell_state)
```

## Arguments

- **state_graph**  
  *igraph*  
  An igraph object representing the state graph.

- **cell_state**  
  *character*  
  A vertex in the state graph for which to find the children.

## Value

A character vector of the names of the children of the given cell state. If there are no children, an empty character vector is returned.

## Details

A state's children are the states its own edges point *out to* — the reverse relationship of [`get_parents()`](get_parents).
