# `get_siblings`

Retrieves the sibling states of a given cell state in a state transition graph.

```r
get_siblings(state_graph, cell_state)
```

## Arguments

- **state_graph**  
  *igraph*  
  An igraph object representing the state transition graph.

- **cell_state**  
  *character*  
  The cell state for which to find siblings.

## Value

A character vector of sibling states. If the cell state has no parents, an empty vector is returned.

## Details

Sibling states are defined as states that share the same parent state — computed by finding the cell state's parents via [`get_parents()`](get_parents) and returning their other children.
