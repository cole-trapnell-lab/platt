# `get_descendants`

> **Branch note:** this function is exported from platt's `develop_refactor` branch only; it is not in `develop` (or `main`) yet.

Retrieves every downstream descendant of a given cell state in a state graph.

```r
get_descendants(cell_type, combined_psg)
```

## Arguments

- **cell_type**  
  *character*  
  The cell state for which to find descendants.

- **combined_psg**  
  *igraph or cell_state_graph*  
  A `cell_state_graph` object or an `igraph`.

## Value

A character vector of descendant node names, excluding `cell_type` itself. Returns an empty character vector if `cell_type` is not in the graph or has no descendants.

## Details

`get_descendants()` is the downstream counterpart of [`get_all_parents()`](get_all_parents): where `get_all_parents()` walks up the graph to return every upstream ancestor, `get_descendants()` walks down to return every node reachable from `cell_type`.
