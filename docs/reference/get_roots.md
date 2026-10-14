# `get_roots`

> **Branch note:** this function is exported from platt's `develop_refactor` branch only; it is not in `develop` (or `main`) yet.

Finds the root node(s) of a cell state's connected component in a state graph.

```r
get_roots(ct, combined_psg)
```

## Arguments

- **ct**  
  *character*  
  The cell state whose component's roots should be found.

- **combined_psg**  
  *igraph or cell_state_graph*  
  A `cell_state_graph` object or an `igraph`.

## Value

A character vector of root node names in `ct`'s connected component. Returns an empty character vector if `ct` is not in the graph.

## Details

A root is any node with no incoming edges. `get_roots()` treats the graph as undirected when finding `ct`'s connected component, so it returns every origin of `ct`'s lineage — not just the single path back through `ct`'s own parents.
