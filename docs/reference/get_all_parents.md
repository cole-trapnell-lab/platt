# `get_all_parents`

Recursively retrieves all upstream ancestor parents for a given cell state in a state graph.

```r
get_all_parents(state_graph, cell_state, visited = character())
```

## Arguments

- **state_graph**  
  *igraph or cell_state_graph*  
  An igraph or `cell_state_graph` object.

- **cell_state**  
  *character*  
  The cell state for which to find all ancestor parents.

- **visited**  
  *character vector*  
  Internal recursion accumulator of visited nodes. Not normally set by the caller.

## Value

A character vector of unique ancestor parent node names.

## Details

`get_all_parents()` is the recursive counterpart of [`get_parents()`](get_parents). Where `get_parents()` returns only the immediate parent(s) of a cell state, `get_all_parents()` walks all the way up the graph to return every upstream ancestor, guarding against cycles by tracking which nodes have already been visited.

## Examples

```r
# Assuming `state_graph` is a pre-defined igraph object and `cell_state` is a valid node in the graph
parents <- get_all_parents(state_graph, cell_state)
print(parents)
```
