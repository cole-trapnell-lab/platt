# `get_parents`

Retrieves the parent nodes of a given cell state in a state graph.

```r
get_parents(state_graph, cell_state)
```

## Arguments

- **state_graph**  
  *igraph*  
  An igraph object representing the state graph.

- **cell_state**  
  *character*  
  The specific cell state for which to find parent nodes.

## Value

A character vector of parent node names. If there are no parent nodes, an empty character vector is returned.

## Details

A state's parents are the states with an edge pointing *into* it. `get_all_parents()` is the recursive counterpart of `get_parents()` — it walks all the way up the graph to return every upstream ancestor, rather than just the immediate parent(s).

## Examples

```r
# Assuming `state_graph` is a pre-defined igraph object and `cell_state` is a valid node in the graph
parents <- get_parents(state_graph, cell_state)
print(parents)
```
