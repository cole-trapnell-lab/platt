# `contract_state_graph`

Simplifies a directed state transition graph by contracting (grouping) nodes according to a specified label. Assumes the graph's nodes correspond to the groups in the associated cell count set (CCS).

```r
contract_state_graph(
  ccs,
  state_graph,
  group_nodes_by,
  edge_attr_policy = list(
    weight = "sum",
    name = "concat",
    num_perturbs_supporting = "sum",
    max_timeseries_path_score_supporting = "sum",
    total_timeseries_path_score_supporting = "sum",
    max_perturb_path_score_supporting = "sum",
    total_perturb_path_score_supporting = "sum",
    max_path_score_supporting = "sum",
    total_path_score_supporting = "sum",
    edge_name = "ignore",
    support_label = "concat",
    supporting_perturbs = "concat",
    "ignore"
  )
)
```

## Arguments

- **ccs**  
  *cell_count_set*  
  The cell count set associated with `state_graph`. Its `cell_group_assignments` metadata and `colData` are used to determine, for each node in `state_graph`, which value of `group_nodes_by` it should be contracted into.

- **state_graph**  
  *igraph*  
  A directed state transition graph whose vertex names correspond to the cell groups in `ccs`.

- **group_nodes_by**  
  *character*  
  Column name (in `colData(ccs@cds)`) specifying the grouping variable to contract nodes by. Each original node is mapped to the most frequent value of this column among its member cells, and nodes sharing a value are merged into one.

- **edge_attr_policy**  
  *list*  
  A named list passed to `igraph::simplify()`'s `edge.attr.comb` describing how to combine edge attributes (e.g. `"sum"`, `"concat"`, `"ignore"`) when parallel edges arise from the contraction. Defaults combine known path-scoring attributes by summation, concatenate label/name attributes, and ignore the rest.

## Value

An igraph object: the contracted state transition graph, with vertex names set to the `group_nodes_by` values and a `cell_group` vertex attribute recording the grouping variable used.

## Details

This function collapses a fine-grained state graph (e.g. per-cluster) into a coarser one (e.g. per-cell-type) by merging nodes that share the same value of `group_nodes_by`. It is commonly used after `run_wildtype_assembly()` (or after mutant/perturbation assembly) to view a state graph at a coarser level of granularity, such as collapsing individual clusters into their parent cell types.

## Examples

```r
contract_graph <- contract_state_graph(ccs, cluster_mt_graph, group_nodes_by = "cell_type")
```
