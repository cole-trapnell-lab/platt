# `categorize_genetic_requirements`

Summarizes perturbation effects to identify cell groups that are directly or indirectly lost for each perturbation, using a state transition graph.

```r
categorize_genetic_requirements(perturb_ccm_tbl, state_graph)
```

## Arguments

- **perturb_ccm_tbl**  
  *tibble*  
  A tibble of fitted genotype models, one row per perturbation, with a `perturb_summary_tbl` list-column (the per-cell-state loss summary).

- **state_graph**  
  *igraph*  
  State transition graph defining lineage relationships among cell groups.

## Value

A tibble with per-perturbation lists of directly and indirectly lost cell groups.

## Details

For each perturbation, cell states flagged as lost (`is_lost_when_present`) are classified using the graph's topology: a lost cell state with no lost parent upstream of it is **direct**; a lost cell state whose loss can be explained by an already-lost parent is **indirect**.

## Examples

```r
genetic_requirements <- categorize_genetic_requirements(perturb_ccm_tbl, muscle_state_graph@graph)
```
