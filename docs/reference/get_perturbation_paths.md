# `get_perturbation_paths`

Finds and scores paths between the cell groups a perturbation gains (or loses with no lost parent) and the cell groups it loses, on a pathfinding graph.

```r
get_perturbation_paths(
  perturbation_ccm,
  perturb_summary_tbl,
  pathfinding_graph,
  delta_log_abund_loss_thresh = 0
)
```

## Arguments

- **perturbation_ccm**  
  A fitted genotype CCM (see [`fit_genotype_ccm()`](fit_genotype_ccm)), passed through to the path-scoring step.

- **perturb_summary_tbl**  
  Per-cell-group loss/gain summary for this perturbation; must have `cell_group`, `loss_when_present`, and `gain_when_present` columns.

- **pathfinding_graph**  
  An `igraph` used for path search.

- **delta_log_abund_loss_thresh**  
  Minimum magnitude of `loss_when_present`/`gain_when_present` for a cell group to count as lost or gained (default `0`).

## Value

A tibble of scored paths between gained/parentless-lost cell groups and lost cell groups (one row per `from`/`to` pair with a path), or `NA` if no lost cell groups are present, no paths are found, or path search fails.

## Details

`get_perturbation_paths()` finds shortest paths on `pathfinding_graph` from cell groups that are gained (or lost with no lost parent already in the loss set) to cell groups that are lost, then scores each path for how well it explains the perturbation's abundance changes — a way to trace a perturbation's lost cell states back to a plausible upstream origin along the lineage.
