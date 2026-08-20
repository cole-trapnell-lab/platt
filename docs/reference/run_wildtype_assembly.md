# `run_wildtype_assembly`

Fits the wild-type model and assembles a state transition graph for a partition of the data.

```r
run_wildtype_assembly(
  cds,
  sample_group,
  cell_group,
  partition_name = NULL,
  main_model_formula_str = NULL,
  start_time = 18,
  stop_time = 72,
  interval_col = "timepoint",
  nuisance_model_formula_str = "~1",
  ctrl_ids = NULL,
  sparsity_factor = 0.01,
  perturbation_col = "perturbation",
  batch_col = "expt",
  verbose = FALSE,
  keep_ccs = TRUE,
  num_threads = 1,
  backend = "nlopt",
  q_val = 0.1,
  vhat_method = "bootstrap",
  num_bootstraps = 10,
  newdata = tibble(),
  edge_allowlist = NULL,
  edge_denylist = NULL,
  links_between_components = c("none", "ctp", "strongest-pcor", "strong-pcor"),
  component_col = "partition",
  embryo_size_factors = NULL,
  log_abund_detection_thresh = -5,
  interval_step = 2,
  min_interval = 4,
  max_interval = 24,
  min_pathfinding_lfc = 0,
  num_time_breaks = 4,
  batches_excluded_from_assembly = c(),
  force_allowlist = FALSE,
  min_penalty = 0.01,
  max_penalty = 1e+06,
  break_cycles = TRUE
)
```

## Arguments

- **cds**  
  *CDS object*  
  A SingleCellExperiment/CellDataSet with expression and metadata.

- **sample_group**  
  *character*  
  Column name in `colData` specifying sample grouping.

- **cell_group**  
  *character*  
  Column name in `colData` specifying cell grouping.

- **partition_name**  
  *character, optional*  
  Optional partition label appended to graph nodes.

- **main_model_formula_str**  
  *character*  
  Model formula string for the main effect.

- **start_time**  
  *numeric*  
  Start time for model fitting. Default is `18`.

- **stop_time**  
  *numeric*  
  Stop time for model fitting. Default is `72`.

- **interval_col**  
  *character*  
  Column name for time intervals. Default is `"timepoint"`.

- **nuisance_model_formula_str**  
  *character*  
  Nuisance model formula string. Default is `"~1"`.

- **ctrl_ids**  
  *character vector, optional*  
  Optional vector of control IDs.

- **sparsity_factor**  
  *numeric*  
  Sparsity factor for model fitting. Default is `0.01`.

- **perturbation_col**  
  *character*  
  Column name for perturbation labels. Default is `"perturbation"`.

- **batch_col**  
  *character*  
  Column name for batch labels. Default is `"expt"`.

- **verbose**  
  *logical*  
  Whether to log progress. Default is `FALSE`.

- **keep_ccs**  
  *logical*  
  Whether to keep `cell_count_set` outputs. Default is `TRUE`.

- **num_threads**  
  *numeric*  
  Number of threads to use. Default is `1`.

- **backend**  
  *character*  
  Optimization backend to use. Default is `"nlopt"`.

- **q_val**  
  *numeric*  
  FDR threshold. Default is `0.1`.

- **vhat_method**  
  *character*  
  Method to estimate `vhat`. Default is `"bootstrap"`.

- **num_bootstraps**  
  *numeric*  
  Number of bootstraps for `vhat`. Default is `10`.

- **newdata**  
  *tibble, optional*  
  Optional data frame of new timepoints for prediction. Default is an empty tibble.

- **edge_allowlist**  
  *optional*  
  Optional allowlist of edges.

- **edge_denylist**  
  *optional*  
  Optional denylist of edges.

- **links_between_components**  
  *character*  
  Strategy for linking graph components. One of `"none"`, `"ctp"`, `"strongest-pcor"`, `"strong-pcor"`.

- **component_col**  
  *character*  
  Column name for component labels. Default is `"partition"`.

- **embryo_size_factors**  
  *optional*  
  Optional size factors for embryo data.

- **log_abund_detection_thresh**  
  *numeric*  
  Log abundance detection threshold. Default is `-5`.

- **interval_step**  
  *numeric*  
  Step size for interval grid. Default is `2`.

- **min_interval**  
  *numeric*  
  Minimum interval length. Default is `4`.

- **max_interval**  
  *numeric*  
  Maximum interval length. Default is `24`.

- **min_pathfinding_lfc**  
  *numeric*  
  Minimum log-fold-change for pathfinding. Default is `0`.

- **num_time_breaks**  
  *numeric*  
  Number of time breaks for fitting. Default is `4`.

- **batches_excluded_from_assembly**  
  *vector*  
  Vector of batches to exclude. Default is empty.

- **force_allowlist**  
  *logical*  
  Whether to force the edge allowlist. Default is `FALSE`.

- **min_penalty**  
  *numeric*  
  Minimum penalty for model fitting. Default is `0.01`.

- **max_penalty**  
  *numeric*  
  Maximum penalty for model fitting. Default is `1e+06`.

- **break_cycles**  
  *logical*  
  Whether to break cycles in the graph. Default is `TRUE`.

## Value

An igraph object for the assembled wild-type graph, or `NA` on failure.

## Details

`run_wildtype_assembly()` is the exported top-level entry point for building a wild-type state transition graph from a CDS: it fits the wild-type cell count model over the specified partition and assembles the resulting graph, handling batch effects, edge allow/deny lists, and cycle breaking along the way. It replaces the wild-type half of what was previously a single `assemble_partition()` call — the perturbation/mutant half now lives in an unexported `run_perturbation_assembly()`.

## Examples

```r
cluster_wt_graph <- run_wildtype_assembly(
  cds,
  sample_group = "sample",
  cell_group = "cell_type",
  partition_name = "partition_1",
  start_time = 18,
  stop_time = 72
)
```
