# `fit_wt_model`

Fits a wild type (WT) model to a cell dataset (CDS) by estimating cell count dynamics over time and accounting for nuisance variables.  
Supports various customization options for model fitting, including user-defined formulas, size factors, and control IDs.

```r
fit_wt_model(
  cds,
  sample_group,
  cell_group,
  main_model_formula_str = NULL,
  num_time_breaks = 4,
  nuisance_model_formula_str = "~1",
  ctrl_ids = NULL,
  sparsity_factor = 1,
  vhat_method = "bootstrap",
  interval_col = "timepoint",
  perturbation_col = "knockout",
  batch_col = "expt",
  start_time = NULL,
  stop_time = NULL,
  interval_step = 2,
  log_abund_detection_thresh = -5,
  keep_ccs = TRUE,
  q_val = 0.1,
  edge_allowlist = NULL,
  edge_denylist = NULL,
  base_penalty = 1,
  keep_cds = TRUE,
  verbose = FALSE,
  num_threads = 1,
  backend = "nlopt",
  penalize_by_distance = TRUE,
  embryo_size_factors = NULL,
  batches_excluded_from_assembly = c(),
  ...
)
```

## Arguments

- **cds**  
  *CDS object*  
  A cell dataset containing cell data.

- **sample_group**  
  *string*  
  Column in `colData(cds)` defining sample groups.

- **cell_group**  
  *string*  
  Column in `colData(cds)` defining cell groups.

- **main_model_formula_str**  
  *string*  
  Main model formula. If `NULL`, generated automatically.

- **num_time_breaks**  
  *integer*  
  Number of time breaks for main model formula. Default is `4`.

- **nuisance_model_formula_str**  
  *string*  
  Nuisance model formula. Default is `"~1"`.

- **ctrl_ids**  
  *vector*  
  Control IDs. If `NULL`, inferred from `perturbation_col`.

- **sparsity_factor**  
  *numeric*  
  Sparsity factor for model selection. Default is `1`.

- **vhat_method**  
  *string*  
  Method for estimating variance. Default is `"bootstrap"`.

- **interval_col**  
  *string*  
  Column in `colData(cds)` defining time intervals. Default is `"timepoint"`.

- **perturbation_col**  
  *string*  
  Column in `colData(cds)` defining perturbation groups. Default is `"knockout"`.

- **batch_col**  
  *string*  
  Column in `colData(cds)` defining batch groups. Default is `"expt"`.

- **start_time**  
  *numeric*  
  Start time for model. If `NULL`, inferred from data.

- **stop_time**  
  *numeric*  
  Stop time for model. If `NULL`, inferred from data.

- **interval_step**  
  *numeric*  
  Step size for time intervals. Default is `2`.

- **log_abund_detection_thresh**  
  *numeric*  
  Threshold for log abundance detection. Default is `-5`.

- **keep_ccs**  
  *logical*  
  Whether to retain the CCS. Default is `TRUE`.

- **q_val**  
  *numeric*  
  Q-value threshold for significance. Default is `0.1`.

- **edge_allowlist**  
  *vector*  
  Edges to allow in the model.

- **edge_denylist**  
  *vector*  
  Edges to deny in the model.

- **base_penalty**  
  *numeric*  
  Base penalty for model selection. Default is `1`.

- **keep_cds**  
  *logical*  
  Whether to retain the CDS. Default is `TRUE`.

- **verbose**  
  *logical*  
  Whether to print verbose messages. Default is `FALSE`.

- **num_threads**  
  *integer*  
  Number of threads to use. Default is `1`.

- **backend**  
  *string*  
  Backend for model fitting. Default is `"nlopt"`.

- **penalize_by_distance**  
  *logical*  
  Whether to penalize by distance. Default is `TRUE`.

- **embryo_size_factors**  
  *named vector*  
  Size factors for embryos. Default is `NULL`.

- **batches_excluded_from_assembly**  
  *vector*  
  Batch IDs to exclude from assembly. Default is an empty vector.

- **...**  
  Additional arguments passed to the underlying model fitting functions.

## Value

A fitted wild type cell count model object, or `NULL` if no control cells are available.

## Details

The function subsets the CDS to include only control cells based on `perturbation_col` and `ctrl_ids`.  
It constructs a cell count set (CCS) and fits a model using the specified formulas and parameters.  
If no control cells or only a single cell group is present, the function returns `NULL`.

## Examples

```r
wt_model <- fit_wt_model(
  cds = my_cds,
  sample_group = "sample",
  cell_group = "cell_type",
  interval_col = "timepoint",
  perturbation_col = "knockout",
  batch_col = "batch"
)
```
