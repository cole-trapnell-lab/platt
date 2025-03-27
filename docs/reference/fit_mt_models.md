# `fit_mt_models`

Fits models for analyzing multi-timepoint perturbation data using a cell dataset (`cds`).  
Allows specification of parameters for model formulae, perturbation settings, batch correction, and time interval handling.

```r
fit_mt_models(
  cds,
  sample_group,
  cell_group,
  main_model_formula_str = NULL,
  num_time_breaks = 3,
  nuisance_model_formula_str = "~1",
  ctrl_ids = NULL,
  mt_ids = NULL,
  sparsity_factor = 1,
  vhat_method = "bootstrap",
  interval_col = "timepoint",
  perturbation_col = "knockout",
  batch_col = "expt",
  newdata = tibble::tibble(),
  start_time = NULL,
  stop_time = NULL,
  interval_step = 2,
  log_abund_detection_thresh = -5,
  q_val = 0.1,
  edge_allowlist = NULL,
  edge_denylist = NULL,
  keep_cds = TRUE,
  keep_ccs = TRUE,
  verbose = FALSE,
  num_threads = 1,
  backend = "nlopt",
  penalize_by_distance = TRUE,
  independent_spline_for_ko = TRUE,
  num_bootstraps = 10,
  embryo_size_factors = NULL,
  batches_excluded_from_assembly = c()
)
```

## Arguments

- **cds**  
  *cell dataset*  
  The input dataset containing cell data.

- **sample_group**  
  *character*  
  Sample grouping variable.

- **cell_group**  
  *character*  
  Cell grouping variable.

- **main_model_formula_str**  
  *character*  
  Main model formula. Default is `NULL`.

- **num_time_breaks**  
  *integer*  
  Number of time breaks. Default is `3`.

- **nuisance_model_formula_str**  
  *character*  
  Nuisance model formula. Default is `"~1"`.

- **ctrl_ids**  
  *vector*  
  Control identifiers. Default is `NULL`.

- **mt_ids**  
  *vector*  
  Multi-timepoint identifiers. Default is `NULL`.

- **sparsity_factor**  
  *numeric*  
  Sparsity adjustment factor. Default is `1`.

- **vhat_method**  
  *character*  
  Method for variance estimation. Default is `"bootstrap"`.

- **interval_col**  
  *character*  
  Column name for time intervals. Default is `"timepoint"`.

- **perturbation_col**  
  *character*  
  Column name for perturbation IDs. Default is `"knockout"`.

- **batch_col**  
  *character*  
  Column name for batch identifiers. Default is `"expt"`.

- **newdata**  
  *tibble*  
  New data for prediction. Default is an empty tibble.

- **start_time**  
  *numeric*  
  Start time. Default is `NULL`.

- **stop_time**  
  *numeric*  
  Stop time. Default is `NULL`.

- **interval_step**  
  *integer*  
  Step size for time intervals. Default is `2`.

- **log_abund_detection_thresh**  
  *numeric*  
  Log abundance detection threshold. Default is `-5`.

- **q_val**  
  *numeric*  
  Q-value threshold for significance. Default is `0.1`.

- **edge_allowlist**  
  *vector*  
  Allowed edges in the analysis. Default is `NULL`.

- **edge_denylist**  
  *vector*  
  Denied edges in the analysis. Default is `NULL`.

- **keep_cds**  
  *logical*  
  Whether to retain the CDS in output. Default is `TRUE`.

- **keep_ccs**  
  *logical*  
  Whether to retain the CCS in output. Default is `TRUE`.

- **verbose**  
  *logical*  
  Whether to print verbose output. Default is `FALSE`.

- **num_threads**  
  *integer*  
  Number of computation threads. Default is `1`.

- **backend**  
  *character*  
  Optimization backend. Default is `"nlopt"`.

- **penalize_by_distance**  
  *logical*  
  Whether to penalize by distance. Default is `TRUE`.

- **independent_spline_for_ko**  
  *logical*  
  Use independent splines for knockouts. Default is `TRUE`.

- **num_bootstraps**  
  *integer*  
  Number of bootstraps for variance estimation. Default is `10`.

- **embryo_size_factors**  
  *named vector*  
  Embryo size factors. Default is `NULL`.

- **batches_excluded_from_assembly**  
  *vector*  
  Batch IDs to exclude. Default is an empty vector.

## Value

A tibble containing fitted perturbation models and associated metadata.

## Details

This function performs multi-timepoint perturbation analysis on cell count data.  
It supports handling of batch effects, customizable formulas, and exclusion of specific data components.  
User-supplied embryo size factors and batch exclusions are supported for fine-tuned control.

## Examples

```r
# Example usage:
fit_mt_models(
  cds,
  sample_group = "sample",
  cell_group = "cell_type"
)
```
