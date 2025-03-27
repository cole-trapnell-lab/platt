# `fit_genotype_ccm`

Fits a cell count model (CCM) for a given genotype using a specified dataset of cell counts (CCS).  
Supports flexible modeling of time intervals, perturbations, batch effects, and various backend optimization methods.

```r
fit_genotype_ccm(
  genotype,
  ccs,
  prior_state_transition_graph = NULL,
  interval_col = "timepoint",
  perturbation_col = "gene_target",
  batch_col = "expt",
  ctrl_ids = default_ctrl_ids,
  contrast_time_start = NULL,
  contrast_time_stop = NULL,
  num_time_breaks = NULL,
  independent_spline_for_ko = TRUE,
  edge_allowlist = NULL,
  edge_denylist = NULL,
  penalize_by_distance = TRUE,
  keep_ccs = TRUE,
  vhat_method = "bootstrap",
  num_threads = 1,
  backend = "nlopt",
  sparsity_factor = 0.01,
  num_bootstraps = 10,
  ftol_rel = 1e-6
)
```

## Arguments

- **genotype**  
  *character*  
  The genotype to fit the model for.

- **ccs**  
  *CCS object*  
  A cell count dataset containing the data to model.

- **prior_state_transition_graph**  
  *optional*  
  A prior state transition graph to use as a model constraint.

- **interval_col**  
  *character*  
  Column name in `ccs` representing time intervals. Default is `"timepoint"`.

- **perturbation_col**  
  *character*  
  Column name in `ccs` representing perturbations (e.g., gene target). Default is `"gene_target"`.

- **batch_col**  
  *character*  
  Column name in `ccs` representing batch info. Default is `"expt"`.

- **ctrl_ids**  
  *character vector*  
  Control IDs to include in the model.

- **contrast_time_start**  
  *numeric*  
  Start time for contrast interval. If `NULL`, uses minimum time in the data.

- **contrast_time_stop**  
  *numeric*  
  Stop time for contrast interval. If `NULL`, uses maximum time in the data.

- **num_time_breaks**  
  *numeric*  
  Number of time breakpoints for spline modeling. Default is `NULL` (interpreted as 3).

- **independent_spline_for_ko**  
  *logical*  
  Use independent spline for knockout modeling. Default is `TRUE`.

- **edge_allowlist**  
  *optional*  
  List of edges allowed in the model.

- **edge_denylist**  
  *optional*  
  List of edges denied in the model.

- **penalize_by_distance**  
  *logical*  
  Whether to penalize edges by distance. Default is `TRUE`.

- **keep_ccs**  
  *logical*  
  Whether to retain the CCS object in the model. Default is `TRUE`.

- **vhat_method**  
  *character*  
  Method to compute `vhat`. Default is `"bootstrap"`.

- **num_threads**  
  *numeric*  
  Number of threads for computation. Default is `1`.

- **backend**  
  *character*  
  Optimization backend. Default is `"nlopt"`.

- **sparsity_factor**  
  *numeric*  
  Sparsity factor for model selection. Default is `0.01`.

- **num_bootstraps**  
  *numeric*  
  Number of bootstraps for `vhat` estimation. Default is `10`.

- **ftol_rel**  
  *numeric*  
  Relative tolerance for optimization. Default is `1e-6`.

## Value

A fitted genotype CCM object containing the model, metadata, and fit diagnostics.

## Details

The function subsets the input CCS based on the genotype, control IDs, and time range.  
It builds a model formula based on time points and parameters, and fits the model using the chosen backend.  
Batch effects are corrected if multiple batches are detected. The model matrix is checked for full rank before fitting.

## Examples

```r
genotype_ccm <- fit_genotype_ccm(
  genotype = "example_genotype",
  ccs = example_ccs,
  interval_col = "timepoint",
  perturbation_col = "gene_target",
  batch_col = "expt"
)
```
