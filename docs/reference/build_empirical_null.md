# `build_empirical_null`

Builds a control-split empirical null for the DEG artifact.

```r
build_empirical_null(
  cds,
  sample_group = "embryo_ID",
  cell_group = "cell_type",
  control_ids = c("ctrl-inj"),
  perturbation_col = "perturbation",
  n_perturb_grid = NULL,
  n_cell_types = 80L,
  min_cells_per_type = 100L,
  gene_subsample = NULL,
  nperm = 10000L,
  cores = 1L,
  seed = 42L,
  max_simultaneous_genes = 2000L,
  write_dir = tempfile("efdr_null_"),
  verbose = TRUE
)
```

## Arguments

- **cds**  
  *cell_data_set*  
  A `cell_data_set` containing the experiment's **control** cells.

- **sample_group**, **cell_group**  
  `colData` columns for replicate and cell type.

- **control_ids**, **perturbation_col**  
  Control values and the perturbation column.

- **n_perturb_grid**  
  Sizes of the pseudo-perturbation (control samples relabelled). If `NULL`, chosen from the control count to give control-heavy log-ratios bracketing real perturbations.

- **n_cell_types**, **min_cells_per_type**  
  Stratified cell-type panel controls.

- **gene_subsample**  
  Optional cap on genes fit per cell type (expression-spanning).

- **nperm**, **cores**, **max_simultaneous_genes**, **seed**, **write_dir**, **verbose**  
  Pass-through / control parameters.

## Value

A tibble with columns `cell_group`, `log_mean_expression`, `z`, `log_ratio`.

## Details

Relabels control samples as a pseudo-perturbation at several split sizes, runs the standard within-state DEG contrast on a stratified subsample of cell types (and, optionally, genes), and returns the pooled null gene-tests tagged with their sampling log-ratio. Training data for [`train_efdr_model()`](train_efdr_model).

_Note: [`train_efdr_model()`](train_efdr_model)'s own docs describe its `null` input as needing `gene_short_name`, `cell_group`, `log_mean_expression`, `z`, and `n_pert_cells` — a `n_pert_cells`/arm-cell-count column not listed here. Worth confirming against the current source before relying on the column list above being exhaustive._

## Examples

```r
null_draws <- build_empirical_null(control_only_cds,
  sample_group = "embryo",
  cell_group = "cell_type",
  control_ids = c("ctrl-inj"),
  perturbation_col = "gene_target",
  cores = 6
)
```
