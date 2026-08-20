# `compare_genes_within_state_graph`

Compares gene expression within a state graph for a given cell count set (CCS), contrasting a perturbation against control(s) within each cell state.

```r
compare_genes_within_state_graph(
  ccs,
  perturbation_col = "perturbation",
  control_ids = c("Control"),
  nuisance_model_formula_str = "0",
  ambient_coeffs = NULL,
  cell_groups = NULL,
  assembly_group = NULL,
  state_graph = NULL,
  perturbations = NULL,
  group_nodes_by = NULL,
  log_fc_thresh = 1,
  abs_expr_thresh = 0.001,
  sig_thresh = 0.05,
  detection_min_samples = 2,
  min_cells_per_pseudobulk = NULL,
  cores = 1,
  write_dir = NULL,
  max_simultaneous_genes = NULL,
  cv_threshold = 100,
  filter_mode = c("global", "by_background_count"),
  filter_only = FALSE,
  min_mean_expr = NULL,
  background_bottom_frac = 0.25,
  background_quantile_p = 0.99,
  background_count_floor = 2,
  alpha = 0.05,
  perf_summary_path = NULL,
  profile = FALSE,
  profile_out = NULL,
  ...
)
```

## Arguments

- **ccs**  
  *cell_count_set*  
  A cell count set object.

- **perturbation_col**  
  *character*  
  The column name in the cell data set that contains perturbation information. Default is `"perturbation"`.

- **control_ids**  
  *character vector*  
  A vector of control IDs. Default is `c("Control")`.

- **nuisance_model_formula_str**  
  *character*  
  A string representing the nuisance model formula. Default is `"0"`.

- **ambient_coeffs**  
  *optional*  
  Coefficients for ambient noise. Default is `NULL`.

- **cell_groups**  
  *character vector*  
  A vector of cell groups to be analyzed. Default is `NULL`.

- **assembly_group**  
  *character*  
  The assembly group to be analyzed. Default is `NULL`.

- **state_graph**  
  *igraph*  
  A state graph object. Default is `NULL`.

- **perturbations**  
  *character vector*  
  A vector of perturbations to be analyzed. Default is `NULL`.

- **group_nodes_by**  
  *character*  
  A parameter to group nodes by. Default is `NULL`.

- **log_fc_thresh**  
  *numeric*  
  Log fold change threshold. Default is `1`.

- **abs_expr_thresh**  
  *numeric*  
  Absolute expression threshold. Default is `0.001`.

- **sig_thresh**  
  *numeric*  
  Significance threshold. Default is `0.05`.

- **detection_min_samples**  
  *numeric*  
  Minimum number of pseudobulk samples required for a gene to be considered detected (interpreted per mode: overall for `"global"`; per-arm above background for `"by_background_count"`). Default is `2`.

- **min_cells_per_pseudobulk**  
  *numeric*  
  Minimum number of cells per pseudobulk. Default is `NULL`.

- **cores**  
  *numeric*  
  Number of cores to use for parallel processing. Default is `1`.

- **write_dir**  
  *character*  
  Directory to write output files. Default is `NULL`.

- **max_simultaneous_genes**  
  *numeric*  
  Maximum number of genes to analyze simultaneously. Default is `NULL`.

- **cv_threshold**  
  *numeric*  
  Coefficient of variation threshold. Default is `100`.

- **filter_mode**  
  *character*  
  Gene filtering strategy: `"global"` uses aggregate detection, `"by_background_count"` uses per-cell-type background rates.

- **filter_only**  
  *logical*  
  Whether to only run the gene filtering step (without fitting models). Default is `FALSE`.

- **min_mean_expr**  
  *numeric*  
  Minimum mean expression required for a gene to pass filtering. Default is `NULL`.

- **background_bottom_frac**  
  *numeric*  
  Fraction of samples treated as the "background" (lowest-expressing) tail when using `"by_background_count"` filtering. Default is `0.25`.

- **background_quantile_p**  
  *numeric*  
  Quantile of the background distribution used to estimate the background count. Default is `0.99`.

- **background_count_floor**  
  *numeric*  
  Minimum floor applied to the estimated background count. Default is `2`.

- **alpha**  
  *numeric*  
  Significance level used for statistical thresholds. Default is `0.05`.

- **perf_summary_path**  
  *character, optional*  
  Optional path to a TSV where runtime and filtering diagnostics will be appended.

- **profile**  
  *logical*  
  Whether to enable `Rprof` profiling for the DEG run. Default is `FALSE`.

- **profile_out**  
  *character, optional*  
  Path for the `Rprof` output (default derived from `write_dir` when `profile = TRUE`).

- **...**  
  Additional arguments.

- **gene_ids**  
  *character vector*  
  A vector of gene IDs to be analyzed. Default is `NULL`.

## Value

A data frame containing the comparison results.

## Details

This function compares gene expression within a state graph for a given cell count set, contrasting each perturbation against the specified control(s) within each cell state — in contrast to `compare_genes_over_graph()`, which compares expression across states of the graph rather than between perturbation and control.

## Examples

```r
genes_within_cell_state <- compare_genes_within_state_graph(
  ccs,
  perturbation_col = "perturbation",
  control_ids = c("Control")
)
```
