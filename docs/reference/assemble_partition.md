# `assemble_partition`

Assembles a partition for a given cell dataset (CDS) by fitting wild-type (WT) and mutant (MT) models, constructing state transition graphs, and assessing perturbation effects.

```r
assemble_partition(
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
  mt_ids = NULL,
  sparsity_factor = 0.01,
  perturbation_col = "gene_target",
  batch_col = "expt",
  max_num_cells = NULL,
  verbose = FALSE,
  keep_ccs = TRUE,
  num_threads = 1,
  backend = "nlopt",
  q_val = 0.1,
  vhat_method = "bootstrap",
  num_bootstraps = 10,
  newdata = tibble::tibble(),
  edge_allowlist = NULL,
  min_lfc = 0,
  links_between_components = c("ctp", "none", "strongest-pcor", "strong-pcor"),
  log_abund_detection_thresh = -5,
  batches_excluded_from_assembly = c(),
  component_col = "partition",
  embryo_size_factors = NULL
)
```

## Arguments

- **cds**  
  *CDS object*  
  Cell dataset to analyze.

- **sample_group**  
  *string*  
  Column name in CDS for sample groups.

- **cell_group**  
  *string*  
  Column name in CDS for cell groups.

- **partition_name**  
  *string*  
  Name of the partition. Default is `NULL`.

- **main_model_formula_str**  
  *string*  
  Main model formula. Default is `NULL`.

- **start_time**  
  *numeric*  
  Start time for analysis. Default is `18`.

- **stop_time**  
  *numeric*  
  Stop time for analysis. Default is `72`.

- **interval_col**  
  *string*  
  Column name for time intervals. Default is `"timepoint"`.

- **nuisance_model_formula_str**  
  *string*  
  Nuisance model formula. Default is `"~1"`.

- **ctrl_ids**  
  *vector*  
  Control IDs. Default is `NULL`.

- **mt_ids**  
  *vector*  
  Mutant IDs. Default is `NULL`.

- **sparsity_factor**  
  *numeric*  
  Sparsity factor. Default is `0.01`.

- **perturbation_col**  
  *string*  
  Perturbation column name. Default is `"gene_target"`.

- **batch_col**  
  *string*  
  Batch column name. Default is `"expt"`.

- **max_num_cells**  
  *numeric*  
  Max number of cells to include. Default is `NULL`.

- **verbose**  
  *logical*  
  Print verbose output. Default is `FALSE`.

- **keep_ccs**  
  *logical*  
  Keep connected components. Default is `TRUE`.

- **num_threads**  
  *numeric*  
  Number of threads to use. Default is `1`.

- **backend**  
  *string*  
  Optimization backend. Default is `"nlopt"`.

- **q_val**  
  *numeric*  
  Q-value threshold for perturbation significance. Default is `0.1`.

- **vhat_method**  
  *string*  
  Method for variance estimation. Default is `"bootstrap"`.

- **num_bootstraps**  
  *numeric*  
  Number of bootstraps. Default is `10`.

- **newdata**  
  *tibble*  
  New data for predictions. Default is an empty tibble.

- **edge_allowlist**  
  *list*  
  Allowed edges in graph. Default is `NULL`.

- **min_lfc**  
  *numeric*  
  Minimum log fold change. Default is `0`.

- **links_between_components**  
  *character vector*  
  Types of links between components. Default is `c("ctp", "none", "strongest-pcor", "strong-pcor")`.

- **log_abund_detection_thresh**  
  *numeric*  
  Log abundance detection threshold. Default is `-5`.

- **batches_excluded_from_assembly**  
  *vector*  
  Batches to exclude. Default is empty.

- **component_col**  
  *string*  
  Column name for components. Default is `"partition"`.

- **embryo_size_factors**  
  *vector*  
  Size factors for embryos. Default is `NULL`.

## Value

A tibble containing the results of the partition assembly, including WT and MT graphs, perturbation effects, and state graph plots.

## Details

The function performs the following steps:

- Prepares the CDS by adding subassembly group and cell state information.  
- Fits a wild-type model and assembles a WT state transition graph.  
- Fits mutant models and assembles MT state transition graphs.  
- Assesses perturbation effects and constructs annotated graphs.  
- Handles errors gracefully and returns `NA` for failed steps.

## Examples

```r
results <- assemble_partition(
  cds = my_cds,
  sample_group = "sample",
  cell_group = "cell_type",
  partition_name = "partition_1",
  main_model_formula_str = "~timepoint",
  start_time = 18,
  stop_time = 72
)
```
