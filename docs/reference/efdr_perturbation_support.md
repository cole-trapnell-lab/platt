# `efdr_perturbation_support`

Per-cell-type perturbation-arm replicate support, computed directly from a `cell_data_set`.

```r
efdr_perturbation_support(
  cds,
  sample_group = "embryo_ID",
  cell_group = "cell_type",
  control_ids = c("ctrl-inj"),
  perturbation_col = "perturbation"
)
```

## Arguments

- **cds**  
  *cell_data_set*  
  A contrast `cell_data_set`.

- **sample_group**, **cell_group**, **control_ids**, **perturbation_col**  
  `colData` columns/values identifying replicate, cell type, control labels, and the perturbation column.

## Value

`tibble(cell_group, n_pert_pb, n_ctrl_pb, n_eff_pb, n_pert_cells, n_ctrl_cells)`.

## Details

Counts, per cell type, the number of **perturbation** embryos (pseudobulks) contributing cells to the within-state contrast — this is the replicate count the pseudobulk GLM actually has on the perturbation side, and it's the honest degrees-of-freedom for the contrast. When a cell type is depleted in the perturbation, the deep control arm supplies a small standard error as if the contrast were well-powered, inflating `|z|` and manufacturing spurious (predominantly down) calls. [`annotate_empirical_fdr_model()`](annotate_empirical_fdr_model) uses this to floor the ordinary p-value by a t-distribution with `n_pert_pb - 1` degrees of freedom, and to gate cell types with fewer than 2 perturbation pseudobulks (no valid two-group contrast).

_This is the `cell_data_set`-driven convenience wrapper around [`efdr_support_from_coldata()`](efdr_support_from_coldata)._

## Examples

```r
support <- efdr_perturbation_support(tbx16_msgn1_cds,
  sample_group = "embryo",
  cell_group = "cell_type",
  control_ids = c("ctrl-inj"),
  perturbation_col = "gene_target"
)
```
