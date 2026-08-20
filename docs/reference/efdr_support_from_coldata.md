# `efdr_support_from_coldata`

Per-cell-type perturbation-arm support, from a plain coldata table.

```r
efdr_support_from_coldata(
  coldata,
  sample_group = "embryo_ID",
  cell_group = "cell_type",
  control_ids = c("ctrl-inj"),
  perturbation_col = "perturbation"
)
```

## Arguments

- **coldata**  
  *data.frame*  
  A table with embryo, perturbation, and cell-type columns.

- **sample_group**, **cell_group**, **control_ids**, **perturbation_col**  
  Column names / control labels.

## Details

The single source of truth for every support quantity the empirical-FDR model needs, computed from a `(embryo, perturbation, cell_type)` table so it can be driven from either a `cell_data_set`'s `colData` ([`efdr_perturbation_support()`](efdr_perturbation_support)) or a lightweight coldata TSV (the mcclintock eFDR pipeline stage) without duplicating the definitions across repos:

- `n_pert_pb` — perturbation embryos with the cell type (the gate's replicate count).
- `n_ctrl_pb` — control embryos with the cell type.
- `n_eff_pb` — cell-weighted effective perturbation replicates; an embryo is a full replicate only once it carries enough cells, which drives the degrees-of-freedom correction.
- `n_pert_cells` / `n_ctrl_cells` — total cells per arm, feeding the theoretical SE floor.

## Examples

```r
support <- efdr_support_from_coldata(coldata_tbl,
  sample_group = "embryo",
  cell_group = "cell_type",
  control_ids = c("ctrl-inj"),
  perturbation_col = "gene_target"
)
```
