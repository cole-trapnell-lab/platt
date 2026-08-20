# `annotate_empirical_fdr_model`

Annotates a DEG table with per-gene empirical p and BH FDR, using a trained `efdr_model`.

```r
annotate_empirical_fdr_model(
  model,
  deg_tbl,
  log_ratio = NULL,
  detection = NULL,
  support = NULL,
  counts = NULL,
  arm_cells = NULL,
  min_pseudobulks = .EFDR_MIN_PB,
  group_col = "cell_group"
)
```

## Arguments

- **model**  
  *efdr_model*  
  An `efdr_model` from [`train_efdr_model()`](train_efdr_model).

- **deg_tbl**  
  *data.frame*  
  DEG table with `gene_short_name`, `cell_group`, `log_mean_expression`, and either `z` or `perturb_to_ctrl_shrunken_lfc` + `perturb_to_ctrl_shrunken_lfc_se`.

- **log_ratio**  
  Deprecated / ignored (the tail is conditioned on the actual per-cell-type arm cell count, not a single global sampling ratio).

- **detection**  
  Deprecated / unused (the %embryos covariate is retired).

- **support**  
  Per-cell-type embryo counts (`cell_group`, `n_pert_pb`) for the <2-pseudobulk gate — see [`efdr_perturbation_support()`](efdr_perturbation_support).

- **counts**  
  Per-(gene, cell_type) arm UMI totals (`id`/`gene_short_name`, `cell_group`, `K_ctrl`, `K_pert`) for the gene-absent and control-absent gates.

- **arm_cells**  
  Per-cell-type arm cell counts (data.frame with `cell_group`, `n_ctrl_cells`, `n_pert_cells`). Supplies the arm-size covariate the tail is conditioned on. Genes in cell types with no match fall back to ashr only.

- **min_pseudobulks**  
  Minimum perturbation pseudobulks required for a valid two-group contrast.

- **group_col**  
  Column to BH-adjust within. Default `"cell_group"`.

## Value

`deg_tbl` with added `empirical_p` and `empirical_fdr`.

## Details

The tail null is conditioned on **expression x thin-arm cell count** (not %embryos / global ratio), so a real depletion (huge `|z|`, deep arm) is kept while a thin-arm sampling artifact of the same `|z|` is demoted. Each call is matched on `min(n_ctrl_cells, n_pert_cells)` (the thinner arm) and the sign of the effect *in that thin arm*, which generalizes to perturbation-heavy designs (where the control arm is the thin one).

`empirical_p = max(p_ashr, p_tail)`: the empirical step may only **demote** a call the ordinary ashr test made — it can never turn a non-significant call significant. Filtering on `empirical_p < 0.05` is a drop-in, artifact-aware replacement for filtering on `perturb_to_ctrl_p_value < 0.05`.

## Examples

```r
support <- efdr_perturbation_support(tbx16_msgn1_cds, sample_group = "embryo", perturbation_col = "gene_target")
null_draws <- build_empirical_null(ctrl_only_cds, sample_group = "embryo", perturbation_col = "gene_target")
efdr_model <- train_efdr_model(null_draws)

genes_within_cell_state <- annotate_empirical_fdr_model(
  model = efdr_model,
  deg_tbl = genes_within_cell_state,
  support = support,
  arm_cells = support
)
```
