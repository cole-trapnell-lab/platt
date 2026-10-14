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
  gates = .EFDR_DEFAULT_GATES,
  min_replicate_embryos = .EFDR_MIN_REP_EMB,
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
  *integer*  
  Minimum number of perturbation pseudobulks (embryos) a cell type needs for a valid two-group contrast; below it every gene in that cell type is gated (`too_few_embryos`). Default `.EFDR_MIN_PB` = `2`.

- **gates**  
  *character vector*  
  Which uncallability gates to apply, from `.EFDR_GATES`. Default `.EFDR_DEFAULT_GATES` = `too_few_embryos`, `not_replicated`, `too_few_cells`. The previous policy's other gates (`not_detected_anywhere`, `gain_from_zero`, `gain_no_counts`, `loss_from_trace`, `flat_cell_floor`) are off by default because each is subsumed; pass them explicitly to reproduce that policy. Gates not listed are still computed and counted in the `efdr_stats` attribute, they just stop pinning `empirical_p` to `1`. Unknown names error.

- **min_replicate_embryos**  
  *integer*  
  Distinct embryos that must detect the gene in at least one arm for the `not_replicated` gate to pass it. Default `.EFDR_MIN_REP_EMB` = `2`. Needs `counts` to carry `n_emb_ctrl` / `n_emb_pert`; without them nothing is gated on this axis.

- **group_col**  
  Column to BH-adjust within. Default `"cell_group"`.

## Value

`deg_tbl` with added `empirical_p` and `empirical_fdr`.

## Details

The tail null is conditioned on **expression x thin-arm cell count** (not %embryos / global ratio), so a real depletion (huge `|z|`, deep arm) is kept while a thin-arm sampling artifact of the same `|z|` is demoted. Each call is matched on `min(n_ctrl_cells, n_pert_cells)` (the thinner arm) and the sign of the effect *in that thin arm*, which generalizes to perturbation-heavy designs (where the control arm is the thin one).

`empirical_p = max(p_ashr, p_tail)`: the empirical step may only **demote** a call the ordinary ashr test made — it can never turn a non-significant call significant. Filtering on `empirical_p < 0.05` is a drop-in, artifact-aware replacement for filtering on `perturb_to_ctrl_p_value < 0.05`.

`empirical_fdr` is **not** `BH(empirical_p)`. BH-adjusting the already-combined per-gene max would floor every gene's combined p near the permutation resolution, compressing the p-distribution so badly that BH crushes even genuine hits. Instead, `p_ashr` and `p_tail` are each BH-adjusted separately within `group_col`, and the max of the two resulting q-values is taken: `empirical_fdr = max(BH(p_ashr), BH(p_tail))`. A gene must still clear both criteria at the FDR level — the intersection of the two FDR-controlled call sets — but the ordinary-significance q, which is genuinely small for real hits, is no longer discarded by an up-front per-gene maximum. Genes gated to `empirical_p = 1` (see `gates` above) are forced to `empirical_fdr = 1` as well.

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
