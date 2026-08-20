# `train_efdr_model`

Trains the empirical-FDR null model from control-split null draws.

```r
train_efdr_model(
  null,
  detection = NULL,
  taus = .EFDR_TAUS,
  egrid = .EFDR_EGRID,
  lcgrid = .EFDR_LCGRID,
  we = 0.75,
  wc = 0.4,
  min_n = 40L,
  max_we = 2.5
)
```

## Arguments

- **null**  
  *tibble*  
  Output of [`build_empirical_null()`](build_empirical_null) (`gene_short_name`, `cell_group`, `log_mean_expression`, `z`, `n_pert_cells`).

- **detection**  
  Unused (kept for back-compat with callers).

- **taus**  
  Quantile levels (tail-concentrated) the surface is evaluated at.

- **egrid**, **lcgrid**  
  Expression and `log10(thin-arm cells + 1)` grids evaluated on.

- **we**, **wc**  
  Half-widths of the local window (expression, log10-cells) the tail quantiles are estimated over; widened adaptively until `min_n` draws are in.

- **min_n**  
  Minimum null draws required in a window before its quantiles are used.

## Value

An `efdr_model`: per-direction 2-D tail surfaces on `(egrid x lcgrid)`.

## Details

Models the null upper tail of `|z|` as a local-empirical surface over **expression and the thin-arm cell count**, per direction:

- **expression**: the dispersion mis-pricing that inflates `|z|` worsens as expression → 0, so the tail stays heavy in the sparse extreme-low region rather than drooping back toward zero.
- **thin-arm cell count**: a thin arm hits exact/near-zero counts which, combined with the trend-dispersion-shrunk SE, manufactures large `|z|` under no real effect — this is the axis that separates a real depletion (huge `|z|`, deep arm) from a thin-arm sampling artifact.

Draws are pooled across the sampling grid: the grid exists only to span the arm-size axis, read directly off `n_pert_cells`. Direction is kept separate (`dn` = the gene is lower in the thin arm, `z < 0`; `up` = higher).

## Examples

```r
efdr_model <- train_efdr_model(null_draws)
```
