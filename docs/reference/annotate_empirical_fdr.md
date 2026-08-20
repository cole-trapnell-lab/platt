# `annotate_empirical_fdr`

Annotates a within-state DEG table with an empirical FDR for the low-expression control-sampling artifact.

```r
annotate_empirical_fdr(
  deg_tbl,
  abundances = NULL,
  unaffected_cell_types = NULL,
  sig_thresh = 0.05,
  n_abund_bins = 4L,
  expr_breaks = c(-Inf, -5, -4, -3.5, -3, -2.5, -2, -1, Inf),
  lower_envelope_frac = 1 / 3
)
```

## Arguments

- **deg_tbl**  
  *data.frame / tibble*  
  Within-state DEG results. Must contain `cell_group`, `log_mean_expression`, `perturb_to_ctrl_p_value` (the ashr local false sign rate used as the significance quantity), and `perturb_to_ctrl_shrunken_lfc`.

- **abundances**  
  *named numeric vector or data.frame*  
  Optional per-cell-type abundance used to stratify the null. Either a named numeric vector (`cell_group` -> total cells) or a data.frame with columns `cell_group` and `abund`. If `NULL`, a depth proxy is derived from `mean_log_sf` (mean cells per pseudobulk), which must then be present in `deg_tbl`. Supplying real cell counts is more accurate and recommended.

- **unaffected_cell_types**  
  *character vector*  
  Optional `cell_group` names known to have no real perturbation effect (e.g. abundance "A0" and no fitness/identity phenotype), used as the self-null training set. If `NULL`, a lower-envelope fallback is used: within each abundance bin the cell types in the bottom tertile of overall call rate are treated as unaffected.

- **sig_thresh**  
  *numeric*  
  Significance threshold defining a call. Default `0.05`.

- **n_abund_bins**  
  *integer*  
  Number of abundance strata. Default `4`.

- **expr_breaks**  
  *numeric vector*  
  Breakpoints used to bin `log_mean_expression`.

- **lower_envelope_frac**  
  *numeric*  
  Fraction of each abundance bin's cell types treated as the self-null when `unaffected_cell_types` is `NULL`. Default `1/3`.

## Value

`deg_tbl` with three added columns joined per (`cell_group`, expression-bin, direction): `empirical_null_rate` (percent of admitted genes called in the matched unaffected cell types), `observed_rate` (percent called in this cell type), and `empirical_fdr` (`min(1, empirical_null_rate / observed_rate)`). Rows in strata with no unaffected estimate receive `empirical_fdr = 0` (i.e. not flagged), the conservative-for-discovery default.

## Details

In control-heavy perturbation designs the deep control arm detects lowly expressed genes that the shallow perturbation arm misses, producing a large excess of low-expression, "on-in-control" (down) DEG calls that are sampling artifacts rather than biology. This function estimates, per (cell_group x expression-bin x direction), the rate of calls attributable to that artifact and returns the input rows annotated with the null/observed rates and the resulting FDR.

The null is read directly off the perturbation's own unaffected cell types (perturbation biology only *adds* calls, so cell types with no real effect are self-null), matched to each cell type by abundance so that shallow cell types are compared against shallow ones. No permutation, model fitting, or re-run is required, so existing DEG tables can be decorated in place.

Estimation is direction-specific because the artifact is directional: in control-heavy designs it inflates down (on-in-control) calls, while up calls at high expression are largely genuine. Consumers can threshold `empirical_fdr` (e.g. drop calls above some cutoff) to suppress the artifact while retaining real signal.

_Pipeline-computed DEG tables (as in the [Perturbation DEGs page](../perturbation_degs/)) are typically already annotated with `empirical_p`/`empirical_fdr` via a more elaborate, arm-size-conditioned model — see [`annotate_empirical_fdr_model()`](../annotate_empirical_fdr_model). `annotate_empirical_fdr()` is the simpler, self-contained version of the same idea for decorating your own DEG tables without training a model first._

## Examples

```r
genes_within_cell_state <- annotate_empirical_fdr(genes_within_cell_state)

genes_within_cell_state %>% dplyr::filter(empirical_fdr < 0.1)
```
