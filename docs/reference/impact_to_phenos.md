# `impact_to_phenos`

Transforms raw impact calls into the standardized fields the phenotype-plotting functions expect.

```r
impact_to_phenos(
  impact_table,
  sev_map = c(none = 0.25, mild = 0.6, moderate = 1.2, severe = 2),
  abundance_code_map = c(
    "A0 No change" = 0,
    "A1 Expansion" = 1,
    "A2 Depletion" = -1,
    "A3 Ablation/Loss" = -1.5,
    "A4 Ectopic/extra state" = 1.2
  ),
  identity_glyph_map = c(
    "I0 Identity intact" = "",
    "I1 Maturation delay" = "<<",
    "I2 Precocious maturation" = ">>",
    "I3 Program failure within identity" = "!!",
    "I4 Fate switch / misspecification" = "!!",
    "I5 Identity fragmentation" = ""
  )
)
```

## Arguments

- **impact_table**  
  *data frame*  
  Must have columns `cell_type`, `abundance_code`, `abundance_severity`, `identity_label`, and `fitness_label`.

- **sev_map**  
  *named numeric vector*  
  Maps abundance severity labels to a magnitude, used to build the abundance proxy.

- **abundance_code_map**  
  *named numeric vector*  
  Maps abundance class labels to a signed direction/magnitude.

- **identity_glyph_map**  
  *named character vector*  
  Maps identity labels to the glyph string drawn on the plot.

## Value

A tibble ready for `plot_phenotypes_glyphs()`, including abundance proxy values, identity glyph fields, fitness axes, and optional dysregulated-gene summaries.

## Examples

```r
tbx16_phenos <- impact_to_phenos(tbx16_impact_table)
```
