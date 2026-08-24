# `plot_phenotype_counts`

Summarizes the number of distinct cell types with abundance, fitness, or identity phenotypes within each `projection_group`, and displays the counts as stacked bars.

```r
plot_phenotype_counts(impact_table, facet_by = NULL)
```

## Arguments

- **impact_table**  
  *data frame*  
  Must have at least `cell_type`, `projection_group`, `effect_type`, `abundance_code`, `fitness_label`, and `identity_label`. Any additional columns referenced in `facet_by` must also be present.

- **facet_by**  
  *character vector, optional*  
  Column names to facet across horizontally. `effect_type` is always used as the vertical facet.

## Value

A `ggplot` object showing stacked phenotype counts by `projection_group`.

## Details

This is a quick way to sanity-check an [impact table](https://cole-trapnell-lab.github.io/platt/phenotyping/) at a glance — e.g. before plotting it on the graph with `plot_phenotypes_from_impact()`, `plot_phenotype_counts()` can show whether a perturbation's abundance/fitness/identity calls are concentrated in one tissue or spread across many, and whether cell-autonomous and non-autonomous effects (`effect_type`) look different in scale.
