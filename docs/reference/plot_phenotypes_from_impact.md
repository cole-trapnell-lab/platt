# `plot_phenotypes_from_impact`

Converts an impact table to phenotype glyph fields with `impact_to_phenos()` and renders the annotated lineage graph, in one call.

```r
plot_phenotypes_from_impact(
  cell_state_graph,
  impact_table,
  filter_by_group = FALSE,
  cell_types = NULL,
  show_node_labels = FALSE,
  label_font_size = 3,
  label_font_size_pt = NULL,
  show_group_labels = FALSE,
  group_label_font_size = 2,
  ...
)
```

## Arguments

- **cell_state_graph**  
  *cell_state_graph*  
  A platt `cell_state_graph` object with node coordinates and layout metadata.

- **impact_table**  
  *data frame*  
  Per-cell-type phenotype calls (see `impact_to_phenos()`).

- **filter_by_group**  
  *logical*  
  If `TRUE` and `cell_types` is set, restrict to the same grouping box.

- **cell_types**  
  *character vector*  
  Optional subset of cell types to keep.

- **show_node_labels**  
  *logical*  
  Label all nodes, unless a subset is given via `label_cell_types` in `...`.

- **label_font_size**, **label_font_size_pt**  
  *numeric*  
  Node-label size in ggplot/ggrepel units, or in points if `label_font_size_pt` is set (overrides the former).

- **show_group_labels**, **group_label_font_size**  
  Draw and size text labels for grouping boxes.

- **...**  
  Additional arguments forwarded to `plot_phenotypes_glyphs()`.

## Value

A `ggplot` object.

## Examples

```r
plot_phenotypes_from_impact(muscle_state_graph,
  tbx16_impact_table,
  show_node_labels = TRUE
)
```
