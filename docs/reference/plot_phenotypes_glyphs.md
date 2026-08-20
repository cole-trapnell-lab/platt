# `plot_phenotypes_glyphs`

Draws a lineage graph with phenotype-driven node fills, optional identity glyphs, fitness badges, and optional node/group labels.

```r
plot_phenotypes_glyphs(
  cell_state_graph,
  phenos_df,
  map = list(
    id = "cell_group",
    lfc = "abundance_log2fc",
    q = "abundance_q",
    ident = "identity_label",
    glyph = "identity_glyph",
    f1 = "F1_dir",
    f2 = "F2_apoptosis",
    f3 = "F3_stress_score",
    f4 = "F4_senescence"
  ),
  lfc_cap = 2,
  arrow_unit = 3,
  node_size = 2.2,
  con_colour = "darkgrey",
  legend_position = "none",
  label_cell_types = NULL,
  show_node_labels = FALSE,
  label_font_size = 3,
  label_font_size_pt = NULL,
  stress_cap = 2.5,
  filter_by_group = FALSE,
  cell_types = NULL,
  draw_group_boxes = TRUE,
  show_group_labels = FALSE,
  group_label_font_size = 2,
  node_overlay = c("none", "glyphs", "badges", "both"),
  glyph_color = NULL,
  badge_color = "black",
  badge_outline_color = "white"
)
```

## Arguments

- **cell_state_graph**  
  *cell_state_graph*  
  A platt `cell_state_graph` object with node coordinates and layout metadata.

- **phenos_df**  
  *data frame*  
  Phenotype table, typically from `impact_to_phenos()`.

- **map**  
  *named list*  
  Maps the expected phenotype roles (`lfc`, `q`, `ident`, `glyph`, `f1`-`f4`) to columns in `phenos_df`.

- **lfc_cap**  
  *numeric*  
  Cap applied to abundance proxy values for plotting.

- **arrow_unit**, **node_size**  
  *numeric*  
  Edge-arrow size (points) and base node-size scalar.

- **con_colour**  
  *character*  
  Connection/outline color.

- **label_cell_types**  
  *character vector*  
  Optional node names to label, or `"all"`.

- **show_node_labels**, **label_font_size**, **label_font_size_pt**  
  Label all nodes (when `label_cell_types` is `NULL`) and set label size, in ggplot units or points.

- **stress_cap**  
  *numeric*  
  Cap for stress-score alpha scaling.

- **filter_by_group**, **cell_types**  
  Restrict nodes/edges to a subset, optionally scoped to a whole grouping box.

- **draw_group_boxes**, **show_group_labels**, **group_label_font_size**  
  Draw grouping boxes and their labels.

- **node_overlay**  
  *character*  
  One of `"none"`, `"glyphs"`, `"badges"`, or `"both"` — which identity/fitness overlays to draw on top of each node. _As of the current `develop` source, only `"glyphs"` and `"none"` actually work (`match.arg()` errors on `"badges"`/`"both"`) — stick to those two until the badge overlay is finished._

- **glyph_color**, **badge_color**, **badge_outline_color**  
  *character*  
  Colors for the identity glyph text and fitness badge fill/outline.

## Value

A `ggplot` object.

## Examples

```r
plot_phenotypes_glyphs(muscle_state_graph,
  tbx16_phenos,
  cell_types = c("paraxial mesoderm (tbx16+)", "fast-committed myocyte, pre-fusion"),
  show_node_labels = TRUE
)
```
