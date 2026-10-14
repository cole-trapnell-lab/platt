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
  abundance_q_cut = 0.05,
  lfc_cap = 2,
  arrow_unit = 3,
  arrow_gap = 0,
  node_size = 2.2,
  con_colour = "darkgrey",
  legend_position = "right",
  legend_scale = 1,
  non_autonomous_alpha = 0.6,
  label_cell_types = NULL,
  show_node_labels = FALSE,
  label_font_size = 3,
  label_font_size_pt = NULL,
  stress_cap = 2.5,
  filter_by_group = FALSE,
  cell_types = NULL,
  draw_group_boxes = TRUE,
  group_box_linewidth = 0.25,
  show_group_labels = FALSE,
  group_label_font_size = 2,
  node_overlay = c("glyphs", "none"),
  glyph_color = "white",
  badge_color = "black",
  badge_outline_color = "white",
  node_outline_color = "black",
  render_mode = c("tissue", "global"),
  tooltip_builder = NULL,
  interactive = FALSE,
  width = NULL,
  height = NULL
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

- **abundance_q_cut**  
  *numeric*  
  Numeric q-value cutoff used to decide whether an abundance call survives into the glyph. Must match the `abundance_q_cut` that produced the codes in `phenos_df` (see `assign_abundance_code()`): a looser value re-admits rows the cascade rejected, a stricter one silently rewrites called rows to `"A0 No change"`. Default `0.05`.

- **lfc_cap**  
  *numeric*  
  Cap applied to abundance proxy values for plotting.

- **arrow_unit**, **node_size**  
  *numeric*  
  Edge-arrow size (points) and base node-size scalar.

- **arrow_gap**  
  *numeric*  
  Distance (in layout units) to pull each arrowhead back from its target node so it does not overlap the node fill. Default `0`.

- **con_colour**  
  *character*  
  Connection/outline color. Default `"darkgrey"`.

- **legend_position**  
  *character*  
  ggplot legend position. Default `"right"`.

- **legend_scale**  
  *numeric*  
  Multiplier applied to legend symbol, key, and text sizing; above `1` enlarges the legend, below `1` shrinks it. Default `1`.

- **non_autonomous_alpha**  
  *numeric*  
  Alpha used for nodes with `effect_type == "Non-autonomous"`; cell-autonomous nodes stay fully opaque. Default `0.6`.

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

- **draw_group_boxes**, **group_box_linewidth**, **show_group_labels**, **group_label_font_size**  
  Draw grouping boxes (with the given line width, default `0.25`) and their labels.

- **node_overlay**  
  *character*  
  One of `"glyphs"` (default) or `"none"`: whether to draw the identity glyph overlay on top of each node.

- **glyph_color**, **badge_color**, **badge_outline_color**, **node_outline_color**  
  *character*  
  Colors for the identity glyph text (default `"white"`), the fitness badge fill/outline, and the node outline (default `"black"`).

- **render_mode**  
  *character*  
  One of `"tissue"` (default) or `"global"`. Global mode hides the glyph/badge overlays and can draw a centre marker for identity phenotypes.

- **interactive**  
  *logical*  
  If `TRUE`, render with ggiraph so nodes and labels carry hover tooltips, and return a `girafe` object instead of a `ggplot`. Default `FALSE`.

- **tooltip_builder**  
  *function*  
  Optional `function(g, render_mode)` returning one tooltip string per plotted node, used when `interactive = TRUE`. Default `NULL`, which uses the built-in phenotype tooltip.

- **width**  
  *numeric*  
  Output width in inches for the interactive rendering. Default `NULL`.

- **height**  
  *numeric*  
  Output height in inches for the interactive rendering. When `width` and `height` are both `NULL` (the default) the girafe widget is sized automatically with equal coordinates; when both are given the widget is sized to them plus room for the legend.

## Value

A `ggplot` object, or a `girafe` htmlwidget when `interactive = TRUE`.

## Examples

```r
plot_phenotypes_glyphs(muscle_state_graph,
  tbx16_phenos,
  cell_types = c("paraxial mesoderm (tbx16+)", "fast-committed myocyte, pre-fusion"),
  show_node_labels = TRUE
)
```
