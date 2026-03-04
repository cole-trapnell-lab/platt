% Manual documentation file
\name{plot_phenotypes_glyphs}
\alias{plot_phenotypes_glyphs}
\title{Plot lineage graph with phenotype overlays}
\usage{
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
  legend_position = "bottom",
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
}
\arguments{
\item{cell_state_graph}{Cell-state graph object with node coordinates and layout metadata.}

\item{phenos_df}{Phenotype table, typically produced by \code{impact_to_phenos()}.}

\item{map}{Named list mapping expected phenotype roles to columns in \code{phenos_df}.}

\item{lfc_cap}{Numeric cap applied to abundance proxy values for plotting.}

\item{arrow_unit}{Numeric edge-arrow size in points.}

\item{node_size}{Base node-size scalar.}

\item{con_colour}{Connection/outline color.}

\item{legend_position}{Legend position setting.}

\item{label_cell_types}{Optional character vector of node names to label, or \code{"all"}.}

\item{show_node_labels}{Logical; if \code{TRUE} and \code{label_cell_types} is \code{NULL}, labels all nodes.}

\item{label_font_size}{Numeric node-label size in ggplot/ggrepel units.}

\item{label_font_size_pt}{Optional node-label size in points; overrides \code{label_font_size}.}

\item{stress_cap}{Numeric cap for stress-score alpha scaling.}

\item{filter_by_group}{Logical; restrict nodes/edges to one grouping box when \code{cell_types} is provided.}

\item{cell_types}{Optional character vector of node names to retain.}

\item{draw_group_boxes}{Logical; draw grouping boxes.}

\item{show_group_labels}{Logical; draw grouping-box labels.}

\item{group_label_font_size}{Numeric size for grouping-box labels.}

\item{node_overlay}{One of \code{"none"}, \code{"glyphs"}, \code{"badges"}, or \code{"both"}.}

\item{glyph_color}{Optional fixed text color for identity glyphs.}

\item{badge_color}{Color used for fitness badge fills/strokes.}

\item{badge_outline_color}{Outline color for filled badge symbols.}
}
\value{
A \code{ggplot} object.
}
\description{
Draws a lineage graph with phenotype-driven node fills, optional identity glyphs,
fitness badges, and optional node/group labels.
}
