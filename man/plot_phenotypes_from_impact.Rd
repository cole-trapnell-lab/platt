% Manual documentation file
\name{plot_phenotypes_from_impact}
\alias{plot_phenotypes_from_impact}
\title{Plot phenotype glyphs from an impact table}
\usage{
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
}
\arguments{
\item{cell_state_graph}{Cell-state graph object with node coordinates and layout metadata.}

\item{impact_table}{Data frame/tibble with per-cell-type phenotype calls.}

\item{filter_by_group}{Logical; if \code{TRUE} and \code{cell_types} is set, restrict to the same grouping box.}

\item{cell_types}{Optional character vector of cell types to keep.}

\item{show_node_labels}{Logical; if \code{TRUE}, labels all nodes unless a subset is supplied via \code{label_cell_types} in \code{...}.}

\item{label_font_size}{Numeric node-label size in ggplot/ggrepel size units (mm-based).}

\item{label_font_size_pt}{Optional node-label size in points. If provided, overrides \code{label_font_size}.}

\item{show_group_labels}{Logical; draw text labels for grouping boxes.}

\item{group_label_font_size}{Numeric text size for grouping-box labels.}

\item{...}{Additional arguments forwarded to \code{plot_phenotypes_glyphs()}.}
}
\value{
A \code{ggplot} object.
}
\description{
Converts an impact table to phenotype glyph fields and renders the annotated lineage graph.
}
\examples{
\dontrun{
p <- plot_phenotypes_from_impact(
  cell_state_graph = cell_state_graph,
  impact_table = impact_table,
  show_node_labels = TRUE,
  label_font_size_pt = 14,
  show_group_labels = TRUE
)
}
}
