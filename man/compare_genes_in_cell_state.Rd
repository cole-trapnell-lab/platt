% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deg.R
\name{compare_genes_in_cell_state}
\alias{compare_genes_in_cell_state}
\title{Compare Genes in Cell State}
\usage{
compare_genes_in_cell_state(
  cell_state,
  state_graph,
  estimate_matrix,
  stderr_matrix,
  state_term = "cell_group",
  log_fc_thresh = 1,
  abs_expr_thresh = 0.001,
  sig_thresh = 0.05,
  cores = 1,
  expected_effect_mode_interval = c(-10, 10),
  cv_threshold = 10
)
}
\arguments{
\item{cell_state}{A character string representing the cell state to be analyzed.}

\item{state_graph}{An igraph object representing the state graph.}

\item{estimate_matrix}{A matrix of estimated gene expression values.}

\item{stderr_matrix}{A matrix of standard errors for the estimated gene expression values.}

\item{state_term}{A character string representing the state term. Default is "cell_group".}

\item{log_fc_thresh}{A numeric value representing the log fold change threshold. Default is 1.}

\item{abs_expr_thresh}{A numeric value representing the absolute expression threshold. Default is 1e-3.}

\item{sig_thresh}{A numeric value representing the significance threshold. Default is 0.05.}

\item{cores}{An integer representing the number of cores to use for parallel processing. Default is 1.}

\item{expected_effect_mode_interval}{A numeric vector representing the expected effect mode interval. Default is c(-10, 10).}

\item{cv_threshold}{A numeric value representing the coefficient of variation threshold. Default is 10.}
}
\value{
A tibble containing gene expression comparisons and interpretations.
}
\description{
This function compares gene expression levels in a given cell state with its parent, sibling, and child states.
}
\details{
The function performs the following steps:
\enumerate{
\item Identifies parent, sibling, and child states of the given cell state.
\item Computes gene expression levels and significance for the cell state and its related states.
\item Compares gene expression levels between the cell state and its parent, sibling, and child states.
\item Interprets the gene expression patterns based on the comparisons.
}
}
\examples{
\dontrun{
compare_genes_in_cell_state(
  cell_state = "state1",
  state_graph = my_state_graph,
  estimate_matrix = my_estimate_matrix,
  stderr_matrix = my_stderr_matrix
)
}

}
