% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assembly_utils.R
\name{fit_wt_model}
\alias{fit_wt_model}
\title{Fit Wild Type Model}
\usage{
fit_wt_model(
  cds,
  sample_group,
  cell_group,
  main_model_formula_str = NULL,
  num_time_breaks = 4,
  nuisance_model_formula_str = "~1",
  ctrl_ids = NULL,
  sparsity_factor = 1,
  vhat_method = "bootstrap",
  interval_col = "timepoint",
  perturbation_col = "knockout",
  batch_col = "expt",
  start_time = NULL,
  stop_time = NULL,
  interval_step = 2,
  log_abund_detection_thresh = -5,
  keep_ccs = TRUE,
  q_val = 0.1,
  edge_allowlist = NULL,
  edge_denylist = NULL,
  base_penalty = 1,
  keep_cds = TRUE,
  verbose = FALSE,
  num_threads = 1,
  backend = "nlopt",
  penalize_by_distance = TRUE,
  embryo_size_factors = NULL,
  batches_excluded_from_assembly = c(),
  ...
)
}
\arguments{
\item{cds}{A cell dataset (CDS) object containing cell data.}

\item{sample_group}{A string specifying the column in \code{colData(cds)} that defines sample groups.}

\item{cell_group}{A string specifying the column in \code{colData(cds)} that defines cell groups.}

\item{main_model_formula_str}{A string specifying the main model formula. If \code{NULL}, it will be generated automatically.}

\item{num_time_breaks}{An integer specifying the number of time breaks for the main model formula. Default is 4.}

\item{nuisance_model_formula_str}{A string specifying the nuisance model formula. Default is \code{"~1"}.}

\item{ctrl_ids}{A vector of control IDs. If \code{NULL}, control IDs are inferred from \code{perturbation_col}.}

\item{sparsity_factor}{A numeric value for sparsity factor used in model selection. Default is 1.}

\item{vhat_method}{A string specifying the method for estimating variance (\code{"bootstrap"} by default).}

\item{interval_col}{A string specifying the column in \code{colData(cds)} that defines time intervals. Default is \code{"timepoint"}.}

\item{perturbation_col}{A string specifying the column in \code{colData(cds)} that defines perturbation groups. Default is \code{"knockout"}.}

\item{batch_col}{A string specifying the column in \code{colData(cds)} that defines batch groups. Default is \code{"expt"}.}

\item{start_time}{A numeric value specifying the start time for the model. If \code{NULL}, it is inferred from the data.}

\item{stop_time}{A numeric value specifying the stop time for the model. If \code{NULL}, it is inferred from the data.}

\item{interval_step}{A numeric value specifying the step size for time intervals. Default is 2.}

\item{log_abund_detection_thresh}{A numeric threshold for log abundance detection. Default is -5.}

\item{keep_ccs}{A logical value indicating whether to retain the cell count set (CCS). Default is \code{TRUE}.}

\item{q_val}{A numeric value specifying the q-value threshold for significance. Default is 0.1.}

\item{edge_allowlist}{A vector of edges to allow in the model. Default is \code{NULL}.}

\item{edge_denylist}{A vector of edges to deny in the model. Default is \code{NULL}.}

\item{base_penalty}{A numeric value specifying the base penalty for model selection. Default is 1.}

\item{keep_cds}{A logical value indicating whether to retain the CDS. Default is \code{TRUE}.}

\item{verbose}{A logical value indicating whether to print verbose messages. Default is \code{FALSE}.}

\item{num_threads}{An integer specifying the number of threads to use. Default is 1.}

\item{backend}{A string specifying the backend for model fitting. Default is \code{"nlopt"}.}

\item{penalize_by_distance}{A logical value indicating whether to penalize by distance. Default is \code{TRUE}.}

\item{embryo_size_factors}{A named vector of size factors for embryos. Default is \code{NULL}.}

\item{batches_excluded_from_assembly}{A vector of batch IDs to exclude from assembly. Default is an empty vector.}

\item{...}{Additional arguments passed to the underlying model fitting functions.}
}
\value{
A fitted wild type cell count model object, or \code{NULL} if no control cells are available.
}
\description{
This function fits a wild type (WT) model to a cell dataset (CDS) by estimating
cell count dynamics over time and accounting for nuisance variables. It supports
various customization options for model fitting, including user-defined formulas,
size factors, and control IDs.
}
\details{
The function first subsets the CDS to include only control cells based on the \code{perturbation_col}
and \code{ctrl_ids}. It then constructs a cell count set (CCS) and fits a model using the specified
formulas and parameters. If no control cells are available or only a single cell group is present,
the function returns \code{NULL}.
}
\examples{
# Example usage:
wt_model <- fit_wt_model(
  cds = my_cds,
  sample_group = "sample",
  cell_group = "cell_type",
  interval_col = "timepoint",
  perturbation_col = "knockout",
  batch_col = "batch"
)
}
