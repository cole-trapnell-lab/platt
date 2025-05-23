% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assembly_utils.R
\name{assemble_mt_graph}
\alias{assemble_mt_graph}
\title{assembles a graph using the perturbation data}
\usage{
assemble_mt_graph(
  wt_ccm,
  perturb_models_tbl,
  interval_col = "timepoint",
  start_time = NULL,
  stop_time = NULL,
  interval_step = 2,
  links_between_components = "none",
  log_abund_detection_thresh = -5,
  q_val = 0.1,
  newdata = tibble(),
  break_cycles = TRUE,
  component_col = "partition",
  edge_allowlist = NULL,
  edge_denylist = NULL,
  verbose = FALSE
)
}
\description{
assembles a graph using the perturbation data
}
