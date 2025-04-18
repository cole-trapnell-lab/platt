% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{get_max_abundance_contrast}
\alias{get_max_abundance_contrast}
\title{Get Maximum Abundance Contrast}
\usage{
get_max_abundance_contrast(
  ccm,
  start_time = NULL,
  stop_time = NULL,
  ctrl_abundances = NULL,
  interval_col = "timepoint",
  newdata = tibble()
)
}
\arguments{
\item{ccm}{A cell clustering model object containing cell clustering data.}

\item{start_time}{(Optional) The starting timepoint for the analysis. Defaults
to the minimum timepoint in the data if not provided.}

\item{stop_time}{(Optional) The ending timepoint for the analysis. Defaults
to the maximum timepoint in the data if not provided.}

\item{ctrl_abundances}{(Optional) A data frame of control abundances. If not
provided, control abundances will be calculated based on the specified time
interval.}

\item{interval_col}{A string specifying the column name in the cell clustering
model that contains the timepoint information. Defaults to "timepoint".}

\item{newdata}{A tibble containing additional data to be used in the analysis.
Defaults to an empty tibble.}
}
\value{
A data frame or tibble containing the maximum abundance contrast
results.
}
\description{
This function calculates the maximum abundance contrast for a given
cell clustering model (ccm) over a specified time interval. It allows for
customization of the start and stop times, control abundances, and additional
data inputs.
}
\details{
The function first determines the range of timepoints to analyze
based on the provided or default start and stop times. It then calculates
perturbation effects and filters control abundances to match the timepoints
present in the perturbation effects.
}
\examples{
# Example usage:
# result <- get_max_abundance_contrast(ccm, start_time = 0, stop_time = 10)

}
