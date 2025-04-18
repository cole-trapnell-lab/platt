% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assembly_utils.R
\name{fit_subset_genotype_ccm}
\alias{fit_subset_genotype_ccm}
\title{Fit a Subset of Genotype Cell Cycle Models (CCM)}
\usage{
fit_subset_genotype_ccm(ccm, umap_space = NULL, ...)
}
\arguments{
\item{ccm}{A CCM object containing the cell cycle model and associated metadata.}

\item{umap_space}{A character string specifying the UMAP space to switch to.
If \code{NULL}, the function attempts to retrieve the UMAP space from the CCM metadata.}

\item{...}{Additional arguments passed to the \code{subset_ccs} function for subsetting CCS.}
}
\value{
A CCM object fitted to the subset of genotype data.
}
\description{
This function fits a subset of genotype cell cycle models (CCM) based on the
provided parameters. It allows for switching to a specified UMAP space and
subsetting the cell cycle states (CCS) before fitting the genotype CCM.
}
\details{
\itemize{
\item If \code{umap_space} is not provided, the function tries to retrieve it from
\code{ccm@ccs@cds@metadata$umap_space}.
\item If a valid UMAP space is found, the function switches the CCM to that space
using \code{switch_ccm_space}.
\item The function subsets the CCS using \code{subset_ccs} and fits the genotype CCM
using \code{fit_genotype_ccm}.
}
}
\examples{
# Example usage:
# Assuming `ccm` is a valid CCM object:
sub_ccm <- fit_subset_genotype_ccm(ccm, umap_space = "UMAP_1", some_filter = TRUE)

}
\seealso{
\itemize{
\item \code{switch_ccm_space} for switching the CCM to a specified UMAP space.
\item \code{subset_ccs} for subsetting cell cycle states.
\item \code{fit_genotype_ccm} for fitting genotype cell cycle models.
}
}
