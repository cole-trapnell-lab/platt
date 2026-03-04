% Manual documentation file
\name{impact_to_phenos}
\alias{impact_to_phenos}
\title{Convert impact-table fields to phenotype glyph fields}
\usage{
impact_to_phenos(
  impact_table,
  sev_map = c(none = 0.25, mild = 0.6, moderate = 1.2, severe = 2),
  abundance_code_map = c(
    "A0 No change" = 0,
    "A1 Expansion" = 1,
    "A2 Depletion" = -1,
    "A3 Ablation/Loss" = -1.5,
    "A4 Ectopic/extra state" = 1.2
  ),
  identity_glyph_map = c(
    "I0 Identity intact" = "",
    "I1 Maturation delay" = "<<",
    "I2 Precocious maturation" = ">>",
    "I3 Program failure within identity" = "!!",
    "I4 Fate switch / misspecification" = "!!",
    "I5 Identity fragmentation" = ""
  )
)
}
\arguments{
\item{impact_table}{Data frame with required columns \code{cell_type}, \code{abundance_code}, \code{abundance_severity}, \code{identity_label}, and \code{fitness_label}.}

\item{sev_map}{Named numeric vector mapping abundance severity labels to magnitude.}

\item{abundance_code_map}{Named numeric vector mapping abundance class labels to signed direction/magnitude.}

\item{identity_glyph_map}{Named character vector mapping identity labels to glyph strings.}
}
\value{
A tibble ready for \code{plot_phenotypes_glyphs()}, including abundance proxy values,
identity glyph fields, fitness axes, and optional dysregulated-gene summaries.
}
\description{
Transforms raw impact calls into standardized phenotype fields used by the plotting layer.
}
