% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui_predictions.R
\name{tidy_data}
\alias{tidy_data}
\title{Extract data from a networkModel object into a tidy tibble.}
\usage{
tidy_data(x)
}
\arguments{
\item{x}{A networkModel object.}
}
\value{
A tibble (note: row ordering is not the same as in the input).
}
\description{
Extract data from a networkModel object into a tidy tibble.
}
\examples{
tidy_data(aquarium_mod)
tidy_data(trini_mod)

}
