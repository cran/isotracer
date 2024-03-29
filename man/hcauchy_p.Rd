% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui_priors.R
\name{hcauchy_p}
\alias{hcauchy_p}
\title{Define a half-Cauchy prior (on [0;+Inf])}
\usage{
hcauchy_p(scale)
}
\arguments{
\item{scale}{Median of the half-Cauchy distribution.}
}
\value{
A list defining the prior.
}
\description{
Define a half-Cauchy prior (on [0;+Inf])
}
\examples{
hcauchy_p(scale = 0.5)

}
