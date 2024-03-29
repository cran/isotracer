% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui_model-building.R
\name{set_params}
\alias{set_params}
\title{Set the parameters in a network model}
\usage{
set_params(nm, params, force = TRUE, quick = FALSE)
}
\arguments{
\item{nm}{A \code{networkModel} object.}

\item{params}{A named vector or a tibble with columns c("parameter",
"value") containing the (global) parameter values.}

\item{force}{Boolean, if FALSE will not overwrite already set parameters.}

\item{quick}{Boolean, if TRUE take some shortcuts for faster parameter
settings when called by another function. This should usually be left to
the default (FALSE) by a regular package user.}
}
\value{
A \code{networkModel} object.
}
\description{
Set the parameters in a network model
}
\examples{
m <- aquarium_mod
p <- sample_params(m)
m2 <- set_params(m, p)
m2$parameters

}
