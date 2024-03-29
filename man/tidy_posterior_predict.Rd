% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui_predictions.R
\name{tidy_posterior_predict}
\alias{tidy_posterior_predict}
\title{Draw from the posterior predictive distribution of the model outcome}
\usage{
tidy_posterior_predict(object, newdata, draw = NULL, cores = NULL, ...)
}
\arguments{
\item{object}{A networkModelStanfit object.}

\item{newdata}{The original model used to fit the networkStanfit object.}

\item{draw}{Integer, number of draws to sample from the posterior. Default
is 100.}

\item{cores}{Number of cores to use for parallel calculations. Default is
\code{NULL}, which means to use the value stored in
\code{options()[["mc.cores"]]} (or 1 if this value is not set).}

\item{...}{Not used for now.}
}
\value{
A tidy table.
}
\description{
Draw from the posterior predictive distribution of the model outcome
}
