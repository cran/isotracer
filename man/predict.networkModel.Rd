% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui_predictions.R
\name{predict.networkModel}
\alias{predict.networkModel}
\title{Add a column with predictions from a fit}
\usage{
\method{predict}{networkModel}(
  object,
  fit,
  draws = NULL,
  error.draws = 5,
  probs = 0.95,
  cores = NULL,
  dt = NULL,
  grid_size = NULL,
  at = NULL,
  end = NULL,
  ...
)
}
\arguments{
\item{object}{Network model}

\item{fit}{Model fit (mcmc.list object)}

\item{draws}{Integer, number of draws from the posteriors}

\item{error.draws}{Integer, number of draws from the error distribution, for
a given posterior draw.}

\item{probs}{Credible interval (default 0.95).}

\item{cores}{Number of cores to use for parallel calculations. Default is
\code{NULL}, which means to use the value stored in
\code{options()[["mc.cores"]]} (or 1 if this value is not set).}

\item{dt, grid_size}{Time step size or grid points, respectively.}

\item{at}{Timepoints at which the predictions should be returned.}

\item{end}{Final timepoint used in the projections.}

\item{...}{Not used.}
}
\value{
A network model object with an added column \code{"prediction"}.
}
\description{
Add a column with predictions from a fit
}
