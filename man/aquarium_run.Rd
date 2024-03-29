% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/datasets_aquarium_mod.R
\docType{data}
\name{aquarium_run}
\alias{aquarium_run}
\title{An MCMC run from a simple aquarium network model}
\format{
An object of class \code{networkModelStanfit} (inherits from \code{mcmc.list}) of length 4.
}
\usage{
aquarium_run
}
\description{
This is an MCMC run on \code{\link{aquarium_mod}}. The code used to run the
MCMC is: \code{aquarium_run <- run_mcmc(aquarium_mod, thin = 4)} (note that
\code{thin = 4} was only used here to reduce the size of the data file
shipped with the package, but for a real-life analysis keeping the default
\code{thin = 1} is usually recommended). The code used to build the model
itself is given in the help page for \code{\link{aquarium_mod}}.
}
\examples{
\dontrun{
plot(aquarium_run)
summary(aquarium_run)
}
}
\keyword{datasets}
