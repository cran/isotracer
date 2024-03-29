% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui_priors.R
\name{gamma_p}
\alias{gamma_p}
\title{Define a gamma prior}
\usage{
gamma_p(alpha, beta)
}
\arguments{
\item{alpha}{Shape parameter (equivalent to the \code{shape} parameter of
R's \code{rgamma}).}

\item{beta}{Rate parameter (equivalent to the \code{rate} parameter of R's
\code{rgamma}).}
}
\value{
A list defining the prior.
}
\description{
Note the name of the function to define a prior (\code{gamma_p}), in order
to avoid confusion with the R mathematical function \code{gamma}.
}
\examples{
gamma_p(9, 2)
hist(sample_from_prior(gamma_p(9, 2), 1e3))

}
