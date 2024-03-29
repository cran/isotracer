% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ui_getters.R
\name{comps}
\alias{comps}
\title{Return the compartments of a network model}
\usage{
comps(nm)
}
\arguments{
\item{nm}{A \code{networkModel} object.}
}
\value{
A list of character vectors, with one list element per row of the
    input network model (list elements are in the same order as the input
    network model rows). Each list element containing the names of the
    compartments in the topology defined in the corresponding row of the
    input network model.
}
\description{
Return the compartments of a network model
}
\examples{
aquarium_mod
comps(aquarium_mod)

trini_mod
comps(trini_mod)

}
