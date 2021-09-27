### * All functions in this file are exported

### * run_mcmc()

#' Run a MCMC sampler on a network model using Stan
#'
#' @param model A \code{networkModel}.
## #' @param dt,grid_size Either the time step size for trajectory calculations
## #'     (\code{dt}) or the number of points for the calculation
## #'     (\code{grid_size}) can be provided. If none is provided, then a default
## #'     grid size of 256 steps is used. Not used for now.
#' @param iter A positive integer specifying the number of iterations for each
#'     chain (including warmup). The default is 2000.
#' @param chains A positive integer specifying the number of Markov chains.
#'     The default is 4.
#' @param cores Number of cores to sue for parallel use. Default is
#'     \code{NULL}, which means to use the value stored in
#'     \code{options()[["mc.cores"]]} (or 1 if this value is not set).
#' @param stanfit If TRUE, returns a `stanfit` object instead of the more
#'     classical `mcmc.list` object.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#'
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#'
#' @export

run_mcmc <- function(model, iter = 2000, chains = 4, cores = NULL,
                     stanfit = FALSE, ...) {
    fit <- mugen_stan(nm = model, iter = iter, chains = chains,
                      cores = cores, stanfit = stanfit, ...)
    if (stanfit) {
        outClass <- "stanfit"
    } else {
        outClass <- c("networkModelStanfit", class(fit))
    }
    return(structure(fit, class = outClass))
}
