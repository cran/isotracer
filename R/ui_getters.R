### * All functions in this file are exported

### * topo()

#' Return the list of topologies, or a unique topology if all identical
#'
#' @param nm A \code{networkModel} object.
#' @param simplify Boolean, return only a unique topology if all topologies are
#'     identical or if there is only one? Default is TRUE.
#'
#' @return A list of the \code{networkModel} topologies or, if all topologies
#'     are identical (or if there is only one) and \code{simplify} is TRUE, a
#'     single topology (not wrapped into a single-element list).
#' 
#' @export

topo <- function(nm, simplify = TRUE) {
    out <- nm[["topology"]]
    if (simplify && length(unique(out)) == 1) {
        return(unique(out)[[1]])
    }
    return(out)
}

### * prop_family()

#' Return the distribution family for observed proportions
#'
#' @param nm A \code{networkModel} object.
#'
#' @return A character string describing the distribution family used to model
#'     observed proportions
#'
#' @examples
#' prop_family(aquarium_mod)
#' prop_family(trini_mod)
#' 
#' @export

prop_family <- function(nm) {
    attr(nm, "prop_family")
}

### * groups() method for networkModel

#' Get the grouping for a \code{networkModel} object
#'
#' @param x A \code{networkModel} object.
#'
#' @return A tibble giving the grouping variable(s) for the input network
#'     model. This tibble is in the same order as the rows of the input network
#'     model. If the input network model did not have any grouping variable,
#'     returns \code{NULL}.
#' 
#' @importFrom dplyr groups
#' @method groups networkModel
#'
#' @examples
#' groups(aquarium_mod)
#' groups(trini_mod)
#' 
#' @export

groups.networkModel <- function(x) {
    nm_get_groups(x, error = FALSE)
}

### * priors()

#' Return the tibble containing the priors of a networkModel
#'
#' @param nm A \code{networkModel} object.
#' @param fix_set_params If TRUE, parameters for which a value is set are given a
#'     fixed value (i.e. their prior is equivalent to a point value).
#' @param quiet Boolean to control verbosity.
#'
#' @return A tibble giving the current priors defined for the input network
#'     model.
#'
#' @examples
#' priors(aquarium_mod)
#' priors(trini_mod)
#' 
#' @export

priors <- function(nm, fix_set_params = FALSE, quiet = FALSE) {
    if (fix_set_params) {
        set_params <- attr(nm, "parameterValues")
        if (!is.null(set_params)) {
            for (i in seq_len(nrow(set_params))) {
                myPrior <- constant(value = set_params[["value"]][i])
                nm <- set_prior(nm, myPrior, param = set_params[["parameter"]][i],
                                use_regexp = FALSE, quiet = quiet)
            }
        }
    }
    out <- attr(nm, "priors")
    return(out)
}

### * params()

#' Return the parameters of a network model
#'
#' @param nm A \code{networkModel} object.
#' @param tibble If TRUE, return a tidy tibble with the parameter mapping.
#' @param expand If TRUE, return a detailed tibble. The "tibble" argument is
#'     disregarded in this case.
#'
#' @return A sorted character vector containing the parameters names.
#'
#' @examples
#' params(aquarium_mod)
#' params(trini_mod)
#' 
#' @export

params <- function(nm, tibble = FALSE, expand = FALSE) {
    stopifnot("parameters" %in% colnames(nm))
    params <- dplyr::bind_rows(nm[["parameters"]])
    if (expand) {
        if (is.null(nm[["group"]])) {
            stopifnot(nrow(nm) == 1)
            return(nm[["parameters"]])
        }
        z <- nm[, c("group", "parameters")]
        z <- tidyr::unnest(z, "parameters")
        groups <- tibble::as_tibble(do.call(rbind, z[["group"]]))
        # TODO Check if it would be ok to have a group column named "in_replicate",
        # "in_model" or "value".
        stopifnot(!(any(c("in_replicate", "in_model", "value") %in% colnames(groups))))
        if ("value" %in% colnames(z)) {
            params <- z[, c("in_replicate", "in_model", "value")]
        } else {
            params <- z[, c("in_replicate", "in_model")]
        }
        return(cbind(groups, params))
    }
    if (tibble) {
        if (!is.null(nm[["group"]])) {
            return(nm[, c("group", "parameters")])
        } else {
            return(nm[, "parameters"])
        }
    }
    return(sort(unique(params[["in_model"]])))
}

### * comps()

#' Return the compartments of a network model
#'
#' @param nm A \code{networkModel} object.
#'
#' @return A list of character vectors, with one list element per row of the
#'     input network model (list elements are in the same order as the input
#'     network model rows). Each list element containing the names of the
#'     compartments in the topology defined in the corresponding row of the
#'     input network model.
#'
#' @examples
#' comps(aquarium_mod)
#' comps(trini_mod)
#' 
#' @export

comps <- function(nm) {
    comps <- nm[, "topology"]
    comps$compartments <- lapply(comps$topology, colnames)
    return(comps[["compartments"]])
}

