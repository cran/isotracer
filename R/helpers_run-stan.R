### * None of the functions in this file is exported

### * mugen_stan()

#' Run a stan model from a network model (temporary name)
#'
#' Incorporate a loglik trace: https://mc-stan.org/loo/reference/extract_log_lik.html
#'
#' @param nm A \code{networkModel} object.
#' @param iter A positive integer specifying the number of iterations for each
#'     chain (including warmup). The default is 2000.
#' @param chains A positive integer specifying the number of Markov chains.
#'     The default is 4.
#' @param cores Number of cores to use for parallel run. Default is
#'     \code{NULL}, which means to use the value stored in
#'     \code{options()[["mc.cores"]]} (or 1 if this value is not set).
#' @param stanfit If TRUE, returns a `stanfit` object instead of the more
#'     classical `mcmc.list` object.
#' @param use_fixed_values Boolean, if TRUE any parameter value set with
#'     \code{set_params()} will be taken as fixed during the MCMC run. Default
#'     is FALSE.
#' @param ... Passed to \code{rstan::sampling}.
#'
#' @keywords internal
#' @noRd

mugen_stan <- function(nm, iter = 2000, chains = 4, cores = NULL,
                       stanfit = FALSE, use_fixed_values = FALSE, ...) {
    # Detect cores
    cores <- get_n_cores(cores = cores)
    # Convert network model to stan data
    stan.data <- prep_stan_data(nm, use_fixed_values = use_fixed_values)
    # Fit the model
    fit <- rstan::sampling(stanmodels[["networkModelMugen"]],
                           data = stan.data,
                           iter = iter,
                           chains = chains,
                           cores = cores,
                           pars = c("nonConstantParams", "log_lik"), ...)
    # Return stanfit object if requested
    if (stanfit) {
        return(fit)
    }
    # Get mcpars
    start <- fit@sim[["warmup"]] + 1
    end <- fit@sim[["iter"]]
    thin <- 1
    n_kept <- end - start + 1
    mcpars <- c(start, end, thin)
    # Prepare the mcmc.list object
    out <- rstan::As.mcmc.list(fit)
    for (i in seq_along(out)) {
        stopifnot(nrow(out[[i]]) == n_kept)
    }
    attr(out, "mcpar") <- mcpars
    rawNames <- coda::varnames(out)
    nonConstantParamNames <- rawNames[startsWith(rawNames, "nonConstantParams[")]
    loglikNames <- rawNames[startsWith(rawNames, "log_lik[")]
    ll <- out[, loglikNames]
    out <- out[, nonConstantParamNames]
    coda::varnames(out) <- stan.data[["allParams"]][stan.data[["mappingParamPriorType"]] != 0]
    llTrace <- lapply(ll, function(x) {
        out <- coda::as.mcmc(apply(as.matrix(x), 1, sum))
        attr(out, "mcpar") <- attr(x, "mcpar")
        return(out)
    })
    llTrace <- coda::as.mcmc.list(llTrace)
    attr(out, "loglik") <- llTrace
    attr(out, "mcpar") <- mcpars
    # Add values of constant parameters (if any)
    n_constant_params <- sum(stan.data[["mappingParamPriorType"]] == 0)
    if (n_constant_params > 0) {
        constant_params <- stan.data[["allParams"]][stan.data[["mappingParamPriorType"]] == 0]
        constant_values <- stan.data[["constantParams"]][stan.data[["mappingParamPriorType"]] == 0]
        constant_params <- setNames(constant_values, nm = constant_params)
        attr(out, "constant_params") <- constant_params
    }
    # Return the mcmc.list object
    return(out)
}

### * prep_stan_data()

#' Prepare stan data from a network model
#'
#' @param nm A \code{networkModel} object.
#' @param dt,grid_size Either the time step size for trajectory calculations
#'     (\code{dt}) or the number of points for the calculation
#'     (\code{grid_size}) can be provided. If none is provided, then a default
#'     grid size of 256 steps is used.
#' @param use_fixed_values Boolean, if TRUE any parameter value set with
#'     \code{set_params()} will be taken as fixed during the MCMC run. Default
#'     is FALSE.
#' 
#' @keywords internal
#' @noRd

prep_stan_data <- function(nm, dt = NULL, grid_size = NULL,
                           use_fixed_values = FALSE) {
    d <- list()
    end <- NULL
    params_nm <- params(nm)
    priors_nm <- priors(nm, fix_set_params = use_fixed_values,
                        quiet = TRUE)
    priors_nm <- priors_nm[match(params_nm, priors_nm[["in_model"]]), ]
    stopifnot(all(params_nm == priors_nm[["in_model"]]))
    # For now the stan model is only implemented for network models with the
    # same number of compartments on each row (a more general case where rows
    # can have different numbers of compartments is easily converted to this
    # case, by adding compartments without connections to fill the topology in
    # each row).
    stopifnot(length(unique(sapply(comps(nm), length))) == 1)
    # Counts
    d[["nComps"]] <- length(comps(nm)[[1]])
    d[["nGroups"]] <- nrow(nm)
    d[["nParams"]] <- length(params_nm)
    # Encode steady state compartments
    dSteady <- encode_steady(nm)
    d <- c(d, dSteady)
    # Encode split compartments
    dSplit <- encode_split(nm, params_nm)
    d <- c(d, dSplit)
    # Encode initial conditions
    dInit <- encode_init(nm)
    d <- c(d, dInit)
    # Encode events
    dEvents <- encode_events(nm, dt = dt, grid_size = grid_size, end = end)
    d <- c(d, dEvents)
    # Encode observations (including eta and zeta parameter indices)
    dObs <- encode_obs(nm, params_nm, dt = dt, grid_size = grid_size, end = end)
    d <- c(d, dObs)
    # Encode parameter priors
    dPriors <- encode_priors(params_nm, priors_nm)
    d <- c(d, dPriors)
    # Encode time schemes
    dTimeSchemes <- encode_time_schemes(nm, dt = dt, grid_size = grid_size,
                                        end = end)
    d <- c(d, dTimeSchemes)
    # Encode uptake rates (upsilons)
    dUpsilons <- encode_upsilons(nm, params_nm)
    d <- c(d, dUpsilons)
    # Encode losses (lambdas)
    dLambdas <- encode_lambdas(nm, params_nm)
    d <- c(d, dLambdas)
    # Encode decay rate for radioactive tracers
    lambda_decay <- attr(nm, "lambda_hl")
    if (is.null(lambda_decay)) {
        lambda_decay <- 0
    }
    d[["lambda_decay"]] <- lambda_decay
    # Encode distribution family for proportions
    known_families <- c("gamma_cv" = 1, "normal_cv" = 2, "normal_sd" = 3,
                        "beta_phi" = 4)
    prop_family <- attr(nm, "prop_family")
    if (!prop_family %in% names(known_families)) {
        stop("Unknown distribution family for proportions. Got value: ",
             prop_family, "\n",
             "Allowed values are: ", names(known_families))
    }
    d[["propFamily"]] <- known_families[prop_family]
    # Encode distribution family for sizes
    size_known_families <- c("normal_cv" = 1, "normal_sd" = 2)
    size_family <- attr(nm, "size_family")
    if (!size_family %in% names(size_known_families)) {
        stop("Unknown distribution family for sizes. Got value: ",
             size_family, "\n",
             "Allowed values are: ", names(size_known_families))
    }
    d[["sizeFamily"]] <- size_known_families[size_family]
    # Return
    d[["allParams"]] <- params_nm
    return(d)
}

### * encode_steady()

#' Encode steady compartments for stan data
#'
#' @param nm A \code{networkModel} object.
#' 
#' @keywords internal
#' @noRd

encode_steady <- function(nm) {
    d <- list()
    steady <- lapply(nm[["topology"]], function(x) {
        match(attr(x, "steadyState"), colnames(x))
    })
    d[["maxNsteady"]] <- max(sapply(steady, length))
    d[["nSteady"]] <- setNames(c(sapply(steady, length), 0), # Padded
                               nm = c(paste0("grp", seq_len(nrow(nm))), "padding"))
    d[["steadyIndices"]] <- array(0, dim = c(d[["maxNsteady"]], nrow(nm)),
                                  dimnames = list(seq_len(d[["maxNsteady"]]),
                                                  paste0("grp", seq_len(nrow(nm)))))
    for (i in seq_len(nrow(nm))) {
        d[["steadyIndices"]][seq_along(steady[[i]]),i] <- steady[[i]]
    }
    # Return
    return(d)
}

### * encode_split()

#' Encode split compartments for stan data
#'
#' @param nm A \code{networkModel} object.
#' @param allParams Parameters of the network model.
#'
#' @keywords internal
#' @noRd

encode_split <- function(nm, allParams) {
    d <- list()
    split <- lapply(nm[["topology"]], function(x) {
        match(attr(x, "split"), colnames(x))
    })
    d[["splitPresent"]] <- as.numeric(max(sapply(split, length)) > 0)
    nGroups <- nrow(nm)
    splitComps <- array(0, dim = c(length(comps(nm)[[1]]), nGroups),
                        dimnames = list(seq_along(comps(nm)[[1]]),
                                        paste0("grp", seq_len(nrow(nm)))))
    for (i in seq_along(split)) {
        splitComps[split[[i]], i] <- 1
    }
    d[["splitComps"]] <- splitComps
    # Parameter mapping
    piMapping <- array(0, dim = c(length(comps(nm)[[1]]), nGroups),
                       dimnames = list(seq_along(comps(nm)[[1]]),
                                       paste0("grp", seq_len(nrow(nm)))))
    for (g in seq_len(nGroups)) {
        compNames <- colnames(nm$topology[[g]])
        stopifnot(length(compNames) == length(comps(nm)[[1]]))
        for (k in seq_along(compNames)) {
            if (splitComps[k, g] > 0) {
                paramName <- paste0("portion.act_", compNames[k])
                paramGlobal <- nm$parameters[[g]]$in_model[nm$parameters[[g]]$in_replicate == paramName]
                stopifnot(length(paramGlobal) == 1)
                piMapping[k, g] <- match(paramGlobal, allParams)
            }
        }
    }
   d[["piMapping"]] <- piMapping
    # Return
    return(d)
}

### * encode_init()

#' Encode initial conditions for a network model
#'
#' For now this function assumes that each row has the same number of
#' compartments.
#'
#' @param nm A \code{networkModel} object.
#'
#' @keywords internal
#' @noRd

encode_init <- function(nm) {
    nComps <- sapply(comps(nm), length)
    stopifnot(all(nComps == nComps[1]))
    nGroups <- nrow(nm)
    d <- array(0, dim = c(nComps[1], 2, nGroups),
               dimnames = list(1:nComps[1], c("unmarked", "marked"),
                               c(paste0("grp", seq_len(nGroups)))))
    for (i in seq_len(nGroups)) {
        comps <- colnames(nm$topology[[i]])
        stopifnot(nrow(nm$initial[[i]]) == length(comps))
        stopifnot(setequal(nm$initial[[i]][["compartment"]], comps))
        comps <- nm$initial[[i]][match(comps, nm$initial[[i]][["compartment"]]), ]
        stopifnot(all(comps$compartment == colnames(nm$topology[[i]])))
        unmarked <- comps$size * (1 - comps$proportion)
        marked <- comps$size * comps$proportion
        d[,,i] <- cbind(unmarked, marked)
    }
    return(list(initialQuantities = d))
}

### * encode_events()

#' Encode events
#'
#' For now, only pulse events are encoded.
#'
#' @param nm A \code{networkModel} object.
#' @param dt,grid_size Either the time step size for trajectory calculations
#'     (\code{dt}) or the number of points for the calculation
#'     (\code{grid_size}) can be provided. If none is provided, then a default
#'     grid size of 256 steps is used.
#' @param end Time value for end point. If not provided, the last observation
#'     or event is used.
#'
#' @examples
#' encode_events <- isotracer:::encode_events
#' encode_events(aquarium_mod)
#' encode_events(trini_mod)
#' 
#' @keywords internal
#' @noRd

encode_events <- function(nm, dt = NULL, grid_size = NULL, end = NULL) {
    d <- nm_get_time_schemes(nm, dt = dt, grid_size = grid_size, end = end)
    nGroups <- nrow(nm)
    o <- list()
    events <- nm[["events"]]
    stopifnot(all(dplyr::bind_rows(events)[["event"]] == "pulse"))
    # Encode compartments and timepoints
    for (i in seq_len(nrow(nm))) {
        comps <- colnames(nm$topology[[i]])
        comps <- setNames(seq_along(comps), nm = comps)
        if (!is.null(events[[i]])) {
            events[[i]]$compartment <- comps[events[[i]]$compartment]
            events[[i]]$timepoints <- match(events[[i]]$time, d$timepoints[[i]])
        }
    }
    # Encode the number of events
    if (length(events) == 0) {
        nEvents <- rep(0, nrow(nm))
    } else {
        nEvents <- sapply(events, function(x) {
            if (is.null(x)) return(0)
            return(nrow(x))
        })
    }
    o[["maxNpulseEvents"]] <- max(nEvents)
    o[["nPulseEvents"]] <- setNames(c(nEvents, 0), # Padded
                                    nm = c(paste0("grp", seq_len(nrow(nm))), "padding"))
    # Encode pulses
    o[["pulseEventsIndices"]] <- array(0, dim = c(max(nEvents), 2, nGroups),
                                       dimnames = list(seq_len(max(nEvents)),
                                                       c("timepoint", "comp"),
                                                       paste0("grp", seq_len(nGroups))))
    o[["pulseEventsQuantities"]] <- array(0, dim = c(max(nEvents), 2, nGroups),
                                          dimnames = list(seq_len(max(nEvents)),
                                                          c("unmarked", "marked"),
                                                          paste0("grp", seq_len(nGroups))))
    for (i in seq_len(nGroups)) {
        if (!is.null(events[[i]])) {
            o[["pulseEventsIndices"]][1:nEvents[i], 1:2, i] <-
                as.matrix(events[[i]][, c("timepoints", "compartment")])
            chars <- lapply(events[[i]]$characteristics, tibble::as_tibble)
            chars <- dplyr::bind_rows(chars)[, c("unmarked", "marked")]
            o[["pulseEventsQuantities"]][1:nEvents[i], 1:2, i] <-
                as.matrix(chars)
        }
    }
    # Return
    return(o)
}

### ** nm_get_time_schemes()

#' Build the time schemes for numerical solving of the system of differential equations
#'
#' This function processes each row of a networkModel. It always assumes that
#' the starting time point is t=0.
#'
#' @param nm A networkModel
#' @param dt,grid_size Either the time step size for trajectory calculations
#'     (\code{dt}) or the number of points for the calculation
#'     (\code{grid_size}) can be provided. If none is provided, then a default
#'     grid size of 256 steps is used.
#' @param end Time value for end point. If not provided, the last observation
#'     or event is used.
#' @param at Optional, vector of time values at which the trajectory must be
#'     evaluated
#'
#' @keywords internal
#' @noRd

nm_get_time_schemes <- function(nm, dt = NULL, grid_size = NULL, end = NULL,
                                at = NULL) {
    # Get the time schemes for each row of nm
    ts <- lapply(seq_len(nrow(nm)), function(i) {
        z <- nm_row_get_time_scheme(nm[i, ], dt = dt, grid_size = grid_size,
                                    end = end, at = at)
        z <- tibble::as_tibble(lapply(z, list))
    })
    ts <- dplyr::bind_rows(ts)
    return(ts)    
}

### ** nm_row_get_time_scheme()

#' Build the time scheme for numerical solving of the system of differential equations
#'
#' This function is applied to one row of a networkModel. It always assumes
#' that the starting time point is t=0.
#'
#' @param nm_row A row from a \code{networkModel} object.
#' @param dt,grid_size Either the time step size for trajectory calculations
#'     (\code{dt}) or the number of points for the calculation
#'     (\code{grid_size}) can be provided. If none is provided, then a default
#'     grid size of 256 steps is used.
#' @param end Time value for end point. If not provided, the last observation
#'     or event is used.
#' @param at Optional, vector of time values at which the trajectory must be
#'     evaluated
#'
#' @examples
#' isotracer:::nm_row_get_time_scheme(aquarium_mod)
#' 
#' @keywords internal
#' @noRd

nm_row_get_time_scheme <- function(nm_row, dt = NULL, grid_size = NULL, end = NULL,
                                   at = NULL) {
    nmRow <- nm_row
    stopifnot(nrow(nmRow) == 1)
    # Parse dt and gridsize
    if (!(is.null(dt) | is.null(grid_size))) {
        stop("Only \"dt\" or \"grid_size\" can be specified, not both.")
    }
    if (is.null(dt) & is.null(grid_size)) {
        grid_size <- 256
    }
    # Get observation times
    observations <- nmRow[["observations"]][[1]]
    obsTimes <- observations[["time"]]
    obsTimes <- sort(unique(obsTimes))
    # Get events time
    eventTimes <- c()
    if (!is.null(nmRow[["events"]][[1]])) {
        eventTimes <- unique(nmRow[["events"]][[1]]$time)
    }
    # Get "at" times
    atTimes <- c()
    if (!is.null(at)) {
        atTimes <- unique(at)
    }
    # Get end time
    if (is.null(end)) {
        maxTime <- max(c(obsTimes, eventTimes, atTimes))
    } else {
        maxTime <- end
    }
    # Calculate dt
    if (is.null(dt)) {
        dt <- maxTime / grid_size
    }
    # Generate a time scheme
    timepoints <- c(seq(0, maxTime, by = dt), obsTimes, eventTimes, atTimes)
    timepoints <- sort(unique(timepoints))
    timepoints <- timepoints[timepoints <= maxTime]
    # (if end = NULL, timepoints is a timeline with timesteps at most dt wide
    # and containing all the observations sampling times and the event times,
    # even if not a multiple of dt.)
    # (else, timepoints is truncated at end.)
    # Annotate the observations tibble with the timepoints ids
    observations[["timepoint"]] <- match(observations[["time"]], timepoints)
    if (is.null(end)) {
        stopifnot(!any(is.na(observations[["timepoint"]])))
    }
    # Get unique dts
    dts <- timepoints_to_dt(timepoints)
    # Return
    return(list(timepoints = timepoints,
                unique_dt = dts[["unique_dt"]],
                dt_i = dts[["dt_i"]],
                observations = observations))
}

### ** timepoints_to_dt()

#' Prepare the set of unique dt from a vector of timepoints
#'
#' Useful to identify how many different transfer matrices have to be calculated.
#'
#' @param timepoints Numeric vector of timepoints, sorted.
#'
#' @keywords internal
#' @noRd

timepoints_to_dt <- function(timepoints) {
    dts <- diff(timepoints)
    uniqueDts <- sort(unique(dts))
    dt_i <- match(dts, uniqueDts)
    stopifnot(!any(is.na(dt_i)))
    return(list(unique_dt = uniqueDts,
                dt_i = dt_i))
}

### * encode_obs()

#' Encode observations for a network model
#'
#' @param nm A \code{networkModel} object.
#' @param allParams Parameters of the network model.
#' @param dt,grid_size Either the time step size for trajectory calculations
#'     (\code{dt}) or the number of points for the calculation
#'     (\code{grid_size}) can be provided. If none is provided, then a default
#'     grid size of 256 steps is used.
#' @param end Time value for end point. If not provided, the last observation
#'     or event is used.
#' 
#' @return NULL if the observations column only contain NULLs.
#'
#' @importFrom stats na.omit
#' 
#' @keywords internal
#' @noRd

encode_obs <- function(nm, allParams, dt = NULL, grid_size = NULL, end = NULL) {
    d <- nm_get_time_schemes(nm, dt = dt, grid_size = grid_size, end = end)
    nGroups <- nrow(nm)
    zeta_by_comp <- attr(nm, "size_zeta_per_compartment")
    if (is.null(zeta_by_comp)) {
        zeta_by_comp<- FALSE
    }
    # TODO Add filtering to keep only compartments present in topo
    # TODO Handle gracefully the case without observations
    if (all(sapply(nm$observations, is.null))) {
        return(NULL)
    }
    # Get sizes
    sizes <- purrr::map(d$observations, function(x) {
        na.omit(x[, c("compartment", "size", "timepoint")])
    })
    # Get proportions
    props <- purrr::map(d$observations, function(x) {
        na.omit(x[, c("compartment", "proportion", "timepoint")])
    })
    # Encode compartments
    for (i in seq_len(nrow(nm))) {
        comps <- colnames(nm$topology[[i]])
        comps <- setNames(seq_along(comps), nm = comps)
        sizes[[i]]$compartment <- comps[sizes[[i]]$compartment]
        props[[i]]$compartment <- comps[props[[i]]$compartment]
    }
    # Encode sizes and props
    o <- list()
    o[["nSizesObs"]] <- c(purrr::map_dbl(sizes, nrow), 0) # Padded
    names(o[["nSizesObs"]]) <- c(paste0("grp", seq_len(nrow(nm))), "padding")
    o[["nPropsObs"]] <- c(purrr::map_dbl(props, nrow), 0) # Padded
    names(o[["nPropsObs"]]) <- c(paste0("grp", seq_len(nrow(nm))), "padding")
    o[["maxNsizesObs"]] <- max(o[["nSizesObs"]])
    o[["maxNpropsObs"]] <- max(o[["nPropsObs"]])
    # Prepare containers
    o[["sizesObsIndices"]] <- array(0, dim = c(o[["maxNsizesObs"]], 3, nGroups),
                                    dimnames = list(seq_len(o[["maxNsizesObs"]]),
                                                    c("comp", "timepoint", "zeta"),
                                                    paste0("grp", seq_len(nGroups))))
    o[["sizesObs"]] <- array(0, dim = c(o[["maxNsizesObs"]], nGroups),
                             dimnames = list(seq_len(o[["maxNsizesObs"]]),
                                             paste0("grp", seq_len(nGroups))))
    o[["propsObsIndices"]] <- array(0, dim = c(o[["maxNpropsObs"]], 3, nGroups),
                                    dimnames = list(seq_len(o[["maxNpropsObs"]]),
                                                    c("comp", "timepoint", "eta"),
                                                    paste0("grp", seq_len(nGroups))))
    o[["propsObs"]] <- array(0, dim = c(o[["maxNpropsObs"]], nGroups),
                             dimnames = list(seq_len(o[["maxNpropsObs"]]),
                                             paste0("grp", seq_len(nGroups))))
    # Fill containers
    for (g in seq_len(nGroups)) {
        if (o[["nSizesObs"]][g] > 0) {
            o[["sizesObsIndices"]][1:o[["nSizesObs"]][g], 1:2, g] <- as.matrix(sizes[[g]][, c("compartment", "timepoint")])
            o[["sizesObs"]][1:o[["nSizesObs"]][g], g] <- as.matrix(sizes[[g]][, c("size")])
            if (!zeta_by_comp) {
                zeta_global <- nm[["parameters"]][[g]]$in_model[nm[["parameters"]][[g]]$in_replicate == "zeta"]
                zeta_index <- match(zeta_global, allParams)
                o[["sizesObsIndices"]][1:o[["nSizesObs"]][g], 3, g] <- zeta_index
            } else {
                comps <- colnames(nm$topology[[g]])
                comps <- setNames(seq_along(comps), nm = comps)
                zeta_replicate <- paste0("zeta_", names(comps)[sizes[[g]][["compartment"]]])
                zeta_global <- nm[["parameters"]][[g]]$in_model[match(zeta_replicate, nm[["parameters"]][[g]]$in_replicate)]
                zeta_index <- match(zeta_global, allParams)
                o[["sizesObsIndices"]][1:o[["nSizesObs"]][g], 3, g] <- zeta_index
            }
        }
    }
    for (g in seq_len(nGroups)) {
        if (o[["nPropsObs"]][g] > 0) {
            o[["propsObsIndices"]][1:o[["nPropsObs"]][g], 1:2, g] <- as.matrix(props[[g]][, c("compartment", "timepoint")])
            o[["propsObs"]][1:o[["nPropsObs"]][g], g] <- as.matrix(props[[g]][, c("proportion")])
            eta_global <- nm[["parameters"]][[g]]$in_model[nm[["parameters"]][[g]]$in_replicate == "eta"]
            eta_index <- match(eta_global, allParams)
            o[["propsObsIndices"]][1:o[["nPropsObs"]][g], 3, g] <- eta_index
        }
    }
    # Return
    return(o)
}

### * encode_priors()

#' Prepare the prior encoding for stan data
#'
#' @param params_nm Output from params(nm).
#' @param priors_nm Output from priors(nm) (must have the same parameter order
#'     as params_nm).
#'
#' @keywords internal
#' @noRd

encode_priors <- function(params_nm, priors_nm) {
    d <- list()
    d[["nPriorConstant_code0"]] <- 0
    d[["nPriorUniform_code1"]] <- 0
    d[["nPriorHcauchy_code2"]] <- 0
    d[["nPriorBeta_code3"]] <- 0
    d[["nPriorTrNormal_code4"]] <- 0
    ### Prior types (used to map the correct prior in the stan model)
    priorTypes <- c("constant" = 0, "uniform" = 1, "hcauchy" = 2, "scaled_beta" = 3,
                    "trun_normal" = 4)
    ### Parameter priors
    d[["mappingParamPriorType"]] <- rep(NA, length(params_nm))
    d[["mappingParamPriorID"]] <- rep(NA, length(params_nm))
    prior_i <- c("constant" = 1, "uniform" = 1, "hcauchy" = 1,
                 "scaled_beta" = 1, "trun_normal" = 1) # Counts for each prior type
    # Set some defaults
    for (priorDistParam in c("constantParams", "lowerParams", "upperParams", "hcauchyScaleParams",
                             "rawBetaAlpha", "rawBetaBeta", "betaScaleParams",
                             "trNormMeanParams", "trNormSdParams")) {
        # Set zero as default for all prior distribution parameters
        d[[priorDistParam]] <- rep(0, length(params_nm))
    }
    # Go through each param prior
    for (i in seq_along(params_nm)) {
        p <- priors_nm[["prior"]][[i]]
        stopifnot(p[["type"]] %in% names(priorTypes))
        d[["mappingParamPriorType"]][i] <- priorTypes[p[["type"]]]
        d[["mappingParamPriorID"]][i] <- prior_i[p[["type"]]]
        prior_i[p[["type"]]] <- prior_i[p[["type"]]] + 1
        if (p[["type"]] == "constant") {
            d[["constantParams"]][i] <- p[["parameters"]][["value"]]
            d[["nPriorConstant_code0"]] <- d[["nPriorConstant_code0"]] + 1
        }
        if (p[["type"]] == "uniform") {
            d[["lowerParams"]][i] <- p[["parameters"]]["min"]
            d[["upperParams"]][i] <- p[["parameters"]]["max"]
            d[["nPriorUniform_code1"]] <- d[["nPriorUniform_code1"]] + 1
        }
        if (p[["type"]] == "hcauchy") {
            d[["hcauchyScaleParams"]][i] <- p[["parameters"]]["scale"]
            d[["nPriorHcauchy_code2"]] <- d[["nPriorHcauchy_code2"]] + 1
        }
        if (p[["type"]] == "scaled_beta") {
            d[["rawBetaAlpha"]][i] <- p[["parameters"]]["alpha"]
            d[["rawBetaBeta"]][i] <- p[["parameters"]]["beta"]
            d[["betaScaleParams"]][i] <- p[["parameters"]]["scale"]
            d[["nPriorBeta_code3"]] <- d[["nPriorBeta_code3"]] + 1
        }
        if (p[["type"]] == "trun_normal") {
            d[["trNormMeanParams"]][i] <- p[["parameters"]][["mean"]]
            d[["trNormSdParams"]][i] <- p[["parameters"]][["sd"]]
            d[["nPriorTrNormal_code4"]] <- d[["nPriorTrNormal_code4"]] + 1
        }
    }
    # Count non-constant priors
    d[["nNonConstantPriors"]] <- length(params_nm) - d[["nPriorConstant_code0"]]
    # Return
    return(d)
}

### * encode_time_schemes()

#' Encode the time schemes for stan data
#'
#' @param nm A \code{networkModel} object.
#' @param dt,grid_size Either the time step size for trajectory calculations
#'     (\code{dt}) or the number of points for the calculation
#'     (\code{grid_size}) can be provided. If none is provided, then a default
#'     grid size of 256 steps is used.
#' @param end Time value for end point. If not provided, the last observation
#'     or event is used.
#'
#' @keywords internal
#' @noRd

encode_time_schemes <- function(nm, dt = NULL, grid_size = NULL, end = NULL) {
    # (timestep and dt are used interchangeably)
    # Get the time schemes
    ts <- nm_get_time_schemes(nm, dt = dt, grid_size = grid_size, end = end)
    # Build the arrays
    nGroups <- nrow(nm)
    n_unique_dts <- purrr::map_dbl(ts$unique_dt, length)
    maxN_unique_dts <- max(n_unique_dts)
    n_timesteps <- purrr::map_dbl(ts$dt_i, length)
    maxN_timesteps <- max(n_timesteps)
    unique_dts <- array(0, dim = c(maxN_unique_dts, nGroups),
                        dimnames = list(c(1:maxN_unique_dts),
                                        paste0("grp", 1:nGroups)))
    timesteps <- array(0, c(maxN_timesteps, nGroups),
                       dimnames = list(c(1:maxN_timesteps),
                                        paste0("grp", 1:nGroups)))
    # Fill the arrays
    for (i in seq_len(nrow(nm))) {
        unique_dts[1:n_unique_dts[i], i] <- ts[["unique_dt"]][[i]]
        timesteps[1:n_timesteps[i],  i] <- ts[["dt_i"]][[i]]
    }
    # Prepare data
    d <- list()
    d[["nTimesteps"]] <- setNames(c(n_timesteps, 0), # Padded
                                  nm = c(paste0("grp", 1:nGroups), "padding"))
    d[["nUniqueDts"]] <- setNames(c(n_unique_dts, 0), # Padded
                                  nm = c(paste0("grp", 1:nGroups), "padding"))
    d[["timesteps"]] <- timesteps
    d[["unique_dts"]] <- unique_dts
    d[["maxNtimesteps"]] <- max(d[["nTimesteps"]])
    d[["maxNuniqueDts"]] <- max(d[["nUniqueDts"]])
    return(d)
}

### * encode_upsilons()

#' Encode the uptake rates for stan data
#'
#' @param nm A \code{networkModel} object.
#' @param allParams Parameters of the network model.
#' 
#' @keywords internal
#' @noRd

encode_upsilons <- function(nm, allParams) {
    nGroups <- nrow(nm)
    upsilons <- nm_get_upsilons(nm, allParams)
    nUpsilons <- setNames(c(upsilons[["nUpsilons"]], 0), # Padding
                          nm = c(paste0("grp", seq_len(nGroups)), "padding"))
    maxNupsilons <- max(nUpsilons)
    upsilonMapping <- array(0, dim = c(maxNupsilons, 3, nGroups),
                            dimnames = list(1:maxNupsilons,
                                            c("from", "to", "param"),
                                            paste0("grp", seq_len(nGroups))))
    for (i in seq_len(nGroups)) {
        upsilonMapping[1:nUpsilons[i], 1:3, i] <- as.matrix(upsilons[["upsilons"]][[i]])
    }
    return(list(nUpsilons = nUpsilons,
                maxNupsilons = maxNupsilons,
                upsilonMapping = upsilonMapping))
}

### ** nm_get_upsilons()

#' Get upsilons for a network model
#' 
#' @param nm A \code{networkModel} object.
#' @param allParams Parameters of the network model.
#' 
#' @keywords internal
#' @noRd

nm_get_upsilons <- function(nm, allParams) {
    upsilons <- lapply(seq_len(nrow(nm)), function(i) {
        z <- nm_row_get_upsilons(nm[i, ], allParams = allParams)
        tibble::tibble(upsilons = list(z), nUpsilons = nrow(z))
    })
    return(dplyr::bind_rows(upsilons))
}

### ** nm_row_get_upsilons()

#' Get upsilons for one row of a network model
#' 
#' @param nm_row A row from a \code{networkModel} object.
#' @param allParams Parameters of the network model.
#' 
#' @keywords internal
#' @noRd

nm_row_get_upsilons <- function(nm_row, allParams) {
    nmRow <- nm_row
    stopifnot(nrow(nmRow) == 1)
    # Get data
    topo <- nmRow[["topology"]][[1]]
    mapping <- nmRow[["parameters"]][[1]][, c("in_replicate", "in_model")]
    mapping <- tibble::deframe(mapping)
    # Build output
    from <- vector()
    to <- vector()
    param <- vector()
    compartments <- colnames(topo)
    stopifnot(all(compartments == rownames(topo)))
    for (j in seq_len(ncol(topo))) {
        for (i in seq_len(nrow(topo))) {
            if (topo[i,j] == 1) {
                from <- c(from, j)
                to <- c(to, i)
                paramName <- paste0("upsilon_", compartments[j], "_to_",
                                    compartments[i])
                param <- c(param, match(mapping[paramName], allParams))
            }
        }
    }
    return(tibble::tibble(from = from, to = to, param = param))
}

### * encode_lambdas()

#' Encode the loss rates for stan data
#'
#' @param nm A \code{networkModel} object.
#' @param allParams Parameters of the network model.
#'
#' @keywords internal
#' @noRd

encode_lambdas <- function(nm, allParams) {
    nGroups <- nrow(nm)
    lambdas <- nm_get_lambdas(nm, allParams)
    nLambdas <- setNames(c(lambdas[["nLambdas"]], 0), # Padding
                         nm = c(paste0("grp", seq_len(nGroups)), "padding"))
    maxNlambdas <- max(nLambdas)
    lambdaMapping <- array(0, dim = c(maxNlambdas, 2, nGroups),
                           dimnames = list(1:maxNlambdas,
                                           c("from", "param"),
                                           paste0("grp", seq_len(nGroups))))
    for (i in seq_len(nGroups)) {
        lambdaMapping[1:nLambdas[i], 1:2, i] <- as.matrix(lambdas[["lambdas"]][[i]])
    }
    return(list(nLambdas = nLambdas,
                maxNlambdas = maxNlambdas,
                lambdaMapping = lambdaMapping))
}

### ** nm_get_lambdas()

#' Get lambdas for a network model
#' 
#' @param nm A \code{networkModel} object.
#' @param allParams Parameters of the network model.
#' 
#' @keywords internal
#' @noRd

nm_get_lambdas <- function(nm, allParams) {
    lambdas <- lapply(seq_len(nrow(nm)), function(i) {
        z <- nm_row_get_lambdas(nm[i, ], allParams = allParams)
        tibble::tibble(lambdas = list(z), nLambdas = nrow(z))
    })
    return(dplyr::bind_rows(lambdas))
}

### ** nm_row_get_lambdas()

#' Get lambdas for one row of a network model
#' 
#' @param nm_row A row from a \code{networkModel} object.
#' @param allParams Parameters of the network model.
#' 
#' @keywords internal
#' @noRd

nm_row_get_lambdas <- function(nm_row, allParams) {
    nmRow <- nm_row
    stopifnot(nrow(nmRow) == 1)
    # Get data
    topo <- nmRow[["topology"]][[1]]
    mapping <- nmRow[["parameters"]][[1]][, c("in_replicate", "in_model")]
    mapping <- tibble::deframe(mapping)
    # Build output
    from <- vector()
    param <- vector()
    compartments <- colnames(topo)
    stopifnot(all(compartments == rownames(topo)))
    for (j in seq_len(ncol(topo))) {
        paramName <- paste0("lambda_", compartments[j])
        if (paramName %in% names(mapping)) {
            from <- c(from, j)
            param <- c(param, match(mapping[paramName], allParams))
        }
    }
    return(tibble::tibble(from = from, param = param))
}
