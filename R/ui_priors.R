### * All functions in this file are exported

### * constant(value)

#' Define a fixed-value prior
#'
#' This is equivalent to having a fixed parameter.
#'
#' @param value The constant value of the parameter.
#'
#' @return A list defining the prior.
#'
#' @examples
#' constant(2)
#'
#' @export

constant <- function(value) {
    x <- list(type = "constant",
              parameters = c(value = value))
    x <- structure(x, class = "prior")
    return(x)
} 

### * hcauchy(scale)

#' Define a half-Cauchy prior (on [0;+Inf])
#'
#' @param scale Median of the half-Cauchy distribution.
#'
#' @return A list defining the prior.
#'
#' @importFrom stats rcauchy
#'
#' @examples
#' hcauchy(scale = 0.5)
#'
#' @export

hcauchy <- function(scale) {
    x <- list(type = "hcauchy",
              parameters = c(scale = scale))
    x <- structure(x, class = "prior")
    return(x)
}

### * normal(mean, sd)

#' Define a truncated normal prior (on [0;+Inf])
#'
#' @param mean Mean of the untruncated normal.
#' @param sd Standard deviation of the untruncated normal.
#'
#' @return A list defining the prior.
#'
#' @importFrom stats rnorm
#'
#' @examples
#' normal(mean = 0, sd = 4)
#'
#' @export

normal <- function(mean, sd) {
    x <- list(type = "trun_normal",
              parameters = c(mean = mean,
                             sd = sd))
    x <- structure(x, class = "prior")
    return(x)
}

### * uniform(min, max)

#' Define a uniform prior
#'
#' @param min,max Minimum and maximum boundaries for the uniform prior.
#'
#' @return A list defining the prior.
#'
#' @importFrom stats runif
#' 
#' @examples
#' uniform(min = 0, max= 1)
#'
#' @export

uniform <- function(min, max) {
    x <- list(type = "uniform",
              parameters = c(min = min, max = max))
    x <- structure(x, class = "prior")
    return(x)
}


### * scaled_beta(alpha, beta, scale=1)

#' Define a beta prior (on [0;scale])
#'
#' If a random variable X follows a scaled beta distribution with parameters
#' (alpha, beta, scale), then X/scale follows a beta distribution with
#' parameters (alpha, beta).
#' 
#' @param alpha Alpha parameter of the unscaled beta distribution.
#' @param beta Beta parameter of the unscaled beta distribution.
#' @param scale The upper boundary of the prior.
#'
#' @return A list defining the prior.
#'
#' @examples
#' scaled_beta(0.8, 20, scale = 10)
#'
#' @export

scaled_beta <- function(alpha, beta, scale = 1) {
    x <- list(type = "scaled_beta",
              parameters = c(alpha = alpha, beta = beta, scale = scale))
    x <- structure(x, class = "prior")
    return(x)
}

### * Methods for nice display of priors

### ** format.prior()

#' Pretty formatting of a \code{prior} object
#'
#' @param x An object of class \code{prior}.
#' @param ... Not used.
#'
#' @return A character string for pretty printing of a prior.
#' 
#' @export

format.prior <- function(x, ...) {
    params <- paste0("(",
                     paste(paste(names(x[["parameters"]]), x[["parameters"]],
                                 sep = "="), collapse = ","),
                     ")")
    out <- paste0(x[["type"]], " ", params)
    return(out)
}

### ** print.prior()

#' Pretty printing of a \code{prior} object
#'
#' @param x An object of class \code{prior}.
#' @param ... Not used.
#'
#' @return Mostly called for its side effect of printing, but also returns its
#'     input invisibly.
#' 
#' @export

print.prior <- function(x, ...) {
    cat(format(x), sep = "\n")
    invisible(x)
}

### * Extending tibbles

# https://cran.r-project.org/web/packages/tibble/vignettes/extending.html

#' Function used for displaying \code{prior} object in tibbles
#'
#' @param x An object of class \code{prior}.
#'
#' @return Input formatted with \code{format(x)}.
#' 
#' @importFrom pillar type_sum
#' @export
type_sum.prior <- function(x) {
    format(x)
}

#' Function used for displaying \code{prior} object in tibbles
#'
#' @param x An object of class \code{prior}.
#' 
#' @return Input formatted with \code{format(x)}.
#' 
#' @importFrom pillar obj_sum
#' @export
obj_sum.prior <- function(x) {
    format(x)
}

#' Function used for displaying \code{prior} object in tibbles
#'
#' @param x An object of class \code{prior}.
#' @param ... Not used.
#'
#' @return An object prepared with pillar::new_pillar_shaft_simple.
#' 
#' @importFrom pillar pillar_shaft
#' @export
pillar_shaft.prior <- function(x, ...) {
    out <- format(x)
    out[is.na(x)] <- NA
    pillar::new_pillar_shaft_simple(out, align = "right")
}

### * Methods for Ops on priors (implementing '==' operator)

# https://stackoverflow.com/a/35902710

#' Implementation of the '==' operator for priors
#'
#' @param e1,e2 Objects of class "prior".
#'
#' @return Boolean (or throws an error for unsupported operators).
#' 
#' @examples
#' p <- constant(0)
#' q <- constant(4)
#' p == q
#'
#' p <- hcauchy(2)
#' q <- hcauchy(2)
#' p == q
#'
#' @method Ops prior
#' 
#' @export

Ops.prior <- function(e1, e2) {
    op <- .Generic[[1]]
    switch(op,
           `==` = {
               if (e1$type != e2$type) {
                   return(FALSE)
               }
               if (!all(e1$parameters == e2$parameters)) {
                   return(FALSE)
               }
               return(TRUE)
           },
           stop("Undefined operation for objects of class \"priors\".")
           )
}

### * sample_from_prior()

#' Sample from a prior object
#'
#' @param x A \code{prior} object.
#' @param n Integer, number of samples to draw.
#'
#' @return A numeric vector of length \code{n}.
#' 
#' @examples
#' sample_from_prior(constant(1))
#' sample_from_prior(constant(1), 10)
#' sample_from_prior(hcauchy(0.5), 1)
#' hist(sample_from_prior(hcauchy(0.5), 20))
#' hist(sample_from_prior(uniform(0, 3), 1000))
#' hist(sample_from_prior(scaled_beta(3, 7, 2), 1e4))
#' 
#' @export
#'

sample_from_prior <- function(x, n = 1) {
    switch(x$type,
           "constant" = {
               rep(x$parameters[["value"]], n)
           },
           "hcauchy" = {
               scale <- x$parameters[["scale"]]
               replicate(n,
               {
                   o <- rcauchy(1, location = 0, scale = scale)
                   while (o < 0) {
                       o <- rcauchy(1, location = 0, scale = scale)
                   }
                   return(o)
               })
           },
           "trun_normal" = {
               mean <- x$parameters[["mean"]]
               sd <- x$parameters[["sd"]]
               replicate(n, {
                   o <- rnorm(1, mean = mean, sd = sd)
                   while (o < 0) {
                       o <- rnorm(1, mean = mean, sd = sd)
                   }
                   return(o)
               })
           },
           "uniform" = {
               runif(n, min = x$parameters[["min"]], max = x$parameters[["max"]])
           },
           "scaled_beta" = {
               p <- x$parameters
               p[["scale"]] * rbeta(n, shape1 = p[["alpha"]], shape2 = p[["beta"]])
           },
           stop("Unknown prior type."))
}
