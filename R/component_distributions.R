# Generics -----

#' Convolved logpdf
#'
#' Generic function for `component_distribution` objects, computes the loglikelihood of an observation from
#' the component distribution, corrupted by gaussian measurement error
#'
#' @param betahat observations
#' @param se standard errors of observations
#' @return the convolved log density log p(betahat | se) = log \int N(betahat | beta, se) p(beta) dbeta
#' @export
convolved_logpdf <- function(x, ...) {
  UseMethod("convolved_logpdf", x)
}

#' Update Params
#'
#' Generic function for updating the parameters of a `component_distribution`
#'
#' @param betahat observations
#' @param se standard errors of observations
#' @param weights weights to weigh each observation by when optimizing
#' @return a new `component_distribution` object with update parameters
#' @export
update_params <- function(x, ...) {
  UseMethod("update_params", x)
}

# Defaults ----

#' @export
convolved_logpdf.component_distribution <- function(dist, betahat, se) {
  # TODO: compute convolved pdf by numerical integration?
  # Will probably require logpdf/pdf be implemented for each component
  stop("Generic convolved logpdf not implimented yet")
}

#' @export
update_params.component_distribution <- function(dist, betahat, se) {
  # TODO: compute convolved pdf by numerical integration?
  # Will probably require logpdf/pdf be implemented for each component
  stop("Generic update not implimented")
}

# Point mass Component----

#' @export
point_component <- function(mu = 0) {
  f <- list(mu = mu)
  class(f) <- c("point", "component_distribution")
  return(f)
}

#' @export
is.point <- function(x) {
  inherits(x, "point")
}

#' Just a normal distribution with mean centered at the point mass
#' @export
convolved_logpdf.point <- function(dist, betahat, se) {
  # return(dnorm(betahat, sd=se, log=T))
  sd <- se
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 1e3, -1e3)
  return(logp)
}

# Normal Component----

#' @export
normal_component <- function(mu = 0, var = 1) {
  f <- list(mu = mu, var = var)
  class(f) <- c("normal", "component_distribution")
  return(f)
}

#' @export
is.normal <- function(x) {
  inherits(x, "normal")
}

#' Normal distribution with variance given by sum of component dist. and noise
#' @export
convolved_logpdf.normal <- function(dist, betahat, se) {
  sd <- sqrt(se^2 + dist$var)
  logp <- dnorm(betahat, mean = dist$mu, sd = sd, log = TRUE)
  logp <- .clamp(logp, 100, -100)
  return(logp)
}

#' Note: only updates the variance parameter, these components are assumed
#' to be mean 0
#' @export
update_params.normal <- function(dist, betahat, se, weights) {
  # TODO: right now it's just a grid search, but we can impliment EM update easily
  var.init <- dist$var
  grid_var <- 2^(seq(-1, 1, by = 0.1) + log2(var.init))
  grid_dist <- purrr::map(grid_var, ~ normal_component(mu = dist$mu, var = .x))
  ll <- purrr::map_dbl(grid_dist, ~ sum(convolved_logpdf.normal(.x, betahat, se) * weights, na.rm = T))
  return(grid_dist[[which.max(ll)]])
}



# Beta components -----
#strictly increasing and convex, if alpha \geq 1 and beta \leq1
# strictly decreasing and convex if alpha \leq 1 and beta \geq1
#' @export
beta_component <- function(alpha = 1, beta=1) {
  f <- list(alpha = alpha,
            beta  = beta
            )
  class(f) <- c("beta", "component_distribution")
  return(f)
}

#' @export
convolved_logpdf.beta <- function(dist, p, se = 1) {
  logp <- dbeta(p,
                shape1 = dist$alpha,
                shape2 = dist$beta
                )
  logp <- .clamp(logp, 100, -100)
  return(logp)
}

#' @export
update_params.beta <- function(dist, p, se = 1, weights) {
  # TODO: decide how to update alpha
  # weights come from posterior assignment probabilities.
}

