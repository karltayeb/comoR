# interface for covariate moderated prior weights

#' Compute the prior probability of assignment to each mixture component
#' @export
compute_log_prior_assignment <- function(x, ...){
  UseMethod("compute_log_prior_assignment", x)
}


#' Update the multinomial regression
#' @export
update_prior <- function(x, ...){
  UseMethod("update_prior", x)
}

#' Compute ELBO, but if there is not ELBO
#' @export
compute_elbo <- function(x, ...){
  UseMethod("compute_elbo", x)
}

#' Default no ELBO
#' @export
compute_elbo.default <- function(mnreg, ...){
  return(0)
}


# Constant prior------------
# Prior is the same for all observations
initialize_constant_mnreg <- function(K){
  logpi <- rep(log(1/K), K)
  mnreg <- list(logpi=logpi, K=K)
  class(mnreg) <- 'constant_mnreg'
  return(mnreg)
}

#' @param mnreg multinomial regression object
#' @param data a list at least containing covariates `X`
#' @export
compute_log_prior_assignment.constant_mnreg <- function(mnreg, data){
  n <- nrow(data$X)
  logpi <- matrix(rep(mnreg$logpi, n), ncol = mnreg$K, byrow = T)
  return(logpi)
}

#' @param mnreg multinomial regression object
#' @param resp responsibilities, posterior assignment probability
#'     of each observation to each class
#' @param data a list at least containing covariates `X`
#' @export
update_prior.constant_mnreg <- function(mnreg, resps, data){
  new_log_pis <- log(colMeans(resps))
  mnreg$logpi <- new_log_pis
  return(mnreg)
}


#' @param mnreg multinomial regression object
#' @param resp responsibilities, posterior assignment probability
#'     of each observation to each class
#' @param data a list at least containing covariates `X`
#' @export
compute_elbo.constant_mnreg <- function(mnreg, resps, data){
  # E[log p(y | X, \theta)] - KL[q(theta) | p (theta)]
  # in this case theta is a point estimate so just compute
  # E[log p(y | pi)] where expectations are over q(y)
  logpi <- compute_log_prior_assignment(mnreg, data)
  ll <- sum(resps * exp(logpi))
  return(ll)
}


# Another model------------


# SuSiE prior------------


