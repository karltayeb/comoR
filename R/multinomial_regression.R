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


#Multinomial/ NNET model------------
initialize_mnreg  <- function(mnreg_type,K, n,  p,param_nnet=list( size=1, decay=1)){


  tt <-rlang::exec( "nnet",
                    !!! param_nnet ,
                    y = matrix (1/K, nrow=n, ncol = K),
                    x =matrix (rnorm(n*p), nrow=n, ncol = p),
                    softmax=TRUE  ,
                    trace=FALSE )

  logpi <- log(tt$fitted.values )
  coef  <- tt
  mnreg  <- list(logpi=logpi,
                 K=K,
                 coef=coef,
                 param_nnet=param_nnet)
  class(mnreg) <- 'mult_reg'

  return(mnreg)
}
update_prior.mult_reg<- function(mnreg, resps, data   ){
  X  = as.matrix(data$X)
  tt <-rlang::exec( "nnet",
                    !!!mnreg$param_nnet ,
                    y = resps,
                    x = X,
                    softmax=TRUE  ,
                    trace=FALSE )
  mnreg$logpi <- log(tt$fitted.values)
  mnreg$coef  <- tt

  return(mnreg)
}
compute_log_prior_assignment.mult_reg <- function(mnreg, data){

  X <- data$X
  fitted_pi <-predict(mnreg$coef, X)# in case of multinomial cbind( 1,exp(tt))/(1+apply(exp(tt),1,sum))
  logpi <- log( fitted_pi)
  return(logpi)
}

# SuSiE prior------------






