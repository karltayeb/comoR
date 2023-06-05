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

#' @export
compute_log_prior_assignment.constant_mnreg <- function(mnreg, data){
  n <- nrow(data$X)
  logpi <- matrix(rep(mnreg$logpi, n), ncol = mnreg$K, byrow = T)
  return(logpi)
}

#' @export
update_prior.constant_mnreg <- function(mnreg, resps, data){
  new_log_pis <- log(colMeans(resps))
  mnreg$logpi <- new_log_pis
  return(mnreg)
}

#' @export
compute_elbo.constant_mnreg <- function(mnreg, resps, data){
  # E[log p(y | X, \theta)] - KL[q(theta) | p (theta)]
  # in this case theta is a point estimate so just compute
  # E[log p(y | pi)] where expectations are over q(y)

  logpi <- compute_log_prior_assignment(mnreg, data)
  ll <- sum(resps * exp(logpi))
  return(ll)
}



update_prior.mult_reg<- function(mnreg, resps, data){
  new_log_pis <- mult_reg (resps, X=data$X)$logpi
  mnreg$logpi <- new_log_pis
  mnreg$coef <- coef

  return(mnreg)
}

compute_log_prior_assignment.mult_reg <- function(mnreg, data){

  X <- data$X
  tt <-   cbind(rep(1,nrow(X)),X)%*%   coef
  fitted_pi <- exp(tt)/apply(exp(tt),1,sum)
  logpi <- log( fitted_pi)
  return(logpi)
}





mult_reg <-  function(assign_mat,X){

  coef <- do.call(cbind,
                  lapply(1:ncol(assign_mat),
                         function(k) reg_log(X=X, y=assign_mat[,k])$coef
                  )
  )

  tt <-   cbind(rep(1,nrow(X)),X)%*%   coef
  fitted_pi <- exp(tt)/apply(exp(tt),1,sum)
  out <- list( logpi=log(fitted_pi),
               coef=coef)
  return(out)
}
