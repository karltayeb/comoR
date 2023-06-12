# interface for covariate moderated prior weights

#' Compute the prior probability of assignment to each mixture component
#' @export
compute_log_prior_assignment <- function(x, ...){
  UseMethod("compute_log_prior_assignment", x)
}


#' Update the multinomial regression
#' @export
update_prior <- function(mnreg, ...){
  UseMethod("update_prior", mnreg )
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


#' @export

initialize_mnreg  <- function(mnreg_type,K, n,  p){
  logpi <- rep(log(1/K), K)
  mnreg  <- list(logpi=logpi, K=K)
  if(mnreg_type == 'constant'){
    class(mnreg ) <- 'constant_mnreg'
  }
  if(mnreg_type== 'mult_reg'){
    class(mnreg) <- 'mult_reg'
    mnreg$coef   <- matrix(0, ncol=K-1, nrow=(p+1))
  }
  if(mnreg_type== 'susie'){
    class(mnreg) <- 'mult_susie'
    mnreg$coef   <- matrix(0, ncol=K-1, nrow=(p+1))
  }


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


#' @export
update_prior.mult_reg<- function(mnreg, resps, data){
   X  = as.matrix(data$X)
   tt <- nnet (y = resps,
              x = X,
              size=1, decay=1,softmax=TRUE     )
  mnreg$logpi <- log(tt$fitted.values)
  mnreg$coef  <- tt

  return(mnreg)
}
#' @export
compute_log_prior_assignment.mult_reg <- function(mnreg, data){

  X <- as.matrix(data$X)


  if( is.matrix(mnreg$coef )){
    tt <-   cbind(rep(1,nrow(X)),X)%*%   mnreg$coef
    fitted_pi <- cbind( 1,exp(tt))/(1+apply(exp(tt),1,sum))
  }else{
    fitted_pi <-  predict(mnreg$coef,X)
  }

  logpi <- log( fitted_pi)
  return(logpi)
}




mult_reg_susie <-  function(assign_mat,X){



  #### TODO correct for susie calll ----
  coef <- do.call(cbind,
                  lapply(2:(ncol(assign_mat) ),
                         function(k) lm( log(assign_mat[,k]/assign_mat[,1])~X)$coefficients
                  )
  )


  tt <-   cbind(rep(1,nrow(X)),X)%*%   coef
  fitted_pi <- cbind( 1,exp(tt))/(1+apply(exp(tt),1,sum))

  out <- list( logpi=log(fitted_pi),
               coef=coef )
  return(out)
}

