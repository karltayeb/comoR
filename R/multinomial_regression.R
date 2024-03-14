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
update_prior.constant_mnreg <- function(mnreg, resps,loglik, data){
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
update_prior.mult_reg<- function(mnreg, resps, loglik, data   ){
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


#' @rdname compute_log_prior_assignment
#'
#' @method compute_log_prior_assignment
#'
#' @export compute_log_prior_assignment.mult_reg
#' @export
#' @keywords internal

compute_log_prior_assignment.mult_reg <- function(mnreg, data){
  # in case of multinomial cbind( 1,exp(tt))/(1+apply(exp(tt),1,sum))
  logpi <- mnreg$logpi
  return(logpi)
}

## Keras object
initialize_nnet_keras <-function (mnreg_type,K ,n , param_nnet,epoch, verbose=0 ){
  logpi <-  matrix (log(1/K), nrow=n, ncol = K)
  mnreg <- list(logpi=logpi,
                K=K,
                param_nnet= param_nnet,
                model= NULL,
                epoch=epoch,
                verbose=verbose)
  class(mnreg)="keras_obj"
  return(mnreg)

}
#' @rdname compute_log_prior_assignment
#'
#' @method compute_log_prior_assignment
#'
#' @export compute_log_prior_assignment.keras_obj
#' @export
#' @keywords internal

compute_log_prior_assignment.keras_obj <- function(mnreg, data){
 # in case of multinomial cbind( 1,exp(tt))/(1+apply(exp(tt),1,sum))
  logpi <- mnreg$logpi
  return(logpi)
}

#' @rdname update_prior
#'
#' @method update_prior
#'
#' @export update_prior.keras_obj
#' @importFrom keras clone_model
#' @export
#' @keywords internal
update_prior.keras_obj<- function(mnreg, resps,loglik, data   ){
  X  = as.matrix(data$X)


  model1 <- keras::clone_model(mnreg$param_nnet )
  if (!is.null(mnreg$model)){

    set_weights(model1, get_weights(mnreg$model))
  }


  model1 <-  model1 %>% compile(
    loss = custom_loss,
    optimizer = 'adam',
    metrics = c('accuracy')
  )

  x_train=X
  y_train= loglik
  candidate_batch_size =divisors(nrow(y_train))

  idx = which.min(abs( divisors(nrow(y_train))-100))
  custom_batch_size <- candidate_batch_size[idx]

  history <-model1 %>% fit(
    x_train, y_train,
    epochs = mnreg$epoch,
    batch_size = custom_batch_size,
      verbose  = mnreg$verbose
  )



  mnreg$logpi <- log(predict(model1, x_train))
  mnreg$model <- model1

  return(mnreg)
}


custom_loss <- function(y_true, y_pred) {
  case_wise_tt <- tf$vectorized_map(
    elems = c(y_true, y_pred),
    fn = function(x) {
      c(y_true1, y_pred1) %<-% x
      log(sum(exp(y_true1) * y_pred1))
    }
  )
   out <- -sum(case_wise_tt)

  #out <- out -tf$math$maximum(0,out) + tf$math$sqrt(out -tf$math$minimum(out,0)+1 )+1# doing the same thing as below
  # but tensorflow do not like logical operation in cus tom loss

  return(out)
}



