

#' Compute Posterior Assignment Probabilities
#' For each data point return posterior assignment probabilities
#' @param fit a MoCoCoMo fit object
#' @return an n x K matrix of log posterior probabilities for each data point
compute_posterior_assignment <- function(fit, data, log = FALSE) {
  data_loglik <- fit$data_loglik

  # TODO: generalize to other models
  assignment_loglik <- compute_log_prior_assignment(fit$mnreg, data)
  assignment_loglik[, 1] <- assignment_loglik[, 1] + fit$nullweight

  # normalize
  res <- do.call(
    rbind,
    apply(data_loglik + assignment_loglik, 1, function(x) x - logSumExp(x), simplify = F)
  )

  # exponentiate if log=FALSE
  if (!log) {
    res <- exp(res)
  }
  return(res)
}


#' @title Compute individual posterior variance from marginal normal mean model
#' @description internal function to compute posterior mean and sds
t_ind_var.mococomo <- function(fit, i) {
  do.call(
    c,
    lapply(
      1:length(fit$f_list),
      function(k) {
        1 / ((1 / fit$data$se[i]^2) + (1 / fit$f_list[[k]]$var))
      }
    )
  )
}


#' @title Compute individual posterior first and second moment
#' @description Compute individual posterior first and second moment
#'
#' # TODO: currently only for center prior
#' @param fit a mococomo object
#' @param  t_ind_var output of t_ind_var (using same mocomo object!)
#' @param i individual of interest
#' @exemple
#' t_post_var <-   do.call(rbind,
#'                        lapply( 1:length(fit$data$betahat),
#'                                function(i)t_ind_var.mococomo(fit, i)
#'                        )
#' )
#'
#'
#' post_beta <-     do.call(c, lapply( 1: length(fit$data$betahat), function(i)cal_ind_postmean(fit, t_post_var,i,) ))

cal_ind_moment12 <- function(fit, t_post_var, i) {
  temp <- do.call(
    c,
    lapply(
      1:length(fit$f_list),
      function(k) {
        (t_post_var[i, k] / (fit$data$se[i]^2) )*
          (fit$data$betahat[i])
      }
    )
  )

  ind_mean <- sum(fit$post_assignment[i, ] * temp)
  ind_second_moment <- sum(fit$post_assignment[i, ] * (t_post_var[i, ] + temp^2))

  ind_var <- ind_second_moment - ind_mean^2

  return(list(
    mean = ind_mean,
    sd = sqrt(ind_var)
  ))
}





#' @title Compute individual posterior mean and sd under mococomo model
#' @description Compute individual posterior mean and sd under mococomo model
#'
#' # TODO: allow using new observation from another data set (e.g. testing)
#' @param fit a mococomo object
#' @export
#' @example
#' see \link{\code{fit.mococomo}}
post_mean_sd.mococomo <- function(fit) {
  t_post_var <- do.call(
    rbind,
    lapply(
      1:length(fit$data$betahat),
      function(i) t_ind_var.mococomo(fit, i)
    )
  )
  out <- do.call(
    rbind,
    lapply(
      1:length(fit$data$betahat),
      function(i) cal_ind_moment12(fit, t_post_var, i)
    )
  )
  out <- data.frame(
    mean = do.call(c, out[, 1]),
    sd = do.call(c, out[, 2])
  ) # could be skip for speed



  return(out)
}


#' @title Compute individual fdr value  mococomo model with centered normal mixture
#' @descriptionCompute individual fdr value  mococomo model with centered normal mixture
#'
#' @param fit a mococomo object
#' @export
get_fdr <- function(fit) {
  tt1 <- fit$post_assignment[, 1] * dnorm(fit$data$betahat, mean = 0, sd = fit$data$se)
  tt2 <- Reduce("+", lapply(
    2:ncol(fit$post_assignment),
    function(k) {
      fit$post_assignment[, k] * dnorm(fit$data$betahat,
                                       mean = 0,
                                       sd = sqrt(fit$data$se^2 + fit$f_list[[k]]$var)
      )
    }
  ))
  out <- tt1 / (tt1 + tt2)

  return(out)
}





back_fit_component  <- function( fit, k, new_L){


  new_L <- max(c(1,new_L))

  if( new_L>=fit$logreg_list[[k]]$hypers$L){
    return(fit$logreg_list[[k]])
  }
  tt <- fit$logreg_list[[k]]
  tt$hypers$L  <- new_L

  if( new_L==1){
    tt$hypers$pi <-  matrix(tt$hypers$pi[1:new_L,]  , nrow = 1)
    tt$hypers$prior_mean <- tt$hypers$prior_mean[1:new_L]
    tt$hypers$prior_variance <- tt$hypers$prior_variance[1:new_L]

    tt$params$alpha <- matrix(tt$params$alpha[1:new_L,], nrow = 1)
    tt$params$mu    <- matrix(fit$logreg_list[[1]]$params$mu   [1:new_L,], nrow = 1)
    tt$params$var   <- matrix(tt$params$var  [1:new_L,], nrow = 1)
    tt$params$delta <- matrix(tt$params$delta[1:new_L,],ncol = 1)
  }else{
    tt$hypers$pi <-  tt$hypers$pi[1:new_L,]
    tt$hypers$prior_mean <- tt$hypers$prior_mean[1:new_L]
    tt$hypers$prior_variance <- tt$hypers$prior_variance[1:new_L]


    tt$params$alpha <- tt$params$alpha[1:new_L,]
    tt$params$mu    <- fit$logreg_list[[1]]$params$mu   [1:new_L,]
    tt$params$var   <- tt$params$var  [1:new_L,]
    tt$params$delta <- matrix(tt$params$delta[1:new_L,],ncol = 1)
  }

  return( tt)

}


backfit_effect_mococomo <- function(fit, l_dummy_cs){

  new_L  <- lapply(1:length(fit$logreg_list), function(k )
    max(c(1, fit$logreg_list[[k]]$hypers$L- length(l_dummy_cs[[k]])))
  )

 out <-  lapply(1:length(fit$logreg_list),
         function( k)
           back_fit_component(fit, k=k, new_L = new_L[[k]] )
  )
 return( out)
}




