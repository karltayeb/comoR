#' Compute Data Log Likelihood
#' Compute p(betahat_i| z=k, se_i) i = 1..., n, k=1... K
#' @title Compute data likelihood
#' @description Compute data likelihood
#' @param fit an MoCoCoMo fit object
#' @return n x K matrix of  log likelihood of each data point for each component
compute_data_loglikelihood <- function(fit,...){
  UseMethod("compute_data_loglikelihood")
}


compute_data_loglikelihood.default <- function(fit, data){
  data_loglik <- do.call(cbind,
                         purrr::map(
                           fit$f_list, ~ convolved_logpdf(.x,  data$betahat,  data$se)
                         )
  )
  return(data_loglik)
}
