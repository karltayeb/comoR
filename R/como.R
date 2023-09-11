# Implements covariate moderated ASH "MOre COmponents COvariate MOderated"

#' @title Function implementation the como mode
#' @details Function implementation the como mode
#'
#' @param data an object of class data_como  see \link{\code{set_data_como}}
#' @param modeltype of model currently supported (normal and beta )
#' @param maxiter numeric, maximum numerous of iteration set to 100 by defaults
#' @param tol tolerance in term of change in ELBO value for stopping criterion
#' @param upper, logical, set to FALSE by default. Specific to beta distribution.
#'  If true use a to set of mixture for fitting both end of the of the distribution as in the ZAP paper by Leung and Sunn
#' @parma nullweight  numeric value for penalizing likelihood at point mass 0/null component (should be larger than  1, 1 corresponds to no penalty , 2 corresponds to considering 1 individual "being null" and so on)
#' (usefull in small sample size)
#' @export
#' @example
#' #Simulate data under the como model
#' sim  <- sim_twococomo()
#' #preparing the data
#' data <- set_data_como(betahat = sim$betahat,
#'                                se = sim$se ,
#'                                 X = sim$X)
#' #fit como model
#' fit <- fit.como(data, maxiter=20)
#' plot(fit$elbo)
#' .monotone(fit$elbo)
#'
#' #get posterior quantities
#' est<- post_mean_sd.como (fit)
#' head(est)
#'  plot( est$mean, data$betahat)
#'
#' #comparison with ash
#'
#' t_ash <- ash(sim $betahat, sim $se, mixcompdist = "normal")
#' post_mean_ash <- t_ash$result$PosteriorMean
#' plot(est$mean, post_mean_ash)
#' # TODO make a more convincing example
#'
#'  sim  <- logisticsusie:::sim_como_beta(n=100)
#'#preparing the data
#'data <- set_data_como(p = sim$p,
#'                          X = sim$X)
#'
#' fit <- fit.como(data, maxiter=20)

initialize_como <- function(scales,
                            n,
                            p,
                            p2,
                            mu0=0,
                            var0=1,
                            nullweight=0,
                            mnreg='constant',
                            param_nnet =list( size=1, decay=1)
                            ){
  # initialize multinomial susie-- but could be any multinomial regression
  K <- length(scales)

  if(mnreg == 'constant'){
    mnreg <- initialize_constant_mnreg(K)
  }
  if( mnreg== "mult_reg"){

    mnreg <- initialize_mnreg (mnreg_type = mnreg_type,
                               K          = K,
                               n          = n,
                               p          = p,
                               param_nnet = param_nnet)

  }

  #mn_reg <- logisticsusie:::initialize_sbmn_susie(K, n, p, p2, L, mu0, var0)

  # initialize_scales
  f_list <- purrr::map(scales, ~ normal_component(mu = 0, var = .x^2))

  fit <- list(
    mnreg = mnreg, # multinomial regression function. takes X, returns pi
    f_list = f_list, # component distributions
    nullweight = nullweight, # penalty promoting the first component,
    K = K,
    elbo = -Inf
  )
  class(fit) <- c('como')
  return(fit)
}

#' Use data to autoselect scales
data_initialize_como <- function(data, max_class,
                                 scales=NULL,
                                 mu0=0,
                                 var0=1,
                                 nullweight=0,
                                 mnreg='constant',
                                 param_nnet =list( size=1, decay=1)) {
  como_check_data(data)

  if(is.null(scales)){
    scales <- autoselect_scales(data$betahat, data$se, max_class)
  }

  K <- length(scales) # K <= max_class
  p <- ncol(data$X)
  n <- nrow(data$X)
  p2 <- ncol(data$Z)

  fit <- initialize_como(scales=scales,
                         n=n,
                         p=p,
                         p2=p2,
                         mu0=mu0,
                         var0=var0,
                         nullweight,
                         mnreg=mnreg,
                         param_nnet=param_nnet)
  return(fit)
}

#' @export
update_model.como <- function(fit, data, update_assignment = T, update_logreg=T, fit_prior_variance=F, track_elbo=T){
  K <- fit$K

  # pre-compute data likelihood, if we haven't already
  # TODO: use digest::digest to hash the scales and recompute if scales change?
  if(is.null(fit$data_loglik)){
    fit$data_loglik <- compute_data_loglikelihood(fit, data)
  }

  # updates posterior assignments probabilities
  # these are the response variable for mn_regression
  if (update_assignment) {
    fit$post_assignment <- compute_posterior_assignment(fit, data)
    #data$Y <- fit$post_assignment
    #data$Nk <- logisticsusie:::stick_breaking_N(data$Y)
  }

  if (update_logreg) {
    fit$mnreg <- update_prior(fit$mnreg, resps=fit$post_assignment, data=data)
  }

  if (track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }
  return(fit)
}

#' @export
compute_elbo.como <- function(fit, data) {
  # E[log p(beta | y)] -- expected data likelihood
  ll <- sum(fit$post_assignment * fit$data_loglik)

  # Entropy term # E[-log q(y)]
  assignment_entropy <- sum(apply(fit$post_assignment, 1, logisticsusie:::categorical_entropy))

  # E[log p(y | X, theta)] - KL[q(theta) || p(theta)] SuSiE ELBO
  elbo <- compute_elbo(fit$mnreg, fit$post_assignment, data)

  # put it all together
  elbo <- ll - assignment_entropy + elbo
  return(elbo)
}



#' @export
fit.como <- function(x, ...){
  return(fit_model(x, ...))
}


get_KL.como <- function(fit,data){
  ll <- sum(fit$post_assignment * fit$data_loglik)

  # Entropy term # E[-log q(y)]
  assignment_entropy <- sum(apply(fit$post_assignment, 1, logisticsusie:::categorical_entropy))

  # E[log p(y | X, theta)] - KL[q(theta) || p(theta)] SuSiE ELBO

  # put it all together
  KL <- -ll + assignment_entropy
  return(KL)
}
