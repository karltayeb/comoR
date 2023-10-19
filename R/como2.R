# Covariate moderated two groups model

#' @export
prep_data_como2 <- function(betahat, se, X, Z){
  data <- logisticsusie:::binsusie_prep_data(X, rep(0, length(betahat)), 1, Z)
  data$y <- NULL
  data$betahat <- betahat
  data$se <- se
  return(data)
}

compute_post_assignment <- function(fit, data) {
  # ln p(z=1)/p(z=0)
  prior_log_odds <- compute_prior_log_odds(fit$logreg, data)
  f0_loglik <- convolved_logpdf(fit$f0, data$betahat, data$se)+ fit$penalty
  f1_loglik <- convolved_logpdf(fit$f1, data$betahat, data$se)
  logits <- (f1_loglik - f0_loglik) + prior_log_odds
  post_assignment <- sigmoid(logits)
  return(post_assignment)
}

compute_assignment_entropy <- function(p) {
  q <- 1 - p
  h <- -1 * ((p * log(p)  + q * log(q)))
  return(sum(h, na.rm = T))
}

#' E[logp(y | f0, f1, z)]
#' @export
loglik.como2 <- function(fit, data) {
  f0_loglik <- fit$f0_loglik # convolved_logpdf(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- fit$f1_loglik # convolved_logpdf(fit$f1, fit$data$betahat, fit$data$se)
  q <- compute_post_assignment(fit, data)
  loglik <- (1 - q) * f0_loglik + q * f1_loglik
  return(loglik)
}

#' @export
compute_elbo.como2 <- function(fit, data) {
  post_assignment <- compute_post_assignment(fit, data)
  data_loglik <- sum(loglik.como2(fit, data))
  assignment_entropy <- compute_assignment_entropy(post_assignment)
  logreg_elbo <- compute_elbo(fit$logreg, post_assignment, data)
  if(is.null(logreg_elbo)){
    logreg_elbo <- -Inf
  }
  elbo <- (data_loglik + assignment_entropy) + logreg_elbo
  return(elbo)
}

#' @export
update_model.como2 <- function(fit, data, estimate_f1=FALSE, track_elbo=T){
  # compute posterior assignment
  prior_log_odds <- compute_prior_log_odds(fit$logreg, data)
  post_log_odds <- (fit$f1_loglik - fit$f0_loglik) + prior_log_odds
  qz <- sigmoid(post_log_odds)
  fit$logreg <- update_prior(fit$logreg, qz, data)

  if (estimate_f1) {
    fit$f1 <- update_params(
      fit$f1, data$betahat, data$se,
      weights = qz
    )
    fit$f1_loglik <- convolved_logpdf(fit$f1, data$betahat, data$se)
  }

  if(track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }
  return(fit)
}

#' @export
data_initialize_como2 <- function(data,
                                  f1_dist='normal',
                                  f1_params = list(),
                                  logreg='constant',
                                  logreg_params = list(),
                                  penalty=0.2){

  # initialize component distribution
  f0 <- point_component(mu = 0)

  f1 <- rlang::exec(
    paste0(f1_dist, "_component"), # constructor name
    !!!f1_params # unpack list of params
  )

  # fit f1 with full data to initialize
  f1 <- update_params(f1, data$betahat, data$se, rep(1, length(data$betahat)))

  # initialize posterior assignment
  f0_loglik <- convolved_logpdf( f0,
                                 data$betahat,
                                 data$se)
  f1_loglik <- convolved_logpdf( f1,
                                 data$betahat,
                                 data$se)
  logits <- f1_loglik - f0_loglik

  # initialize logreg
  if(logreg == 'constant'){
    logreg <- rlang::exec(initialize_constant_logreg, !!!logreg_params)
  }
  else if(logreg == 'linear_susie'){
    logreg_params$n <- length(f0_loglik)
    logreg <- rlang::exec(initialize_linear_susie, !!!logreg_params)
  }
  else if(logreg == 'logistic_ibss'){
    logreg_params$n <- length(f0_loglik)
    logreg <- rlang::exec(initialize_logistic_ibss, !!!logreg_params)
  }

  fit <-list(
    logreg = logreg,
    f0 = f0,
    f0_loglik = f0_loglik,
    f1 = f1,
    f1_loglik = f1_loglik,
    elbo = -Inf,
    penalty =penalty
  )
  class(fit) <- 'como2'
  return(fit)
}
