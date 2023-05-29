# Covariate moderated two groups model

#' @export
prep_data_como2 <- function(betahat, se, X, Z){
  data <- logisticsusie:::binsusie_prep_data(X, rep(0, length(betahat)), 1, Z)
  data$y <- NULL
  data$betahat <- betahat
  data$se <- se
  return(data)
}

compute_post_assignment <- function(fit) {
  # ln p(z=1)/p(z=0)
  logit_pi <- compute_Xb.binsusie(fit)
  f0_loglik <- convolved_logpdf(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- convolved_logpdf(fit$f1, fit$data$betahat, fit$data$se)
  logits <- f1_loglik - f0_loglik + logit_pi
  post_assignment <- sigmoid(logits)
  return(post_assignment)
}

compute_assignment_entropy <- function(fit, data) {
  p <- data$y
  entropy <- -1 * (p * log(p) + (1 - p) * log(1 - p))
  return(entropy)
}

#' E[logp(y | f0, f1, z)]
#'
loglik.como2 <- function(fit, data) {
  f0_loglik <- fit$f0_loglik # convolved_logpdf(fit$f0, fit$data$betahat, fit$data$se)
  f1_loglik <- fit$f1_loglik # convolved_logpdf(fit$f1, fit$data$betahat, fit$data$se)
  y <- data$y
  loglik <- (1 - y) * f0_loglik + y * f1_loglik
  return(loglik)
}

#' @export
compute_elbo.como2 <- function(fit, data) {
  data_loglik <- sum(loglik.como2(fit, data))
  assignment_entropy <- sum(compute_assignment_entropy(fit, data), na.rm = TRUE)
  logreg_elbo <- tail(fit$logreg$elbo, 1)
  if(is.null(logreg_elbo)){
    logreg_elbo <- -Inf
  }
  elbo <- data_loglik + assignment_entropy + logreg_elbo
  return(elbo)
}

#' @export
update_model.como2 <- function(fit, data, estimate_f1=FALSE, track_elbo=T){
  # compute posterior assignment
  fit$logits <- fit$f1_loglik - fit$f0_loglik + fit$psi$mu
  data$y <- sigmoid(fit$logits)

  # update logreg, store expected predictions in fit$psi$mu
  fit$logreg <- logisticsusie::update_model(fit$logreg, data, track_elbo=T)
  psi <- logisticsusie::compute_psi(fit$logreg, data)
  fit$psi$mu <- psi

  if (estimate_f1) {
    fit$f1 <- update_params(
      fit$f1, data$betahat, data$se,
      weights = data$y
    )
    fit$f1_loglik <- convolved_logpdf(fit$f1, data$betahat, data$se)
  }

  if(track_elbo){
    fit$elbo <- c(fit$elbo, compute_elbo(fit, data))
  }
  return(fit)
}

#' @export
data_initialize_como2 <- function(data, L, mu0=0, var0=1, f1_dist='normal', f1_params = list()){

  # initialize component distribution
  f0 <- point_component(mu = 0)

  f1 <- rlang::exec(
    paste0(f1_dist, "_component"), # constructor name
    !!!f1_params # unpack list of params
  )

  # initialize posterior assignment
  f0_loglik <- convolved_logpdf(f0, data$betahat, data$se)
  f1_loglik <- convolved_logpdf(f1, data$betahat, data$se)
  logits <- f1_loglik - f0_loglik

  # initialize logreg
  logreg <- logisticsusie:::data_initialize_binsusie(data, L, mu0, var0)

  fit <-list(
    logreg = logreg,
    f0 = f0,
    f0_loglik = f0_loglik,
    f1 = f1,
    f1_loglik = f1_loglik,
    psi = list(mu=0, mu2=0),
    elbo = -Inf
  )
  class(fit) <- 'como2'
  return(fit)
}
